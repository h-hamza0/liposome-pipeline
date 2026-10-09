from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .discovery import SystemEntry

log = logging.getLogger(__name__)

DEFAULT_BOX_SIDE_NM = 35.0
DEFAULT_LIPID_SELECTION = "resname DOPC PEL CHEMS and not (name NH O1)"
OUTER_RADIUS_PERCENTILE = 90.0
ANGSTROM_PER_NM = 10.0
PROGRESS_EVERY = 50
Results = dict[tuple[str, str, str], pd.DataFrame]


def analyze_partition(
    topology: str | Path,
    trajectory: str | Path,
    drug_resname: str,
    lipid_selection: str = DEFAULT_LIPID_SELECTION,
    box_side_nm: float = DEFAULT_BOX_SIDE_NM,
    stride: int = 1,
) -> pd.DataFrame:
    """Per-frame drug partition ratio ``C_lipid / C_bulk`` for one liposome trajectory.

    The vesicle radius is the 90th percentile of lipid centre-of-geometry distances from the vesicle
    centre. Drugs within that radius count as encapsulated; the remainder are in the bulk volume
    ``box_side_nm**3 - 4/3 pi r^3``.
    """
    import MDAnalysis as mda

    universe = mda.Universe(str(topology), str(trajectory))
    lipids = universe.select_atoms(lipid_selection)
    drugs = universe.select_atoms(f"resname {drug_resname}")
    if len(lipids) == 0:
        raise ValueError(f"no lipids selected in {topology}")
    if len(drugs) == 0:
        raise ValueError(f"no atoms with resname {drug_resname} in {topology}")

    n_drugs = len(drugs.residues)
    box_volume = box_side_nm**3
    frames = universe.trajectory[::stride]
    rows = []
    for i, ts in enumerate(frames):
        if i % PROGRESS_EVERY == 0 or i == len(frames) - 1:
            log.info("frame %d/%d (t = %.1f ns)", i + 1, len(frames), ts.time / 1000.0)

        centre = lipids.center_of_geometry()
        lipid_radii = _radii_nm(lipids.center_of_geometry(compound="residues"), centre)
        outer_radius = float(np.percentile(lipid_radii, OUTER_RADIUS_PERCENTILE))
        drug_radii = _radii_nm(drugs.center_of_geometry(compound="residues"), centre)

        n_in = int(np.sum(drug_radii <= outer_radius))
        n_out = n_drugs - n_in
        vesicle_volume = 4.0 / 3.0 * np.pi * outer_radius**3
        bulk_volume = box_volume - vesicle_volume
        c_in = n_in / vesicle_volume
        c_out = n_out / bulk_volume if bulk_volume > 0 else np.nan
        rows.append(
            {
                "time_ps": ts.time,
                "outer_r_nm": outer_radius,
                "n_in": n_in,
                "n_out": n_out,
                "c_in_per_nm3": c_in,
                "c_out_per_nm3": c_out,
                "partition_ratio": c_in / c_out if c_out and c_out > 0 else np.nan,
                "encapsulation_fraction": n_in / n_drugs,
            }
        )
    return pd.DataFrame(rows)


def _radii_nm(positions: np.ndarray, centre: np.ndarray) -> np.ndarray:
    return np.linalg.norm(positions - centre, axis=1) / ANGSTROM_PER_NM


def cache_path(cache_dir: str | Path, entry: SystemEntry) -> Path:
    return Path(cache_dir) / f"{entry.ratio}__{entry.drug}__{entry.replicate}.csv"


def compute_one(
    entry: SystemEntry,
    cache_dir: str | Path,
    force: bool = False,
    **kwargs: object,
) -> Path:
    out = cache_path(cache_dir, entry)
    if out.exists() and not force:
        log.info("cached result exists, skipping: %s", out)
        return out
    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    frame = analyze_partition(entry.topology, entry.trajectory, entry.drug_resname, **kwargs)  # type: ignore[arg-type]
    frame.insert(0, "replicate", entry.replicate)
    frame.insert(0, "drug", entry.drug)
    frame.insert(0, "ratio", entry.ratio)
    frame.to_csv(out, index=False)
    return out


def compute_all(entries: list[SystemEntry], **kwargs: object) -> Results:
    results: Results = {}
    for i, entry in enumerate(entries, start=1):
        log.info("[%d/%d] ratio=%s drug=%s replicate=%s", i, len(entries), *entry.key)
        results[entry.key] = analyze_partition(  # type: ignore[index, arg-type]
            entry.topology, entry.trajectory, entry.drug_resname, **kwargs
        )
    return results


def load_cache(cache_dir: str | Path) -> Results:
    paths = sorted(Path(cache_dir).glob("*.csv"))
    if not paths:
        raise FileNotFoundError(f"no cached CSVs found in {cache_dir}")
    results: Results = {}
    for path in paths:
        frame = pd.read_csv(path, dtype={"ratio": str, "drug": str, "replicate": str})
        key = (frame["ratio"].iloc[0], frame["drug"].iloc[0], frame["replicate"].iloc[0])
        results[key] = frame.drop(columns=["ratio", "drug", "replicate"])
    return results


def smooth(series: pd.Series, window: int) -> pd.Series:
    return series.rolling(window, center=True, min_periods=1).mean()


def timeseries_table(results: Results, window: int) -> pd.DataFrame:
    frames = []
    for (ratio, drug, replicate), df in sorted(results.items()):
        frames.append(
            pd.DataFrame(
                {
                    "ratio": ratio,
                    "drug": drug,
                    "replicate": replicate,
                    "time_ns": df["time_ps"].to_numpy() / 1000.0,
                    "partition_ratio_raw": df["partition_ratio"].to_numpy(),
                    "partition_ratio_smoothed": smooth(df["partition_ratio"], window).to_numpy(),
                }
            )
        )
    return pd.concat(frames, ignore_index=True)


def summarize(results: Results, tail_fraction: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Mean partition ratio over the trailing ``tail_fraction`` of each trajectory, then across replicates."""
    rows = []
    for (ratio, drug, replicate), df in sorted(results.items()):
        n_keep = max(1, int(len(df) * tail_fraction))
        rows.append(
            {
                "ratio": ratio,
                "drug": drug,
                "replicate": replicate,
                "mean_partition_ratio": df["partition_ratio"].iloc[-n_keep:].mean(),
            }
        )
    per_replicate = pd.DataFrame(rows)
    grouped = (
        per_replicate.groupby(["ratio", "drug"])["mean_partition_ratio"].agg(["mean", "std"]).reset_index()
    )
    return per_replicate, grouped
