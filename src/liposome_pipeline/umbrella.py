from __future__ import annotations

import math
from pathlib import Path
from typing import Sequence, TypeVar

import numpy as np

from .gromacs import read_index_groups

T = TypeVar("T")

LIPID_RESNAMES = "DOPC PEL CHEMS"
DRUG_SEARCH_RADIUS_NM = 9.0
MAX_WINDOW_DISTANCE_NM = 11.0
DENSE_REGIONS = ((2.0, 4.0, 0.05),)
SPARSE_SPACING_NM = 0.05


def chunk_evenly(items: Sequence[T], parts: int) -> list[list[T]]:
    size, extra = divmod(len(items), parts)
    chunks = [
        list(items[i * size + min(i, extra) : (i + 1) * size + min(i + 1, extra)]) for i in range(parts)
    ]
    return [c for c in chunks if c]


def target_distances(
    min_dist: float,
    max_dist: float,
    dense_regions: Sequence[tuple[float, float, float]] = DENSE_REGIONS,
    sparse_spacing: float = SPARSE_SPACING_NM,
) -> list[float]:
    targets: list[float] = []
    for start, end, spacing in dense_regions:
        targets += list(np.arange(start, end + spacing, spacing))

    current = math.ceil(min_dist * 10) / 10
    ceiling = math.floor(max_dist * 10) / 10
    for start, end in sorted((s, e) for s, e, _ in dense_regions):
        while current < start:
            targets.append(current)
            current += sparse_spacing
        current = end + sparse_spacing
    while current <= ceiling:
        targets.append(current)
        current += sparse_spacing
    return sorted(set(np.round(targets, 4).tolist()))


def select_frames(
    times: np.ndarray,
    distances: np.ndarray,
    max_distance: float = MAX_WINDOW_DISTANCE_NM,
) -> list[tuple[float, float]]:
    selected: list[tuple[float, float]] = []
    used: set[int] = set()
    for target in target_distances(float(distances.min()), float(distances.max())):
        index = int(np.abs(distances - target).argmin())
        if index not in used:
            used.add(index)
            selected.append((target, float(times[index])))
    return [(d, t) for d, t in selected if d < max_distance]


def read_xvg(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = [
        line.split()[:2]
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.startswith(("#", "@"))
    ]
    data = np.array(rows, dtype=float)
    return data[:, 0], data[:, 1]


def format_distance(value: float) -> str:
    return f"{value:g}"


def frame_extraction_script(
    trajectory: Path, tpr: Path, index: Path, frames: Sequence[tuple[float, float]], out_dir: Path, gmx: str = "gmx"
) -> str:
    lines = []
    for distance, time in frames:
        target = out_dir / f"{format_distance(distance)}.gro"
        lines.append(
            f'echo 0 | {gmx} trjconv -s {tpr} -f {trajectory} -dump {time:g} -o {target} -n {index}'
        )
    return "\n".join(lines)


def window_script(
    frames: Sequence[Path], mdp: Path, topology: Path, index: Path, minimization: bool, gmx: str = "gmx"
) -> str:
    mdrun_flags = "-ntmpi 1 -nb gpu -pin off -v" + ("" if minimization else " -bonded gpu")
    quoted = " ".join(f'"{f}"' for f in frames)
    return f"""
for frame in {quoted}; do
    name=$(basename "$frame" .gro)
    {gmx} grompp -f {mdp} -c "$frame" -r "$frame" -p {topology} -n {index} -o "$name.tpr" -maxwarn 10 || exit 1
    {gmx} mdrun -deffnm "$name" {mdrun_flags} || exit 1
done
"""


def nearest_to_center(structure: Path, index_file: Path, group: str = "LIPD") -> int:
    import MDAnalysis as mda

    members = read_index_groups(index_file)[group]
    universe = mda.Universe(str(structure))
    atoms = universe.atoms[np.array(members) - 1]
    centre = atoms.center_of_mass()
    nearest = int(np.linalg.norm(atoms.positions - centre, axis=1).argmin())
    return int(atoms[nearest].index) + 1


def select_drug_command(structure: Path, drug: str, percentile: float, group_name: str = "SELECT") -> str:
    import MDAnalysis as mda

    universe = mda.Universe(str(structure))
    lipids = universe.select_atoms(f"resname {LIPID_RESNAMES}")
    atoms = universe.select_atoms(f"resname {drug}")
    if len(lipids) == 0 or len(atoms) == 0:
        raise ValueError(f"cannot select a {drug} molecule: lipids or drug missing from {structure}")
    centre = lipids.center_of_mass()
    distances = np.linalg.norm(atoms.positions - centre, axis=1) / 10.0
    near = np.flatnonzero(distances < DRUG_SEARCH_RADIUS_NM)
    if near.size == 0:
        raise ValueError(f"no {drug} atoms within {DRUG_SEARCH_RADIUS_NM} nm of the lipid centre")
    ranked = near[np.argsort(distances[near])]
    position = min(math.ceil(percentile / 100 * ranked.size), ranked.size - 1)
    resid = int(atoms[int(ranked[position])].resid)
    return f"ri {resid}\nname {{next_id}} {group_name}\n"
