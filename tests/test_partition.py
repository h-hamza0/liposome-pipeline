import MDAnalysis as mda
import numpy as np
import pytest

from liposome_pipeline.analysis import discovery, partition
from liposome_pipeline.analysis.cli import main

CENTRE = 175.0
OUTER_A = 80.0


def sphere(n, radius, rng=None):
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * k / n)
    theta = np.pi * (1 + 5**0.5) * k
    unit = np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)
    return CENTRE + radius * unit


def write_system(directory, n_in=30, n_out=10, n_frames=3):
    rng = np.random.default_rng(0)
    n_lipid = 200
    resnames = ["DOPC"] * n_lipid + ["DEXT"] * (n_in + n_out)
    n = len(resnames)
    u = mda.Universe.empty(n, n_residues=n, atom_resindex=np.arange(n), trajectory=True)
    u.add_TopologyAttr("name", ["NC3"] * n_lipid + ["D1"] * (n_in + n_out))
    u.add_TopologyAttr("resname", resnames)
    u.add_TopologyAttr("resid", np.arange(1, n + 1))
    lipids = sphere(n_lipid, OUTER_A, rng)
    inside = CENTRE + rng.uniform(-30, 30, size=(n_in, 3))
    outside = sphere(n_out, 150.0, rng)
    positions = np.vstack([lipids, inside, outside]).astype("float32")
    u.dimensions = [350, 350, 350, 90, 90, 90]

    directory.mkdir(parents=True)
    gro = directory / discovery.STRUCTURE_NAME
    xtc = directory / discovery.TRAJECTORY_NAME
    u.atoms.positions = positions
    u.atoms.write(str(gro))
    with mda.Writer(str(xtc), n) as writer:
        for _ in range(n_frames):
            u.atoms.positions = positions
            writer.write(u.atoms)
    return gro, xtc


def expected_ratio(n_in, n_out, radius_nm=8.0, box=35.0):
    v = 4 / 3 * np.pi * radius_nm**3
    return (n_in / v) / (n_out / (box**3 - v))


def test_analyze_partition_matches_analytic_ratio(tmp_path):
    gro, xtc = write_system(tmp_path / "52_47_1_298_LoadedDEXT" / "rep1")
    df = partition.analyze_partition(gro, xtc, "DEXT")
    assert len(df) == 3
    assert (df["n_in"] == 30).all() and (df["n_out"] == 10).all()
    assert df["outer_r_nm"].iloc[0] == pytest.approx(8.0, abs=0.05)
    assert df["partition_ratio"].iloc[0] == pytest.approx(expected_ratio(30, 10), rel=0.03)
    assert df["encapsulation_fraction"].iloc[0] == pytest.approx(0.75)


def test_missing_drug_raises(tmp_path):
    gro, xtc = write_system(tmp_path / "a")
    with pytest.raises(ValueError, match="resname IBU"):
        partition.analyze_partition(gro, xtc, "IBU")


def test_discovery_and_filtering(tmp_path):
    write_system(tmp_path / "52_47_1_298_LoadedDEXT" / "rep1", n_frames=1)
    write_system(tmp_path / "52_47_1_298" / "rep1", n_frames=1)
    write_system(tmp_path / "bad_name" / "rep1", n_frames=1)
    found = discovery.discover_systems(tmp_path)
    assert len(found) == 2
    loaded = discovery.drug_loaded(found)
    assert [(s.ratio, s.drug, s.replicate) for s in loaded] == [("52_47_1", "DEXT", "rep1")]


def test_parse_system_dir_rejects_garbage():
    with pytest.raises(ValueError):
        discovery.parse_system_dir("nonsense")


def test_summary_statistics():
    import pandas as pd

    def frame(values):
        return pd.DataFrame({"time_ps": range(len(values)), "partition_ratio": values})

    results = {("a", "X", "r1"): frame([0, 0, 10, 10]), ("a", "X", "r2"): frame([0, 0, 20, 20])}
    per_replicate, grouped = partition.summarize(results, 0.5)
    assert per_replicate["mean_partition_ratio"].tolist() == [10, 20]
    assert grouped["mean"].iloc[0] == 15


def test_cli_end_to_end(tmp_path, capsys):
    write_system(tmp_path / "data" / "52_47_1_298_LoadedDEXT" / "rep1")
    write_system(tmp_path / "data" / "52_47_1_298_LoadedDEXT" / "rep2")
    root = str(tmp_path / "data")
    cache, out = str(tmp_path / "cache"), str(tmp_path / "out")

    assert main(["count", root]) == 0
    assert capsys.readouterr().out.strip() == "2"
    for i in (0, 1):
        assert main(["compute", root, "--index", str(i), "--cache-dir", cache]) == 0
    assert main(["compute", root, "--index", "5", "--cache-dir", cache]) == 1
    assert main(["aggregate", "--cache-dir", cache, "--output-dir", out]) == 0

    produced = {p.name for p in (tmp_path / "out").iterdir()}
    assert {
        "partition_ratio_timeseries.pdf",
        "partition_ratio_bar.pdf",
        "partition_ratio_summary.csv",
        "partition_ratio_grouped.csv",
        "partition_ratio_timeseries.csv",
    } <= produced


def test_umbrella_geometry_helpers(tmp_path):
    from liposome_pipeline import umbrella

    gro, _ = write_system(tmp_path / "s")
    command = umbrella.select_drug_command(gro, "DEXT", 50)
    assert command.startswith("ri ") and "{next_id}" in command
    resid = int(command.split()[1])
    assert 201 <= resid <= 240

    ndx = tmp_path / "i.ndx"
    ndx.write_text("[ LIPD ]\n" + " ".join(str(i) for i in range(1, 201)) + "\n")
    atom = umbrella.nearest_to_center(gro, ndx)
    assert 1 <= atom <= 200
