import shutil

import pytest

from liposome_pipeline.errors import ConfigError
from liposome_pipeline.gromacs import index_group_names
from liposome_pipeline.pipeline import Pipeline

CONFIG = """
system
    input_dir = {root}/input
    output_dir = {out}
    forcefield = {root}/input/martini
    index_groups = LIPD:['DOPC','CHOL'] SOLV:['W']
    box_size = 10 10 10
    water_radius = 0.21
end
molecule_dopc
    name = DOPC
    nmol = 30
    insertion_radius = 0.3
end
molecule_chol
    name = CHOL
    nmol = 10
    insertion_radius = 0.3
end
stage_em
    type = EM
    nsteps = 50
end
"""

pytestmark = pytest.mark.skipif(shutil.which("gmx") is None, reason="GROMACS not installed")


@pytest.mark.gromacs
def test_build_and_minimize(tmp_path, repo_root):
    out = tmp_path / "run"
    cfg = tmp_path / "c.txt"
    cfg.write_text(CONFIG.format(root=repo_root, out=out))

    final = Pipeline.from_file(cfg).run().structure
    assert final == out / "EM" / "EM.gro" and final.is_file()
    assert {"LIPD", "SOLV", "CORE"} <= set(index_group_names(out / "setup" / "index.ndx"))
    topology = (out / "setup" / "topol.top").read_text()
    assert "DOPC   30" in topology and "CHOL   10" in topology and "W " in topology


@pytest.mark.gromacs
def test_existing_output_requires_overwrite(tmp_path, repo_root):
    out = tmp_path / "run"
    out.mkdir()
    cfg = tmp_path / "c.txt"
    cfg.write_text(CONFIG.format(root=repo_root, out=out))
    with pytest.raises(ConfigError, match="overwrite"):
        Pipeline.from_file(cfg).run()


@pytest.mark.gromacs
def test_unknown_index_member_is_reported(tmp_path, repo_root):
    cfg = tmp_path / "c.txt"
    cfg.write_text(CONFIG.format(root=repo_root, out=tmp_path / "r").replace("['W']", "['W','XYZ']"))
    with pytest.raises(ConfigError, match="XYZ"):
        Pipeline.from_file(cfg).run()


def test_duplicate_stage_types_rejected(tmp_path, repo_root):
    cfg = tmp_path / "c.txt"
    cfg.write_text(CONFIG.format(root=repo_root, out=tmp_path / "r") + "stage_again\n type = EM\nend\n")
    with pytest.raises(ConfigError, match="unique"):
        Pipeline.from_file(cfg)


def test_example_configs_validate(repo_root, monkeypatch):
    monkeypatch.chdir(repo_root)
    for name in ("quickstart", "drug_loaded_liposome", "umbrella_sampling"):
        assert "stages" in Pipeline.from_file(repo_root / "examples" / f"{name}.txt").describe()
