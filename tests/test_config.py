import pytest

from liposome_pipeline.config import load_config
from liposome_pipeline.errors import ConfigError


def write(tmp_path, text):
    path = tmp_path / "cfg.txt"
    path.write_text(text)
    return path


SYSTEM = """
system
    input_dir = {inp}
    output_dir = out
    forcefield = {inp}
    index_groups = LIPD:['DOPC', 'CHOL'] SOLV:['W']
    box_size = 10
    water_radius = 0.21
end
molecule_a
    name = DOPC
    nmol = 5
    insertion_radius = 0.3
end
"""


def test_parses_all_blocks(tmp_path):
    cfg = write(
        tmp_path,
        SYSTEM.format(inp=tmp_path)
        + "stage_em\n type = EM\n nsteps = 10\n tc-grps = A B\nend\n"
        + "slurm\n partition = gpu\n setup_commands = module purge; source rc\nend\n",
    )
    config = load_config(cfg)
    assert config.system.box_size == (10.0, 10.0, 10.0)
    assert config.system.index_groups == {"LIPD": ["DOPC", "CHOL"], "SOLV": ["W"]}
    assert config.molecules[0].count == 5
    assert config.stages[0].mdp_overrides == {"nsteps": "10", "tc-grps": ["A", "B"]}
    assert config.slurm.setup_commands == ("module purge", "source rc")


def test_legacy_keys_are_accepted(tmp_path):
    text = SYSTEM.format(inp=tmp_path).replace("index_groups", "idxGRPS").replace("box_size", "boxSize")
    assert load_config(write(tmp_path, text)).system.box_size[0] == 10.0


def test_missing_parameter_is_reported(tmp_path):
    text = SYSTEM.format(inp=tmp_path).replace("    water_radius = 0.21\n", "")
    with pytest.raises(ConfigError, match="water_radius"):
        load_config(write(tmp_path, text))


def test_unterminated_block(tmp_path):
    with pytest.raises(ConfigError, match="missing its 'end'"):
        load_config(write(tmp_path, "system\n input_dir = x\n"))


def test_unknown_block(tmp_path):
    with pytest.raises(ConfigError, match="unknown block"):
        load_config(write(tmp_path, "bogus\n a = 1\nend\n"))


def test_missing_file(tmp_path):
    with pytest.raises(ConfigError, match="not found"):
        load_config(tmp_path / "nope.txt")
