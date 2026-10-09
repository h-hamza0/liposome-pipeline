import pytest

from liposome_pipeline import mdp
from liposome_pipeline.errors import ConfigError


def test_all_bundled_templates_parse():
    names = mdp.available_templates()
    assert {"EM", "NVT", "MD", "PULL"} <= set(names)
    for name in names:
        assert mdp.read_mdp(mdp.template_path(name))


def test_unknown_template():
    with pytest.raises(ConfigError):
        mdp.template_path("DOES_NOT_EXIST")


def test_inline_comments_do_not_drop_parameters(tmp_path):
    src = tmp_path / "a.mdp"
    src.write_text("; header\nnsteps = 5 ; steps\ntc-grps = A B\n")
    assert mdp.read_mdp(src) == {"nsteps": "5", "tc-grps": "A B"}


def test_overrides_match_hyphen_and_underscore_spellings():
    merged = mdp.apply_overrides({"tc-grps": "A B", "nsteps": "1"}, {"tc_grps": ["X", "Y"], "dt": 0.02})
    assert merged == {"tc-grps": "X Y", "nsteps": "1", "dt": "0.02"}


def test_materialize_roundtrip(tmp_path):
    out = mdp.materialize(mdp.template_path("EM"), tmp_path / "EM.mdp", {"nsteps": 7})
    assert mdp.read_mdp(out)["nsteps"] == "7"
