from __future__ import annotations

import shutil
from importlib import resources
from pathlib import Path

from .errors import ConfigError

Mdp = dict[str, str]


def _canonical(key: str) -> str:
    return key.strip().lower().replace("-", "_")


def available_templates() -> list[str]:
    root = resources.files("liposome_pipeline") / "data" / "mdp"
    return sorted(p.name[: -len(".mdp")] for p in root.iterdir() if p.name.endswith(".mdp"))


def template_path(name: str) -> Path:
    root = resources.files("liposome_pipeline") / "data" / "mdp"
    candidate = root / f"{name}.mdp"
    if not candidate.is_file():
        raise ConfigError(
            f"no bundled MDP template named {name!r}; available: {', '.join(available_templates())}"
        )
    return Path(str(candidate))


def read_mdp(path: Path) -> Mdp:
    params: Mdp = {}
    for raw in Path(path).read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        params[key.strip()] = " ".join(value.split())
    return params


def apply_overrides(params: Mdp, overrides: dict[str, object]) -> Mdp:
    merged = dict(params)
    lookup = {_canonical(k): k for k in merged}
    for key, value in overrides.items():
        text = " ".join(str(v) for v in value) if isinstance(value, (list, tuple)) else str(value)
        merged[lookup.get(_canonical(key), key)] = text
    return merged


def write_mdp(params: Mdp, path: Path) -> None:
    width = max((len(k) for k in params), default=0)
    lines = [f"{key:<{width}} = {value}" for key, value in params.items()]
    Path(path).write_text("\n".join(lines) + "\n")


def materialize(source: Path, destination: Path, overrides: dict[str, object]) -> Path:
    if overrides:
        write_mdp(apply_overrides(read_mdp(source), overrides), destination)
    else:
        shutil.copy(source, destination)
    return destination
