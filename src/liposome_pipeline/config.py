from __future__ import annotations

import ast
import re
from dataclasses import dataclass, field
from pathlib import Path

from .errors import ConfigError

SYSTEM_ALIASES = {
    "idxGRPS": "index_groups",
    "boxSize": "box_size",
    "boxDist": "center_box",
    "Wradius": "water_radius",
    "US": "restart_from",
}
STAGE_ALIASES = {"US": "umbrella", "mdp_fname": "mdp_file"}
STAGE_RESERVED = {"type", "umbrella", "mdp_file", "lipid", "drug"}
TRUE_VALUES = {"1", "true", "yes", "on"}
FALSE_VALUES = {"0", "false", "no", "off"}
GROUP_PATTERN = re.compile(r"([\w-]+)\s*:\s*(\[[^\]]*\])")


@dataclass
class SystemConfig:
    input_dir: Path
    output_dir: Path
    forcefield: Path
    index_groups: dict[str, list[str]]
    box_size: tuple[float, float, float]
    water_radius: float
    center_box: float = 0.0
    restart_from: Path | None = None
    drug: str | None = None
    percentile: float | None = None
    gmx_path: str | None = None


@dataclass
class MoleculeConfig:
    name: str
    count: int
    insertion_radius: float
    insert: bool | None = None


@dataclass
class StageConfig:
    type: str
    umbrella: bool = False
    mdp_file: Path | None = None
    lipid: str | None = None
    drug: str | None = None
    mdp_overrides: dict[str, object] = field(default_factory=dict)


@dataclass
class SlurmConfig:
    partition: str | None = None
    cpus_per_task: int = 12
    mem: str = "32gb"
    gres: str | None = "gpu:1"
    mail_user: str | None = None
    mail_type: str = "END,FAIL"
    setup_commands: tuple[str, ...] = ()
    frame_jobs: int = 5
    window_jobs: int = 15
    poll_seconds: float = 30.0
    timeout_seconds: float | None = None


@dataclass
class PipelineConfig:
    system: SystemConfig
    molecules: list[MoleculeConfig]
    stages: list[StageConfig]
    slurm: SlurmConfig = field(default_factory=SlurmConfig)


def _is_none(value: str) -> bool:
    return value.strip().lower() in {"", "none", "null"}


def _scalar(value: str) -> str | list[str]:
    tokens = value.split()
    return tokens if len(tokens) > 1 else value.strip()


def _as_bool(key: str, value: object) -> bool:
    text = str(value).strip().lower()
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES or _is_none(text):
        return False
    raise ConfigError(f"{key}: expected a boolean, got {value!r}")


def _as_float(key: str, value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ConfigError(f"{key}: expected a number, got {value!r}") from None


def _as_int(key: str, value: object) -> int:
    number = _as_float(key, value)
    if number != int(number):
        raise ConfigError(f"{key}: expected an integer, got {value!r}")
    return int(number)


def _parse_index_groups(raw: str) -> dict[str, list[str]]:
    groups: dict[str, list[str]] = {}
    for name, members in GROUP_PATTERN.findall(raw):
        try:
            parsed = ast.literal_eval(members)
        except (ValueError, SyntaxError):
            raise ConfigError(f"index_groups: cannot parse members of {name!r}: {members}") from None
        groups[name] = [str(m) for m in parsed]
    if not groups:
        raise ConfigError(f"index_groups: no groups found in {raw!r} (expected NAME:['A','B'] ...)")
    return groups


def _require(block: dict[str, object], keys: list[str], label: str) -> None:
    missing = [k for k in keys if k not in block]
    if missing:
        raise ConfigError(f"{label}: missing required parameter(s): {', '.join(missing)}")


def _read_blocks(text: str) -> list[tuple[str, dict[str, str], int]]:
    blocks: list[tuple[str, dict[str, str], int]] = []
    header: str | None = None
    body: dict[str, str] = {}
    start = 0
    for number, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if header is None:
            header, body, start = line, {}, number
        elif line == "end":
            blocks.append((header, body, start))
            header = None
        elif "=" in line:
            key, value = line.split("=", 1)
            body[key.strip()] = value.strip()
        else:
            raise ConfigError(f"line {number}: expected 'key = value' or 'end', got {line!r}")
    if header is not None:
        raise ConfigError(f"block {header!r} starting at line {start} is missing its 'end'")
    return blocks


def _build_system(raw: dict[str, str], base: Path) -> SystemConfig:
    block = {SYSTEM_ALIASES.get(k, k): v for k, v in raw.items()}
    _require(block, ["input_dir", "output_dir", "forcefield", "index_groups", "box_size", "water_radius"], "system")

    box = block["box_size"].split()
    if len(box) == 1:
        box = box * 3
    if len(box) != 3:
        raise ConfigError(f"box_size: expected 1 or 3 values, got {block['box_size']!r}")

    def resolve(value: str) -> Path:
        path = Path(value).expanduser()
        return path if path.is_absolute() else base / path

    input_dir = resolve(block["input_dir"])
    if not input_dir.is_dir():
        raise ConfigError(f"input_dir does not exist: {input_dir}")

    restart = block.get("restart_from", "None")
    drug = block.get("drug", "None")
    percentile = block.get("percentile", "None")
    gmx_path = block.get("gmx_path", "None")
    return SystemConfig(
        input_dir=input_dir,
        output_dir=resolve(block["output_dir"]),
        forcefield=resolve(block["forcefield"]),
        index_groups=_parse_index_groups(block["index_groups"]),
        box_size=tuple(_as_float("box_size", v) for v in box),  # type: ignore[arg-type]
        water_radius=_as_float("water_radius", block["water_radius"]),
        center_box=_as_float("center_box", block.get("center_box", "0")),
        restart_from=None if _is_none(restart) else resolve(restart),
        drug=None if _is_none(drug) else drug,
        percentile=None if _is_none(percentile) else _as_float("percentile", percentile),
        gmx_path=None if _is_none(gmx_path) else gmx_path,
    )


def _build_molecule(raw: dict[str, str], label: str) -> MoleculeConfig:
    _require(raw, ["name", "nmol", "insertion_radius"], label)
    return MoleculeConfig(
        name=raw["name"],
        count=_as_int("nmol", raw["nmol"]),
        insertion_radius=_as_float("insertion_radius", raw["insertion_radius"]),
        insert=_as_bool("insert", raw["insert"]) if "insert" in raw else None,
    )


def _build_stage(raw: dict[str, str], base: Path, label: str) -> StageConfig:
    block = {STAGE_ALIASES.get(k, k): v for k, v in raw.items()}
    _require(block, ["type"], label)
    mdp_file = block.get("mdp_file")
    umbrella = block.get("umbrella", "no")
    return StageConfig(
        type=block["type"],
        umbrella=False if _is_none(umbrella) else _as_bool("umbrella", umbrella),
        mdp_file=None if mdp_file is None else (Path(mdp_file) if Path(mdp_file).is_absolute() else base / mdp_file),
        lipid=block.get("lipid"),
        drug=block.get("drug"),
        mdp_overrides={k: _scalar(v) for k, v in block.items() if k not in STAGE_RESERVED},
    )


def _build_slurm(raw: dict[str, str]) -> SlurmConfig:
    config = SlurmConfig()
    for key, value in raw.items():
        if key == "setup_commands":
            config.setup_commands = tuple(c.strip() for c in value.split(";") if c.strip())
        elif key in {"cpus_per_task", "frame_jobs", "window_jobs"}:
            setattr(config, key, _as_int(key, value))
        elif key in {"poll_seconds", "timeout_seconds"}:
            setattr(config, key, None if _is_none(value) else _as_float(key, value))
        elif key in {"partition", "mem", "gres", "mail_user", "mail_type"}:
            setattr(config, key, None if _is_none(value) else value)
        else:
            raise ConfigError(f"slurm: unknown parameter {key!r}")
    return config


def load_config(path: str | Path) -> PipelineConfig:
    path = Path(path)
    if not path.is_file():
        raise ConfigError(f"configuration file not found: {path}")
    base = Path.cwd()

    system: SystemConfig | None = None
    molecules: list[MoleculeConfig] = []
    stages: list[StageConfig] = []
    slurm = SlurmConfig()

    for header, body, line in _read_blocks(path.read_text()):
        label = f"{header} (line {line})"
        kind = header.split("_", 1)[0]
        if header == "system":
            if system is not None:
                raise ConfigError("more than one 'system' block")
            system = _build_system(body, base)
        elif header == "slurm":
            slurm = _build_slurm(body)
        elif kind == "molecule":
            molecules.append(_build_molecule(body, label))
        elif kind == "stage":
            stages.append(_build_stage(body, base, label))
        else:
            raise ConfigError(f"unknown block {header!r} at line {line}")

    if system is None:
        raise ConfigError("no 'system' block found")
    if not molecules:
        raise ConfigError("no 'molecule_*' blocks found")
    return PipelineConfig(system=system, molecules=molecules, stages=stages, slurm=slurm)
