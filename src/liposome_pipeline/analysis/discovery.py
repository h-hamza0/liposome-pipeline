from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

STRUCTURE_NAME = "aligned_centered.gro"
TRAJECTORY_NAME = "MD_CENTERED_NOJUMP.xtc"
DRUG_RESNAME = {"DEXT": "DEXT", "IBU": "IBU", "nfd": "NFD"}

_SYSTEM_DIR = re.compile(r"^(?P<ratio>[\d_]+)_(?P<temperature>\d+)(?:_Loaded(?P<drug>\w+))?$")


@dataclass(frozen=True)
class SystemEntry:
    ratio: str
    temperature: str
    drug: str | None
    replicate: str
    stage: str | None
    topology: Path
    trajectory: Path

    @property
    def key(self) -> tuple[str, str | None, str]:
        return (self.ratio, self.drug, self.replicate)

    @property
    def drug_resname(self) -> str:
        if self.drug is None:
            raise ValueError("drug-free system has no drug residue name")
        return DRUG_RESNAME[self.drug]


def parse_system_dir(name: str) -> tuple[str, str, str | None]:
    match = _SYSTEM_DIR.match(name)
    if match is None:
        raise ValueError(f"unrecognized system directory name: {name!r}")
    return match["ratio"], match["temperature"], match["drug"]


def discover_systems(roots: Iterable[str | Path] | str | Path = ".") -> list[SystemEntry]:
    """Find ``<ratio>_<T>[_Loaded<drug>]/<replicate>[/<stage>]/aligned_centered.gro`` under each root."""
    if isinstance(roots, (str, Path)):
        roots = [roots]

    systems = []
    for root in map(Path, roots):
        for structure in sorted(root.glob(f"**/{STRUCTURE_NAME}")):
            parts = structure.relative_to(root).parts[:-1]
            if len(parts) < 2:
                continue
            try:
                ratio, temperature, drug = parse_system_dir(parts[0])
            except ValueError:
                continue
            if drug is not None and drug not in DRUG_RESNAME:
                continue
            systems.append(
                SystemEntry(
                    ratio=ratio,
                    temperature=temperature,
                    drug=drug,
                    replicate=parts[1],
                    stage=parts[2] if len(parts) > 2 else None,
                    topology=structure,
                    trajectory=structure.with_name(TRAJECTORY_NAME),
                )
            )
    return systems


def drug_loaded(systems: Iterable[SystemEntry]) -> list[SystemEntry]:
    return [s for s in systems if s.drug is not None]
