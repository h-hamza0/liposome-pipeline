from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .config import MoleculeConfig


@dataclass
class Molecule:
    name: str
    count: int
    insertion_radius: float
    insert: bool | None = None
    itp_file: Path | None = None
    gro_file: Path | None = None

    @classmethod
    def from_config(cls, config: MoleculeConfig) -> "Molecule":
        return cls(config.name, config.count, config.insertion_radius, config.insert)

    @property
    def bead_count(self) -> int:
        if self.gro_file is None:
            raise ValueError(f"{self.name}: structure file has not been staged yet")
        return int(self.gro_file.read_text().splitlines()[1])
