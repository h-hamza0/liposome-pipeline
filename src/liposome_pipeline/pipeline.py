from __future__ import annotations

import logging
from pathlib import Path

from .config import PipelineConfig, load_config
from .errors import ConfigError
from .gromacs import Gromacs
from .molecules import Molecule
from .stages import Artifacts, Stage
from .system import System

log = logging.getLogger(__name__)


class Pipeline:
    def __init__(self, config: PipelineConfig, gmx: Gromacs) -> None:
        names = [s.type for s in config.stages]
        duplicates = sorted({n for n in names if names.count(n) > 1})
        if duplicates:
            raise ConfigError(f"stage types must be unique (each gets its own directory): {', '.join(duplicates)}")
        self.config = config
        self.gmx = gmx
        self.system = System(config.system, [Molecule.from_config(m) for m in config.molecules], gmx)
        self.stages = [Stage(s, gmx, config.slurm) for s in config.stages]

    @classmethod
    def from_file(cls, path: str | Path, gmx: Gromacs | None = None) -> "Pipeline":
        config = load_config(path)
        return cls(config, gmx or Gromacs(config.system.gmx_path or "gmx"))

    def use_gromacs(self, gmx: Gromacs) -> None:
        self.gmx = gmx
        self.system.gmx = gmx
        for stage in self.stages:
            stage.gmx = gmx

    def describe(self) -> str:
        c = self.config.system
        lines = [
            f"output   : {c.output_dir}",
            f"box (nm) : {' x '.join(f'{v:g}' for v in c.box_size)}",
            "molecules: " + ", ".join(f"{m.name} x{m.count}" for m in self.config.molecules),
            "stages   : " + (" -> ".join(s.type for s in self.config.stages) or "(setup only)"),
        ]
        return "\n".join(lines)

    def run(self, overwrite: bool = False) -> Artifacts:
        self.gmx.ensure_available()
        structure = self.system.setup(overwrite=overwrite)
        for stage in self.stages:
            stage.prepare(self.system.output_dir / stage.name, self.system.topology, self.system.index)

        state = Artifacts(structure=structure)
        for stage in self.stages:
            state = stage.run(state)
        return state
