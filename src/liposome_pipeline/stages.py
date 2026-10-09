from __future__ import annotations

import logging
from dataclasses import dataclass, field, replace
from pathlib import Path

from . import mdp, slurm, umbrella
from .config import SlurmConfig, StageConfig
from .errors import ConfigError, GromacsError
from .gromacs import Gromacs

log = logging.getLogger(__name__)

CONFIG_STAGE = "CONFIG"
MINIMIZATION_STAGES = frozenset({"EM", "EM_US"})
MAX_RESUMES = 3


@dataclass
class Artifacts:
    structure: Path | None = None
    trajectory: Path | None = None
    tpr: Path | None = None
    windows: list[Path] = field(default_factory=list)


class Stage:
    def __init__(self, config: StageConfig, gmx: Gromacs, slurm_config: SlurmConfig) -> None:
        self.config = config
        self.gmx = gmx
        self.slurm = slurm_config
        self.out_dir: Path | None = None
        self.topology: Path | None = None
        self.index: Path | None = None
        self.template: Path | None = None

    @property
    def name(self) -> str:
        return self.config.type

    def prepare(self, out_dir: Path, topology: Path, index: Path) -> None:
        self.out_dir = out_dir
        self.topology = topology
        self.index = index
        out_dir.mkdir(parents=True)
        if self.name == CONFIG_STAGE:
            return
        self.template = self.config.mdp_file or mdp.template_path(self.name)
        if not self.template.is_file():
            raise ConfigError(f"MDP file for stage {self.name} not found: {self.template}")

    def run(self, state: Artifacts) -> Artifacts:
        log.info("stage %s", self.name)
        if self.name == CONFIG_STAGE:
            return self._generate_windows(state)
        if self.config.umbrella:
            return self._run_windows(state)
        return self._run_single(state)

    def _write_mdp(self, extra: dict[str, object] | None = None) -> Path:
        overrides = {**self.config.mdp_overrides, **(extra or {})}
        return mdp.materialize(self.template, self.out_dir / f"{self.name}.mdp", overrides)

    def _run_single(self, state: Artifacts) -> Artifacts:
        if state.structure is None:
            raise GromacsError(f"stage {self.name} has no input structure")
        tpr = self.out_dir / f"{self.name}.tpr"
        self.gmx.grompp(self._write_mdp(), state.structure, self.topology, tpr, self.index)

        gro = self.out_dir / f"{self.name}.gro"
        self.gmx.mdrun(self.name, cwd=self.out_dir)
        for attempt in range(MAX_RESUMES):
            if gro.exists():
                break
            log.warning("%s did not finish; resuming from checkpoint (%d/%d)", self.name, attempt + 1, MAX_RESUMES)
            self.gmx.mdrun(self.name, cwd=self.out_dir, resume=True)
        if not gro.exists():
            raise GromacsError(f"stage {self.name} produced no final structure after {MAX_RESUMES} resume attempts")
        return Artifacts(structure=gro, trajectory=self.out_dir / f"{self.name}.xtc", tpr=tpr)

    def _generate_windows(self, state: Artifacts) -> Artifacts:
        if state.trajectory is None or state.tpr is None:
            raise ConfigError("CONFIG stage needs the trajectory and run input of the preceding stage")
        if not (self.config.lipid and self.config.drug):
            raise ConfigError("CONFIG stage requires 'lipid' and 'drug' index group names")

        xvg = self.out_dir / "dist.xvg"
        self.gmx.distance(state.trajectory, state.tpr, self.index, self.config.lipid, self.config.drug, xvg)
        times, distances = umbrella.read_xvg(xvg)
        frames = umbrella.select_frames(times, distances)

        windows = [self.out_dir / f"{umbrella.format_distance(d)}.gro" for d, _ in frames]
        for chunk in umbrella.chunk_evenly(frames, self.slurm.frame_jobs):
            body = umbrella.frame_extraction_script(
                state.trajectory, state.tpr, self.index, chunk, self.out_dir, self.gmx.executable
            )
            slurm.submit(body, "FRAMEGEN", self.slurm, self.out_dir, gpus=False)
        slurm.wait_for_files(windows, self.slurm.poll_seconds, self.slurm.timeout_seconds)
        return replace(state, windows=windows)

    def _run_windows(self, state: Artifacts) -> Artifacts:
        if not state.windows:
            raise ConfigError(f"umbrella stage {self.name} has no window structures; run a CONFIG stage first")
        centre_atom = umbrella.nearest_to_center(state.windows[0], self.index)
        mdp_path = self._write_mdp({"pull_group2_pbcatom": centre_atom})

        for chunk in umbrella.chunk_evenly(state.windows, self.slurm.window_jobs):
            body = umbrella.window_script(
                chunk, mdp_path, self.topology, self.index, self.name in MINIMIZATION_STAGES, self.gmx.executable
            )
            slurm.submit(body, "WINDOWS", self.slurm, self.out_dir)
        outputs = [self.out_dir / w.name for w in state.windows]
        slurm.wait_for_files(outputs, self.slurm.poll_seconds, self.slurm.timeout_seconds)
        return replace(state, windows=outputs)
