from __future__ import annotations

import logging
import shutil
from pathlib import Path

from .config import SystemConfig
from .errors import ConfigError
from .gromacs import Gromacs
from .index import build_index
from .molecules import Molecule

log = logging.getLogger(__name__)

RESTART_INSERTABLE = frozenset({"NA", "DEXT", "IBU"})
WATER_FILE = "water.gro"
RESTART_STRUCTURE = Path("MD") / "aligned.gro"
RESTART_TOPOLOGY = Path("setup") / "topol.top"


def atom_count(gro: Path) -> int:
    return int(gro.read_text().splitlines()[1])


class System:
    def __init__(self, config: SystemConfig, molecules: list[Molecule], gmx: Gromacs) -> None:
        self.config = config
        self.molecules = molecules
        self.gmx = gmx
        self.output_dir = config.output_dir.absolute()
        self.setup_dir = self.output_dir / "setup"
        self.forcefield = self.setup_dir / "forcefield"
        self.topology = self.setup_dir / "topol.top"
        self.index = self.setup_dir / "index.ndx"
        self.structure: Path | None = None
        self._restart_topology: str | None = None

    def setup(self, overwrite: bool = False) -> Path:
        self._prepare_directories(overwrite)
        self._stage_molecules()
        self._write_topology()
        if self.structure is not None and self.config.center_box:
            self._center()
        self._insert_molecules()
        if self.structure is not None and self.config.center_box and not self.config.restart_from:
            self._center()
        self._solvate()
        self._write_index()
        assert self.structure is not None
        return self.structure

    def _prepare_directories(self, overwrite: bool) -> None:
        if self.output_dir.exists():
            if not overwrite:
                raise ConfigError(f"output directory already exists: {self.output_dir} (pass --overwrite to replace it)")
            shutil.rmtree(self.output_dir)
        for sub in ("molecules", "itp"):
            (self.setup_dir / sub).mkdir(parents=True)
        shutil.copytree(self.config.forcefield, self.forcefield)

        water = self.config.input_dir / WATER_FILE
        if not water.is_file():
            raise ConfigError(f"solvent structure not found: {water}")
        shutil.copy(water, self.setup_dir / "molecules" / WATER_FILE)

        if self.config.restart_from is not None:
            source = self.config.restart_from / RESTART_STRUCTURE
            if not source.is_file():
                raise ConfigError(f"restart structure not found: {source}")
            topology = self.config.restart_from / RESTART_TOPOLOGY
            if not topology.is_file():
                raise ConfigError(f"restart topology not found: {topology}")
            self.structure = self.setup_dir / "clust.gro"
            shutil.copy(source, self.structure)
            self._restart_topology = topology.read_text()

    def _stage_molecules(self) -> None:
        for molecule in self.molecules:
            for suffix, subdir in ((".itp", "itp"), (".gro", "molecules")):
                source = self.config.input_dir / f"{molecule.name}{suffix}"
                if not source.is_file():
                    raise ConfigError(f"missing input file for molecule {molecule.name}: {source}")
                destination = self.setup_dir / subdir / source.name
                shutil.copy(source, destination)
                setattr(molecule, "itp_file" if suffix == ".itp" else "gro_file", destination)

    def _write_topology(self) -> None:
        if self._restart_topology is not None:
            self._write_restart_topology()
            return
        lines = [f"#include <{itp.absolute()}>" for itp in sorted(self.forcefield.glob("*.itp"))]
        lines += [f"#include <{m.itp_file}>" for m in self.molecules]
        lines += ["", "[ system ]", "System", "", "[ molecules ]"]
        lines += [f"{m.name}   {m.count}" for m in self.molecules]
        self.topology.write_text("\n".join(lines) + "\n")

    def _write_restart_topology(self) -> None:
        text = self._restart_topology
        inserted = [m for m in self.molecules if self._should_insert(m)]
        includes = [f"#include <{m.itp_file}>" for m in inserted if f"{m.name}.itp" not in text]
        if includes:
            head, marker, tail = text.partition("[ system ]")
            text = head + "\n".join(includes) + "\n\n" + marker + tail
        text = text.rstrip("\n") + "\n" + "".join(f"{m.name}   {m.count}\n" for m in inserted)
        self.topology.write_text(text)

    def _center(self) -> None:
        output = self.setup_dir / "centered.gro"
        self.gmx.center_in_box(self.structure, output, self.config.center_box)
        self.structure = output

    def _should_insert(self, molecule: Molecule) -> bool:
        if molecule.insert is not None:
            return molecule.insert
        if self.config.restart_from is not None:
            return molecule.name in RESTART_INSERTABLE
        return True

    def _insert_molecules(self) -> None:
        for molecule in self.molecules:
            if not self._should_insert(molecule):
                continue
            output = self.setup_dir / f"{molecule.name}_boxed.gro"
            self.gmx.insert_molecules(
                molecule.gro_file,
                output,
                molecule.count,
                molecule.insertion_radius,
                base=self.structure,
                box=self.config.box_size if self.structure is None else None,
            )
            self._verify_insertion(molecule, output)
            self.structure = output

    def _verify_insertion(self, molecule: Molecule, output: Path) -> None:
        before = atom_count(self.structure) if self.structure is not None else 0
        added = atom_count(output) - before
        expected = molecule.count * molecule.bead_count
        if added != expected:
            raise ConfigError(
                f"inserted {added // molecule.bead_count} of {molecule.count} {molecule.name} molecules; "
                f"the box is too crowded for insertion_radius={molecule.insertion_radius} nm"
            )

    def _solvate(self) -> None:
        output = self.setup_dir / "solvated.gro"
        self.gmx.solvate(
            self.structure, self.setup_dir / "molecules" / WATER_FILE, self.topology, output, self.config.water_radius
        )
        self.structure = output

    def _write_index(self) -> None:
        extra = ""
        if self.config.drug is not None:
            from .umbrella import select_drug_command

            if self.config.percentile is None:
                raise ConfigError("'drug' requires 'percentile' to be set as well")
            extra = select_drug_command(self.structure, self.config.drug, self.config.percentile)
        build_index(self.gmx, self.structure, self.index, self.config.index_groups, extra)
