from __future__ import annotations

import logging
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from .errors import GromacsError

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class MdrunOptions:
    ntomp: int | None = None
    ntmpi: int | None = None
    nb: str | None = None
    pme: str | None = None
    bonded: str | None = None
    pin: str | None = None
    checkpoint_minutes: int = 5

    def args(self) -> list[str]:
        flags = {
            "-ntomp": self.ntomp,
            "-ntmpi": self.ntmpi,
            "-nb": self.nb,
            "-pme": self.pme,
            "-bonded": self.bonded,
            "-pin": self.pin,
        }
        out = ["-cpt", str(self.checkpoint_minutes)]
        for flag, value in flags.items():
            if value is not None:
                out += [flag, str(value)]
        return out


class Gromacs:
    def __init__(
        self,
        executable: str = "gmx",
        launcher: Sequence[str] = (),
        mdrun_options: MdrunOptions | None = None,
        max_warnings: int = 10,
    ) -> None:
        self.executable = executable
        self.launcher = list(launcher)
        self.mdrun_options = mdrun_options or MdrunOptions()
        self.max_warnings = max_warnings

    def ensure_available(self) -> None:
        if shutil.which(self.executable) is None:
            raise GromacsError(f"GROMACS executable not found on PATH: {self.executable!r}")

    def run(self, *args: object, stdin: str | None = None, cwd: Path | None = None) -> None:
        command = [*self.launcher, self.executable, *map(str, args)]
        log.info("$ %s", shlex.join(command))
        result = subprocess.run(command, input=stdin, text=True, cwd=cwd)
        if result.returncode != 0:
            raise GromacsError(f"command failed with exit code {result.returncode}: {shlex.join(command)}")

    def insert_molecules(
        self,
        molecule: Path,
        output: Path,
        count: int,
        radius: float,
        base: Path | None = None,
        box: Sequence[float] | None = None,
        tries: int = 500,
    ) -> None:
        if base is None and box is None:
            raise ValueError("insert_molecules needs either an existing structure or a box size")
        args: list[object] = ["insert-molecules", "-ci", molecule, "-o", output]
        args += ["-f", base] if base is not None else ["-box", *box]  # type: ignore[misc]
        args += ["-nmol", count, "-radius", radius, "-try", tries]
        self.run(*args)

    def center_in_box(self, structure: Path, output: Path, edge: float) -> None:
        self.run("editconf", "-f", structure, "-o", output, "-c", "-box", edge, edge, edge, "-pbc", "yes", "-bt", "cubic")

    def solvate(self, structure: Path, solvent: Path, topology: Path, output: Path, radius: float) -> None:
        self.run("solvate", "-cp", structure, "-cs", solvent, "-p", topology, "-o", output, "-radius", radius)

    def grompp(
        self,
        mdp: Path,
        structure: Path,
        topology: Path,
        output: Path,
        index: Path | None = None,
        restraints: Path | None = None,
    ) -> None:
        args: list[object] = ["grompp", "-f", mdp, "-c", structure, "-p", topology, "-o", output, "-po", output.with_name("mdout.mdp")]
        if index is not None:
            args += ["-n", index]
        args += ["-r", restraints or structure, "-maxwarn", self.max_warnings]
        self.run(*args)

    def mdrun(self, deffnm: str, cwd: Path | None = None, resume: bool = False) -> None:
        args: list[object] = ["mdrun", "-deffnm", deffnm, "-v", *self.mdrun_options.args()]
        if resume:
            args += ["-cpi", f"{deffnm}.cpt"]
        self.run(*args, cwd=cwd)

    def make_index(self, structure: Path, output: Path, commands: str = "q\n") -> None:
        self.run("make_ndx", "-f", structure, "-o", output, stdin=commands)

    def distance(self, trajectory: Path, tpr: Path, index: Path, group_a: str, group_b: str, output: Path) -> None:
        selection = f'com of group "{group_a}" plus com of group "{group_b}"'
        self.run("distance", "-s", tpr, "-f", trajectory, "-n", index, "-select", selection, "-oall", output)


def index_group_names(index_file: Path) -> list[str]:
    names = []
    for line in Path(index_file).read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            names.append(stripped[1:-1].strip())
    return names


def read_index_groups(index_file: Path) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = {}
    current: list[int] | None = None
    for line in Path(index_file).read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            current = groups.setdefault(stripped[1:-1].strip(), [])
        elif stripped and current is not None:
            current.extend(int(token) for token in stripped.split())
    return groups
