from __future__ import annotations

from pathlib import Path

from .errors import ConfigError
from .gromacs import Gromacs, index_group_names

CORE_GROUP = "CORE"
CORE_SOURCE_GROUP = "LIPD"
CORE_BEADS = ("NC3", "NH3", "PO4")


def build_index(
    gmx: Gromacs,
    structure: Path,
    output: Path,
    groups: dict[str, list[str]],
    extra_commands: str = "",
) -> Path:
    gmx.make_index(structure, output)
    available = {name: position for position, name in enumerate(index_group_names(output))}

    unknown = sorted({member for members in groups.values() for member in members} - available.keys())
    if unknown:
        raise ConfigError(
            f"index group member(s) not present in {structure.name}: {', '.join(unknown)} "
            f"(available: {', '.join(available)})"
        )

    next_id = len(available)
    commands = []
    ids: dict[str, int] = {}
    for name, members in groups.items():
        commands.append("|".join(str(available[m]) for m in members))
        commands.append(f"name {next_id} {name}")
        ids[name] = next_id
        next_id += 1

    if CORE_SOURCE_GROUP in ids:
        commands.append(f"{ids[CORE_SOURCE_GROUP]} & a {' '.join(CORE_BEADS)}")
        commands.append(f"name {next_id} {CORE_GROUP}")
        next_id += 1

    script = "\n".join(commands) + "\n" + extra_commands.format(next_id=next_id) + "q\n"
    gmx.make_index(structure, output, script)
    return output
