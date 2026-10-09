from __future__ import annotations

import argparse
import logging
import os
import shlex
import sys

from .errors import PipelineError
from .gromacs import Gromacs, MdrunOptions
from .pipeline import Pipeline


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="liposome-build", description="Build and equilibrate a Martini 3 liposome system.")
    parser.add_argument("-v", "--verbose", action="store_true")
    commands = parser.add_subparsers(dest="command", required=True)

    check = commands.add_parser("check", help="validate a configuration file and print the planned pipeline")
    check.add_argument("config")

    run = commands.add_parser("run", help="execute the pipeline described by a configuration file")
    run.add_argument("config")
    run.add_argument("--gmx", help="GROMACS executable (default: gmx_path from the config, else 'gmx')")
    run.add_argument("--launcher", default=os.environ.get("LIPOSOME_LAUNCHER", ""), help="command prefix such as 'srun'")
    run.add_argument("--overwrite", action="store_true", help="delete an existing output directory first")
    for flag in ("ntomp", "ntmpi"):
        run.add_argument(f"--{flag}", type=int)
    for flag in ("nb", "pme", "bonded", "pin"):
        run.add_argument(f"--{flag}")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(levelname)s %(message)s")
    try:
        pipeline = Pipeline.from_file(args.config)
        if args.command == "check":
            print(pipeline.describe())
            return 0
        options = MdrunOptions(args.ntomp, args.ntmpi, args.nb, args.pme, args.bonded, args.pin)
        pipeline.use_gromacs(Gromacs(args.gmx or pipeline.gmx.executable, shlex.split(args.launcher), options))
        final = pipeline.run(overwrite=args.overwrite)
        print(f"Done. Final structure: {final.structure}")
        return 0
    except PipelineError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
