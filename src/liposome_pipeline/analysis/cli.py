from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

from ..errors import PipelineError
from . import discovery, partition, plots

DEFAULT_SMOOTH_WINDOW = 15
DEFAULT_TAIL_FRACTION = 0.75


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="liposome-partition", description="Drug partition-ratio analysis for liposome trajectories.")
    parser.add_argument("-v", "--verbose", action="store_true")
    sub = parser.add_subparsers(dest="command", required=True)

    def roots(p: argparse.ArgumentParser) -> None:
        p.add_argument("roots", nargs="*", default=["."], help="directories to search (default: .)")

    def compute_opts(p: argparse.ArgumentParser) -> None:
        p.add_argument("--box-nm", type=float, default=partition.DEFAULT_BOX_SIDE_NM, help="cubic box edge in nm")
        p.add_argument("--stride", type=int, default=1, help="analyse every Nth frame")

    def plot_opts(p: argparse.ArgumentParser) -> None:
        p.add_argument("--output-dir", default=".")
        p.add_argument("--smooth-window", type=int, default=DEFAULT_SMOOTH_WINDOW, help="rolling-mean window in frames")
        p.add_argument("--tail-fraction", type=float, default=DEFAULT_TAIL_FRACTION, help="trailing fraction of each trajectory to average")
        p.add_argument("--dpi", type=int, default=300)

    roots(sub.add_parser("count", help="print the number of drug-loaded systems (sizes a SLURM array)"))

    comp = sub.add_parser("compute", help="analyse one system and cache its result")
    roots(comp)
    comp.add_argument("--index", type=int, required=True, help="0-based system index (e.g. $SLURM_ARRAY_TASK_ID)")
    comp.add_argument("--cache-dir", default="results_cache")
    comp.add_argument("--force", action="store_true", default=bool(os.environ.get("FORCE")))
    compute_opts(comp)

    agg = sub.add_parser("aggregate", help="plot and summarise all cached results")
    agg.add_argument("--cache-dir", default="results_cache")
    plot_opts(agg)

    every = sub.add_parser("run", help="analyse every system serially and plot")
    roots(every)
    compute_opts(every)
    plot_opts(every)
    return parser


def _report(results: partition.Results, args: argparse.Namespace) -> None:
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    plots.plot_timeseries(results, args.smooth_window, out / "partition_ratio_timeseries.pdf", args.dpi)
    partition.timeseries_table(results, args.smooth_window).to_csv(out / "partition_ratio_timeseries.csv", index=False)
    per_replicate, grouped = partition.summarize(results, args.tail_fraction)
    per_replicate.to_csv(out / "partition_ratio_summary.csv", index=False)
    grouped.to_csv(out / "partition_ratio_grouped.csv", index=False)
    plots.plot_bars(grouped, args.tail_fraction, out / "partition_ratio_bar.pdf", args.dpi)
    print(grouped.to_string(index=False))


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(levelname)s %(message)s")
    try:
        if args.command == "aggregate":
            _report(partition.load_cache(args.cache_dir), args)
            return 0

        entries = discovery.drug_loaded(discovery.discover_systems(args.roots))
        if args.command == "count":
            print(len(entries))
            return 0
        if not entries:
            raise PipelineError(f"no drug-loaded systems found under: {' '.join(args.roots)}")

        options = {"box_side_nm": args.box_nm, "stride": args.stride}
        if args.command == "compute":
            if not 0 <= args.index < len(entries):
                raise PipelineError(f"--index {args.index} out of range (0-{len(entries) - 1})")
            print(partition.compute_one(entries[args.index], args.cache_dir, args.force, **options))
            return 0

        _report(partition.compute_all(entries, **options), args)
        return 0
    except (PipelineError, FileNotFoundError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
