from __future__ import annotations

import logging
import subprocess
import time
from pathlib import Path
from typing import Iterable

from .config import SlurmConfig
from .errors import SlurmError

log = logging.getLogger(__name__)


def build_script(body: str, config: SlurmConfig, workdir: Path) -> str:
    header = ["#!/bin/bash", *config.setup_commands, f"cd {workdir}"]
    return "\n".join(header) + "\n" + body.strip("\n") + "\n"


def sbatch_command(job_name: str, config: SlurmConfig, gpus: bool) -> list[str]:
    command = [
        "sbatch",
        "--parsable",
        f"--job-name={job_name}",
        "--output=%x_%j.log",
        "--ntasks=1",
        f"--cpus-per-task={config.cpus_per_task}",
        f"--mem={config.mem}",
    ]
    if gpus and config.gres:
        command += [f"--gres={config.gres}"]
    if config.partition:
        command += [f"--partition={config.partition}"]
    if config.mail_user:
        command += [f"--mail-user={config.mail_user}", f"--mail-type={config.mail_type}"]
    return command


def submit(body: str, job_name: str, config: SlurmConfig, workdir: Path, gpus: bool = True) -> str:
    script = build_script(body, config, workdir)
    result = subprocess.run(
        sbatch_command(job_name, config, gpus), input=script, text=True, capture_output=True
    )
    if result.returncode != 0:
        raise SlurmError(f"sbatch failed for {job_name}: {result.stderr.strip()}")
    job_id = result.stdout.strip().split(";")[0]
    log.info("submitted %s as job %s", job_name, job_id)
    return job_id


def wait_for_files(paths: Iterable[Path], poll_seconds: float, timeout_seconds: float | None) -> None:
    pending = [Path(p) for p in paths]
    deadline = None if timeout_seconds is None else time.monotonic() + timeout_seconds
    while pending:
        pending = [p for p in pending if not p.exists()]
        if not pending:
            return
        if deadline is not None and time.monotonic() > deadline:
            raise SlurmError(f"timed out waiting for {len(pending)} file(s), first missing: {pending[0]}")
        log.info("waiting on %d file(s)", len(pending))
        time.sleep(poll_seconds)
