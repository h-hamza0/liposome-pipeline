# Configuration reference

A configuration file is a sequence of blocks. Each block starts with a header line, contains
`key = value` lines, and ends with `end`. Lines starting with `#` are ignored. Relative paths are
resolved against the directory you run the command from.

| Header | Purpose | Count |
|---|---|---|
| `system` | Paths, box, solvent, index groups | exactly one |
| `molecule_<label>` | One molecular species to insert | one or more |
| `stage_<label>` | One simulation stage, run in file order | zero or more |
| `slurm` | Cluster settings for umbrella-sampling stages | optional |

## `system`

| Key | Legacy name | Required | Description |
|---|---|---|---|
| `input_dir` | | yes | Directory with `<name>.itp`, `<name>.gro` per molecule and `water.gro` |
| `output_dir` | | yes | Run directory (must not exist unless `--overwrite`) |
| `forcefield` | | yes | Directory of force-field `.itp` files, included in alphabetical order |
| `index_groups` | `idxGRPS` | yes | `NAME:['RES1','RES2'] ...` custom index groups built from default groups |
| `box_size` | `boxSize` | yes | One value or three values in nm |
| `water_radius` | `Wradius` | yes | `gmx solvate -radius` in nm |
| `center_box` | `boxDist` | no | If non-zero, re-centre in a cubic box of this edge (nm) before solvation |
| `restart_from` | `US` | no | Previous run directory; starts from its `MD/aligned.gro` |
| `drug`, `percentile` | | no | Adds a `SELECT` index group for the drug molecule at that distance percentile from the lipid centre |
| `gmx_path` | | no | GROMACS executable |

If an index group named `LIPD` is defined, a `CORE` group (lipid `NC3`, `NH3` and `PO4` beads) is
added automatically for use as a pull reference.

When restarting, only molecules named `NA`, `DEXT` or `IBU` are inserted unless a molecule sets
`insert = yes|no` explicitly.

## `molecule_<label>`

`name` (matches `<name>.itp`/`.gro`), `nmol`, `insertion_radius` (nm), optional `insert`.

## `stage_<label>`

`type` selects a bundled MDP template (`EM`, `NVT`, `NPT`, `NPT_01`, `NPT_015`, `NPT_018`, `EQUIL`,
`MD`, `PULL`, `PULL_EQUIL`, `PULL_US`, `EM_US`) or the special `CONFIG` stage. Stage types must be
unique because each stage gets its own directory. Other keys:

| Key | Description |
|---|---|
| `mdp_file` | Use this MDP file instead of the template |
| `umbrella` | `yes` runs the stage over all umbrella windows through SLURM |
| `lipid`, `drug` | `CONFIG` only: index groups whose centres of mass define the pull distance |
| anything else | Written into the MDP file, overriding the template (`nsteps = 1000`, `tc-grps = LIPD SOLV`) |

## `slurm`

`partition`, `cpus_per_task`, `mem`, `gres`, `mail_user`, `mail_type`, `setup_commands`
(`;`-separated shell lines such as `module purge; source GMXRC`), `frame_jobs`, `window_jobs`,
`poll_seconds`, `timeout_seconds`.
