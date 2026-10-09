# Changelog

## 0.1.0

Initial public release, refactored from the original research scripts.

- Installable package with `liposome-build` and `liposome-partition` commands.
- Validated configuration parser with descriptive errors (legacy key names still accepted).
- Index groups are resolved by name instead of a hard-coded `make_ndx` group number.
- MDP overrides keep parameters that carry inline comments and match `-`/`_` spellings.
- GROMACS jobs run inside their stage directory and resume from checkpoints a bounded number of times.
- Existing output directories are never deleted unless `--overwrite` is passed.
- Umbrella sampling: the pull reference atom is the lipid atom nearest the vesicle centre,
  identified by atom index (previously its residue number, and the nearest-atom search never updated).
- Partition analysis: system discovery works for any root directory; the per-residue geometry is vectorised;
  single-replicate groups no longer produce NaN error bars.
- Site-specific settings (SLURM partition, e-mail, GROMACS setup scripts) moved out of the code.
