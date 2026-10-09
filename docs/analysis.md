# Analysis notes

- **Vesicle centre:** centre of geometry of all selected lipid atoms
  (`resname DOPC PEL CHEMS and not (name NH O1)`).
- **Outer radius:** 90th percentile of residue-centre distances from that centre.
- **Encapsulated drug:** residue centre of geometry within the outer radius.
- **Box volume:** `box_side³`, supplied with `--box-nm`; the vesicle sphere volume is subtracted
  to obtain the bulk volume.
- **Summary:** per replicate, the mean over the trailing `--tail-fraction` of frames; then mean and
  sample standard deviation across replicates. Check the time series for equilibration before
  trusting the tail average.
