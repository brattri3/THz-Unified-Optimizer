# Specac measurement drop-folder

Put the THz time-domain traces here. `DataManager` picks them up automatically
as dataset **`specac`** when the filename matches:

```
specac_<angle>deg_rep<N>_<sig|bg>.txt
# e.g.  specac_-10deg_rep1_sig.txt  /  specac_-10deg_rep1_bg.txt
```
- `<angle>`: integer degrees, may be signed; one `sig` + one `bg` per (angle, rep).
- Two-column ASCII (time, field), same format as the other datasets.

Micrographs go to `research/validation/micrographs/specac/` using the same
naming as PureWave (`specac__cal__<mag>mag__div0p01mm.*`,
`specac__wgp__<mag>mag[_NN].*`); then run `research/experiments/micrograph_period.py`
and fill the row in `research/validation/MANIFEST.md`.
