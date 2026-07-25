# PureWave measurement drop-folder

Put the THz time-domain traces here. `DataManager` (rglob over `data_pool/`)
picks them up automatically as dataset **`purewave`** as long as the filename
matches:

```
purewave_<angle>deg_rep<N>_<sig|bg>.txt
# e.g.  purewave_-10deg_rep1_sig.txt  /  purewave_-10deg_rep1_bg.txt
```
- `<angle>`: integer degrees, may be signed (`-40`, `0`, `40`); one `sig` + one `bg` per (angle, rep).
- Two-column ASCII (time, field), same format as `356att_*` / `test_grid_40_20_*`.

Geometry ground truth from micrographs (add to `GEOMETRY` in
`research/experiments/fit_lib.py`): **P ≈ 25.5 µm, D ≈ 10.6 µm** (fill ≈ 0.42) —
see `research/validation/MANIFEST.md`.
