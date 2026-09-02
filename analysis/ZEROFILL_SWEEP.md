# Zero-fill sweep, 2026-09-02

Four separate zero-fills turned up one at a time over two days, always in the same place:
where an OPTIONAL layer meets an array that has no "unmeasured" value. Andy asked for a
sweep rather than a fifth one-at-a-time fix.

The defect is not "a zero appears". It is **absence becoming a measurement**: a value a
downstream cut reads as data. Zero is a measurement for a scan angle (nadir), a boolean
(no forest here), a stratum code (other), a timestamp, and a mask (nothing excluded).

Idioms searched across all tracked Python: `else np.zeros`, `else np.full`, `nan_to_num`,
`fillna(0)`, `np.where(isfinite(x), x, 0)`, `.get(..., 0)`, `else False`, `zeros_like`.

## FIXED

### 1. `io.read_tile` read every PF6+ file as all-nadir, at time zero -- LIBRARY, LIVE

```python
scan_angle = np.asarray(f.scan_angle_rank) if "scan_angle_rank" in dims else np.zeros(...)
gps_time   = np.asarray(f.gps_time)        if "gps_time" in dims        else np.zeros(...)
```

Point format 6 and later DROP `scan_angle_rank` in favour of a scaled `scan_angle`. Every
CSF cache this project writes is PF7, so every one took the fallback: scan angle 0 for
every return, which is a claim that the whole tile was flown at nadir. Consumers of
`pc.scan_angle` include the across-track intercept tie in `coreg.coregister_swaths`, which
fits against `tan(scan_angle)`, and `localtie`.

Reads PF6+ `scan_angle` (0.006-deg units) and PF<=5 `scan_angle_rank` (degrees), both
returned in DEGREES; NaN only if a file carries neither. Verified on real files -- the PF7
cache and the PF1 source now agree at +-17.00 deg, the scanner's actual half-angle:

    elba.las        (PF7)  finite 6,771,612/6,771,612   min -17.00  max +17.00 deg
    4342-29-64.laz  (PF1)  finite 7,728,747/7,728,747   min -17.00  max +17.00 deg

**Scope, checked, not assumed:** `difference_dem` is NOT affected. It builds its PointCloud
explicitly from `f.scan_angle * 0.006` (pipeline.py:686, 697, 707) rather than through
`read_tile`. The live callers were `ground_control/our_surface.py` and
`scripts/swath_consistency.py`.

Five tests in `tests/test_io_scan_angle.py`, shown to bite: restoring the PF<=5-only lookup
fails 2 of them.

### 2. `datum_from_mass_balance.py` -- absent floodplain mask masked NOTHING

`flood = ... if fld_p else np.zeros(dod.shape, bool)`. An all-False mask asserts no cell is
floodplain. The standing rule keeps floodplain cells out of hillslope mass balance, so a run
without the mask is a DIFFERENT population. Now refuses, as `refcells.reference_cells`
already did for the same layer.

### 3. `dod_corrections.py` -- unmeasured vegetation received a zero correction

`d = d - m_*np.nan_to_num(veg)`: cells with no `veg_frac` stayed in the product, uncorrected
but indistinguishable from cells that needed no correction. They are set NaN instead --
uncorrectable, and visibly so -- and the count is printed.

## ALREADY CORRECT (checked, not assumed)

- `plot_q2_fit.py:46` fills an absent mask with False, but the panel title already says
  `NO MASK EXISTS -- cut skipped`. The figure is honest; only a comment was added.
- `gen1_save_angles_slope.py:140` uses `np.full(len(x), np.nan)` -- for a float, NaN IS the
  unmeasured value.
- `coherence.py:161` zeroes only where `valid` is already False.
- `detect.py:108-110` is scale-space bookkeeping, not a data substitution.
- `rasterize_change.py:97` fills for DISPLAY only (`nan -> nanmin` for a colour ramp).
- `control_band_window_sweep.py:89` returns a `(n, 0)` array -- a shape, not a value.
- `save_incidence_linked.py:34` writes -9999 as a documented LAS sentinel.

## STILL OPEN

### `beam_offset_table.py:106` -- `overlap` and `scanner_channel`

`_opt(name, dt)` returns zeros for LAS dimensions absent in PF<=5. On gen1 every point then
reads `overlap = 0`, "measured, not in overlap". Unlike the cases above, REFUSING is wrong:
PF<=5 genuinely cannot record it. It needs a representation -- a NaN column, or a
`has_overlap: false` flag in the table's metadata. **Andy's call; not decided.**

### `pipeline.py:818-820` -- a half-measured cell gets a half-sized standard error

```python
stderr = np.sqrt(np.nan_to_num(r08)**2 / max(n08,1) + np.nan_to_num(r21)**2 / max(n21,1))
stderr[~(np.isfinite(r08) | np.isfinite(r21))] = np.nan
```

The guard on the next line means both-NaN correctly yields NaN, so this was considered. But
a cell with roughness in ONE epoch contributes only that epoch's term, understating the
standard error that feeds the LoD. Not a plain zero-fill, and **not sized**: `r08`/`r21` are
internal to `difference_dem` and not persisted, so the affected cell count cannot be read
off the products. Flagged, not changed.
