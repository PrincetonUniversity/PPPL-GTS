# Post-processing

GTS-S writes its main results to `snapshots_hdf5/snapXXXXX.h5`, one file every `msnap` steps. The physical time is in `time_s`.
This folder has scripts for reading and plotting these files, and a list of the main variables.

## Plotting script
`plot_snapshot.py` needs numpy, h5py and matplotlib.
```
python plot_snapshot.py snapshots_hdf5/snap00010.h5
python plot_snapshot.py snapshots_hdf5/snap0019[1-9].h5 snapshots_hdf5/snap00200.h5   # average over 191-200 us
```
It writes two figures:
- `<name>_profiles.png`: flux-surface-averaged n_e, n_i, T_e and T_i versus a = √ψt
- `<name>_2D.png`: Φ, n_e, T_e and T_i on one poloidal plane

With several files, every quantity is averaged over the files.

Options:
- `--instant`: plot profiles at the snapshot step instead of the interval averages.
- `--min-density FRAC`: hide temperatures where the density is below FRAC × its maximum (default 0.05).
- `--outdir DIR`: write the figures to DIR.
- `--show`: also open the figures in a window.

## Radial profiles
Flux-surface averages are in group `DG/ZON`. The radius a = √ψt is in `DG/a_1D`.

| Variable | Quantity | Unit |
|---|---|---|
| `mean_ne_tot_by_fh` | electron density | m^-3 |
| `mean_ni_by_ff` | ion density | m^-3 |
| `mean_Te_by_fh` | electron temperature | eV |
| `mean_Ti_by_ff` | ion temperature | eV |

- `mean_*`: average over the interval since the previous snapshot.
- The same name without `mean_`: value at the snapshot step.

## 2D maps on one poloidal plane
| Variable | Quantity | Unit | Grid |
|---|---|---|---|
| `SG/phi` | potential Φ | V | `SG/a_2D`, `SG/gq_2D` |
| `DG/ne_tot_by_fh` | electron density | m^-3 | `DG/a_2D`, `DG/gq_2D` |
| `DG/ni` | ion density | m^-3 | `DG/a_2D`, `DG/gq_2D` |
| `DG/Te`, `DG/Ti` | electron and ion temperatures | eV | `DG/a_2D`, `DG/gq_2D` |

- The arrays are ordered (a, gq). `gq` is the poloidal angle, from 0 to 2π.
- The maps are values at the snapshot step. For time-averaged maps, average over several snapshots.
- Where the density is almost zero, `DG/Te` and `DG/Ti` can take extreme values.

## Electron density
With `ifelectron=2`, the electron distribution is split into a kinetic part h and an adiabatic part, so several variables have similar names.
- The total electron density is `ne_tot_by_fh` = `ne_by_fh` (full-h part) + `dne_adiab` (adiabatic part).
- `ne_by_fh` alone is not the total. With the Poisson solve it can be far from the total, and even negative.
- `dne`, `dni`, `DG/dne_tot` and `DG/dni` are perturbations (δh or δf), not total densities.
- Without the Poisson solve, `dne_adiab` = 0, so `ne_tot_by_fh` equals `ne_by_fh`.

## Datasets that are not filled
These datasets are always zero in GTS-S:
- some top-level datasets from the original GTS diagnostics, such as `dphi`, `full_dne` and `mean_full_dphi`;
- `SG/mean_phi`.
