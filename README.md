# GTS-S

GTS-S is a version of the gyrokinetic particle-in-cell code GTS, developed by Min-Gu Yoo.
It extends GTS so that it can simulate plasma transport in the tokamak edge, where static magnetic perturbations make the magnetic field lines stochastic and open to the wall.
To do this, GTS-S adds several new algorithms and numerical methods, which are listed below.

This branch is a snapshot of the GTS-S code, placed on top of `master`.
It does not share history with `master`, and later changes to `master` are not included.

## New algorithms and methods
- General δf particle-in-cell algorithm (Ref. 3), which stays efficient with physical boundaries and strong forces:
  - radial particle boundaries with an adaptive particle source;
  - kinetic electrons with full-h weights, with the adiabatic response treated separately.
- Gyrokinetic Poisson solver on a structured flux-coordinate grid (PETSc)
- Diagnostics on a separate grid: flux-surface averages and 2D maps, written to HDF5 snapshots

A field solver using the directional finite difference method (DFDM, Ref. 4) is implemented in `src/mgy_petsc.F90`, but simulations do not use it yet.

## Building
The steps are the same as for `master`:
1. Load the compiler, MPI, HDF5, NetCDF and FFTW modules.
2. Install PETSc, PSPLINE and SPRNG.
3. In `src/`, run:
   ```
   make _ARCH=<arch>
   ```

System settings are in `src/<arch>.mk`, for example:
- `stellar`, `flux`, `perlmutter`
- `perlmuttergpu`, for the OpenACC GPU build (add `ACC=y`)

NAG is not needed.

## Examples
- `examples/POP2022_UniformPhase`: stochastic edge from 8 MP harmonics with the same phase, with the GK Poisson equation solved (32 MPI ranks). Its README explains:
  - the setup and the main input flags;
  - how to switch off the Poisson solve;
  - the outputs.
- `examples/POP2022_RandomPhase`: the same case with random MP phases.

## Output
The main results are written to HDF5 snapshots, `snapshots_hdf5/snapXXXXX.h5`. `postprocessing/` lists the main variables for density, temperature and potential, and has a Python script that plots radial profiles and 2D maps.

## References
1. M.-G. Yoo, W. X. Wang, E. Startsev, C. H. Ma, S. Ethier, J. Chen, and X. Z. Tang, "Collisionless plasma transport mechanisms in stochastic open magnetic field lines in tokamaks," *Nucl. Fusion* **61**, 126036 (2021). https://doi.org/10.1088/1741-4326/ac30c6
2. M.-G. Yoo, W. X. Wang, E. Startsev, C. H. Ma, S. Ethier, J. Chen, and X. Z. Tang, "The 3D magnetic topology and plasma dynamics in open stochastic magnetic field lines," *Phys. Plasmas* **29**, 072502 (2022). https://doi.org/10.1063/5.0085304
3. M.-G. Yoo, W. X. Wang, and E. Startsev, "Innovative delta-f PIC algorithm for efficient homogeneous simulation of fusion plasmas from core to edge," 29th IAEA Fusion Energy Conference (FEC 2023), London, UK (2023). https://conferences.iaea.org/event/316/contributions/28229/
4. M.-G. Yoo, W. X. Wang, E. Startsev, and S. Ethier, "Directional finite difference method for directly solving 3D gyrokinetic field equations with enhanced accuracy," *Comput. Phys. Commun.* **312**, 109597 (2025). https://doi.org/10.1016/j.cpc.2025.109597
