# Example: stochastic edge from magnetic perturbations (uniform phase)

Static magnetic perturbations (MP) make the edge magnetic field lines stochastic. GTS-S follows gyrokinetic ions and kinetic electrons in this field and solves the gyrokinetic (GK) Poisson equation.
`../POP2022_RandomPhase` is the same case with random MP phases.

## Setup
- **Geometry:** cyclone-like equilibrium (`esiP.00_cyclon`), with radial domain a = √ψ_t = 0.1–0.9.
- **Plasma:** flat initial profiles, n = 1.6e19 m^-3 and Te = Ti = 5 keV.
- **MP:** 8 static harmonics, (m,n) = (3,2), (2,1), (5,2), (4,2), (5,3), (6,3), (7,3), (8,3). All have amplitude `ga0=6e-4` and the same phase (`gagf0`).
- **Electrons:** kinetic, with the adiabatic response separated (`ifelectron=2`) and full-h weights (`Flag_Load_and_WeightDef=1`, `Flag_Weight_dfh=1`, `Flag_w_Evol=4`).
- **Ions:** gyrokinetic, with a 4-point gyroaverage.
- **Fields:** GK Poisson equation on the structured grid (`Flag_Poisson=1`, `Flag_Poisson3D_on_SG=1`), solved with PETSc (`PETSc_input`).
- **Radial boundaries** (`Flag_Ptl_BC0_type = Flag_Ptl_BC1_type = 2`): every `Ptl_Src_nstep=10` steps, markers that left the domain are removed.
  - Inner boundary: an adaptive source keeps its initial density.
  - Outer boundary: absorbs particles.
- **Numerics:** dt = 10 ns and 1000 steps (10 µs), with 100 radial grid points and 8 poloidal planes.

## Running
Put the GTS executable next to these files and run it with 32 MPI ranks. This number must equal `mzetamax` (8) × `npartdom` (4).
`submit.sh` is an example Slurm script.

## Without the Poisson solve
Change these lines in `INPUT.d`:
```
Flag_Poisson=0
tstep=20e-9
msnap=50, ndiag=50
mstep=500
```
These values keep the run at 10 µs, with one snapshot every 1 µs.
Without the field solve, Φ = 0 and f_AD = 0, so full-h equals full-f.

## Output
Snapshots are written to `snapshots_hdf5/snapXXXXX.h5` every `msnap` steps (1 µs). `../../postprocessing` lists the main variables and has a plotting script.
