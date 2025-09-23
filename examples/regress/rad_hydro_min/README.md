# rad_hydro_min (regression example)

**Purpose:** Minimal rad-hydro model to flag output changes after code edits.

## Files
- `cv_rad_hydro.pf` — Sirocco parameter file.
- `cv_rad_hydro.pluto` — single PLUTO snapshot (final cycle) referenced by `cv_rad_hydro.pf`.

## How to run (local, optional)
In this example, the number of photons is set to 1e7 (10,000,000); you may change this if needed. The number of ionisation cycles is set to 30, although in our PLUTO–Sirocco rad-hydro calculations we typically use two cycles per call to Sirocco.

Build Sirocco (see the repo README). From the repository root:

- **Quick run (MPI):**
  ```bash
  mpirun -np <N> sirocco examples/regress/rad_hydro_min/cv_rad_hydro.pf > sirocco_log.txt 2>&1
