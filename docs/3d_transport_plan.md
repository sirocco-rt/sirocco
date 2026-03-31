# Plan: Breaking Z-Symmetry as a Path to 3D Photon Transport

## Background

Sirocco currently exploits bi-conical symmetry: the wind is assumed identical above and
below the disk plane (z=0). Photons are implicitly reflected to positive z via `fabs(x[2])`
in the `where_in_grid` functions of each coordinate system. Only one hemisphere of wind
cells is stored in `wmain[]`.

This document describes a staged implementation to break that symmetry as the first step
toward a general one-to-many transport scheme, where many geometric transport cells map
to a single plasma cell.

## Long-term vision

The eventual goal is a transport grid that is denser or more complex than the plasma grid:

    Transport cells (geometry)  →  Plasma cells (physics)
         many                   →       one

For example: 8 azimuthal sectors all mapping to one azimuthally-averaged plasma cell.
This allows genuinely 3D photon paths while keeping the plasma physics tractable in memory.
A 100 x 100 x 8 transport grid (80,000 cells) is within the current ~1e5 cell memory limit.

## Key architectural invariant

The existing design already supports this cleanly:

    p->grid  →  wmain[p->grid].nplasma  →  plasmamain[nplasma]

`p->grid` is always a wind cell index. The wind cell contains `nplasma`, which maps to
the plasma cell. Estimators, ionization, and MPI all follow the `nplasma` indirection and
need no changes. The one-to-many mapping is expressed entirely through `nplasma` in the
new transport-only wind cells.

No new field is needed on the photon struct.

## Where z-symmetry is currently encoded

The symmetry is implemented at several levels — not in one place:

| File | Mechanism |
|---|---|
| `cylindrical.c:466`, `rtheta.c:515`, `wind_util.c:127` | `z = fabs(x[2])` before grid lookup |
| `cylindrical.c:96-100` | Flips z-wall signs in `ds_in_cell` for z < 0 |
| `phot_util.c:231-235` | Reflects photon to northern hemisphere before cone geometry |
| `gradv.c:74-78` | Reflects photon and direction to upper hemisphere for gradient |
| `sv.c:227-228`, `knigge.c:306-307` | Wind models flip `v[2]` when `x[2] < 0` — already correct |
| `wind2d.c:209-210` | Flips `v_z` during interpolation for z < 0 — already correct |
| `estimators_simple.c:281-284, 410-411, 429-430` | Flips z-component of flux estimators |

Note: velocity is already handled correctly for z < 0 in the wind models and interpolation.
The reflection hacks in `gradv.c` and `phot_util.c` are workarounds that can be removed once
proper lower-hemisphere wind cells exist.

Cell volumes currently include both hemispheres (factor of 2 in `cylind_cell_volume:349`).
This remains correct for the 1→2 case since both transport cells share one plasma cell.

## Implementation steps

### Step 0: Establish baseline

Run a standard test model (CV with SV wind). Save output spectra and key plasma diagnostics
(electron temperature, ionization fractions in several cells). Use a fixed random seed.
Every subsequent step is compared against this reference.

---

### Step 1: Extend `wmain` to cover both hemispheres

Create `source/transport_grid.c` with `make_transport_grid(ndom)`, called from
`define_wind.c` after the plasma grid is built.

- Reallocate `wmain` from N to 2N cells for the domain
- For each lower hemisphere cell `k + N`:
  - Copy all fields from `wmain[k]`
  - Negate `x[2]` and `xcen[2]`
  - Set `nplasma = wmain[k].nplasma` (same plasma cell)

No behavior change yet — nothing sets `p->grid` to the new cells.

**Test:** Identical output to Step 0.

---

### Step 2: Modify `where_in_grid` to return the correct hemisphere cell

In `cylindrical.c`, `rtheta.c`, `spherical.c` (and `cylind_var.c`, `import_cylindrical.c`):

- Remove `z = fabs(x[2])` / `theta = acos(fabs(x[2]/r))`
- Compute the upper-hemisphere cell index k as before
- Return `k + N` when `x[2] < 0`, `k` when `x[2] >= 0`

`p->grid` now points to the geometrically correct upper or lower wmain cell.
`wmain[p->grid].nplasma` still resolves to the same plasma cell either way.

**Test:** Identical output to Step 0. Add a diagnostic counting photons with
`p->grid >= nstart + N` — should be ~50% for a symmetric biconical wind.

---

### Step 3: Remove the reflection hacks

These were workarounds for the missing lower hemisphere cells. Remove one at a time:

**3a.** `gradv.c:74-78` — remove reflection of photon position and direction before
gradient calculation. The on-the-fly path uses `model_velocity` which already handles
z < 0 correctly.

**Test:** Statistically identical output to Step 0.

**3b.** `phot_util.c:231-235` — remove reflection before cone intersection geometry.
The cones are defined symmetrically about z=0 (via `dzdr` and `z`), so the geometry
is identical with or without the reflection.

**Test:** Statistically identical output to Step 0.

At this point the code is clean: no fabs reflections, no hemisphere hacks, and output
matches the original.

---

### Step 4: Demonstrate asymmetry capability

Add a test hook that assigns different properties to upper and lower hemisphere wind cells
(e.g., scale density by a factor in lower cells). Verify:

- Spectra seen from inclinations above and below the disk plane differ
- The difference scales correctly with the density factor
- Photon counts per hemisphere in diagnostics are consistent

This confirms the infrastructure is correct and provides a regression test for future work.

---

## What does NOT change

- `p->grid` semantics — still always a wmain index
- `wmain[n].nplasma` usage — the existing indirection to plasma cells
- All estimator code — uses `nplasma` indirection already
- All ionization/temperature solver code — operates on plasma cells only
- MPI communication — `communicate_plasma.c` operates on plasma cells only
- Wind model velocity functions (`sv.c`, `knigge.c`) — already correct for z < 0

## Extension to azimuthal sectors (future work)

Once the 1→2 infrastructure is validated, adding M azimuthal sectors follows the same
pattern:

- `wmain` grows from N to M*N cells
- Each sector cell copies geometry from the base cell, rotated by `k * 2*PI/M`
- `nplasma` in all M sector cells points to the same plasma cell
- `where_in_grid` adds a binary search in phi after the existing (rho, z) or (r, theta)
  search to select the correct sector cell

A grid of 100 x 100 plasma cells with M=8 azimuthal sectors gives 80,000 transport cells,
within the current ~1e5 cell memory budget.
