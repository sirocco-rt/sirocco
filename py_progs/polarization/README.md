# Polarization prototype

A Python prototype of polarized Monte Carlo radiative transfer, written to be
ported into Sirocco. The photon structure mirrors `struct photon` in
`source/sirocco.h`, and routines keep Sirocco's names where there is an
equivalent (`init_dummy_phot`, `stuff_phot`, `lorentz_transform`, ...).

The conventions (Stokes basis, sky frame, scattering matrix, peel-off, frame
transforms) are in the note *Sirocco-2026-Polarization-Conventions*. In
summary:

- `w` is I. Each photon carries `stokes` = (q, u, v) = (Q/I, U/I, V/I) relative
  to `pol_ref` = e1, a unit vector perpendicular to `lmn`. e2 = `lmn` × e1.
- The sky frame follows the astronomical (IAU 1974) convention: North is the
  projected system +z axis, East = n_obs × North, Q > 0 is a North–South
  E-vector, and the position angle runs from North through East.
- Electron scattering only, in the Thomson limit. Klein–Nishina is planned.
  All other processes leave photons unpolarized.

## Modules

| file | contents |
|---|---|
| `photon.py` | `photon_dtype`, Sirocco enum values, `make_photons`, `init_dummy_phot`, `set_unpolarized`, `stuff_phot` |
| `stokes.py` | Stokes rotation, observer direction and sky North/East (as in `spectra.c`), P and PA, the Lorentz transform of `pol_ref` |
| `es_scatter.py` | Thomson matrix, sampling of polarized scattered directions (`scatter_es`), peel-off toward an observer (`peel_es`) |
| `toy_transport.py` | Step 2: vectorized transport through uniform ellipsoids (point source) and a slab lit from below, with peel-off (`Spectrum`) and live-or-die (`LiveOrDie`) tallies, and the Brown & McLean (1977) prediction |
| `milne.py` | Deterministic solution for polarized electron scattering in a slab, the reference for the slab test (limb polarization 11.70%; Chandrasekhar: 11.71%) |
| `step2_validate.py` | Runs the Step 2 validation cases and optionally plots them |

All routines work on whole arrays of photons. Use a slice, `p[i:i+1]`, for a
single photon.

## Step 2 validation

```bash
python py_progs/polarization/step2_validate.py -plot step2.png    # ~1.5 min; -n 4 for 4x the photons
```

The cases:
- A uniform sphere must give no net polarization.
- Optically thin oblate and prolate ellipsoids must match Brown & McLean (1977): Q > 0 for oblate and Q < 0 for prolate.
- In a thicker ellipsoid, peel-off and live-or-die must agree.
- A tau = 10 slab must match `milne.py` in both intensity and polarization.

Peel-off errors come from the scatter between batches of photons. A photon's contributions at successive scatters are correlated, so a simple sum of squares underestimates the errors.

## Tests

```bash
python -m pytest py_progs/polarization/tests
```

Requires numpy, scipy and pytest.
