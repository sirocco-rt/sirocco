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

All routines work on whole arrays of photons. Use a slice, `p[i:i+1]`, for a
single photon.

## Tests

```bash
python -m pytest py_progs/polarization/tests
```

Requires numpy, scipy and pytest.
