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
- Electron scattering in the Thomson limit (Klein–Nishina is planned) and
  resonance-line scattering. All other processes leave photons unpolarized.

## Modules

| file | contents |
|---|---|
| `photon.py` | `photon_dtype`, Sirocco enum values, `make_photons`, `init_dummy_phot`, `set_unpolarized`, `stuff_phot` |
| `stokes.py` | Stokes rotation, observer direction and sky North/East (as in `spectra.c`), P and PA, the Lorentz transform of `pol_ref` |
| `es_scatter.py` | Thomson matrix, sampling of polarized scattered directions (`scatter_es`), peel-off toward an observer (`peel_es`) |
| `toy_transport.py` | Step 2: vectorized transport through uniform ellipsoids (point source) and a slab lit from below, with peel-off (`Spectrum`) and live-or-die (`LiveOrDie`) tallies, and the Brown & McLean (1977) prediction |
| `milne.py` | Deterministic solution for polarized electron scattering in a slab, the reference for the slab test (limb polarization 11.70%; Chandrasekhar: 11.71%) |
| `step2_validate.py` | Runs the Step 2 validation cases and optionally plots them |
| `line_scatter.py` | Step 3: W2(J_l, J_u) from 6j symbols, the resonance-line phase matrix, `scatter_line`, `peel_line` |
| `sobolev_flow.py` | Step 3: transport through an expanding ellipsoid (homologous, or a radial power-law flow with direction-dependent Sobolev escape) with electron scattering and one Sobolev line. It has three treatments of trapping in the resonance region (`explicit`, `single`, `effective`), and spectra resolved in frequency |
| `step3_validate.py` | Runs the Step 3 validation cases and the line-effect spectra, and optionally plots them |
| `step3b_validate.py` | Step 3b: trapping when the Sobolev escape probability depends on direction (radial power-law flow) |

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

## Step 3 validation

```bash
python py_progs/polarization/step3_validate.py -plot step3.png    # ~2.5 min
```

Main result: trapping in a Sobolev resonance region depolarizes. Each local
scatter keeps (7/10) W2 of the polarized (K = 2) part of the radiation. So the
line-scattered polarization is F = beta / (1 - (1 - beta)(7/10) W2) times that
of a single scatter, where beta is the escape probability. A single scatter with
W2_eff = W2 F (the `effective` mode) reproduces exact trapping (`explicit`).
A plain single scatter (`single`, what Sirocco's line scattering would do now)
overestimates the polarization of thick lines by up to 1/F.

## Step 3b: escape probability that depends on direction

```bash
python py_progs/polarization/step3b_validate.py -plot step3b.png   # several minutes
```

In a radial flow v ~ r^alpha, the Sobolev optical depth varies with
direction as 1/(1 - (1 - alpha) mu^2). The results:

- For scattered continuum, `effective` (W2 F with the direction-averaged beta)
  still tracks `explicit`, and `single` is again much too high.
- Explicit trapping also *polarizes unpolarized line emission*: q ~ 2% at
  tau ~ 1-3, falling roughly as 1/tau at large tau, with opposite signs for
  radial (alpha < 1) and tangential (alpha > 1) trapping. Neither `single`
  nor `effective` can produce this, because both emit line photons unpolarized.

## Tests

```bash
python -m pytest py_progs/polarization/tests
```

Requires numpy, scipy and pytest.
