# Polarization tests

These tests check Sirocco's polarization option (`@Spectrum.polarization yes`,
which needs `sirocco -d`) against calculations made independently of
Sirocco. Each test is a full Sirocco run, compared with an analytic result, a
published result or an independent calculation, not with an earlier version
of Sirocco.

- The unit tests in `source/tests` check individual routines.
- The regression tests in `examples/regress` check that results don't change
  between versions. They include `pol_cone`, a short version of the cone test
  below.
- These tests check that the polarization is right.

Run them after changing the electron scattering, the frame transformations,
the extraction (peel-off) or the spectra.

## Running the tests

Run the tests in a working directory **outside the repository**. The driver
copies everything it needs there, and writes nothing in `examples/`.

```
export SIROCCO=/path/to/Sirocco
python $SIROCCO/examples/polarization/pol_tests.py ~/pol_test_runs                   # quick tests, ~8 min
python $SIROCCO/examples/polarization/pol_tests.py -full ~/pol_test_runs_full        # full tests, several hours
```

Options:

| option | meaning |
|---|---|
| `-np N` | number of MPI processes (default 8; 1 runs without `mpirun`) |
| `-exe NAME` | the Sirocco executable to test (default `sirocco`) |
| `-w2t NAME` | a `windsave2table` built from the same code as `-exe` (default `windsave2table`). Only the cvcone test needs it |
| `-tests T ...` | run only some of `brown_mclean hillier sr_shell cvcone` |
| `-check_only` | don't run Sirocco; redo the analysis and checks on the results already in the working directory |

The driver needs Python with numpy and matplotlib, and `$SIROCCO`.

For each test, the driver:
- copies the test's files to `workdir/<test>/` and runs `Setup_Sirocco_Dir`
  there;
- writes the models, runs them with `-d`, and runs the analysis;
- checks the results against the tolerances in the test's `testdef.py`.

It writes `workdir/Summary.txt`. For each check, the summary gives the value,
the range the value must lie in, and PASS or FAIL. The summary starts with the
executable and its version and git hash, so that runs of two versions can be
compared. The model outputs and the comparison plots stay in
`workdir/<test>/`. The driver exits with status 0 only if every check passes.

The plots each test writes:

| test | plots |
|---|---|
| brown_mclean | `compare_bm.png`: q against inclination for Sirocco, single scattering and Brown & McLean, and the ratio that is checked |
| hillier | `compare_hillier.png` (full tier only): q against inclination for each optical depth, and against albedo, with Hillier's values |
| sr_shell | `compare_sr.png`: Sirocco / prediction for the scattered I and Q |
| cvcone | `es_vs_prediction.png`: q of the scattered and total light against wavelength, with the predictions; `root_pol.png` (`plot_pol.py`): the spectrum, P and PA at each inclination; `root_geompol.png`: the prediction's breakdown by source |

### Run times on 8 processes

| test | quick tier | full tier |
|---|---|---|
| brown_mclean | 1 model, ~2 min | 3 models, ~15 min |
| hillier | 1 model, ~3 min | 10 models, ~3 h (one of them takes 85 min) |
| sr_shell | 1 model, ~1 min | 5 models, ~15 min |
| cvcone | 1 model + prediction, ~2.5 min | 2 models + predictions, ~45 min (the wedge prediction alone takes ~25 min) |

## The tests

All four use electron scattering only. They test polarized Thomson
scattering, transport, peel-off and the frame transformations. Line
scattering is unpolarized in Sirocco and isn't tested here.

Each test directory holds:
- the script that writes the models;
- the independent prediction;
- scripts to compare and plot the results;
- `testdef.py`, which tells the driver what to run and what to check.

The docstring of each `testdef.py` explains its tolerances and the runs they
were set from. `h_only_master.dat` and `h_elem_ions.dat` are a hydrogen-only
atomic data set, used where electron scattering must be the only opacity.

### brown_mclean: thin envelopes around a point source

- **Model:** `make_bm.py root tau_eq` writes an envelope with
  n_e ∝ r⁻² sin²θ (10¹¹–10¹² cm) of pure H at 10⁶ K around a small star.
- **Prediction:** `single_scatter.py` computes single scattering exactly,
  with attenuation. In the thin limit this is Brown & McLean (1977),
  q = τ₀ (1 − 3γ) sin² i.
- **Check:** q (4500–6500 Å) at 50°, 70° and 90°, divided by the prediction,
  must lie between 0.96 and 1 + τ_eq + 0.04. Multiple scattering raises it
  above 1, by about 9% at τ_eq = 0.1. |u| must be below 2% of q(90°).
- **Tiers:** quick, τ_eq = 0.1; full, τ_eq = 0.1, 0.03 and 0.01.
- **Note:** τ_eq = 0.01 has given ratios of 0.97–0.985, slightly below 1.
  This is close to the limit, and not understood.

### hillier: Hillier (1994) prolate shell

The test of Bulla et al. (2015, §3.1).

- **Model:** a point-like star inside a detached shell with
  χ ∝ r⁻⁴ (1 + 10 cos²θ). τ_ave runs from 0.03 to 2.9 (`hx*`). Two runs
  (`hal90`, `hal56`) select photons by number of scatters, which gives the
  polarization as a function of albedo.
- **Prediction:** Hillier's values, read by eye from his Figs 2–5
  (`hillier_ref.py`).
- **Check:** q within 0.2% of Hillier at every inclination and optical
  depth. Against albedo: within 0.1% at 90° and 0.3% at 56.25°, where
  Hillier's own methods disagree. |u| < 0.1%.
- **Tiers:** quick, χ₀ = 0.1 only; full, everything, with plots in
  `compare_hillier.png`.

### sr_shell: special relativity

- **Model:** the Hillier shell, made optically thin and moving radially:
  - at β = 3×10⁻⁴, 0.1 and 0.3;
  - at 0.3 with `-nonrel`;
  - as a spherical shell at 0.3.

  The static case uses β = 3×10⁻⁴ because Sirocco fails at exactly zero
  velocity.
- **Prediction:** `sr_single.py`, an exact relativistic single-scattering
  calculation.
- **Check:** for the light scattered once, the ratios of Sirocco's I and Q
  to the prediction must lie between 0.96 and 1.01 (0.95 and 1.03 in the
  quick tier, which has fewer photons). The prediction leaves
  out the attenuation of the scattered light, so Sirocco is 1–2% low. The
  sphere must have |q| < 0.3%.
- **Tiers:** quick, 0.3c only; full, all five, with plots in
  `compare_sr.png`.

### cvcone: a CV-like wind in a cone or a wedge

- **Model:** a central source lights a wind of constant density with radial
  τ_es = 1. The wind is either in a polar cone at 0–20° (`pbl_es`) or in an
  equatorial wedge at 60–88° (`eq_es`). It has simple atoms and no wind
  emission, and Sirocco calculates its ionization. `make_cvcone.py` wrote
  both `.pf` files and their import files, and its defaults reproduce them.
- **Prediction:** `geom_pol.py` computes single scattering of the central
  source, with attenuation, using the n_e of the model. It reads n_e from
  `windsave2table`'s output, which is why this test needs `-w2t`.
- **Check:** q of the light scattered once, in the windows between the
  lines, must be within 1 point of the prediction (1.5 points in the quick
  tier).
- **Tiers:** quick, `pol_cone` (a small version of `pbl_es`); full, the cone
  at 10°, 45° and 80°, and the wedge at 80°.
- **Prediction noise:** the full tier makes one 30,000-ray prediction for each
  model, which is uncertain by about 0.2–0.3 points. The validation in the
  plan note used the mean of 5 seeds for the wedge, so single runs can differ
  from it by that much.
- **Note:** for the wedge at 10° and 45°, Sirocco is about 2 points more
  positive than the prediction. This is not understood, so these angles are
  reported (as INFO) but not checked.
- **Other scripts:** `analyze_split.py root` breaks the light down by number
  of scatters, and calculates τ with `sirocco_optd`.

## Cartoons of the geometries

Two scripts draw the test models, to scale, in a cross-section through the
symmetry axis. They show the density, the star, the observers, and the
direction of the E-vector of the scattered light:
- `shell_cartoons.py` draws brown_mclean, hillier and sr_shell, and writes
  `shell_cartoons.png`.
- `cvcone/cvcone_cartoon.py` draws the cone and the wedge, and writes
  `cvcone_cartoon.png`.

They need no runs and take a second each.

## Running a test by hand

Every script in the test directories can also be run on its own, in a copy
of the directory. The usage is in each docstring. For example:

```
mkdir ~/bm && cp $SIROCCO/examples/polarization/brown_mclean/* $SIROCCO/examples/polarization/h_*.dat ~/bm && cd ~/bm
Setup_Sirocco_Dir
python make_bm.py bm01 0.1
mpirun -np 8 sirocco -d bm01
python compare_bm.py bm01:0.1
```
