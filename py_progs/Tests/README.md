# Tests

Tests that run Sirocco and compare the results with an answer known
independently of Sirocco: an analytic result, a published one, or a separate
calculation. Each kind of test has its own subdirectory, with a README that
says how to run it.

| directory | what it tests |
|---|---|
| `polarization/` | the polarization option (`@Spectrum.polarization`) |

These are not the only tests:

- the unit tests in `source/tests` check individual routines;
- the regression tests, models in `examples/regress` run with
  `py_progs/regression.py`, show whether results change between versions,
  not whether they are right.

The models in `examples/` are for users to run and adapt.
