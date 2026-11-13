# Polarization examples

Models that use Sirocco's polarization option (`@Spectrum.polarization yes`,
which needs `sirocco -d`). The tests of the polarization calculation itself,
against known answers, are in `py_progs/Tests/polarization`; these are
examples of models to run, without a known answer.

## torus_wedge: an AGN torus in X-rays

The "wedge torus" of Vander Meulen et al. (2024, A&A 689; arXiv:2409.02986),
who used it for their X-ray polarization models: gas of constant density
between an inner and an outer radius, within θ_t of the equatorial plane,
around a central point source.

| file | contents |
|---|---|
| `make_torus.py` | writes the density structure, a polar grid that Sirocco imports |
| `torus_wedge.model.pf` | the parameters of the torus, read by `make_torus.py` |
| `torus_wedge.import.txt` | the grid, written by `make_torus.py torus_wedge.model.pf` |
| `torus_wedge.pf` | the Sirocco parameter file |

**The torus** (`torus_wedge.model.pf`):

- r = 10¹⁷ – 10¹⁸ cm, within 30° of the equator, so the covering factor is
  sin 30° = 0.5. Vander Meulen et al. write the covering factor as cos θ,
  with θ = 60° the half-opening angle measured from the pole.
- N_H = 10²⁴ cm⁻² along every line of sight from the source through the
  torus, so n_H = 1.1×10⁶ cm⁻³ and τ_es ≈ 0.8.
- 30 × 30 cells across the torus, logarithmic in r.
- It rotates at the circular speed about a 10⁸ M☉ black hole. Sirocco
  cannot run a wind with exactly zero velocity.

`make_torus.py` creates only the density structure: Sirocco calculates the
ionization of the gas from the radiation field. Vander Meulen et al. treat
the torus as cold, neutral gas, with photo-absorption, fluorescence and
scattering by bound electrons, so this is a Sirocco model of their
geometry, not a reproduction of their spectra.

**The Sirocco run** (`torus_wedge.pf`):

- A 10⁸ M☉ black hole with a power-law source, L(2–10 keV) =
  3×10⁴² erg/s, L_ν ∝ ν^−0.8, i.e. a photon index of 1.8, the reference
  case of Vander Meulen et al. (`Central_object.power_law_index` is α in
  L_ν ∝ ν^α, so the photon index is 1 − α). No disk.
- Macro atoms for H and He (`h20_hetop_standard80.dat`), with the ionization
  from `matrix_pow`.
- 10⁶ photons per cycle, 20 ionization and 10 spectral cycles. Spectra
  from 1 to 100 Å (0.12–12 keV), at inclinations of 10° and 55° (above
  the torus) and 65° and 80° (through it).

On 8 processes the run takes about 40 minutes. About 87% of the cells
converge, by cycle 15; more ionization cycles do not raise this. The
cells that do not converge are almost all in the outer 40% of the torus
(in radius), where heating and cooling have not come into balance. They
have little effect on the spectra: between 10 and 20 ionization cycles,
P changed by 0.2 percentage points or less at every inclination.

The polarization over 1–100 Å is about 0.07% at 10°, 1.0% at 55°, 3.1%
at 65° and 3.8% at 80°, with the electric vector along the projected
axis of the torus. It rises with energy at every inclination: at 80°,
for example, from 0.7% at 0.5–2 keV to 3.3% at 2–5 keV and 4.5% at
8–12 keV. Through the torus (65°, 80°) the spectrum is absorbed below
about 2 keV.

**To run it**, in a directory of your own (the run writes many files):

```
cp $SIROCCO/examples/polarization/torus_wedge.* .
Setup_Sirocco_Dir
mpirun -np 8 sirocco -d torus_wedge
python $SIROCCO/py_progs/run_check.py torus_wedge
```

`run_check.py` shows the spectra with their polarization. Because the
spectra are X-ray spectra, it plots them against energy, as E F_E, from
`torus_wedge.log_spec` (`plot_pol.py -xray`).

**To change the torus**, edit a copy of `torus_wedge.model.pf` (or give
`make_torus.py` a new name, and it will ask for each parameter), run
`python $SIROCCO/examples/polarization/make_torus.py newname.model.pf`,
and point `Wind.model2import` in the Sirocco .pf at `newname.import.txt`.
