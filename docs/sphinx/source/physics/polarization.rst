Polarization
------------

SIROCCO can follow the linear polarization of photon packets and produce Stokes Q and U spectra for each
observer, alongside the usual (intensity) spectra.  This is an option of the spectral cycles only; it has no
effect on the ionization cycles, and when it is off the code behaves exactly as it does without it.

At present the only process that polarizes light is **electron scattering**.  All other emission and
scattering processes, including line scattering, are treated as producing unpolarized light (see
:ref:`pol_not_included`).  Even so, the polarization changes across line features, because unpolarized line
photons dilute the polarized continuum, and line opacity removes some of the polarized, scattered light.

Polarization is switched on with the advanced parameter :ref:`Spectrum.polarization`, so SIROCCO must be run
with the ``-d`` flag.


Physics
=======

Stokes parameters and conventions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The polarization of a beam is described by the Stokes parameters :math:`(I, Q, U, V)`.  :math:`Q` and
:math:`U` describe linear polarization and are measured with respect to a reference axis perpendicular to the
direction of propagation: :math:`Q > 0` means the electric vector lies along the reference axis, :math:`Q < 0`
perpendicular to it, and :math:`\pm U` at :math:`\pm 45^\circ`.  :math:`V` describes circular polarization.
The degree and position angle of the linear polarization are

.. math::

    P = \frac{\sqrt{Q^2 + U^2}}{I}, \qquad \chi = \frac{1}{2} \arctan \frac{U}{Q}.

If the reference axis is rotated by an angle :math:`\alpha` (towards the second axis of the basis), the
Stokes parameters transform as

.. math::

    Q' = Q \cos 2\alpha + U \sin 2\alpha, \qquad U' = -Q \sin 2\alpha + U \cos 2\alpha ,

and :math:`I`, :math:`V` and :math:`P` are unchanged.

SIROCCO reports :math:`Q` and :math:`U` on the sky in the **IAU convention**.  The reference direction is
North, defined as the projection of the system's symmetry (+z) axis onto the sky, and position angles increase
from North through East.  So :math:`Q > 0` means the electric vector is **parallel to the projected symmetry
axis**, and :math:`Q < 0` perpendicular to it.  For a system that is symmetric about the z axis, :math:`U`
should be zero, apart from noise.

Electron scattering
~~~~~~~~~~~~~~~~~~~

For scattering by free electrons at rest, the Stokes parameters of the scattered light are related to those of
the incident light by a scattering matrix.  Taking the reference axis to be the normal to the scattering plane
(the plane that contains the incoming and outgoing directions), and with :math:`\Theta` the scattering
angle, the matrix for linear polarization is (Fano 1949)

.. math::

    \begin{pmatrix} I \\ Q \\ U \end{pmatrix}_{\rm out}
    \propto \epsilon^2
    \begin{pmatrix}
      \epsilon + \epsilon^{-1} - \sin^2\Theta & \sin^2\Theta & 0 \\
      \sin^2\Theta & 1 + \cos^2\Theta & 0 \\
      0 & 0 & 2\cos\Theta
    \end{pmatrix}
    \begin{pmatrix} I \\ Q \\ U \end{pmatrix}_{\rm in}

where :math:`\epsilon = \nu_{\rm out}/\nu_{\rm in} = 1 / [1 + x(1-\cos\Theta)]` and :math:`x = h\nu_{\rm
in}/m_e c^2`.  In the Thomson limit, :math:`\epsilon = 1` and this is the familiar Rayleigh/Thomson matrix
(Chandrasekhar 1960), with :math:`P_{11} = 1 + \cos^2\Theta`.

Some consequences:

* Unpolarized light scattered through :math:`\Theta` acquires a polarization
  :math:`P = \sin^2\Theta / (1 + \cos^2\Theta)`, perpendicular to the scattering plane.  Light scattered
  through :math:`90^\circ` is completely polarized.

* For polarized incident light, the scattered intensity depends on the azimuth of the scattering plane: light
  is preferentially scattered perpendicular to the incident electric vector.  This modulation is what the
  polarized phase function adds to the unpolarized one; the distribution of :math:`\Theta` itself, averaged over
  azimuth, is the same as for unpolarized light.

* The Klein-Nishina (Compton) case differs from Thomson only in :math:`P_{11}`, by an extra unpolarized term
  :math:`\epsilon + \epsilon^{-1} - 2 \ge 0`.  This reduces the polarization, to
  :math:`P = \sin^2\Theta / (\epsilon + \epsilon^{-1} - \sin^2\Theta)`, and moves the maximum forward of
  :math:`90^\circ`.  At :math:`\Theta = 90^\circ` the fractional reduction is about :math:`x^2/(1+x)`: less than
  0.1% of :math:`P` below about 30 keV, a few per cent at 100 keV, and a third at 511 keV.

* Unpolarized electrons do not create circular polarization from unpolarized or linearly polarized light, so
  :math:`V = 0` throughout, and SIROCCO does not follow it.

For an optically thin, axisymmetric distribution of electrons illuminated by a point source, the net
polarization is given analytically by Brown & McLean (1977).  Optically thick, multiply scattering envelopes
are treated by, e.g., Hillier (1994), whose models are a standard test (see :ref:`pol_validation`).

Moving electrons and special relativity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Scattering takes place in the rest frame of the scattering electron, and the light is then transformed back
to the observer frame.  Under a Lorentz boost a photon's direction is aberrated, so its polarization reference
axis must be transformed as well.  The degree of polarization is a Lorentz invariant, but the plane of
polarization rotates.  The new polarization direction follows from the transformation of the electric
field of the wave,

.. math::

    \mathbf{E}' = \gamma \left( \mathbf{E} + \boldsymbol{\beta} \times \mathbf{B} \right)
                 - \frac{\gamma^2}{\gamma + 1} \boldsymbol{\beta} (\boldsymbol{\beta} \cdot \mathbf{E}),
    \qquad \mathbf{B} = \hat{\mathbf{n}} \times \mathbf{E} .

A useful property, which SIROCCO's unit tests check: the polarization angle measured from the normal to the
plane that contains the photon direction and the velocity, :math:`\hat{\mathbf{n}} \times \mathbf{v}`, is
unchanged by the boost.

This matters both for the bulk motion of the wind and for the thermal motions of the electrons: in a hot
plasma, scattering in each electron's rest frame and transforming back reduces the polarization.

Line features
~~~~~~~~~~~~~

Even when the line emission is itself unpolarized, line formation changes the polarized spectrum
(e.g. Bulla, Sim & Kromer 2015):

* Unpolarized line photons dilute the polarized continuum, so :math:`P` falls across emission lines.

* Electron-scattered continuum light that passes through the line-forming region on its way out can be
  absorbed in the line, which reduces :math:`Q` and :math:`U` as well.

* In absorption troughs, the direct, unpolarized light is removed, while scattered light from other parts of
  the wind partly fills in the trough, so :math:`P` rises.  Together these give the "inverted P Cygni" profile
  often seen in polarization spectra.

SIROCCO's transport includes all three effects.


Implementation
==============

When polarization is followed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With :ref:`Spectrum.polarization` set to ``yes``, the flag ``modes.polarization`` is set (``setup.c``).
Polarization is followed only in the spectral cycles (``pol_active()`` in ``polarization.c``).  In the
ionization cycles electron scattering is done exactly as before, so the ionization and temperature structure of
the wind does not depend on the option.

The photon packet
~~~~~~~~~~~~~~~~~

Each photon packet carries two extra numbers, ``q`` :math:`= Q/I` and ``u`` :math:`= U/I`, where the packet's
weight ``w`` plays the role of :math:`I`.  The reference axis is **not stored**.  It is always a default axis
computed from the photon's direction :math:`\hat{\mathbf{n}}` (``pol_basis``):

.. math::

    \hat{\mathbf{e}}_1 = \frac{\hat{\mathbf{z}} \times \hat{\mathbf{n}}}{|\hat{\mathbf{z}} \times \hat{\mathbf{n}}|},
    \qquad \hat{\mathbf{e}}_2 = \hat{\mathbf{n}} \times \hat{\mathbf{e}}_1 ,

with :math:`\hat{\mathbf{x}}` in place of :math:`\hat{\mathbf{z}}` within :math:`10^{-6}` of the poles, so that
:math:`(\hat{\mathbf{e}}_1, \hat{\mathbf{e}}_2, \hat{\mathbf{n}})` is right-handed.  For a packet, ``q`` > 0 means
the electric vector is along :math:`\hat{\mathbf{e}}_1`.  Whenever the direction of a packet changes, ``q`` and
``u`` are re-expressed against the default axis for the new direction (``pol_to_axis``, ``pol_from_axis``,
``pol_rotate``).  Keeping the axis implicit keeps the photon structure small (two doubles) and avoids having to
keep a stored axis perpendicular to the direction.

All sources (star, disk, boundary layer, AGN and wind emission) create unpolarized packets
(``photon_gen.c``).

Electron scattering
~~~~~~~~~~~~~~~~~~~

Electron scattering in SIROCCO is done in ``compton_scatter`` (``compton.c``), which draws a velocity for the
scattering electron from the local temperature and transforms the packet into the electron's rest frame.
``compton_dir`` then chooses the scattering angle :math:`\Theta` as before, from the dipole distribution
for :math:`x < 10^{-4}` (about 51 eV) and from the Klein-Nishina distribution above that.  When polarization
is active:

1. :math:`\cos\Theta` is drawn from the unpolarized distribution, which is the exact marginal distribution of
   :math:`\Theta` for polarized light as well.  In the Thomson case the cumulative distribution is inverted
   exactly (``pol_thomson_cos``).

2. The azimuth :math:`\phi` of the normal to the scattering plane, :math:`\hat{\mathbf{e}}_\perp`, measured from
   :math:`\hat{\mathbf{e}}_1`, is drawn from its conditional distribution,

   .. math::

       p(\phi \mid \Theta) \propto 1 + a \left( q \cos 2\phi + u \sin 2\phi \right),
       \qquad a = \frac{\sin^2\Theta}{\epsilon + \epsilon^{-1} - \sin^2\Theta},

   by rejection (``pol_scatter_es``).  The acceptance rate is at least 50%.

3. The new direction is :math:`\cos\Theta\, \hat{\mathbf{n}} + \sin\Theta\, (\hat{\mathbf{n}} \times
   \hat{\mathbf{e}}_\perp)`.  The incoming ``q``, ``u`` are rotated to the :math:`\hat{\mathbf{e}}_\perp` axis,
   the matrix above is applied (``pol_es_matrix``), and the result is re-expressed against the default axis of
   the new direction.

The Thomson and Klein-Nishina cases use the same code, with :math:`\epsilon = 1` below the switch and
:math:`\epsilon` from the Klein-Nishina energy change above it.  The cost of including Klein-Nishina in the
polarization is one division per scatter.

Frame transformations
~~~~~~~~~~~~~~~~~~~~~

Every Lorentz transformation of a packet goes through ``lorentz_transform`` (``frame.c``), which calls
``pol_lorentz`` (``polarization.c``).  This transforms the old reference axis as the electric field (equation
above), projects it perpendicular to the new direction, and re-expresses ``q`` and ``u`` against the default
axis of the new direction.  This covers the transformations between the observer frame and the frame of the
wind or disk, and into and out of the rest frame of a thermal electron.

With ``-nonrel`` (``REL_MODE_LINEAR``) directions are not aberrated, so neither are the polarization axes; the
intensity and the polarization are then both treated to first order in :math:`v/c`.

Note that SIROCCO treats model densities as observer-frame values, and uses :math:`\rho_{\rm cmf} =
\rho/\gamma` (see :doc:`relativity_comoving`).

Extracted (peel-off) spectra
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the usual extract mode, at each electron scatter a copy of the packet is sent towards each observer, with
its weight multiplied by the phase function for that direction (``extract.c``).  With polarization on, this
weight is multiplied by the ratio of the polarized to the unpolarized phase function,
:math:`(P_{11} + \sin^2\Theta \, q_\perp) / P_{11}`, and the extracted packet's ``q``, ``u`` are set from the
scattering matrix (``pol_peel_es``).  The directions and weights are handled in the electron's rest frame and
then transformed, as described above.

The extracted packet's contribution, :math:`w e^{-\tau}`, is added to the I spectrum as usual, and
:math:`w e^{-\tau} q_s` and :math:`w e^{-\tau} u_s` to the Q and U spectra, where :math:`q_s`, :math:`u_s` are the
packet's polarization re-expressed against North on the sky for that observer (``pol_tally``).  Live-or-die
spectra are handled the same way, using North for each escaping packet's direction.

Processes treated as unpolarized
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following set ``q = u = 0`` for the packet that leaves them:

* resonance-line scattering, in both the simple-atom and macro-atom schemes, including thermal trapping
  (``resonate.c``);
* free-free and bound-free emission, k-packets, and all macro-atom emission;
* reflection from the disk or star (``trans_phot.c``).

Peel-offs from these events are therefore unpolarized as well, but the transport of all packets, polarized or
not, includes the line opacity.

Output
~~~~~~

When polarization is on, two columns are added for each observer to the ``.spec``, ``.spec_wind``,
``.log_spec`` and ``.log_spec_wind`` files, after all of the existing columns, e.g. ``A80P0.50_Q`` and
``A80P0.50_U``.  Q and U are in the same units as the intensity columns, and are referred to North on the sky
(IAU convention, above).  So, for each wavelength bin,

.. math::

    q = Q/I, \quad u = U/I, \quad P = \sqrt{q^2 + u^2}, \quad
    \chi = \tfrac{1}{2}\,\mathrm{atan2}(U, Q) \;\; \text{(from North through East)} .

The Q and U spectra are combined across MPI processes with the other spectra, and are saved in the
``.spec_save`` file, so that a restarted run continues to accumulate them.  If Q and U are to be integrated over a
band, integrate Q, U and I separately and then form the ratios.

Because the polarization of most systems is small, Q and U are noisy; integrating over bands, or combining
bins, is usually necessary.  For an axisymmetric model the U spectrum is a useful estimate of the noise.  The
advanced option :ref:`Spectrum.select_specific_no_of_scatters_in_spectra` can be used to separate the
contributions of unscattered, singly scattered and multiply scattered light.

``py_progs/polarization/c_spectra.py`` reads these columns, prints band-integrated :math:`q`, :math:`u`,
:math:`P` and :math:`\chi` for each observer, and plots I, q and u against wavelength.

Cost
~~~~

With the option off, the only change is two extra doubles in each photon.  With it on, each electron scatter
and each peel-off of a scattered packet does a few rotations of the Stokes parameters and, for scatters, a
short rejection loop.  Memory for the Q and U spectra is allocated only when polarization is on.

.. _pol_validation:

Validation
==========

* **Unit tests** (``source/tests/tests/test_polarization.c``): the basis and rotations, North on the sky,
  the scattering matrix (Thomson and Klein-Nishina), the peel-off weights, and the Lorentz transformation of
  the polarization, including the invariance described above.

* **Polarization off**: the regression spectra are identical to those of a build without the option,
  including the random-number sequence.

* **Point source and thin envelopes**: the net polarization of optically thin, axisymmetric envelopes agrees
  with Brown & McLean (1977) to about 1% of the maximum.  Spherical models give zero polarization, and
  axisymmetric models give :math:`U \approx 0`.

* **Hillier (1994) prolate shell**, the test used by Bulla et al. (2015): pure electron scattering for
  angle-averaged optical depths from 0.03 to 3, and as a function of albedo.  SIROCCO agrees with Hillier to
  about 0.1% in :math:`P`, about the accuracy to which his figures can be read, including the change of sign of
  the polarization at large optical depth.

* **Special relativity**: an optically thin shell moving radially at :math:`0.1c` and :math:`0.3c`
  agrees with an exact single-scattering calculation to about 1-2% in the scattered I and Q, although at
  :math:`0.3c` the relativistic effects halve the scattered Q.  Runs with ``-nonrel`` reproduce the static
  result.

.. _pol_not_included:

What is not included
====================

* **Intrinsic line polarization.** Resonance scattering in a line is not, in general, depolarizing: the
  scattering matrix is a mixture of a Rayleigh-like part, with weight :math:`W_2` that depends on the angular
  momenta of the levels involved (Hamilton 1947), and an isotropic, unpolarized part.  SIROCCO currently treats
  all line scattering as unpolarized (the :math:`W_2 = 0` limit, as in Bulla et al. 2015).  How best to include
  it, in particular for macro-atoms and for line photons that scatter many times locally before escaping, is
  still under consideration.

* Rayleigh scattering by atoms and molecules, scattering by dust, and polarization of the light from the
  disk or star themselves.

* Circular polarization, and magnetic effects such as the Hanle effect.

* General-relativistic rotation of the plane of polarization near compact objects.

References
==========

* Brown, J. C. & McLean, I. S. 1977, A&A, 57, 141
* Bulla, M., Sim, S. A. & Kromer, M. 2015, MNRAS, 450, 967
* Chandrasekhar, S. 1960, *Radiative Transfer* (Dover)
* Fano, U. 1949, J. Opt. Soc. Am., 39, 859
* Hamilton, D. R. 1947, ApJ, 106, 457
* Hillier, D. J. 1994, A&A, 289, 492
