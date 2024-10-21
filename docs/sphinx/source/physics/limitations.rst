Limitations and Caveats
-------------------------

There are a number of limitations of SIROCCO as a code. This page is a non-exhaustive list, partly taken from the release paper.

Specific issues and code problems
======================================

**Pecular emission features in high optical depth models:** As described in `issue 659 <https://github.com/sirocco-rt/sirocco/issues/659>`_

**Free-free Gaunt factors:** As described in `issue 33 <https://github.com/sirocco-rt/sirocco/issues/33>`_, we use a free-free integrated Gaunt factor for both free-free opacity and cooling calculations, but one should really use a frequency-specific (photon by photon) Gaunt factor for the opacity calculation.

**Compton heating energy leak:** As described in `issue 295 <https://github.com/sirocco-rt/sirocco/issues/295>`_,
there is the potential for a slight energy leak (or poor statistics) for Compton heating in non-macro atom mode.

**Iron Fluorescence:** As described in `issue 499 <https://github.com/sirocco-rt/sirocco/issues/499>`_,
Fluorescence data is not dealt with properly in non-macro-atom mode. The user should use a Fe macro-atom data set for
situations where this is important.


Conceptual and Practical difficulties
======================================

SIROCCO is a complex code that includes a wide range of physical processes. As a result, it is capable of synthesizing the electromagnetic signatures of a diverse set of astrophysical systems on all scales. However, there are, of course, important limitations that the user should be aware of. Perhaps most importantly, SIROCCO is *not* a stellar atmosphere code: it is not designed for hugely optically thick and/or static media. Below, we explore these and some other constraints in more detail.

- **Extreme optical depths**

  Since SIROCCO is a Monte Carlo code that tracks the discrete interactions of simulated photon packets, the computational cost of a wind model scales with the characteristic optical depth, :math:`\tau_{ch}`, it presents to these packets. In most environments of interest, at least hydrogen is nearly fully ionized, so it is convenient to define :math:`\tau_{ch}` as the characteristic *electron scattering* optical depth. More specifically, the computational cost depends on the mean number of scatterings, :math:`\overline{N}_{scat}`, packets have to undergo before escaping. As a practical guide to what is feasible, we note that in some of our TDE model models :math:`\overline{N}_{scat}` is as high as `100 - 1000`, corresponding to :math:`\tau_{ch} \sim 10 - 100`. Such models take :math:`O(1000)` core-hours to converge and produce synthetic spectra.

  It is important to understand that optically thick media can give rise to a variety of symptoms. Moreover, geometry matters. In 2-D models, optical depths can be highly direction dependent, with :math:`\tau_{min} << \tau_{max}`. In this case, photon packets can and will escape preferentially along the most transparent direction(s), i.e. :math:`\tau_{ch} \sim \tau_{min}`. This is both good and bad. Good, because :math:`\tau_{ch} << \tau_{max}` implies that even models with very high :math:`\tau_{max}` are, in principle, computationally feasible. Bad, because such models often encounter serious convergence problems, as there are not enough photon packets crossing cells embedded deep within (or located behind) high-:math:`\tau` structures. In practice, unconverged cells in such regions can often be safely ignored, as they are bound to contribute little to the synthesized spectra. However, this decision has to be made by the user.

  In 1-D (spherical) models, photon packets can more easily get trapped, as there are no preferred directions along which they can escape. In practice, computational resources and/or convergence issues tend to be the limiting factor for optically thick 1-D models. However, in "simple atom" mode, there is an additional failure mode for such models. In the limit where most photons packets are (re)absorbed and (re)emitted many times in a single cell, it becomes extremely difficult for the code to maintain energy conservation. The issue is that the radiation field locally becomes nearly isotropic. The *net (outward) flux* can then be a tiny fraction of the *mean (angle-averaged) intensity*. In Monte Carlo terms, the luminosity (number of photon packets) that has to be created within a cell is far greater than the net number of packets that actually leaves the cell. The luminosity carried by the escaping population can then become highly uncertain simply due to Poisson scatter. Yet global energy conservation across a simulation relies on the correct flow rates of photon packets across cells. This is not a concern when SIROCCO is run in macro-atom mode, as this automatically enforces radiative equilibrium and hence guarantees global energy conservation.

- **The Sobolev and related approximations**

  SIROCCO is designed to synthesize spectra formed in *moving* media and treats line transfer in the Sobolev approximation. This means that a photon packet is only allowed to interact with a bound-bound transition at the exact location where its co-moving frequency coincides with the rest frequency of the transition. In reality, interactions take place across a spatial scale defined by the *Sobolev length*, :math:`\Delta s \simeq v_{th} \left(dv_{s}/ds\right)`, where :math:`v_{th}` is the local thermal speed, :math:`s` is distance along the photon packet's line of flight, and :math:`v_{s}` is the projected velocity along this direction. The Sobolev approximation therefore amounts to the assumption that the physical properties of the flow do not change significantly on scales as small as :math:`\Delta s`.

  The optical depth associated with the interaction is proportional to the *Sobolev length*, :math:`\Delta s \simeq v_{th} \left(dv_{s}/ds\right)`, where :math:`v_{th}` is the local thermal speed, :math:`s` is distance along the photon packet's line of flight, and :math:`v_{s}` is projected velocity along this direction. Since :math:`v_{th}` is usually close to the sound speed, the Sobolev approximation is sometimes also called the supersonic approximation.

  In line with this treatment, SIROCCO neglects any thermal or microturbulent broadening associated with bound-bound transitions. It is the responsibility of the user to ensure that the Sobolev approximation is valid for the models they run. Static media and flows in which thermal speeds are comparable to bulk velocities should not be simulated with SIROCCO.

- **The size of the line list**

  An important challenge for SIROCCO is to identify spectral lines that might be in resonance with a particular photon packet as it passes through a given cell. The computational cost of this task scales directly with the size of the line list. It is therefore not feasible to simply import the latest Kurucz line list (http://kurucz.harvard.edu/linelists/gfnew/), for example, which contains :math:`\simeq 10^9` transitions. Again, SIROCCO is *not* a stellar atmosphere code. The default line list we use contains :math:`\simeq 10^4` transitions, which is sufficient for many applications.

  One use case for SIROCCO where it is important to have as complete a line list as possible is its deployment as the radiation module in radiation-hydrodynamics simulations where line-driving forces are important (Higginbottom et al., 2020; Higginbottom et al., 2014). Our approach here is to use SIROCCO to estimate the ionization state and SED in each cell. This is then be passed to a stand-alone code that estimates the total line forces by summing over the complete Kurucz line list.

- **General Relativity**

SIROCCO self-consistently carries out the special relativistic frame transformations required as photon packets travel through the grid and interact with the moving material in each cell. However, it does not account for any purely *general* relativistic effects. In particular, photon packets in SIROCCO always travel in straight lines, rather than along geodesics. This will primarily affect the angular distributions of photon packets that are emitted or travel within, say, :math:`\simeq 10` gravitational radii of a compact object. This caveat should be kept in mind when modelling AGN and XRBs with SIROCCO. However, we console ourselves with the thought that the physical and radiative properties of accretion flows in this regime remain quite uncertain.

- **Polarization**

In principle, polarization can be included quite naturally in Monte-Carlo radiative transfer, but the current version of SIROCCO does not include a treatment of polarization.

- **Thermal and statistical equilibrium**

SIROCCO assumes that the flow is always and everywhere in thermal and statistical equilibrium. That is, the code iterates towards a temperature and ionization state for each cell in which the heating and cooling rates in each cell balance and the net transition rate *into* any given atomic/ionic level matches the net transition rate *out of* that level. This implies that there is no concept of time in SIROCCO -- the code is not designed to deal with non-equilibrium and/or time-dependent conditions.

This limitation can be important even if the input radiation field is steady. For example, if the flow velocity in a grid cell with characteristic size :math:`\Delta x` is given by :math:`v`, matter will flow through the cell on a time-scale :math:`t_{flow} \sim \Delta x / v`. However, ionization equilibrium can only be established on a time-scale of :math:`t_{rec} \sim \alpha N_e`, where :math:`\alpha` is the relevant recombination coefficient, and :math:`N_e` is the local electron density. Thus if :math:`t_{flow} < t_{rec}`, the cell cannot be in ionization equilibrium. In sufficiently fast-moving flows, the ionization state can then become "frozen-in", i.e. fixed to approximately the state at the last point where equilibrium could be established. Since SIROCCO currently has no concept of these time scales, it does not check for such non-equilibrium conditions. It is up to the user to carry out the relevant sanity checks on their models.

- **Dimensionality and resolution limits**

At present, SIROCCO is (at most) a 2.5-dimensional code. That is, the coordinate grid is restricted to 2D and assumed to be symmetric about the y-axis. However, photon packet transport takes place in 2D and allows for a rotational component of motion around the y-axis. In principle, upgrading SIROCCO to "full" 3D is fairly straightforward, but running such models would require significantly more computing resources. Memory requirements, in particular, scale directly with the number of grid cells used in a simulation. This is actually already a limitation for high-resolution models in 2D (or 2.5-D) and the reason we have set an upper limit of 500x500 cells as a default.

It is worth noting that SIROCCO's memory requirements for computationally demanding simulations are driven by the approach to parallelization that has been adopted. Currently, parallelization relies exclusively on ``MPI`` and requires the computational grid to be copied in full to each ``MPI`` process (i.e. each core). Memory requirements therefore increase rapidly as the number of processors is increased. This situation could be improved by adopting a hybrid ``OpenMP`` and ``MPI`` approach.
