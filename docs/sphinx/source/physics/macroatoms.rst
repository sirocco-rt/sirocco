Macro Atoms
-------------------------

.. todo:: This page shoud contain a description of how macro atoms work. The below is copied from JM's thesis.

.. todo:: Add description of accelerated macro-atom scheme

The macro-atom scheme was created by Leon Lucy and is outlined in his 2002/03 papers. It was implemented in SIROCCO by Stuart Sim, initially for the study of recombination lines in YSOs (Sim et al. 2005).

Lucy (2002,2004) hereafter L02, L03 has shown that it is possible to calculate the emissivity of a gas in statistical equilibrium without approximation for problems with large departures from LTE. His `macro-atom` scheme allows for all possible transition paths from a given level, dispensing with the two-level approximation, and provides a full non-LTE solution for the level populations based on Monte Carlo estimators. The macro-atom technique has already been used to model Wolf-Rayet star winds (Sim 2004), AGN disc winds (Sim et al. 2008), supernovae (Kromer and Sim 2009, Kerzendorf and Sim 2014) and YSOs (Sim et al. 2005). A full description of the approach can be found in L02 and L03.

Following L02, let us consider an atomic species interacting with a radiation field. If the quantity :math:`\epsilon_j` represents the ionization plus excitation energy of a level :math:`j` then the rates at which the level absorbs and emits radiant energy are given by

.. math::

    \dot{A}_{j}^{R} = R_{l j} \epsilon_{j l}

and

.. math::

    \dot{E}_{j}^{R} = R_{j l} \epsilon_{j l}


where :math:`\epsilon_{j l} = \epsilon_j - \epsilon_l`.
Here, we have adopted Lucy's convention, in which the subscript
:math:`l` denotes a summation over all lower states (:math:`<j`), and
(:math:`u`) a summation over all upper states (:math:`>j`).
Similarly, the rates corresponding to _kinetic_ (collisional)
energy transport can then be written as

.. math::

    \dot{A}_{j}^{C} = C_{l j} \epsilon_{j l}

and


.. math::

    \dot{E}_{j}^{C} = C_{j l} \epsilon_{j},

Let us define :math:`{\cal R}` as a total rate, such that
:math:`{\cal R}_{l j}  = R_{l j} + C_{l j}`.
If we now impose statistical equilibrium

.. math::

    ({\cal R}_{l j}-{\cal R}_{j l})+({\cal R}_{u j}-{\cal R}_{ju})=0 \;\;\;,

we obtain

.. math::

    \dot{E}_{j}^{R}+\dot{E}_{j}^{C}+{\cal R}_{ju}\epsilon_{j}+ {\cal R}_{j \ell}\epsilon_{l}  \nonumber \\  = \dot{A}_{j}^{R}+\dot{A}_{j}^{C}+{\cal R}_{u j} \epsilon_{j} +{\cal R}_{l j} \epsilon_{l}.

This equation is the starting point for the macro-atom scheme. It shows that, when assuming radiative equilibrium, the energy flows through a system depend only on the transition probabilities and atomic physics associated with the levels the energy flow interacts with. By quantising this energy flow into radiant (:math:`r-`) and kinetic (:math:`k-`) packets, we can simulate the energy transport through a plasma discretised into volume elements (macro-atoms),whose associated transition probabilities govern the interaction of radiant and kinetic energy with the ionization and excitation energy associated with the ions of the plasma.

Although the equation above assumes strict radiative equilbrium,it is trivial to adjust it to include non-radiative source and sink terms. For example, in an expanding parcel of plasma, adiabatic cooling may be included with a simple modification to the RHS. Currently, we include adiabatic cooling by destroying packets with a probability
:math:`p_{i,\mathrm{destruct}} = {\cal C}_{\rm adiabatic} / {\cal C}_{\rm tot}`.


A Hybrid Scheme
=============================

A pure macro-atom approach with only H and He can be easily used for some situations -- for example, in the YSO application described by, which uses a H-only model. However, in accretion disc winds, the densities can be very high, and higher :math:`Z` elements must be  included. Including all these elements as macro-atoms is not currently computationally feasible in the code for anything but the simplest models. We thus often use a `hybrid scheme`, which treats H and He with the macro-atom approach, but models all other atoms as `simple-atoms`.

Simple-atoms still interact with :math:`r`- and :math:`k`-packets but do not possess internal transition probabilities. As a result, they are analogous to the two-level atom treatment, as any excitation is immediately followed by a deactivation into an :math:`r`- or :math:`k`-packet. The choice of radiative or kinetic deactivation is made according to the relative rates in the two-level atom formalism. For a bound-bound transition :math:`u\to j`, these two probabilities

are then

.. math::
    p_{uj}^{S,R} = \frac{ A_{uj} \beta_{uj} } { A_{uj} \beta_{uj} + C_{uj} \exp(-h\nu_{uj} / k T_e) } = 1 - q

and

.. math::
    p_{uj}^{S,C} = \frac{ C_{uj} \exp(-h\nu_{uj} / k T_e) } { A_{uj} \beta_{uj} + C_{uj} \exp(-h\nu_{uj} / k T_e) } = q.


For a bound-free transition, the code assumes radiative recombination, and thus any bound-free simple-atom activation is immediately followed by the creation of an :math:`r`-packet. This approximates the bound-free continuunm, even when compared to other two-level atom radiative transfer schemes. This is discussed further and tested in section~\ref{sec:line_test}.

This hybrid approach preserves the fast treatment of, for example, UV resonance lines, while accurately modelling the recombination cascades that populate the levels responsible for, e.g., H and He line emission. As a result of this hybrid scheme, a separate set of estimators must be recorded for simple-atoms,  and the ionization and excitation of these elements is calculated with a different, approximate approach. In order to include simple-atoms, we must add in a few extra pathways, so that energy packets can also activate simple-atoms, through either bound-free or bound-bound processes. The relative probabilities of these channels are set in proportion with the simple-atom opacities.

Macro-atom Emissivity Calculation
========================================

In order to preserve the philosophy that a detailed spectrum is calculated in a limited wavelength regime, SIROCCO carries out a macro-atom emissivity calculation before the spectral cycles. The aim of this step is to calculate the luminosity contributed by macro-atoms -- equivalent to the total amount of reprocessed emission -- in the wavelength range being considered.

This process can be very computationally intensive, especially if the wavelength regime being simulated has very little emission from bound-free and line processes in the wind, but the overall broad-band emissivity is high. During the ionization cycles, the amount of energy absorbed into :math:`k`-packets and every macro-atom level is recorded using MC estimators. Once  the ionization cycles are finished, and the model has converged, these absorption energies are split into a certain number of packets and tracked through the macro-atom machinery until a deactivation occurs. When this happens, the emissivity of the level the macro-atom de-activated from is incremented if the packet lies in the requested wavelength range. If it does not, then  the packet is thrown away. It is easy to see how what is essentially a MC rejection method can be an inefficient way of sampling this parameter space. Fortunately, this problem is parallelised in the code.

Once the emissivities have been calculated, the spectral synthesis can proceed. This is done in a different way to the ionization cycles. Photons are generated from the specified photon sources over the required wavelength range, but are now also generated according to the calculated macro-atom and :math:`k`-packet emissivities in each cell. These photons are "extracted" as with normal photon packets. In order to ensure that radiative equilibrium still holds, any photon that interacts with a macro-atom or :math:`k`-packet is immediately destroyed. The photons are tracked and extracted as normal until they escape the simulation; resonant scatters are dealt with by a combination of macro-atom photon production and destruction.

.. admonition :: Developer note: Emissivities

    We are a little lax in terms of what we actually call an emissivity in the code. The quantities stored in variables like ``kpkt_emiss`` and ``matom_emiss`` in the plasma and macro-atom structures are actually *comoving-frame energies* in erg, which are sampled when generating :math:`r`-packets in each cell. Roughly speaking, these are luminosities given that the code assumes a time unit of 1s. Similarly, when the code prints out level *emissivities* to screen and to the diag file, these are really a sum over all these quantities (and can approximately be thought of as level *luminosities*).

Bound-free Continua of Simple Atoms
=============================================
.. todo:: this section is not yet completely accurate.

Historically, when using the indivisible packet form of radiative transfer (`macro_atoms_thermal_trapping`, for example), the bound-free continua of simple atoms were treated in a simplified two-level framework. In this case, simple atoms are those `without` a full macro-atom model atom, usually the metals. In this two-level scheme, whenever a simple atom undergoes a bound-free interaction, it is excited into the continuum state, and this is immediately followed by recombination, and an :math:`r`-packet or :math:`k`-packet is created immediately. As a result, the scheme does not capture the physical situation whereby a recombination cascade can occur from an initial recombination to excited levels, leading to a gradual reddening of the photon if there are many interactions. This situation **is** modelled well by a full macro-atom treatment. As of 2024, this is once again the default behaviour.

To try and slightly improve this scheme, we implemented a "total emissivity" upweighting scheme around 2018. The basic idea is that we pay attention to only the heating and cooling. In particular, the rates of all simple atom bound-free emission are governed by the `emissivity` of the bound-free process.
**Currently, this mode is turned off by default**, due to `various issues <https://github.com/sirocco-rt/sirocco/issues?q=is%3Aissue+label%3Aupweighting%3F+>`_ associated with energy conservation, as also described below.

This result in two changes to the code for ionization cycles:
   * whenever a k-packet is eliminated via a bound-free channel of a simple macro atom (simulating energy flow from the :math:`k`-packet pool to the radiation pool, :math:`k \to r`), we have that packet carry additional energy corresponding to the required ionization energy for that particular bf process. This means we upweight the energy of the packet by a factor :math:`f_{\rm up} = \nu / (\nu - \nu_0)`, where :math:`\nu` is the frequency of the new bound-free photon and :math:`\nu_0` is the threshold frequency. This quantity is the ratio of the total energy carried by photons in the packet to the energy supplied to photons in the packet from the thermal pool.
   * whenever an r-packet is “absorbed” by a simple macro atom bound-free process we track explicitly only the flow of energy to the thermal pool. This means we force the creation of a :math:`k`-packet, whereas before there woud be a choice, but we only take the contribution of the absorption to heating only: i.e. we downweight the packet energy by a factor :math:`f_{\rm down} = (\nu - \nu_0) / \nu`.

In the spectral cycles, interactions with simple bound-free continua now kill the photon, and :math:`k \to r` follow the same behaviour as above, because in these cycles we introduce a precalculated band-limited :math:`k`-packet emissivity.

**It is possible for some numerical problems to occur.** For example, there is nothing to stop the value of :math:`f_{\rm up}` being quite large, if the photon is being emitted close to the edge. This is most likely to happen when the electron temperature :math:`T_e` is quite low, but there is nothing to stop it happening anywhere. This is most likely to lead to problems when the factor :math:`f_{\rm up}` is comparable to the typical number of photon passages per cell, since then a single photon can dominate the heating or ionization estimators in a given cell and lead to convergence problems by dramatically exacerbating shot noise.

.. admonition :: Activating the scheme

    This mode can be turned on using the :ref:`Diag.use_upweighting_of_simple_macro_atoms`.
    In this case the code will go back to using the two-level framework for simple atom bound free continua.
