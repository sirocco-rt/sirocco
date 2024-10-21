Spectral Cycles
###############

The purpose of the ionization cycles is to establish the ionization state
of the plasma. The purpose of the spectral cycles is to created simulated
spectra given a defined ionization state in a wavelength range that is
(usually) less wide than required to establish the ionization state.  For
a cataclysminc variable, One might, for example, want to create a simulated
spectrum to compare with an observed spectrum of that object in the ultraviolet.
which is observed with a specific inclination angle with respect to the disk.

For the simple atom case, the process is relatively straightforward.  One
begins by generating photon packets that cover a range that slightly boarder
than the spectral range of interest, slightly broader because one needs to
allow for Doppler effects associated with scatters that can occur.  One then
follows these photons through the wind, where the number of photons carred
by the packet diminishes as it moves through the wind.

Then one could simply create a spectrum from photon packets that exit the
simulation volume at a particular incination angle (plus or minus some
delta) to construct the spectrum.  This so-called live-or-die  option
is implemented in SIROCCO, but only as a diagnostic option, because it is
inefficient since most photon packets exit the system at inclination
angles that are not of interest.

Instead, the standard way of constructing detailed spectra is to use the method,
termed the viewpoint technique, as described by `Knigge, Woods, & Drew 1995
<https://ui.adsabs.harvard.edu/abs/1995MNRAS.273..225K/abstract>`_, also
known as the peel-off method (`Yusef-Zedeh, Morris & White 1984 <https://ui.adsabs.harvard.edu/abs/1984ApJ...278..186Y/abstract>`_.
In this method, one follows a photon through the grid as before, but at point
where the photon changes direction (including the inital point of photon generation),
one creates a dummy photon headed in the desired direction.  One adjust the
weight of the dummy photon accoring to
the relative probability that a photon  will escape in the desired
direction compared to the angle averaged probability, and adjusts the number
of photons by that fraction, that is

.. math::

    w_{\rm out}=\frac{P(\theta)}{\langle P \rangle} w_{\rm in}.

For isotropic scattering the :math:`w_{\rm out}==w_{\rm in}` but for resonant scattering the
weight will increase if the desired photon direction is in the direction of maximum
velocity gradient and decrease if it is along the direction of minumum velocity gradient (see :doc:`/physics/aniso`).
For photons generated at the surface of a star of disk, the weight of the dummy photon
is determined by the limb darkening law assumed. One then extracts the dummy photon along
the line of sight reducing the weight of the photon by the total optical depth along that
lien of sight.  Evidently, one can repeat this process at every interaction when one
wishes to construct a spectrum along multiple lines of sight.

Detailed Spectral Calculation when Macro-atoms are used
-------------------------------------------------------

When a macro-atom is excited, photon packets can emerge at very different frequencies than
the frequency of the photon packet before an interaction.  This requires a modification of
the methods used during ionization cycles, where, in the macro-atom case, no photons or r-packets
originate in the wind and a strict radiative equilibrium constraint is enforced
(with a few exceptions, e.g. adiabatic cooling).

During the ionization cycles, the amounts of energy flowing into each macro-atom level,
and into the thermal k-packet pool, are recorded in the matom_abs and kpkt_abs quantities.
In the spectral cycles, one needs to know where this energy comes out - if energy flows into
a given state, what proportion of that energy comes out via the various possible transitions?
This issue is dealt with in the "macro-atom emissivity calculation", which is carried out
at the start of the spectral cycles. The current procedure is to do a Monte Carlo sampling of
the macro-atom machinery -- a large number of packets are generated with initial macro-atom
states in proportion to the estimators matom_abs and kpkt_abs. The fraction of times these packets
de-activate from given states is then recorded, and the corresponding r-packet frequency is
calculated. If the frequency falls within the requested range, the relevant macro-atom or k-packet
emissivity is incremented by the appropriate fraction of matom_abs. If the frequency falls outside
the range, the contribution is ignored. This procedure can be speeded up by using an implicit/matrix
scheme where the matrix contains the mapping between the absorbed and emergent radiation; This
method is currently in the development stage in the code.

In the actual photon transport stage, r-packets are generated in the wind in proportion with
these frequency-limited emissivities, in a broadly similar to wind photon generation in the non-macro atoms scheme.
In the process, we also ensure that the photons are only generated over the correct frequency range.
The photon transport is then carried out as normal, except that whenever a macro-atom is activated, or a k-packet is created,
that photon / energy packet is immediately thrown away to avoid double-counting -- the emission resulting from the interaction
has already been (statistically speaking) accounted for in the emissivity calculation. The main difference to the approach in the
ionization cycles is that now the radiative equilibrium condition is enforced by the fact that the
emissivities should be consistent with the "absorbed" radiation, instead of being explicitly enforced by never
destroying packets.

When using indivisible packet line transfer (commonly referred to as macro-atom mode), the code is
often used in a hybrid mode where some elements are treated as macro-atoms and some as simple-atoms.
Simple atoms are treated differently to macro-atoms; there is no detailed model atom and internal
transitions to other states are not possible. Instead, for both bound-free and bound-bound interactions,
a fake two-level atom is excited, and the excited state either radiatively or collisionally decays.
If it collisionally decays, this would normally excite a k-packet so the packet is destroyed for the
same reasons as above (the emissivity from this is already accounted for). The k-packets generated from
the emissivity are allowed to create simple-atom emission for consistency. If it radiatively decays, we treat
the interaction like a resonant scatter and proceed, still tracking the packet.
This is a reasonable approximation for resonant lines, but less so for bound-free continua, since the possible
recombination cascade and potential reddening of the photons is not dealt with.
Partially addressing this was the aim of the bound-free simple emissivity approach.

Other Issues
------------
Spectral cycles are executed after the ionization and temperature state of the wind is computed 
by the ionization cycles. It is possible to modify the requested spectra (both wavelength
range and observer angles) by making the changes to the spectral cycle parameters, setting the
number if ionization cycles to zero, and then restarting the simulation.
