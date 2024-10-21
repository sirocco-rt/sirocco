Radiative Transfer Modes
########################################################

SIROCCO has a number of different radiative transfer modes, which affect the treatment of lines and scattering, and also whether we use indivisible packet constraints or allow photon weights to be attenuated by continuum absorption. These modes are selected with the parameter :ref:`Line_transfer`. The different modes are briefly described on that parameter page. This page is designed to give an overview of the assumptions and concepts behind, as well as the basic operation of, the different techniques. The aim is that, in partnership with the parameter page and the atomic data documentation, all the information regarding the different radiative transfer modes should be present.

For introductions and references regarding Monte Carlo radiative transfer techniques generally, we recommend reading `Noebauer & Sim 2019 <https://ui.adsabs.harvard.edu/abs/2019LRCA....5....1N/abstract>`_. For specifics regarding SIROCCO, we recommend reading `Long & Knigge 2002 <https://ui.adsabs.harvard.edu/abs/2002ApJ...579..725L/abstract>`_ as well as  PhD theses by `Higginbottom <https://eprints.soton.ac.uk/368584/1/Higginbottom.pdf>`_ and `Matthews <https://ui.adsabs.harvard.edu/abs/2016PhDT.......348M/abstract>`_. 

Sobolev Approximation
======================
SIROCCO always uses the Sobolev approximation to treat line transfer. In this approximation, it is assumed that the thermal line width is small compared to the velocity gradient. The Sobolev approximation is described extensively in astrophysics literature, and so we do not describe it in detail here. We refer the users to section 8.2 of `Noebauer & Sim 2019 <https://ui.adsabs.harvard.edu/abs/2019LRCA....5....1N/abstract>`_ and references there in for a discussion of the Sobolev escape probabilities approach.

Weight Reduction v Indivisible Packets 
=======================================
SIROCCO was originally written in such a way that photon packet weights were not indivisible and allowed to be attenuated. This is the way the code is described in the original `Long & Knigge 2002 <https://ui.adsabs.harvard.edu/abs/2002ApJ...579..725L/abstract>`_ paper. In the standard, weight reduction mode, photon weights are attenuated by continuum opacities (free-free, bound-free). Conservation of energy is then hopefully achieved by calculating the emission from the wind .

In indivisible packet mode, there is a fundamental shift in philosophy. All energy packets are strictly indivisible and conserve energy whenever they undergo radiative processes (the only exception is adiabatic cooling). Thus, even bound-free absorption is dealt with at a single interaction point.

Indivisible packet mode is activated by setting the :ref:`Line_transfer` parameter to either ``macro_atoms`` or ``macro_atoms_thermal_trapping``. The terminology adopted here is slightly confusing, since the line transfer mode does not explicitly include a macro-atom treatment of atomic species (see next subsection).

.. admonition:: Developer Note
  
  The radiative transfer mode is stored using the code variable geo.rt_mode.

Macro-atoms and 2-level-atoms 
==============================
The macro-atom scheme was devised by Leon Lucy in the early 2000s (`Lucy 2002 <https://ui.adsabs.harvard.edu/abs/2002A%26A...384..725L/abstract>`_, `Lucy 2003 <https://ui.adsabs.harvard.edu/abs/2003A%26A...403..261L/abstract>`_). 
It involves the reformulation of the process of radiation transport through a plasma in radiative equilibrium into a traffic-flow problem. Lucy showed that, when in radiative equilibrium, the energy flows through a system depend only on the transition probabilities and atomic physics associated with the levels the energy flow interacts with. By quantising this energy flow into radiant (r-) and kinetic (k-) packets, we can simulate the energy transport through a plasma discretised into volume elements (*macro-atoms*), whose associated transition probabilities govern the interaction of radiant and kinetic energy with the ionization and excitation energy associated with the ions of the plasma.

.. todo:: add refs, describe properly.

.. admonition:: Developer Note

  Macro-atoms are identified using their atomic data, in particular by providing data with identifiers
  LevMacro, LinMacro, PhotMacro. 

Simple-atoms still interact with r- and k-packets, but do not possess internal transition probabilities. As a result, they are analogous to the two-level atom treatment, as any excitation is immediately followed by a deactivation into an r- or k-packet. The choice of radiative or kinetic deactivation is made according  to the relative rates in the two-level atom formalism. 

Isotropic v Anisotropic Line Scattering 
============================================
SIROCCO always treats electron scattering as an isotropic process, and continuum emission processes are also treated as isotropic, except for Compton scattering. For Compton scattering, the direction and energy change is calculated self-consistently according to the energy change formula :math:`E/E'=1+(h \nu/mc^2)(1+\cos\theta)`. We first draw a random cross section that our photon packet will see. This cross section represents an energy change and hence a direction. The distribution of angles is taken care of by using a differential cross section vs energy change function. 

.. admonition:: Caution

  Compton scattering is currently not accounted for when using indivisible packet mode. 

Line emission and scattering is isotropic unless one of the  ``thermal_trapping`` line transfer modes is selected. In the thermal trapping mode, any line interaction or emission results in an anisotropic direction being generated. This direction is generated by a rejection method which samples the Sobolev escape probability in each direction from the line interaction region. Unless you specifically want to consider isotropic line emission, we recommend always using the anisotropic thermal trapping mode. 

.. todo:: move the below to where we describe photon sources and generation?

In the case of isotropic emission, the direction of a photon packet is chosen so that the probability of emission in each bin of solid angle is the same. It follows that 

.. math::
    p(\Omega)d\Omega \propto \cos \theta \sin \theta d\theta d\phi

where the angles are in polar coordinates and relative to the local outward normal. For a spherical emitting source, such as a star, one must first generate a location on the star's surface and then calculate the photon direction relative to the normal at the point. For emission from optically thick surfaces the above equation can be modified to include linear limb darkening, :math:`\eta(\theta)`, such that

.. math::
    p(\theta, \phi) d\theta d\phi = \eta(\theta) \cos \theta \sin \theta d\theta d\phi.

The Eddington approximation is usually adopted in the code, so that $\eta(\theta)$
is given by

.. math::
    \eta(\theta) = a (1 - \frac{3}{2} \cos \theta).

The constant :math:`a` is normalised such that the total probability sums to 1. Whenever a radiation packet undergoes an electron scatter, the new direction is chosen to be isotropic. However, when the photon is a line photon, the new direction is chosen according to a line trapping model, which samples a probability distribution according to the Sobolev escape probability in different directions. 

Doppler Shifts and The Comoving Frame  
============================================
When calculating opacities, the photon frequency must be shifted from the rest frame of the photon into the rest frame of the plasma. This shift depends on the before and after directions of the photon. Let us denote these two directions with unit vectors :math:`\vec{n}_i` and :math:`\vec{n}_f`, respectively, and consider a situation when a photon scatters off an electron in a region of the wind moving at velocity :math:`\vec{v}`. The final frequency of the photon with initial frequency is 

.. math::
    \nu_f = \nu_i ~\frac{1 - (\vec{v} \cdot \vec{n}_i) / c}{1 - (\vec{v} \cdot \vec{n}_f) / c}.

In the case of a resonance scatter with line transition u to j, the new frequency is

.. math::
    \nu_f = \frac{\nu_{uj}}{1 - (\vec{v} \cdot \vec{n}_f) / c}.

The above formulae are the non-relativistic case, which is currently used in the code. However, this should in general be improved to use the special relativistic formula. This would produce more accurate Doppler shifts for the fastest regions of an outflow, as the current treatment introduces errors of order 5 Angstroms at the blue edges of the highest velocity absorption lines in quasar and CV wind models.

When real photons resonantly (or electron) scatter off real plasma in a flow, they conserve energy and frequency in the co-moving frame of the plasma. In the case of an outflow, doing the frame transformation from system->flow->system over the course of an interaction results in a redshifting of a photon, and as a result an energy loss - in other words, the photon does work on the flow even though energy is conserved in the co-moving frame. Indivisible packet schemes (such as macro-atoms) often enforce strict energy conservation in the frame of a given cell (physically, but see also `Lucy 2002 <https://ui.adsabs.harvard.edu/abs/2002A%26A...384..725L/abstract>`_). This means that, when keeping track of packets in the observer frame, one needs to correct the energies (not just the frequencies) using a Doppler shift. SIROCCO does **not** currently conserve energy in the co-moving frame.

.. todo:: test whether this is an issue.
