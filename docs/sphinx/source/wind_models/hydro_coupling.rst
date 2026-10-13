.. _hydro_coupling:

Coupling SIROCCO to a Hydrodynamics Code
########################################

SIROCCO can be coupled to a hydrodynamics code, such as PLUTO, to carry out
radiation-hydrodynamics calculations (see the `PLUTO-Sirocco
<https://github.com/sirocco-rt/pluto-sirocco>`_ framework). The hydrodynamics
code advances the density, velocity and temperature of the flow; SIROCCO
calculates the ionization state and the radiative heating, cooling and forces
for the current state of the flow, and these are passed back to the
hydrodynamics code.

This page defines the interface between the two codes: the sequence of runs,
the quantities that are exchanged, and the conventions they must follow.

The coupled sequence
====================

Each exchange between the two codes follows the same sequence:

1. The first snapshot of the flow is run with SIROCCO as an imported model
   (``Wind.type imported``, see :doc:`importing_models`)::

       sirocco -f snapshot1.pf

2. The hydrodynamics code advances the flow and writes the next snapshot in the
   same imported-model format.

3. ``hydro2sirocco`` paints the new snapshot onto the wind save from the
   previous SIROCCO run, keeping the ionization state as a starting point::

       hydro2sirocco -model_file snapshot2.txt snapshot1

   This reads ``snapshot1.wind_save`` and writes ``new.wind_save``.

4. SIROCCO continues from the updated wind save, carrying out more ionization
   cycles. For a restart, ``Ionization_cycles`` in the parameter file is the
   *total* number of cycles, so it must be larger than the number already
   completed::

       cp new.wind_save input.wind_save
       sirocco -f -r input.pf

5. ``rad_hydro_files`` writes the heating, cooling, forces and ionization state
   for the hydrodynamics code (see :ref:`rad_hydro_files`)::

       rad_hydro_files input

The ``-f`` switch holds the electron temperature fixed at the value supplied by
the hydrodynamics code; SIROCCO then finds the ionization balance at that
temperature rather than solving for thermal equilibrium.

.. note::

   The ``hydro`` wind type and the ``-z`` switch provide an older route for
   coupling to Zeus, in which SIROCCO reads the hydrodynamics model itself on
   restart. It is deprecated; the sequence above should be used.

Quantities passed to SIROCCO
============================

Each snapshot is a model file in the format described in
:doc:`importing_models`. For the coupling, it must satisfy the following.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Quantity
     - Convention
   * - Grid
     - Identical in every snapshot: the same coordinate system, dimensions and
       cell positions as the first snapshot. The grid cannot change during a
       coupled calculation.
   * - ``inwind``
     - Identical in every snapshot. A wind save has a fixed number of plasma
       cells, so cells cannot enter or leave the wind.
   * - Velocity
     - Cartesian components :math:`v_x, v_y, v_z` in cm/s, in the observer
       frame, at the grid points.
   * - Density
     - The **volume-averaged** mass density in g/cm\ :sup:`3`, in the observer
       frame. This is the density a hydrodynamics code normally works with; it
       is not the density inside clumps.
   * - Temperature
     - Optional. One column is taken as :math:`T_e`, with
       :math:`T_r = 1.1\,T_e`; two columns are :math:`T_e` and :math:`T_r`.
       If there are none, ``hydro2sirocco`` keeps the temperatures in the wind
       save.

SIROCCO stores the density inside clumps. For a filling factor :math:`f`
(``Wind.filling_factor``) and Lorentz factor :math:`\gamma` at the cell centre,
the density SIROCCO uses is

.. math::

   \rho_{\rm clump} = \frac{\rho_{\rm model}}{f\,\gamma},

and the volume filled with material is :math:`f\,V_{\rm cell}`. The same
conversion is made when a model is first imported and when a wind save is
updated by ``hydro2sirocco``. With the ``-nonrel`` or ``-sr_doppler_only``
switches, :math:`\gamma = 1`.

.. important::

   The relativity switches (``-nonrel``, ``-sr_doppler_only``) change
   :math:`\gamma`, so the same switches must be given to ``sirocco`` and
   ``hydro2sirocco`` in a coupled calculation.

What hydro2sirocco updates
==========================

``hydro2sirocco`` makes the wind save what SIROCCO would have produced from the
new snapshot, except that it keeps the ionization state of the previous run as
the starting point for the next ionization cycles. It updates the following
quantities in every cell of the imported domain.

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Quantity
     - How it is set
     - Why
   * - Velocity at the grid points
     - From the snapshot
     - The flow has moved on.
   * - Velocity gradient and divergence
     - Recalculated from the new velocities
     - They depend on the velocity, and set the Sobolev optical depths.
   * - Average and maximum :math:`dv/ds`
     - Recalculated, after all the gradients have been updated
     - They depend on the velocity gradients, and are used in line transfer.
   * - Lorentz factors, :math:`\gamma`
     - Recalculated from the new velocities
     - They enter the density and volume.
   * - Cell volume
     - Rescaled by the change in :math:`\gamma` at the cell centre
     - The volume includes :math:`\gamma`, as it does when the wind is created.
   * - Density, :math:`\rho`
     - :math:`\rho_{\rm model}/(f\,\gamma)`
     - SIROCCO stores the density inside clumps.
   * - Ion densities
     - Multiplied by :math:`\rho_{\rm new}/\rho_{\rm old}`
     - Keeps the ion fractions from the previous run, as a starting guess.
   * - Electron density, :math:`n_e`
     - Recalculated from the ion densities
     - Otherwise it would belong to the old density.
   * - :math:`T_e` and :math:`T_r` (if in the snapshot)
     - From the snapshot
     - The hydrodynamics code sets the temperature.
   * - Partition functions and level populations (if the temperature changed)
     - Set to LTE ground-state values at the new temperature
     - Otherwise they would belong to the old temperature. This follows the
       older ``-z`` route; the level populations are recalculated in the next
       ionization cycle.

Everything else in the wind save is kept, in particular the ion fractions, the
radiation field estimators and the number of cycles already completed.

``hydro2sirocco`` stops without writing a wind save if the wind save does not
have exactly one imported domain, if the snapshot has a different grid or
different ``inwind`` values, or if a cell has zero or negative density.

.. note::

   The average and maximum :math:`dv/ds` are estimated using random
   directions, so they differ slightly from one run to another, as they do
   when SIROCCO creates a wind.

.. admonition:: History

   Until September 2026 this step was carried out with ``modify_wind
   -model_file``. That did not divide the density by the filling factor,
   recalculate :math:`n_e`, or recalculate the velocity gradients, so clumped
   models were run with a density that was too low by a factor :math:`1/f`
   (issue #1190). ``modify_wind`` no longer accepts ``-model_file``.

Quantities passed back to the hydrodynamics code
================================================

``rad_hydro_files`` writes the heating and cooling rates, the radiative
driving forces and the ionization state of each cell; the files are described
in :ref:`rad_hydro_files`. Which of these quantities should be volume averaged
and which should refer to the gas inside clumps, when :math:`f<1`, is being
reviewed.
