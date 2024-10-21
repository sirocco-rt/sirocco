Wind.type
=========

Sets the type of wind model the user wants to use. More information of model specifics can be found in :ref:`Wind Model Parameters` and :ref:`Wind Models`.

Type
  Enumerator

Values
  SV
    The Shlosman and Vitello parameterization of a bi-conical flow. (:ref:`SV`)

  corona
    A simple model for a corona above the disk. (:ref:`Corona`)

  homologous
    A homologous expansion law useful for simulating SNe. (:ref:`Homologous`)

  hydro
    A special purpose mode used by the python collaboration for importing models from Zeus and Pluto. (Depreciated) (:ref:`Hydro`)

  imported
    A general purpose mode for importing a wind from an ascii file (see also :doc:`SIROCCO Script documentation <../wind_models/importing_models>`).

  kwd
    The Knigge-Woods-Drew parameterization of a bi-conical flow. (:ref:`KWD`)

  shell
    A model of a thin shell useful for diagnostic studies. (:ref:`Shell`)

  star
    A stellar-wind model. A fairly standard parameterization of a spherical outflow for a hot star. (:ref:`Stellar_wind`)


File
  `setup_domains.c <https://github.com/agnwinds/python/blob/master/source/setup_domains.c>`_


Parent(s)
  * :ref:`Wind.number_of_components`: Greater than 0. Once per domain.


Child(ren)
  * :ref:`Shell.wind_v_at_rmin`

  * :ref:`Corona.radmax`

  * :ref:`Wind.mdot`

  * :ref:`KWD.mdot_r_exponent`

  * :ref:`Corona.base_den`

  * :ref:`KWD.v_zero`

  * :ref:`Stellar_wind.mdot`

  * :ref:`Homologous.radmin`

  * :ref:`KWD.acceleration_length`

  * :ref:`Corona.radmin`

  * :ref:`KWD.rmax`

  * :ref:`Homologous.radmax`

  * :ref:`SV.thetamax`

  * :ref:`SV.acceleration_exponent`

  * :ref:`Corona.zmax`

  * :ref:`Corona.scale_height`

  * :ref:`Homologous.density_exponent`

  * :ref:`Hydro.thetamax`

  * :ref:`Wind.dim.in.z_or_theta.direction`

  * :ref:`SV.diskmin`

  * :ref:`SV.diskmax`

  * :ref:`SV.acceleration_length`

  * :ref:`Hydro.file`

  * :ref:`KWD.acceleration_exponent`

  * :ref:`Corona.vel_frac`

  * :ref:`Stellar_wind.radmin`

  * :ref:`Shell.wind.radmax`

  * :ref:`Stellar_wind.radmax`

  * :ref:`Wind.model2import`

  * :ref:`Homologous.vbase`

  * :ref:`Homologous.boundary_mdot`

  * :ref:`KWD.rmin`

  * :ref:`Shell.wind_mdot`

  * :ref:`SV.mdot_r_exponent`

  * :ref:`KWD.d`

  * :ref:`Shell.wind.v_at_rmax`

  * :ref:`Stellar_wind.acceleration_exponent`

  * :ref:`Stellar_wind.v_infinity`

  * :ref:`Shell.wind.radmin`

  * :ref:`Shell.wind.acceleration_exponent`

  * :ref:`SV.v_zero_mode`

  * :ref:`SV.v_infinity`

  * :ref:`Stellar_wind.vbase`

  * :ref:`Wind.dim.in.x_or_r.direction`

  * :ref:`KWD.v_infinity`

  * :ref:`SV.thetamin`

