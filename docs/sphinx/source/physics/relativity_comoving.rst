Special Relativity and Co-Moving Frames
---------------------------------------

The current version of SIROCCO incorporates special relativity and takes co-moving frame
effects into account by default.  

Global properties of the wind, such as a densities are defined in the observer
, or global frame, but are immediately converted to co-moving frame values.  

(As an example, if the density of cell and volume (or a cell) in the global frame are
:math:`\rho_{obs}` and :math:`V_{obs}` then 


.. math::

    \rho_{cmf} = \frac{\rho_{obs}}{\gamma}

.. math::

    V_{cmf}=\gamma V_{obs}

the product of the two quantities, being a Lorentz invariant.)

Photon generation takes place in the local, or co-moving, frame (of the disk or wind), 
but photons are immediately converted to the observer, or global, frame for 
photon transport, allowing both for Doppler frequency shifts and directional 
correction due to Doppler abberation.  The number of photons generated is the number of
photons that would be generated in the in one observer frame second. Photons are 
transported in the observer frame, but coverted back to the local frame within i
ndividual wind cells to determine whether the interact with the wind.

Interactions take place in the local frame. Estimators used to calculate, for example, 
photoionization rates also take place in the local frame.  This allows one to calculate
ionization fractions in the local frame as is required, since the numbers of ions in
a region defined by the edges of the cell must also be Lorentz invariant. Allowances are 
made for time dilation effects in calculating the rates in the co-moving frame.


For mainly historical and diagnostic reasons, command line options exist to fall back 
to simple order v/c corrections.  `
