Anisotropic Scattering
--------------------------

SIROCCO has a number of radiative transfer modes, controlled via the :doc:`/input/parameters/rt_ion/Line_transfer` keyword. Included in this mode is the treatment of line anisotropy; whether re-emission of a line photon is isotropic or not. When the scattering is isotropic, a new direction is simply chosen by choosing a random direction in the co-moving frame of the plasma. 

If anisotropic scattering is on, via one of the thermal trapping modes, the new direction is chosen according to a rejection method. The aim is to account for the fact the photon undergoes many interactions in the resonant zone due to the thermal width of the line, and finds it easier to escape along the direction in which the optical depth is lowest (where the velocity gradient is highest). Specifically, the code undergoes the following steps:

* choose a random number between 0 and 1, :math:`z`
* choose random direction, :math:`\hat{n}`
* evalute Sobolev optical depth along this direction, :math:`\tau_\hat{n}`
* calculate the escape probability along this direction :math:`P_{\rm esc} (\tau_\hat{n})`. 
* If :math:`P_{\rm esc} \geq z`, then escape the loop, otherwise increment a counter ``nnscat`` and proceed to the first step.

This process is repeated until the loop is exited or 10,000 iterations are reached. The rejection method is trying to sample the volume bounded in :math:`\theta,\phi` space by the complicated surface :math:`P_{\rm esc} (\theta,\phi)`. 

In highly optically thick regions, escape probabilities in all directions can be small, in which case the above rejection method can be extremely inefficient (the volume bounded by :math:`P_{\rm esc} (\theta,\phi)` is extremely small). Because of this, the code re-normalises the rejection method by first calculating the escape probability along the maximum velocity gradient, which is the maximum escape probability. 

.. admonition :: Developer note: the re-normalisation scheme

	Describe.

Anisotropy within the viewpoint technique
==================================================
Within the viewpoint technique (also called extract or the peel-off method), described under :doc:`../operation/spectral_cycles`, anisotropy has to be accounted for. At each interaction or wind photon generation, the photon packet is forced along a direction :math:`\theta`, with its weight  adjusted according to 

.. math::

    w_{\rm out}=\frac{P(\theta)}{\langle P (\theta) \rangle} w_{\rm in}.

For anisotropic scattering, :math:`P(\theta) \neq \langle P \rangle`. To deal with this, we need to calculate the escape probability along the desired direction, given by 

.. math::

    P(\theta) = \frac{1 - \exp [-\tau(\theta)]}{\tau(\theta)}

where :math:`\tau(\theta)` is the Sobolev optical depth in a given direction. This is a local quantity evaluated at the point of resonance. :math:`\langle P (\theta) \rangle` is calculated using a by-product of the rejection method. For a rejection method that samples a properly normalised probability space -- a probability space that has a (hyper)volume of 1 -- the number of iterations in the rejection method, :math:`N_{\rm it}` tells us (in this case) about the mean escape probability. More correctly, the expectation value of :math:`1/N_{\rm it}` is the mean escape probability. Thus, we multiply by a factor of :math:`1/N_{\rm it}` in the code to account for the :math:`\langle P (\theta) \rangle` factor in the denominator.

.. todo:: check the above statement about the expectation value of :math:`1/N_{\rm it}` is really true -- I think it must be, since it's basically the definition of a probability. Does :math:`N_{\rm it}` also correspond to the actual physical number of scatters? 

.. admonition :: Developer note

    The above calculation is split up within the code. The factor :math:`P(\theta)` is applied in the function ``extract_one``, whereas the division by :math:`\langle P \rangle` is applied using the variable ``nnscat`` in extract, which is :math:`N_{\rm it}` in the above notation. This is because the mean escape probability is (statistically speaking) equal to :math:`1/N_{\rm it}` as described above.

    Note, also, that in practice we have to account for the renormalisation of the rejection method, so rather than multiply by :math:`N_{\rm it}`, we multiply by :math:`N_{\rm it}/P_{\rm norm}` (see pevious developer note).


