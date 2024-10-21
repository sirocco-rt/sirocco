Photon Banding Strategies
#########################

Photon packets are emitted from a number of different radiation sources, such as
the accretion disk or from the wind itself. When a photon is created, it
is defined by its frequency :math:`\nu` and weight :math:`w`. Photons are 
generated at the beginning of each cycle and can either be generated uniformly
over the entire frequency range, or can be generated in pre-defined frequency
bands, where certain frequency bands are biased to have more photons.

Uniform Sampling
================

In the most simple case, the frequency of a photon is sampled uniformly over the
entire frequency range. The total weight of all photons is equal to the 
luminosity of the system and each photon has weight has equal weight given by,

.. math ::
    w_{i} = \frac{\sum_{i = \text{sources}} ~ \int_{\nu_{\text{min}}}^{\nu_{\text{max}}} L_{\nu, i} d\nu}{N},

where :math:`N` is the total number of photons, :math:`\nu_{\text{min}}` and 
:math:`\nu_{\text{max}}` define the frequency range and :math:`L_{\nu, i}` is the
luminosity for a radiation source. Note that a summation is used to find the
luminosity for each radiation source and that :math:`w` has units of 
:math:`\text{ergs s}^{-1}`.

Banded Sampling
===============

In practice, uniform sampling is generally an inefficient approach to generating
photon packets. For example, it is often desirable to produce a sufficient 
number of photons within a specific frequency range, i.e. around a 
photoionisation edge. However, if these frequencies lie on the Wien tail of a
blackbody distribution, it is unlikely that a sufficient number of photons will
be generated as most of the luminosity is generated at lower frequencies.
It is possible to get around this problem by generating an increasingly large
number of photons. But, this is computationally expensive and inefficient.

In order to cope with cases this, SIROCCO implements *importance sampling* which
effectively increases the number of photons which are sampled from specific 
frequency bands considered important. Photons are now generated with the weight,

.. math ::
    w_{j} = \frac{\sum_{j = \text{sources}} ~ \int_{\nu_{i}}^{\nu_{i + 1}} L_{\nu, j} ~ d\nu}{f_{i} N},

where, again, this is a summation over all radiation sources. :math:`N` is the
total number of photons, :math:`f_{i}` is the fraction of photons emerging from
frequency band :math:`i`, :math:`\nu_{i}` and :math:`\nu_{i+1}` are the lower
and upper frequency boundaries for each frequency band and :math:`L_{\nu, j}` is
the luminosity of the radiation source. Hence, more photons from frequency bands
with a larger fraction :math:`f_{i}` will be generated. However, photons from
*important* bands (where more photons are sampled from) will have reduced 
weight, whilst photons from the *less important* frequency bands will have
increased weight. 

This scheme has the benefit of allowing the generation of a lower number of
photons whilst still sufficiently sampling important frequency ranges, 
decreasing the computational expense of a simulation. As such, this is the
preferred sampling method.

Available Sampling Schemes
==========================

SIROCCO currently implements seven pre-defined frequency bands and and two
flexible *run time* banding schemes. The parameter used to define the photon
sampling scheme is,

``Photon_sampling.approach(T_star,cv,yso,AGN,min_max_freq,user_bands,cloudy_test,wide,logarithmic)``

.. admonition :: Minimum and Maximum Wavelengths

    At present, the largest wavelength a photon can be is hardwired to 20,000
    Angstroms. The smallest wavelength a photon can take is defined by the 
    temperature of hottest radiation source, but is at least 115 Angstroms - 
    twice that of the Helium edge.

T_star
------

Create a single frequency band given a temperature T, which is the temperature
of the hottest radiation source in the model. All photons will then be 
generated from this single frequency band.

CV
--

Pre-defined bands which have been tuned for use with CV systems, where a hot
accretion disk (~100,000 K) is assumed to exist. In this scheme, there are four
bands where the majority of photons are generated with a wavelength of 912
Angstroms or less.

YSO
---

Pre-defined bands which have been tuned for use with YSO systems. In this 
scheme, there are four bands.

AGN
---

Pre-defined which have been tuned for use with AGN system. In this scheme, there
are ten bands, with a minimum frequency of :math:`1 \times 10^{14}` Hz and a 
maximum frequency of :math:`1 \times 10^{20}` Hz.

min_max_freq
------------

Create a single band using the minimum and maximum wavelength as described by
the minimum and maximum wavelengths calculated for the current model.

user_bands
----------

This allows a user to create their own frequency bands, defined by photon
energies measured in electron volts. The first band has the lowest photon energy
and each subsequent band must have a larger energy than the previous band. Each
band also requires a minimum fraction of photons to be sampled from this band,
where the sum of the fractions for each band must be equal to or less than one.

.. admonition :: Maximum Number of Bands

    Currently, a maximum of 20 frequency bands can be defined. If a user 
    attemps to specify more than than 20 bands, SIROCCO will create an error
    message and fallback to using 20 bands.

cloudy_test
-----------

This set of bands were created for use in testing against the photoionisation 
and spectral synthesis code Cloudy_.

.. _Cloudy: https://www.nublado.org

wide
----

Pre-defined bands which have very wide frequency range. The purpose of this
band is for testing, hence is best to avoid using this band for a working model. 

logarithmic
-----------

This is the same as ``user_bands``, however the frequency bands are now defined
in log space. This allows one to better sample a frequency range which spans
many orders of magnitude. 

.. admonition :: Maximum Number of Bands

    Currently, a maximum of 20 frequency bands can be defined. If a user 
    attemps to specify more than than 20 bands, SIROCCO will create an error
    message and fallback to using 20 bands.

.. admonition :: Minimum Fraction

    For logarithmic user defined bands, the fraction of each band is set to
    1 / nbands.
