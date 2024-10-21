Ionization Cycles
#################

In order simulate a spectrum from a parameterized model of an outflow, one must first 
determine the ionization state of the wind.  In order to accomplish this, one begins 
with a guess at the ionization structure, usually by setting the temperature of the wind 
at a specific value and assuming that the ionization equilibrium is simple given by the 
Saha equation for that particular temperature.

In SIROCCO, one then generates a set of photon bundles over a wide frequency range, and then 
causes these photons to pass through and interact via various processes with the wind.  
As the photons transit the wind, estimatores for various processes are accumulated, which
characterize the intensity and spectrum  of the radiation field in various parts of 
the wind, the amount of heating and rate at which ions are photoionized, etc. 

Once all of these photon bundles have passed through the wind, one uses the various 
estimators to modify the ionization state and electron temperature in each cell, and then one repeats the process in order to try to find the actual state of the wind, given the 
assumed density and velocity field of the wind.  There are a variety of approaches to 
carrying out this calculation and various limitations placed on the rate at which the 
plasma is is allowed to change between cycles.  As the accuracy of any Monte Carlo simulation depends on numbers of photons bundles one uses to approximate the spectrum there are various
options within SIROCCO to choose the number of photons with various energy/wavelength bins, 
and other options to begin with a smaller number of photons and increase this number in 
later cycles.  
