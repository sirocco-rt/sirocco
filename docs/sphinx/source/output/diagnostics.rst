Diagnostic files
################

SIROCCO logs a considerable amount of information as it runs.
Some of this information is printed to the screen but a much more voluminous version of progress of the program is placed in a sub-directory,
named diag_whatever, where whatever is the root name of the model being run.

In this directory one will find log files, e.g. **whatever_0.diag**, **whatever_1.diag**, and so on,
where the in a multiprocessor run, the number refers to the log of a specific thread.

Inspecting these logs is important for understanding whether a SIROCCO run is successful,
and for understanding why if failed if that is the case. 

Additionally in this directory one will find a file **whatever.disk.diag** if the system contained a disk.  This file contains information
about the structure of the disk as a function of disk radius. It also indicates how many photons hit the disk with how much total energy
as a function of radius.  It provides several estimates of the effective temperature of the radiation that hits the disk, and determines 
how the disk temperature of the disk would so that it would reradiate the sum of the viscous energy and the energy associated with illumination.
