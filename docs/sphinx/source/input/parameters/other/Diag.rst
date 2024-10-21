Diag
====
A series of advanced/diagnostic commands (Accessed using the ``-d`` flag, see :ref:`Running SIROCCO`). 
The commands generally allow access to additional information from the simulation, 
or allow more precise control. Advanced commands have an @ symbol in front of them in the parameter file.

Note that some of these commands will also change the radiative transfer treatment in the code. 
A number of them are also only accessed when the ``Diag.extra`` command is answered with yes. 

Below is all possible diagnostic options. Some diagnostics may not be available depending on options selected.


.. toctree::
   :glob:

   Diag/*