Overview
########

SIROCCO uses a keyword-based parameter file to specify a model. Below is an example portion of a parameter file for a star system type (which must have the extension .pf). The parameters look as follows:

.. code::

   System_type(star,cv,bh,agn,previous)                 star

   ### Parameters for the Central Object
   Central_object.mass(msol)                      52.5
   Central_object.radius(cm)                  1.32e+12
   Central_object.radiation(yes,no)                  yes
   Central_object.rad_type_to_make_wind(bb,models)                   bb
   Central_object.temp                           42000
   ...


Each line begins with a keyword followed optionally by a comment in parentheses, and then a value, e.g

* **Keyword:** :code:`System_type`
* **Comment:** :code:`star,cv,bh,agn,previous`
* **Value:** :code:`star`

The comment generally specifies a set of valid choices or the units in which information is expected.

When a series of choices is presented, one does not need to enter the complete word, just enough to provide unique match to the choice.

The user does not need to create a parameter .pf file before running SIROCCO. Invoking SIROCCO without a parameter file will cause SIROCCO to prompt the user for the necessary information to create a parameter file. The user can specfiy any name for the parameter file. The example below calls the filename 'my_new_model'

.. Instead, assuming one is not working from a template parameter file, one simply invokes SIROCCO.

.. code:: console

   py my_new_model

SIROCCO then queries the user for answers to a series of questions, creating in the process a pf file, my_new_model.pf,
that can be edited and used in future runs.

An example of a line presented to the user in interactive mode is:

.. code::

   Disk.mdot(msol/yr) (1e-08) :

There the number in the second set of parenthesis is a suggested value of the parameter.
The user types in a new value and a carriage return, or, if the the suggested value seems appropriate,
responds with a carriage return, in which case the suggested value will be used.

The user can use the :code:`-i` switch when invoking SIROCCO. This indicates that SIROCCO should accumulate all of the necessary inputs, write out the parameter file, and exit, which is useful if one is not completely sure what one wants.

.. code:: console

   py -i my_new_model


Changes in the input files as the code evolves
----------------------------------------------

Occassionally, new input variables will be introduced into SIROCCO.  In this case, when one tries to run a parameter file 
created with a previous version of SIROCCO in single processor mode, the program will query the user for the parameters that are missing, and then
attempt to run the program as normal. 

If the original name of the parameter file was :code:`test.pf`, the modified version of the parameter file will be written to  :code:`test.out.pf`, so
one normally copies, in this case :code:`test.out.pf` to  :code:`test.pf` to avoid having the reenter the variable by hand if one wishes to run the parameter file a second time.

A better approach, if one is aware a change to the inputs has been made, is to run the old parameter file with :code:`-i` switch, copy the :code:`test.out.pf` to  :code:`test.pf`, and then
run the program normally.

Alternatively, if one heeds to modify a number of input files, once one knows what the change is, one can simply edit the .pf files directly.

(In multiprocessor mode, if the inputs have changed, the program will fail at the outset, requiring one to got through the process of runnning the program with  the  :code:`-i` switch, copying the :code:`test.out.pf` 
to  :code:`test.pf`, and then running normally.)

