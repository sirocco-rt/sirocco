Atomic Data
###########

Any SIROCCO model is only as good as the atomic data which goes into making the model.  
All of the atomic data that SIROCCO accepts is read in by the routine get_atomicdata,
and all of the data is read in from a series of ascii data files and stored in structures
that are defined in atomic.h.


The purpose of documentation is as follows:

* to explain the atomic data formats used by SIROCCO and the relationship of different sets
  of data to one another

* to explain where the data currently used in SIROCCO and to explain how the raw data 
  is translated in to a format the SIROCCO accepts

The routines used to translate raw data format for two-level atoms (as well as much of the raw data)
are contained in a separate github `repository <https://github.com/agnwinds/data-gen>`_
These routines are very "rough-and-ready", and not extensively documented, but so users
should beware.  On the other hand, they are not exceptionally complicated so in most cases
it should be fairly clear from the code what the various routines do. 

The routines used to generate data for MacroAtoms are described in :doc:`Generating Macro Atom data <./py_progs/MakeMacro>`

Choosing a dataset 
-----------------------
The "masterfile" that determines what data will be read into SIROCCO is determined by the
line in the parameter file, which will read something like::

    Atomic_data                                data/standard80.dat
   
where the file `data/standard80.dat` will contain names (one to a line) of files which will
be read in sequentially.  

All of the atomic data that comes as standard with SIROCCO is stored in the `xdata` directory (and its subdirectories) but users are not required to put their data
there. Various experimental or testing dataset masterfiles are stored in the `zdata` directory. Symbolic links to these directories
are setup by running `Setup_Py_Dir`.

.. todo::

    Add table of recommended data sets

Data hierarchy and I/O 
-----------------------
As mentioned above, the masterfile will contain names (one to a line) of files which will
be read in sequentially. Every line in the atomic data files read by SIROCCO consists of a keyword that defines the type
of data and various data values that are required for that particular data type.  Lines that
beging with # or are empty are ignored.

The data from the various files are read as if they were one long file, so how the data is split up into files is a matter of convenience.  

However, the data must be read in a logical order.  As an simple example, information about elements 
must be read in prior to information about ions.  This allows one to remove all data about, say Si,
from a calculation simply by commenting out the line in the atomic data that gives the 
properties of the element Si, without having to removed all the ion and other information about
a data file from the calculation.

The main hierarchy is as follows

 * Elements
 * Ions
 * Energy levels
 * Lines

Once these sets of data have been read in the order in which other information is read in 
is irrelevant, that is one can read the collision data (which is indexed to lines) and 
photoionization cross sections  (which are indexed to energy levels) in either order

*It is important that all of the ions, energy levels, and line for an atom, appear before 
photoionization cross sections and collisional data for that atom this data is only read in if it can be indexed to the corresponding levels and lines*

    

(Note that although the  approach has advantages of allowing one to remove or add atoms or ions from a model easily, it requires
one to be careful when doing this.  The 
program does not give you a good summary of 
what data has been omitted.  If concerned about this one should use the advanced 
command::

 @Diag.write\_atomicdata(yes,no)

which prints out an ascii version of the input data that was used.)

More information on the various types of input data can be found below:

.. toctree::
   :glob:

   atomic/*
