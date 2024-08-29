# Sirocco

Sirocco is currently in the process of being renamed (from Python, for obvious reasons). Soon, this repository URL will point to the full repository with issue tracking and so on, but for now it is a placeholder. The code in the source code and documentation is referred to as Python. 

Python is a (confusingly named) Monte Carlo radiative transfer code which uses the Sobolev approximation. It has been developed by Knox Long, Christian Knigge, Stuart Sim, Nick Higginbottom, James Matthews, Sam Manghamm Edward Parkinson, Mandy Hewitt and Nico Scepi.The code has been used for a variety of research projects invovling the winds of cataclysmic variables, young stellar objects, X-ray binaries and AGN.

The code is under active development, but we are looking for beta users to test the code, and potentially use it for their own research. If you are interested in using Python please contact Knox Long via long[at]stsci[dot]edu. 

Documentation is hosted on [ReadTheDocs](http://agnwinds.readthedocs.io/en/dev/).

## CI \& Docs Build Status

[![C/C++ CI](https://github.com/agnwinds/python/actions/workflows/build.yml/badge.svg)](https://github.com/agnwinds/python/actions/workflows/build.yml)

[![Documentation Status](https://readthedocs.org/projects/agnwinds/badge/?version=latest)](https://agnwinds.readthedocs.io/en/latest/?badge=latest)

## Installation

Python and the various routines associated are set up in a self-contained directory structure. The basic directory structure and the data files that one needs to run Python need to be retrieved and compiled. 

If you want to obtain a stable (!) release, go to the [Releases](https://github.com/agnwinds/python/releases) page.

If you want to download the latest dev version, you can zip up the git repository by clicking on the zip icon to the right of the GitHub page. Aternatively, you can clone the repository using 

    $ git clone https://github.com/agnwinds/python.git 

If you anticipate contributing to development we suggest forking the repository and submitting pull requests with any proposed changes.

Once you have the files, you need to cd to the new directory and set your environment variables
    
    $ export PYTHON = /path/to/python/
    $ cd $PYTHON 
    $ ./configure
    $ make install  (or better make install 2>&1 | tee today.txt)
    $ make clean

If you have any difficulties with the installation, please submit an issue, along with the file today.txt

Note that the export syntax is for bash- for csh use 
  
    $ setenv PYTHON /path/to/python/

Atomic data for the current version of Python stored in the directory xdata which is part of the main repository,

However, if one wishes to use model stellar spectra to simulate the spectra of the star or disk, one my wish to
also to download various model grids of spectra that have been used in conjunction with Python over the years. These
are in a separate [models repository]((https://github.com/agnwinds/xmod).  

These can be downloaded as follows:

    $ cd $PYTHON; git clone https://github.com/agnwinds/xmod xmod 

(Previously, both the atomic data and the model grids were stored in a separate repository.  Users wishing
to run older versions of the code pre-84b may need to download the 
[old data repository](https://github.com/agnwinds/data)  This repository can be downloaded as follows


    $ cd $PYTHON; git clone https://github.com/agnwinds/data data

Those users interested in the current version of Python should not need to do this)

## Running python

To run python you need to add the following to your $PATH variable:

    $PYTHON/bin

You can then setup your symbolic links by running 

    $ Setup_Py_Dir

and run the code by typing, e.g.

    $ py root.pf


Please see the [ReadTheDocs](http://agnwinds.readthedocs.io/en/dev/) or the docs folder for how to use the code. You will need sphinx installed to build the documentation yourself. 

Any comments, email one of the addresses above.
