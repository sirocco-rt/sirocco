### examples

This folder contains a number of standard parameter files, used for testing and porting to the user.

### Basic

Core files which are intended to be examples that new users would use to explore Pyton

* cv_standard.pf -- Standard CV model (Shlosman & Vitello 1993; Long & Knigge 2002).
* fiducial_agn.pf -- QSO model from Higginbottom et al. (2013).
* 1d_sn.pf -- simple Supernovae model as presented in Kerzendorf & Sim (2013).
* star.pf -- simple spherical stellar wind model


### Extended
More performance intensive  models that might serve as the basis for serius reseach with Python:

* cv_macro_benchmark.pf -- Macro-atom CV run from Matthews et al. (2015).
* m16_agn.pf -- Macro-atom AGN run from Matthews et al. (2016).


### beta

Files of current interest but which may not be fully tested.  These may be parameter files
we are working on for scientific or development reaasons

* ngc5548.pf -- designed to represent rough params for NGC 5548
* ulx1.pf -- First attempt at a ULX model
* agn_ss_2010_modela.pf -- AGN model A from Sim et al. (2010).
* lamp_post.pf -- beta model for testing lamp post geometry

### polarization

Models that use the polarization option (`@Spectrum.polarization yes`, run with `sirocco -d`).
See the README there.

* torus_wedge.pf -- an AGN torus in X-rays, the wedge torus of Vander Meulen et al. (2024), with
  make_torus.py, which writes its density structure.

### regress

This is a standard set of files which are used for regression testing.  They are intended
for use with the python script regression.py

### regress_hydro

A special directory for testing that we can still run the type of model needed for working
with the hydro code zeus and being used for studies of X-ray binary thermal winds

### travis

A special directory which contains the parameter files necessary to run Travis.
