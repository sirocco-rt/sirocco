Elements and Ions
#################


The first file that must be read into *SIROCCO* is the file that defines the elements and ions.  The 

Source:
=======
This data comes from `Verner, Barthel & Tytler, 1994, ApJ 108, 287. <http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1994A%26AS..108..287V&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf>`_



Translation to python:
======================

The original data and the translation can be found in py\_verner.  A simple awkscript converts the downloaded data to SIROCCO format.


Data Format
===========


There are two sections to the file, first elements are defined

The first portion of a typical file is as follows::

    # Source Verner, Barthel & Tytler 1994 A&AS 108, 287
    Element    1    H   12.00     1.007940
    Element    2   He   10.99     4.002602
    #Element    3   Li    3.31     6.941000
    #Element    4   Be    1.42     9.012182
    #Element    5    B    2.88    10.811000
    Element    6    C    8.56    12.011000
    Element    7    N    8.05    14.006740
    Element    8    O    8.93    15.999400


And the columns are as follows

* A label for this type of data entry
* The z of the elment  
* The common abbreviation for the elemen
* The atomic weight of the element

Lines beginning with # (and empty lines) are treated as comments.  So in this case, Li, B and B are ignored, because
of their relatively low abundance.

Abundances are generally defined logarithmically 
with respect to H at 12.00.  In principle, there are two choices if one
wished to defien a plasma where, for example, He was the dominant 
element.  One could leave the H abundance at 12 and define the He 
abundance as for example 13.00 Alternatively, one could set the He 
abundance to 12.00 and define all of the other elements with respect
to this.  Either choice should work but none has been tested. It is
unclear whether code will work at all for a plasma with no H.


The ion section (which could be in a separate file) has the following format::

    IonV    H   1   1   2   13.59900 1000  10     1s(2S_{1/2})
    IonV    H   1   2   1 1.0000e+20   1   1     Bare

    IonV   He   2   1   1   24.58800 1000  10     1s^2(1S_0)
    IonV   He   2   2   2   54.41800 1000   10     1s(2S_{1/2})
    IonV   He   2   3   1 1.0000e+20   1   1     Bare


and the columns have the following meaning

* Label for an ion that will be treated as a simple ion
* The common abbreviation for the element
* z of the ion
* ionization state of the ion 
* g-factor for the ground state
* the ionizaton potential in eV
* maximum number of (simple) levels allowed
* maximum number of nlte (macro-atom) levels
* the configurations (which for information only)

The label for the ion entries determines whether an element will be treated as simple atom or as a macro-atom.  For case where H is to be treated as
a macro atom, but He is to be treated as a simple atom, this file would become::


    IonM    H   1   1   2   13.59900 1000  10     1s(2S_{1/2})
    IonM    H   1   2   1 1.0000e+20   1   1     Bare

    IonV   He   2   1   1   24.58800 1000  10     1s^2(1S_0)
    IonV   He   2   2   2   54.41800 1000   10     1s(2S_{1/2})
    IonV   He   2   3   1 1.0000e+20   1   1     Bare

Note that only evident changed is the label, but in this case the number of nlte levels, and not the number of levels  is what is important.  



SIROCCO structure:
=================
This data is held in SIROCCO in various fields in structures **elements** and **ions**.

Comments:
=========

**Maximun numbers of levels**

As indicated the numbers here are maximum values, and the actual numbers of levels for particular ion will depend on the data that follows. 
One can use the numbers here to limit the complexity of, for example, a macro-atom to see whether making a more complicated macro-atom affects
the reusult of a calculation.  One does not need to change the "downstream" data to make this happen, *SIROCCO* will simply ignore the extra
data.

