Python Scripts
################################

SIROCCO can be analysed using the **PySI** package which can be found `on github <https://github.com/sirocco-rt/pysi>`_, with some example usage under :ref:`Plotting and Processing Outputs <plotting>`.

There are also several python scripts written to prepare input for and analyse the output of SIROCCO. Some of the more useful scripts/modules are documented below, organised approximately by application. 

.. admonition :: Warning to user

    The scripts documented here form an incomplete and inhomogenous list, in the sense that they have been developed by different people at different times and do not fit nicely together as a single python package.
    Some of the scripts should still be useful, particularly if you consult example notebooks, but use with caution!

You can also generate documentation for all the scripts by navigating to :code:`docs/pydocs` and running :code:`write_docs.py`.
The resulting output file can be found `here <../../pydocs/doc_index.html>`_.



.. toctree::
    :maxdepth: 2
    :glob:

    py_progs/*
