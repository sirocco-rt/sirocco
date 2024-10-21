Meta-documentation
##################

How to document SIROCCO
======================

This documentation is written in **ReStructured Text**, and parsed by **Sphinx**.
A general guide to **ReStructured Text** can be found `here <http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html>`_.
We're trying to maintain a roughly consistent format for the documentation.

Installing the documentation tools
----------------------------------

This guide is produced using **Sphinx**.
**Sphinx** is written in SIROCCO and available from the **pip** package manager.
We require the **Python 3** version of **Sphinx**. Install it, and the other modules required, as:

.. code :: bash

    cd docs/sphinx
    pip3 install -r requirements.txt

Building the documentation
--------------------------

Once **Sphinx** is installed, you can make the documentation using a **Makefile** as:

.. code :: bash

    make html

You can tell if the documentation was built successfully by looking at the output of ``make html``.
You should see:

.. code ::

    build succeeded.

    The HTML pages are in html.

If the build was successful then the documentation can be viewed by opening ``docs/sphinx/html/index.html``.
Many errors will not stop the build process.
Details on the build errors can be found in the section on :ref:`Common errors & warnings`.

You can make minor changes to the documentation and recompile using :code:`make html` again.
If you add new pages or move existing ones, the table of contents will need to be regenerated.
Do this via:

.. code :: bash

    make clean
    make html

General documentation
=====================

Conventions
-----------

Each file should have a title, with subsections within it created by underlining the titles as:

.. code :: rst

    Title
    #####

    Section
    =======

    Subsection
    ----------

    Subsubsection
    ^^^^^^^^^^^^^

When referring to a parameter, link to the documentation page as:

.. code :: rst

    The number of domains can be set by :ref:`Wind.number_of_components`.

When referring to files, code (e.g. shell script) or values used by the code, render it as ``monospaced`` text as:

.. code :: rst

    Run the program using ``py``.
    Set the parameter :ref:`Wind.type` to ``SV``.
    Outputs can be found in ``filename.rst``.

When referring to a library or program name, render it in **bold** as:

.. code :: rst

    Though this program is called *SIROCCO*, it is written in **C**, using the **GSL** library.

Content of interest to developers but not users should be broken out into a callout as:

.. code :: rst

    .. admonition :: Developer Note

        This value is only stored to single-precision

.. admonition :: Developer Note

    This is a developer note

Documentation that needs expanding should be indicated with a to-do callout as:

.. code :: rst

    .. todo :: Expand this section

.. todo :: This is a to-do note

Content relating to a specific **GitHub** issue/pull request can be linked directly to it as :issue:`1`/:pr:`56`:

.. code :: rst

  This arose due to issue :issue:`1`, which was fixed by :user:`kslong` using :pr:`56`.

When writing a table, use the full form where possible as:

.. code :: rst

    +----+----+
    |Name|X   |
    +----+----+

+----+----+
|Name|X   |
+----+----+


Parameter documentation
=======================

Formatting
----------
Parameters are documented in a consistent way. They have a set of properties.
Not every parameter will have all properties but you should fill them all in where possible.
A full example outline is:

.. literalinclude :: reference_rst.txt
   :language: rst

The sections we expect are entered as a definition list.
A definition list consists of titles followed by a definition block indented by 2 characters.
The headings, in the order we expect, are:

Name
  The parameter name, as used by SIROCCO input files.

Description
  A description of the parameter and its function.
  This can include links to other pages and parameters, using the format

  .. literalinclude :: reference_rst.txt
     :language: rst
     :lines: 4-4

Type
  This is whether the parameter is an integer, float, or enumerator (a list of choices).

Unit
  This is the unit. It can be something like ``cm``, ``m`` or even derived from other parameters
  (e.g. :ref:`Central_object.radius`).

Values
  If the parameter is an integer or float, this should describe the range of values it can take.
  For example, ``Greater than 0`` or ``0-1``.

  If the variable type is ``Enumerator``, then instead it should include a nested definition list of
  the possible choices. Where each choice implies a different set of possible children
  (e.g. :ref:`Wind.type`) then each choice should have its own **Children** definition list, e.g.

  .. code :: rst

    SV
      * :ref:`SV.thetamin`
      * :ref:`SV.thetamax`

File
  The file the parameter is found in. This is a link to the file on the `master` branch.

Child(ren)
  If the parameter implies any others.
  For example, :ref:`Spectrum.no_observers` has child parameters :ref:`Spectrum.angle`.

Parent(s)
  If the parameter depends on another.
  For example, :ref:`KWD.rmax` is only required for a specific choice of :ref:`Wind.type`.


Locations
---------

Parameters are stored in ```docs/sphinx/source/inputs/parameters/``.

If multiple parameters share a root (i.e. ``SV.radmin``, ``SV.radmax``), then they should be stored within a directory with the
same root name as the parameters (i.e. ``SV/SV.radmin.rst``, ``SV/SV.radmax.rst``). In the level above that directory, there should
be a  ``.rst`` file with the same name that serves to link those files to the table of contents, as:


.. code :: rst

    SV
    ==

    Some description of the parameter group.

    .. toctree::
       :glob:

       SV/*

Storing all the parameters in one folder would result in it being unreadably busy. Instead, we sift the parameters into groups.
Where multiple different parameters or parameter folders fall into the same rough category (e.g. central object parameters,
wind types and the like) we create subfolders to group them into. The order that these appear in the sidebar can be set if you
enter the filenames explicitly in the ``docs/sphinx/source/input/parameters.rst`` file.


Common errors & warnings
========================

Undefined Label
  .. code ::

      /path/to/file.rst:line_number:
      WARNING: undefined label: label_name (if the link has no caption the label must precede a section header)

  This warning occurs when the :code:`:ref:'location'` format is used to link to a section that does not exist.
  Check the spelling

Duplicate Label
  .. code ::

      /path/to/file.rst:line_number:
      WARNING: duplicate label label_name, other instance in /path/to/other/file.rst

  This warning occurs when two sections have the same name. The **autosectionlabel** addon automatically creates a label
  for each title and section component. This is generally not a problem unless you *really* need to

Inline emphasis
  .. code ::

      /path/to/file:line_number:
      WARNING: Inline emphasis start-string without end-string.

  This warning occurs when a line contains an un-escaped \* character, as \* is used to denote *emphasis*.
  Either escape it with \\ (i.e. :code:`\*`) or wrap it in a \:code\: tag.


Bullet list ends without a blank line
  .. code ::

      /path/to/file.rst:line_number:
      WARNING: Bullet list ends without a blank line; unexpected unindent.

  This warning occurs when a bullet-list doesn't have a blank line between it and the next bit of text.
  It commonly happens when you accidentally forget to space a bullet and the text following it, e.g.

  .. code ::

      * blah1
      * blah2
      *blah3

Inline substitution_reference
  .. code ::

      /path/to/file:line_number:
      WARNING: Inline substitution_reference start-string without end-string.

  This warning occurs when you have a word immediately followed by a pipe that is not part of a table, e.g. :code:`something|`.
  It tends to occur during typos in table creation e.g.

  .. code :: rst

      +---+---+
      | a||b  |
      +---+---+

Documenting Python Scripts
===========================

The :doc:`py_progs` page is intended to document various python scripts contained within the py_progs folder. The aim is to do this using Sphinx's `autodoc extension <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_, invoked by adding ``sphinx.ext.autodoc`` to extensions list in the conf.py file. py_progs is also added to the path using ``sys.path.insert(0, '../../py_progs/')``.

The above link contains full documentation of the commands. A module in py_progs can be documented by adding the following text to the rst file, where module.py is the name of the module you wish to document.

.. code :: rst

    .. automodule:: py_read_output
        :members:

For this to work properly, docstrings have to be in a reasonable rst format. We might consider using the `napoleon extension <https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html>`_ if this is not to our taste.

