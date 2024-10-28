Diag.invoke_searchlight_option
================================

The searchlight option is helpful for tracking the photon paths from a specific angle and their contribution to the final spectrum. Hence, this option is best used in tandem with saving the photon paths with :ref:`Diag.save_photons`.


Type
  Boolean (yes/no)

File
  `diag.c <https://github.com/sirocco-rt/sirocco/blob/master/source/diag.c>`_


Parent(s)
  * :ref:`Diag.extra`: :code:`yes`

Child(ren)
  * :ref:`Diag.location` :code:`central_object`
      * :ref:`Diag.angle`
  * :ref:`Diag.location` :code:`disk`
      * :ref:`Diag.angle`
      * :ref:`Diag.r`