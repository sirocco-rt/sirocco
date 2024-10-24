Disk.temperature.profile
========================
The choice of disk temperature profile

Type
  Enumerator

Values
  standard
    A Shakura - Sunyaev  disk, with a hard inner boundary

  readin
    Read the profile in from a file; the user will be queried for the name of the file


File
  `setup_disk.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_disk.c>`_


Parent(s)
  * :ref:`Disk.radiation`: ``True``


Child(ren)
  * :ref:`Disk.mdot`

  * :ref:`Disk.T_profile_file`

