Reverb.path_bins
================
Number of bins for photon paths. Reverb modes that record the distribution of
path lengths in every wind cell bin them in this number of bins. Bins are
logarithmically spaced between the minimum scale in the system (the smallest
'minimum radius' in any domain) and the 10 * the maximum scale in the system
(10 * the 'maximum radius' in any domain). Default value is 1000, going much
higher does not lead to qualitative differences in TF, going lower makes the
bin boundaries show up in the TF.

Type
  Integer

Values
  Greater than 0

File
  `setup_reverb.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_reverb.c>`_


Parent(s)
  * :ref:`Reverb.type`: ``wind``, ``matom``


