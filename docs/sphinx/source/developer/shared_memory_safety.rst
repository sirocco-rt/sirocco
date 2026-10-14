Hazards of the Shared Memory Model
########################################

:doc:`mixed_memory_model` describes how the shared memory model is implemented.  
This page describes some pitfalls or hazars to be aware of when writing new code
for SIROCCO. 

Since version 2.0, some arrays that look like ordinary
per-process C arrays are in fact **one copy per node, written by several
processes at once**.  AS a result, innocuous-looking code that is correct 
in a serial or pure-MPI build can silently give wrong answers when run on 
more than one rank. 

.. contents:: On this page
   :local:
   :depth: 1


Who may write what
------------------

The rules that keep shared data correct are:

1. **One writer per element.**  A shared element must be written by exactly
   one rank (per node), or by every rank with the *same* value at a time when
   no rank is reading it.  Normally the writer of cell ``n`` is the rank whose
   range from ``get_parallel_nrange()`` contains ``n``.
2. **Barrier between writing and reading.**  If rank A writes shared data that
   rank B on the same node will read (or overwrite), there must be an
   ``MPI_Barrier (node_comm)`` between the write and the read.
3. **Shared data is read-only during photon transport.**  Anything written
   while photons are in flight belongs in a private ``est`` field.
4. **Cross-node consistency still needs communication.**  Shared memory only
   spans one node.  Values computed on one node reach other nodes through the
   usual broadcast routines described in :doc:`mixed_memory_model`.


Running without shared memory
=============================

Shared memory can be switched off at run time with the ``-no-shared-memory``
command-line switch:

.. code:: bash

    mpirun -np 4 sirocco -no-shared-memory model

The log then reports ``MPI topology: 4 ranks across 4 node(s)`` and ``MPI
shared memory disabled``.

This works by giving every rank a node communicator of its own
(``MPI_Comm_split`` on the rank number, in place of ``MPI_Comm_split_type``
with ``MPI_COMM_TYPE_SHARED``), so each rank behaves as a single-rank node.
The allocation routines still call ``MPI_Win_allocate_shared``, but on a
one-rank communicator, so every rank gets a private copy of every array, and
the grids are kept consistent by the cross-node broadcasts described in
:doc:`mixed_memory_model`.  No other code path changes.

Results are bit-for-bit identical with and without the switch for the same
number of ranks; only memory use changes.  This makes the switch useful as:

- **a fallback** for users if shared memory causes problems on a particular
  system or MPI installation;
- **a diagnostic**: if a problem disappears with ``-no-shared-memory``, it
  is related to shared memory;
- **a test**: see `Comparison with shared memory switched off`_.

The node communicator is created immediately after ``MPI_Init``, before the
command line is parsed in ``parse_command_line()``, so ``main()`` checks for
the switch separately with ``no_shared_memory_requested()``
(``para_update.c``).  ``parse_command_line()`` then accepts the switch and
logs it.  ``rad_hydro_files`` accepts the switch in the same way.

Code that is only correct when *every* rank on a node shares memory (for
example, relying on a node partner's write instead of a broadcast) will give
different results with the switch on.  That is a bug, and the comparison
below is designed to catch it.


Hazard 1: zero everything, then fill your own cells
====================================================

This pattern is common in the older code and is correct when every rank has
its own copy:

.. code:: c

    /* every rank zeroes ALL cells ... */
    for (n = 0; n < NPLASMA; n++)
      for (m = 0; m < nlevels_macro; m++)
        macromain[n].derived.matom_emiss[m] = 0.0;

    /* ... then fills only its own cells */
    get_parallel_nrange (rank_global, NPLASMA, np_mpi_global, &my_nmin, &my_nmax);
    for (n = my_nmin; n < my_nmax; n++)
      for (m = 0; m < nlevels_macro; m++)
        macromain[n].derived.matom_emiss[m] += ...;

With ``matom_emiss`` in shared memory, a slow rank can still be in the zeroing
loop while a fast rank has finished filling its cells, and the slow rank wipes
them out.  The result depends on timing, so repeated runs with identical seeds
give different spectra.  For a while, this was a real bug in the memory development branch
(``get_matom_f()`` and ``get_matom_f_accelerate()`` in ``macro_gen_f.c``) and
lost a small fraction of the macro-atom emissivity on 4 ranks.

**Fix** — either zero only your own cells, or wait for all ranks on the node to
finish zeroing before anyone fills:

.. code:: c

    #ifdef MPI_ON
      if (np_mpi_global > 1)
        MPI_Barrier (node_comm);
    #endif

The same pattern applies to any shared array that is reset and then
accumulated, e.g. per-ion recombination arrays in ``derived``.


Hazard 2: read-modify-write on shared data
==========================================

An update that reads a shared value and writes back something computed from
it is applied **once per rank on the node** if every rank executes it.
Consider this rescaling of ion densities after the mass density changes
(adapted from ``hydro_restart()`` in ``hydro_import.c``):

.. code:: c

    for (n = 0; n < NPLASMA; n++)        /* executed by every rank */
    {
      old_density = plasmamain[n].state.rho;            /* private scalar */
      plasmamain[n].state.rho = new_rho;                /* private scalar */
      for (nion = 0; nion < nions; nion++)              /* SHARED array   */
        plasmamain[n].state.density[nion] *= plasmamain[n].state.rho / old_density;
    }

``state.rho`` is a private scalar, so each rank rescales it once, correctly.
``state.density[]`` is shared, so on a node with *R* ranks it is multiplied by
the factor *R* times.  Interleaving between the ranks can also make the result
non-reproducible.

Typical read-modify-write operations to look out for: ``*=``, ``+=``,
``-=``, ``x = f(x)``, ``vadd (a, b, a)``, ``rescale (a, s, a)`` and damping or
under-relaxation updates such as ``x = (1 - g) * x_old + g * x``.

**Fix** — make sure only one rank per node performs the update, then
synchronise:

.. code:: c

    if (node_rank == 0)                  /* or loop only over this rank's cells */
    {
      for (nion = 0; nion < nions; nion++)
        plasmamain[n].state.density[nion] *= scale;
    }
    #ifdef MPI_ON
    MPI_Barrier (node_comm);
    #endif

If the update is part of the normal per-cell wind update (each rank updating
only its own ``get_parallel_nrange()`` cells), it is already safe: each
element has a single writer.


Hazard 3: writing shared data during photon transport
=====================================================

During transport every rank works on every cell at once.  Any write to
shared memory is a race, even a temporary one that is "restored" afterwards.
In the 2.0 development branch ``sobolev()`` used to modify
``state.density`` temporarily; it now passes an override value instead.

**Fix** — keep temporary copied values in local variables, function
arguments, or a private ``est`` field.  If a value must be cached per cell,
compute it during the wind update instead, before transport starts.

Writing the same value from every rank (e.g., the bound-free
look-up table ``state.kbf_use[]`` built by ``kbf_need()``, see below) 
does not change the
result on current hardware, but it is still a data race in the C and MPI
memory models and should be avoided


Hazard 4: collectives on shared buffers
=======================================

Do not pass a pointer into a shared block as the send or receive buffer of an
MPI collective (``MPI_Allreduce``, ``MPI_Bcast`` to all ranks, etc.).  On a
node with *R* ranks the same memory is passed *R* times.  ``MPI_Allreduce``
with ``MPI_SUM`` then adds the node's single copy *R* times, and with
``MPI_IN_PLACE`` several ranks write the result into the same memory at once.

**Fix** — copy into a private helper buffer, reduce that, and have one writer
copy the result back.  This is what the existing ``reduce_*`` routines do.  


Visibility and memory ordering
==============================

Rule 2 above uses ``MPI_Barrier (node_comm)`` to make one rank's writes visible
to another.  Strictly, the MPI-3 standard asks for ``MPI_Win_sync`` (inside a
passive-target epoch such as ``MPI_Win_lock_all``) in addition to process
synchronisation.  SIROCCO currently relies on ``MPI_Barrier`` alone.  This works
in practice with common MPI implementations on x86 processors, 
and all tests of the shared memory model so far have passed
with it

Apple Silicon and other ARM processors (including ARM HPC nodes) order memory
more weakly.  This makes them a stricter test of the synchronisation than
x86 Linux, so it is worth testing an ARM Mac when e.g., Barriers are changed.


Checklist for code review
=========================

When reviewing a change, you can apply this checklist for each modified line that writes to
``wmain``, ``plasmamain`` or ``macromain``:

#. Is the destination a shared array (or any ``wmain`` field)?  If not, stop
   here.
#. Is it written during photon transport?  It must not be.
#. Do several ranks on a node execute this write?  If so, do they all write
   the same value, and does nobody read it in the meantime?
#. Is it a read-modify-write?  Then exactly one rank per node may execute it.
#. Is there a reset of *all* cells followed by a fill of *own* cells?  Then
   there must be a ``MPI_Barrier (node_comm)`` between them.
#. Will another rank read this value next?  Then there must be a barrier (on
   the node) and a broadcast (across nodes) before it does.
#. Is shared memory passed directly to an MPI collective?  It should not be.

New shared arrays should also be added to the tables in :doc:`mixed_memory_model`
together with the barriers that protect them.


An idealised testing workflow
=============================

.. note::

   This is a pretty exhaustive workflow for testing memory changes, the kind of thing
   we carried out on the memory PR merge. It does not test correctness of the code,
   but is also not necessary to run in many cases. 

All of the hazards above produce either results that change
from run to run, or results that depend on how ranks are laid out across
nodes.  SIROCCO's fixed random seeds make both easy to detect, 
because for a given number of ranks,
results must be bit-for-bit identical regardless. 

The tests below compare the files produced by
``windsave2table`` and the ``.spec`` file.  
Any differing value is a failure, regardless of any tolerance. 

Short runs are enough to test this: e.g.,  for the regression models in ``examples/regress``,
20,000 photons per cycle with two ionization cycles and one spectral cycle
take about a minute each on four ranks.  The macro-atom models
(``agn_macro``, ``agn_macro_matrix``) exercise the most shared memory and
are the best to include. LLMs like Claude Code are quite effective for this 
sort of workflow and for spotting the hazards above.

Repeat test
-----------

Run the same model twice with the same number of ranks (at least 4) and
compare:

.. code:: bash

    mpirun -np 4 sirocco agn_macro          # in directory run1/
    mpirun -np 4 sirocco agn_macro          # in directory run2/
    # then compare windsave2table output and agn_macro.spec

Any difference indicates a race.  This test caught the ``matom_emiss`` race
described above within a one-minute run.

Comparison with the previous version
-------------------------------------

Run the same model and number of ranks on the shared-memory build and on a
build without the change being tested.  
Changes that only affect memory layout
must give identical results.

Comparison with shared memory switched off
------------------------------------------

Run the same build twice with the same number of ranks, once normally and
once with ``-no-shared-memory``, and compare.  The two runs must be
identical.  

As with the multi-node test below, switching shared memory off changes the
communicator layout, but for the MPI implementations tested so far this has
not changed the results.  If you see differences on a new system, also compare
the switched-off run with the pre-memory refactor code before concluding there is a bug.

Multi-node layout on one machine
--------------------------------

Code paths that only run when there is more than one node (broadcasts that
cross nodes, ``num_nodes > 1`` branches) can be tested on a single laptop with
MPICH, which can be told to impersonate several nodes:

.. code:: bash

    # two "nodes" (ranks are split between them round-robin)
    MPIR_CVAR_NUM_CLIQUES=2 mpirun -np 4 sirocco agn_macro

    # every rank on its own "node": no memory is shared at all
    MPIR_CVAR_NOLOCAL=1 mpirun -np 4 sirocco agn_macro

Check that the log reports the intended layout, e.g.
``MPI topology: 4 ranks across 2 node(s), 2 ranks on this node``.  The
setting must reach every rank; with some launchers you need to export it
explicitly (``-x`` for Open MPI's ``mpirun``, ``-genv`` for MPICH's Hydra).

.. important::

   Changing the node layout changes how MPICH performs ``MPI_Allreduce``
   internally, which changes the order of floating-point additions.  Results
   for different layouts are therefore **not** expected to agree with each
   other, even for code with no shared memory at all.  Compare the shared
   memory build against a reference build **with the same layout** (and
   repeat runs with the same layout against each other), not one layout
   against another.

These variables are specific to MPICH.  Open MPI has no equivalent, so a
genuinely multi-node Open MPI run needs two machines or a cluster.

Testing on ARM
--------------

Run the repeat test on an Apple Silicon Mac or an ARM cluster node whenever
barriers or shared allocations change (see `Visibility and memory ordering`_).


Known gaps
==========

These items were identified in the review of the shared memory model
(:pr:`1180`) and do not affect results at present, but do not fully follow
the rules above:

.. list-table::
   :header-rows: 1
   :widths: 30 40 30

   * - Location
     - Issue
     - Status
   * - Per-ion recombination arrays (``derived``)
     - Reset on all cells without a barrier before cells are filled
       (hazard 1)
     - Used for diagnostics only
   * - ``kbf_need()`` in ``resonate.c``
     - Every rank writes the same values to shared ``state.kbf_use[]``
       (hazard 3)
     - Benign on current hardware
   * - ``update_wind_from_model()``
     - Read-modify-write on shared data (hazard 2)
     - Only called serially at present
   * - ``hydro_restart()`` in ``hydro_import.c``
     - Read-modify-write on shared ``state.density[]`` from every rank
       (hazard 2)
     - Deprecated ``-z`` route; fix before re-enabling
