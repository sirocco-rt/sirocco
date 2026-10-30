Mixed Memory Model and MPI Communication
#########################################

Overview
========

Version 2.0 of SIROCCO introduced a **mixed memory model** that significantly
reduces per-rank memory consumption on multi-core nodes.  Understanding this
model is essential for any developer working on parallelisation, data
communication, or memory allocation.

SIROCCO uses MPI to run across multiple ranks (independent processes, each
with their own address space).  There are two distinct inter-rank communication
patterns, applied at different points in the simulation cycle:

**Broadcast** (one rank → all ranks)
    Used *after* each ionization or spectral cycle, once a rank has updated
    its share of the wind or plasma grid.  Each rank in turn packs its updated
    cells into a buffer and calls ``MPI_Bcast`` so every other rank receives
    the result.  This is how ``state`` and ``derived`` sub-structure fields are
    kept consistent across ranks.

**Reduce** (all ranks → one combined result)
    Used *after* photon transport, to combine the partial radiation-field
    estimators that each rank accumulated independently while transporting its
    subset of photons.  All ranks call ``MPI_Allreduce`` (or an equivalent
    loop) to sum their per-rank tallies into a single global value.  This is
    how ``est`` sub-structure fields are aggregated.

Which pattern applies to a given variable is determined entirely by which
sub-structure it belongs to — ``state`` and ``derived`` fields are broadcast;
``est`` fields are reduced.  See `Choosing the right communication pattern`_
for a decision guide before diving into the implementation details.

Memory model: version 1.2  vs version 2.0
======================================================

Version 1.2:  pure MPI, all data duplicated
---------------------------------------------------------

In version 1.2 every MPI rank held its **own independent, complete copy** of
all data structures — the plasma grid (``plasmamain``), the macro-atom grid
(``macromain``), and the wind geometry array (``wmain``).  Data was shared
between ranks exclusively through explicit MPI message-passing calls
(``MPI_Pack`` / ``MPI_Bcast`` / ``MPI_Unpack``).

This meant that on a node with *R* ranks, the same read-only state data (ion
populations, spectral-model parameters, transition-probability matrices, etc.)
was stored *R* times in physical memory.  For large models with macro-atom
atomic data this could amount to tens of gigabytes of duplicated memory per
node.

Version 2.0 : mixed shared + private memory
-----------------------------------------------------------

Version 2.0 replaces the duplicated-copy model with a **mixed memory model**:

- **Shared memory (one copy per node)** — data that is read-only during
  photon transport is allocated using MPI-3 ``MPI_Win_allocate_shared``.
  All ranks on the same physical node map their pointers to the *same*
  physical pages.  There is therefore only one copy per node regardless of
  how many ranks that node hosts.

- **Private memory (one copy per rank)** — data that each rank writes
  independently during photon transport (radiation-field estimators, photon
  counters, scatter tallies) is allocated with ordinary ``calloc``.  Each
  rank has its own copy so writes never race.

- **MPI message passing (unchanged)** — after each ionization or spectral
  cycle the results from all ranks are combined via the same ``MPI_Pack`` /
  ``MPI_Bcast`` / ``MPI_Unpack`` pattern used in v1.2.  The broadcasts now
  write into the shared memory regions and are followed by
  ``MPI_Barrier(node_comm)`` calls to ensure coherence before transport
  resumes.

The split between shared and private follows the three sub-structure
decomposition of ``plasma_dummy`` and ``macro_dummy``:

.. list-table::
   :header-rows: 1
   :widths: 20 35 15 30

   * - Sub-structure
     - Typical contents
     - Memory in v2.0
     - Memory in v1.2
   * - ``state``
     - Ion populations, spectral model params
     - **Shared (for arrays), Private (for scalars)**
     - Private (duplicated)
   * - ``est`` (estimators)
     - Radiation field estimators, heating rates
     - **Private**
     - Private
   * - ``derived``
     - Cooling rates, luminosities, transition matrix
     - **Shared (for arrays), Private (for scalars)**
     - Private (duplicated)

The net effect is that for a model with *N* plasma cells, *I* ions, and *R*
ranks on one node, the dominant variable-length state and derived arrays
(roughly ``N × I × 14 × 8`` bytes) exist only once per node instead of *R*
times.  See `Memory savings`_ below for worked examples.

To be more explicit, with more than one MPI rank, the following live in node-shared memory
(``MPI_Win_allocate_shared``).  Every rank on a node reads and writes the
**same physical memory**:

- the whole wind array ``wmain[]``, including its scalar fields
- every *pointer* (array) field of ``plasmamain[n].state`` and
  ``plasmamain[n].derived``, e.g. ``state.density[]``;
- every *pointer* field of ``macromain[n].state`` and ``macromain[n].derived``,
  e.g. ``state.jbar_old[]``, ``derived.matom_emiss[]``;
- the macro-atom transition matrix block behind ``matom_matrix``.

The following are **private** (one copy per rank, as in version 1.x):

- the ``plasmamain[]`` and ``macromain[]`` arrays of structures themselves,
  so every **scalar** field is private, including those in ``state`` and
  ``derived`` such as ``plasmamain[n].state.t_e``, ``state.ne`` and
  ``state.rho``;
- every field of the ``est`` sub-structures (scalars and arrays);
- the exceptions listed below (``scatters``,
  ``xscatters``, ``n_bf_in``, ``n_bf_out``), photon stores, spectra and all
  other globals.

.. warning::

   The split is by **storage**, not by sub-structure.
   ``plasmamain[n].state.rho`` is private but
   ``plasmamain[n].state.density[i]`` is shared, even though both are "state".
   A loop that updates both behaves differently for each.

With one rank (or a build without ``MPI_ON``), everything is allocated with
``calloc`` and the code behaves exactly as it did before version 2.0.  None of
the problems described in :doc:`shared_memory_safety` can occur in serial, which is why serial tests cannot catch
them.

MPI communication files
=======================

All calls to MPI are isolated from the rest of SIROCCO.  Most, if not all, of
the MPI code is contained within five source files:

- :code:`communicate_macro.c`
- :code:`communicate_plasma.c`
- :code:`communicate_spectra.c`
- :code:`communicate_wind.c`
- :code:`para_update.c`

If you need to extend or implement a new function for MPI, please place it
either in one of the above files or create a new file using an appropriately
similar name.  Any parallel code should be wrapped by :code:`#ifdef MPI_ON`
and :code:`#endif` as shown below:

.. code:: c

    void communication_function(void)
    {
    #ifdef MPI_ON
        /* MPI communication code should go between the #ifdef's here */
    #endif
    }

Don't forget to update the Makefile and :code:`templates.h` if you add a new
file or function.

Choosing the right communication pattern
=========================================

The ``plasma_dummy`` and ``macro_dummy`` structures are each divided into three
sub-structures.  The sub-structure a variable belongs to determines which
communication pattern applies — and therefore which function to modify when
adding or removing a variable.

**The rule is simple:**

- ``est`` (estimator) fields are **reduced** — each rank accumulates
  partial results during photon transport; ``MPI_Allreduce`` sums them across
  all ranks.  Any variable that is incremented or accumulated as photons pass
  through the grid belongs here and must be handled via a reduce.
- ``state`` and ``derived`` fields are **broadcast** — one rank computes
  updated values during the wind update phase; ``MPI_Bcast`` distributes them
  to all ranks.  Any variable that is not written during photon transport
  belongs here and is placed in shared memory.

.. list-table:: Plasma sub-structures and their communication
   :header-rows: 1
   :widths: 20 30 20 30

   * - Sub-structure
     - Contents
     - Pattern
     - Key functions
   * - ``state``
     - Thermodynamic state (``ne``, ``t_e``, ``t_r``, ``w``, ``rho``, ``vol``), ion populations (``density``, ``partition``, ``levden``), spectral model parameters, bound-free data
     - **Broadcast**
     - ``broadcast_updated_plasma_properties()``, ``broadcast_plasma_grid()``
   * - ``est``
     - Radiation field estimators (``j``, ``ave_freq``), heating rates (``heat_tot``, ``heat_lines``, etc.), photon counters, flux estimators, cell spectra, ionization estimators
     - **Reduce**
     - ``reduce_simple_estimators()``
   * - ``derived``
     - Cooling rates, luminosities, convergence diagnostics, scatter counts, persistent flux averages, ionization parameter (``xi``)
     - **Broadcast**
     - ``broadcast_updated_plasma_properties()``, ``broadcast_wind_luminosity()``, ``broadcast_wind_cooling()``

.. list-table:: Macro-atom sub-structures and their communication
   :header-rows: 1
   :widths: 20 30 20 30

   * - Sub-structure
     - Contents
     - Pattern
     - Key functions
   * - ``state``
     - Normalized rate coefficients (``jbar_old``, ``gamma_old``, ``alpha_st_old``, etc.), mode flags
     - **Broadcast**
     - ``broadcast_updated_macro_atom_properties()``
   * - ``est``
     - Raw Sobolev mean intensities (``jbar``), photoionization rates (``gamma``), stimulated recombination rates (``alpha_st``), macro-atom absorption, cooling stores
     - **Reduce**
     - ``reduce_macro_atom_estimators()``
   * - ``derived``
     - Macro-atom emissivities (``matom_emiss``), k-packet rate flags, transition probability matrix (``matom_matrix``)
     - **Broadcast**
     - ``broadcast_macro_atom_emissivities()``, ``broadcast_macro_atom_state_matrix()``

The sections below describe how each pattern is implemented and give step-by-step
instructions for adding or removing a variable from each.

Mode 1: Broadcast (state and derived variables)
================================================

Broadcast is used to distribute updated ``state`` and ``derived`` fields after
the wind update phase of each ionization cycle.  One rank at a time acts as
root: it packs its updated cell range into a buffer and calls ``MPI_Bcast`` so
every other rank receives the result.  The full loop ensures each rank
eventually sends its cell subset to all others.

As the data structures in SIROCCO are fairly complex and use
pointers/dynamic memory allocation, data must be manually packed and unpacked
into a contiguous communication buffer — a fairly manual (and error-prone)
process.

Calculating the broadcast buffer size
--------------------------------------

The size of the communication buffer has to be calculated manually, by
counting the number of variables being copied into it and converting this to
the appropriate number of bytes.  This is done by the
:code:`calculate_comm_buffer_size` function which takes two arguments: 1) the
number of :code:`int`'s and 2) the number of :code:`double`'s.  We have to
*manually* count the number of :code:`int` and :code:`double` variables being
communicated.  Due to the manual nature of this, great care has to be taken to
ensure the correct number are counted otherwise MPI will cause a crash during
communication.

When counting variables, one needs to count the number of *single* variables
of a certain type as well as the number of elements in an array of that same
type.  Consider the example below,

.. code:: c

    int my_int;
    int *my_int_arr = malloc(10 * sizeof(int));
    int num_ints = 11;

In this case there are 11 integer variables which will want to be communicated.
In practise, calculating the communication buffer is usually done as in the
code example below:

.. code:: c

    /* We need to ensure the buffer is large enough, as some ranks may be sending a smaller
       communicating buffer. When communicating the plasma grid for example, some ranks may send
       10 cells whilst others may send 9. Therefore we need the buffer to be big enough to receive
       10 cells of data */
    int n_cells_max = get_max_cells_per_rank(NDIM2);

    /* Count the number of integers which will be copied to the communication buffer. In this
       example (20 + 2 * nphot_total + 1) is the number of ints being sent PER CELL;
       20 corresponds to 20 ints, 2 * nphot_total corresponds to 2 arrays with nphot_total elements
       and the + 1 is an extra int to send the cell number. The extra + 1 at the end is used to
       communicate the size of the buffer in bytes */
    int num_ints = n_cells_max * (20 + nphot_total + 1) + 1;

    /* Count the number of doubles to send, following the same arguments as above */
    int num_doubles = n_cells_max * (71 + 2 * NXBANDS + 6 * nphot_total);

    /* Using the data above, we can calculate the buffer size in bytes and then allocate memory*/
    int comm_buffer_size = calculate_comm_buffer_size(num_ints, num_doubles);
    char * comm_buffer = malloc(comm_buffer_size);

Broadcast implementation
-------------------------

The general pattern for packing data into a communication buffer and then
broadcasting it between ranks is as follows:

- Loop over all the MPI ranks (in MPI_COMM_WORLD).
- If the loop variable is equal to a rank's ID, that rank will broadcast its
  subset of data to the other ranks.  This rank uses :code:`MPI_Pack` to copy
  its data into the communication buffer.
- All ranks call :code:`MPI_Bcast`, which sends data from the root rank (this
  is the rank which has just put its data into the communication buffer) and
  receives it into all non-root ranks.
- Non-root ranks use :code:`MPI_Unpack` to copy data from the communication
  buffer into the appropriate location.
- This is repeated until all MPI ranks have sent their data as root, and have
  therefore received data from all other ranks.

In code, this looks something like this:

.. code:: c

    char *comm_buffer = malloc(comm_buffer_size);

    /* loop over all mpi ranks */
    for (int rank = 0 ; rank < np_mpi_global; ++rank)
    {
        /* if rank == your rank id, then pack data into comm_buffer. This is the root rank */
        if (rank_global == rank)
        {
            /* communicates the number of cells the other ranks have to unpack. n_cells_rank
               is usually provided via a function argument  */
            MPI_Pack(&n_cells_rank, 1, MPI_INT, comm_buffer, ...);
            /* start and stop refer to the first cell and last cell for the subset
               of cells which this rank has updated or is broadcasting. stop and start
               usually are provided via function arguments */
            for (int n_plasma = start; n_plasma < stop; ++n_plasma)
            {
                MPI_Pack(&plasmamain[n_plasma]->nwind, 1, MPI_INT, comm_buffer, ...);
            }
        }

        /* every rank calls MPI_Bcast: the root rank will send data and non-root ranks
           will receive data */
        MPI_Bcast(comm_buffer, comm_buffer_size, ...);

        /* if you aren't the root rank, then unpack data from the comm buffer */
        if (rank_global != rank)
        {
            /* unpack the number of cells communicated, so we know how many cells of data,
               for example, we need to unpack */
            MPI_Unpack(comm_buffer, 1, MPI_INT, ..., &n_cells_communicated, ...);
            /* now we can unpack back into the appropriate data structure */
            for (int n_plasma = 0; n_plasma < n_cells_communicated; ++n_plasma)
            {
                MPI_Unpack(comm_buffer, 1, MPI_INT, ..., &plasmamain[n_plasma]->nwind, ...);
            }
        }
    }

This is the standard method for communicating data in SIROCCO, given the
complexity of the data structures.  Unfortunately there are not many structures
or situations where using a derived data type is viable because none of the
structures are contiguous in memory.

Adding a variable to a broadcast
----------------------------------

1. Identify the correct broadcast function from the table above (e.g.
   ``broadcast_updated_plasma_properties()`` for plasma ``state``/``derived``
   fields).
2. Increment the appropriate variable count in the call to
   :code:`calculate_comm_buffer_size`.  For example, if the new variable is an
   :code:`int` in the plasma grid, update
   :code:`n_cells_max * (20 + 2 * n_phot_total + 1)` to
   :code:`n_cells_max * (21 + 2 * n_phot_total + 1)`.
3. In the block where :code:`rank == rank_global`, add a new call to
   :code:`MPI_Pack` following the existing pattern.
4. In the block where :code:`rank != rank_global`, add a matching call to
   :code:`MPI_Unpack`.
5. If the variable lives in a shared-memory region, ensure an
   ``MPI_Barrier(node_comm)`` follows the broadcast so all node-local ranks
   see the updated value before transport resumes (see `Synchronisation`_).

Mode 2: Reduce (estimator variables)
======================================

Reduce is used to combine the ``est`` (estimator) fields that each rank has
accumulated independently during photon transport.  Because every rank
transports a different subset of photons, its estimators represent only a
partial contribution to the total radiation field.  After transport, all ranks
call ``MPI_Allreduce`` to sum their partial values into a single consistent
result that is then available on every rank.

The reduce step is simpler to implement than the broadcast because the
operation is a straightforward element-wise sum rather than a
rank-by-rank pack/unpack cycle.  The key functions are:

- ``reduce_simple_estimators()`` in ``communicate_plasma.c`` — sums plasma
  ``est`` fields (``j``, ``ave_freq``, ``heat_tot``, etc.) across all ranks
  using ``MPI_Allreduce``.
- ``reduce_macro_atom_estimators()`` in ``communicate_macro.c`` — averages
  the macro-atom ``est`` fields (``jbar``, ``gamma``, ``matom_abs``, etc.)
  across all ranks, each contiguous block in place (see
  `In-place reduction of the macro-atom estimators`_).
- ``communicate_spectra.c`` — handles the output spectra, which are also
  accumulated independently per rank during spectral cycles.  Each rank
  records photon contributions to the extracted spectra (``xxspec``) for
  its own subset of photons; after transport the partial spectra are summed
  across all ranks so that every rank holds the complete spectrum.  The
  same reduce logic applies: spectra are private during transport and
  combined afterwards.

Adding a variable to a reduction
----------------------------------

1. Confirm the variable belongs to an ``est`` sub-structure (it is written
   during photon transport and must be summed across ranks).
2. In the appropriate reduce function, add the variable to the existing
   ``MPI_Allreduce`` call, or add a new ``MPI_Allreduce`` call following the
   existing pattern.
3. Ensure the variable is zeroed (reset) at the start of each transport cycle
   so partial sums from the previous cycle do not accumulate.
4. Because estimator arrays are always **private** (one copy per rank), no
   shared-memory or barrier changes are needed.

MPI-3 shared memory implementation
===================================

When running with more than one MPI rank, SIROCCO uses MPI-3 shared memory
windows to reduce per-node memory consumption.  The key idea is that ranks on
the same physical node share a single copy of data that is read-only during
photon transport, rather than duplicating it across every rank.

During MPI initialisation (in ``sirocco.c``), a *node-local communicator* is
created with ``MPI_Comm_split_type(MPI_COMM_TYPE_SHARED, ...)``.  Three global
variables track the node topology:

- ``node_comm`` — communicator for ranks sharing the same node
- ``node_rank`` — rank index within the node (0 = node leader)
- ``node_size`` — number of ranks on the node

If SIROCCO is run with the ``-no-shared-memory`` switch, each rank is instead
placed in a node communicator of its own, so nothing is shared and
every rank holds a private copy of all grids (see :doc:`shared_memory_safety`).

Contiguous block allocation
---------------------------

All variable-length plasma arrays (density, partition, ioniz, etc.) are
allocated as contiguous blocks in ``calloc_dyn_plasma()`` (in
``gridwind.c``), with each cell's pointer set to the appropriate offset
within the block.  This replaces the earlier pattern of separate ``calloc``
calls per cell and is a prerequisite for shared memory, since
``MPI_Win_allocate_shared`` requires contiguous regions.

The allocation is performed by two helper functions, ``alloc_block_double()``
and ``alloc_block_int()``, which accept a ``use_shared`` flag:

- When ``use_shared`` is TRUE and ``np_mpi_global > 1``, only the node leader
  (``node_rank == 0``) allocates memory via ``MPI_Win_allocate_shared``; other
  ranks on the same node obtain a pointer to the same physical memory via
  ``MPI_Win_shared_query``.
- When ``use_shared`` is FALSE, each rank allocates its own private block with
  regular ``calloc``.

The same contiguous block layout is used in non-MPI builds and with a single
MPI rank; the only difference is that ``calloc`` is used unconditionally.

Which arrays are shared
-----------------------

The allocation strategy mirrors the three sub-structures:

.. list-table::
   :header-rows: 1
   :widths: 15 55 15 15

   * - Category
     - Arrays
     - Allocation
     - Reason
   * - **State (dynamic)**
     - ``density``, ``partition``, ``levden``, ``recomb_simple``, ``recomb_simple_upweight``, ``kbf_use``
     - Shared
     - Read-only during photon transport
   * - **State (fixed-size)**
     - ``f1``, ``f2``, ``spec_mod_type``, ``pl_alpha``, ``pl_log_w``, ``exp_temp``, ``exp_w``, ``fmin_mod``, ``fmax_mod``
     - Shared
     - Spectral model parameters, read-only during transport
   * - **Estimators**
     - ``ioniz``, ``heat_ion``, ``heat_inner_ion``, ``inner_ioniz``
     - Private
     - Each rank accumulates independently
   * - **Derived (dynamic)**
     - ``recomb``, ``cool_rr_ion``, ``lum_rr_ion``, ``cool_dr_ion``, ``inner_recomb``
     - Shared
     - Computed during wind updates, then broadcast
   * - **Derived (fixed-size)**
     - ``F_vis_persistent``, ``F_UV_persistent``, ``F_Xray_persistent``, ``rad_force_es_persist``, ``rad_force_ff_persist``, ``rad_force_bf_persist``, ``F_UV_ang_theta_persist``, ``F_UV_ang_phi_persist``, ``F_UV_ang_r_persist``
     - Shared
     - Persistent radiation field averages, read-only during transport
   * - **Derived (exceptions)**
     - ``scatters``, ``xscatters``, ``n_bf_in``, ``n_bf_out``
     - Private
     - Incremented during photon transport (would race in shared memory).  ``n_bf_in``/``n_bf_out`` are dynamically sized to ``nphot_total`` (formerly fixed at ``N_PHOT_PROC=500``).

The same shared/private split applies to macro-atom dynamic arrays in
``calloc_estimators()`` and ``calloc_matom_matrix()`` (both in
``gridwind.c``).  State and derived arrays (``jbar_old``, ``gamma_old``,
``matom_emiss``, and the transition probability matrix ``matom_matrix``) are
shared, while estimator arrays (``jbar``, ``gamma``, ``matom_abs``, etc.)
are private.

The spontaneous recombination rates (``recomb_sp``, ``recomb_sp_e``) and the
k-packet cooling rates (``cooling_bf``, ``cooling_bf_col``, ``cooling_bb``)
are also shared, in ``macro_derived``.  They are not estimators: they depend
only on the plasma state, and are calculated once for every cell at the start
of each ionization cycle and before the spectral cycles, by
``calc_macro_recomb_and_kpkt_rates()`` in ``wind_updates2d.c``, with the cells
divided among the ranks.  Each rank writes its own cells into the node's
copy.  With one node a barrier then completes the arrays; with several nodes
the copies are zeroed first (``prepare_macro_rate_blocks()``) and the node
leaders sum them over ``leader_comm`` (``share_macro_recomb_and_kpkt_rates()``,
in ``communicate_macro.c``).  The cooling totals, which are scalars in each
rank's ``macromain``, are summed over all ranks in the same way.  Before
#1218 every rank filled these arrays itself, privately, for each cell in
which it created a k-packet.

The ``matom_matrix`` (an *nrows × nrows* transition probability matrix per
cell, where *nrows = nlevels_macro + 1*) is allocated as a single contiguous
shared block in ``calloc_matom_matrix()``.  The flat data
(``NPLASMA × nrows × nrows`` doubles) lives in
``macro_block_ptrs.matom_matrix_block`` (shared), while a private per-rank
array of row-pointers (``matom_matrix_rowptrs``) points into the shared block
to preserve the ``double **`` interface used throughout the code.  The matrix
is computed during wind updates — each rank fills its own cell slice — then
broadcast via ``broadcast_macro_atom_state_matrix()`` so all nodes obtain a
complete copy.  The ``MPI_Barrier(node_comm)`` at the end of that function
ensures node-local ranks see the written data before transport begins.
Because the matrix is strictly read-only during photon transport, no further
synchronisation is required.

Block pointer management
------------------------

Base pointers for all contiguous blocks are stored in global structs
``plasma_block_ptrs`` (type ``plasma_blocks``) and ``macro_block_ptrs``
(type ``macro_blocks``), declared in ``sirocco.h``.  These structs also hold
the ``MPI_Win`` handles needed to free shared windows and a
``shared_memory_active`` flag that records whether the current allocation
used shared memory.

Synchronisation
---------------

After any broadcast that writes to shared dynamic arrays, an
``MPI_Barrier(node_comm)`` ensures all node-local ranks see the new data
before proceeding.  These barriers appear at the end of:

- ``broadcast_wind_grid()``
- ``broadcast_updated_plasma_properties()``
- ``broadcast_plasma_grid()``
- ``broadcast_wind_luminosity()``
- ``broadcast_wind_cooling()``
- ``broadcast_updated_macro_atom_properties()``
- ``broadcast_macro_atom_emissivities()``
- ``reduce_macro_atom_estimators()``

Barriers are also needed *before* filling a shared array that every rank has
just reset.  ``get_matom_f()`` and ``get_matom_f_accelerate()``
(``macro_gen_f.c``) zero ``matom_emiss`` for all cells on every rank and then
fill only their own cells, so an ``MPI_Barrier(node_comm)`` separates the two
loops.  See :doc:`shared_memory_safety` for this and other patterns to avoid.

A barrier is also placed in ``calloc_wind()`` (``gridwind.c``) immediately
after the node leader's ``memset`` that zero-initialises the shared ``wmain``
block, ensuring the zeroed memory is visible to all node-local ranks before
any rank begins writing wind-cell fields.

Two additional barriers appear in ``create_wind_grid()``
(``define_wind.c``):

- Before ``make_coordinate_grid()`` — ensures all ranks have completed the
  serial ``wmain`` field initialisation loop (which writes ``inwind =
  W_NOT_ASSIGNED``) before any rank enters ``make_coordinate_grid()``, which
  for imported models overwrites ``inwind`` with values from the import file.
  Without this, a fast rank's import writes can be overwritten by a slow
  rank's init-loop writes, leaving cells with ``inwind = W_NOT_ASSIGNED``.
- After ``wind_complete()`` — ensures all ranks have finished
  ``make_coordinate_grid()`` and ``wind_complete()`` before any rank enters
  the parallel volume/velocity loop.  Without this, a fast rank can read a
  cell's ``inwind`` value before a slow rank has finished writing it from the
  coordinate grid setup.

During photon transport, state arrays are read-only so no synchronisation is
required.  The ``sobolev()`` function in ``resonate.c`` previously modified
``state.density`` temporarily during transport; it now passes a density
override to ``two_level_atom()`` instead, avoiding a race condition on shared
memory.

Cleanup
-------

At program exit, ``free_plasma_grid()`` and ``free_macro_grid()`` in
``janitor.c`` free the contiguous blocks.  For shared blocks the memory is
owned by the MPI window, so the pointer is simply NULLed (the MPI runtime
frees it at ``MPI_Finalize``).  Private blocks are freed with ``free()`` as
usual.

Memory savings
--------------

For a model with *N* plasma cells, *I* ions, and *R* ranks on one node, the
dominant dynamic arrays total roughly ``N * I * 14 * 8`` bytes per rank.
With shared memory the state and derived arrays exist only once per node,
reducing the per-node footprint by approximately ``(R-1)/R`` of the shared
portion.  Estimator arrays remain duplicated across ranks.

In addition to the variable-length dynamic arrays, fixed-size arrays that
were previously embedded in the ``plasma_state`` and ``plasma_derived``
sub-structures (spectral model parameters, persistent flux averages) have
been moved to shared contiguous blocks.  These arrays are declared as
pointers in the struct and point into combined blocks allocated in
``calloc_dyn_plasma()``.  This reduces ``sizeof(plasma_dummy)`` by
approximately 2.3 KB per cell, yielding additional PSS savings of roughly
``2.3 * N * (R-1)/R`` KB.  The savings scale linearly with NPLASMA: for a
model with 80K cells and 29 ranks, this adds approximately 177 MB of
per-rank savings.

The transition probability matrix ``matom_matrix`` (``nrows × nrows`` doubles
per cell, allocated by ``calloc_matom_matrix()``) is also placed in shared
memory.  For the ``h20_hetop_standard80`` atomic dataset (85 macro-atom
levels, *nrows* = 86) and a 300×300 grid with ~12,000 active plasma cells,
this matrix totals approximately 726 MB.  Without shared memory each of the
*R* ranks holds its own copy; with shared memory there is one copy per node.
On a 24-rank single-node run this saves roughly ``726 × 23 ≈ 16.7 GB`` of
physical memory, making it the single largest shared-memory saving in the
code.

``matom_matrix`` allocation and ``mc_jumps`` mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``calloc_matom_matrix()`` in ``gridwind.c`` allocates the shared
``NPLASMA × nrows × nrows`` block **only** when
``modes.store_matom_matrix == TRUE`` (i.e. the ``matrix`` macro-atom
transition mode is selected).  In ``mc_jumps`` mode
``modes.store_matom_matrix`` is ``FALSE``, so the allocation is skipped
entirely and the block remains ``NULL``.  This avoids reserving a shared
memory window of roughly ``NPLASMA × (nlevels_macro+1)^2 × 8`` bytes —
approximately 300–700 MB per node for a typical macro-atom dataset and a large
grid — in runs that never use the matrix.

The check is performed after freeing any previously allocated block, so a
re-run that switches from ``matrix`` to ``mc_jumps`` mode correctly frees the
old window before returning.

Shared wind structure
---------------------

The wind geometry array ``wmain`` (type ``wind_dummy``, indexed by NDIM2) is
allocated via ``MPI_Win_allocate_shared`` in ``calloc_wind()`` so that all
ranks on the same node share a single copy.  This is safe because ``wmain``
is populated during initialization and is strictly read-only during photon
transport.

The reverb path-tracking data (``paths`` and ``line_paths``) was moved out of
``wind_dummy`` into a separate per-rank array ``wind_paths_main`` (type
``wind_paths_store``), because path histograms are accumulated during photon
transport and must remain private per rank.  Code in ``paths.c`` accesses
these via ``wind_paths_main[cell_index]`` instead of
``wmain[cell_index]``.

For a 300x300 grid (NDIM2 = 90,000, ``sizeof(wind_dummy)`` = 288 bytes),
this saves approximately ``90000 * 288 * (R-1)/R`` bytes, or about 25 MB per
rank with 29 ranks.

Single-node optimisation for matom_matrix broadcast
----------------------------------------------------

``broadcast_macro_atom_state_matrix()`` in ``communicate_macro.c`` normally
packs the full transition-probability matrix for each rank's cell range into a
comm buffer and broadcasts it to all other ranks.  When ``matom_matrix`` lives
in shared memory (the normal MPI build) and all ranks are on the same node
(``num_nodes == 1``), this broadcast is unnecessary: the writing rank's data
is already visible to all node-local ranks through shared memory.  The
function therefore returns early with a ``MPI_Barrier(node_comm)`` to ensure
coherence, skipping the pack/Bcast/unpack cycle entirely.  On a single-node
run this avoids allocating the comm buffer (~100 KB) and removes latency
proportional to the number of ranks.

In-place reduction of the macro-atom estimators
-----------------------------------------------

``reduce_macro_atom_estimators()`` in ``communicate_macro.c`` averages the
macro-atom estimators over the ranks.  Each estimator is held in one
contiguous block of ``NPLASMA × size`` doubles (``macro_block_ptrs``), so the
function divides each block by the number of ranks and sums it with
``MPI_Allreduce(MPI_IN_PLACE, ...)``, without temporary copies.  The sum is
done in chunks of 2\ :sup:`27` doubles, so that the count passed to MPI fits
in an int.

Earlier versions packed the estimators into two temporary buffers each, and
also averaged the k-packet cooling rates, including ``cooling_bb``
(``NPLASMA × nlines`` doubles, reduced in ~50 MB chunks).  Those rates are now
calculated once per cycle in shared memory (see above), and are not reduced
(#1218).

Platform differences: macOS vs Linux
=====================================

On both macOS and Linux, ``wmain`` and the plasma and macro-atom arrays
listed above are allocated with ``MPI_Win_allocate_shared``, so the same code
paths and the same races are possible on both.  (An earlier version of the
development branch gave each rank a private ``wmain`` on macOS, guarded by
``#ifdef __APPLE__``; that fallback has been removed.)

The platforms still differ in ways that matter for testing:

- **Memory ordering.**  Apple Silicon (ARM) orders memory more weakly than
  x86.  A missing or misplaced barrier is *more* likely to show up as
  run-to-run differences on an ARM Mac than on an x86 Linux machine, so ARM
  is the stricter test of synchronisation.
- **Temporary directory.**  Open MPI places its session files (including
  shared-memory backing files and sockets) under ``$TMPDIR``.  The default
  macOS ``$TMPDIR`` path can be too long for this; if ``mpirun`` fails at
  start-up, set ``export TMPDIR=/tmp``.
- **Permission faults.**  Earlier builds on macOS with Open MPI 5 saw a
  ``SEGV_ACCERR`` (signal code 2, "address has wrong permissions") on the
  first write to a shared window from a rank other than the node leader.  If
  this reappears, report it together with the Open MPI version.

**Rule of thumb:** whenever you add, remove, or reorder a ``MPI_Barrier``,
``MPI_Win_allocate_shared``, ``MPI_Win_shared_query``, or a write to a shared
block, run the tests described in :doc:`shared_memory_safety` on x86 Linux
*and* on an ARM machine (for example an Apple Silicon Mac) before merging.

3D coordinate types (CYLIND3D and SPH3D)
=========================================

The 3D coordinate types add a third azimuthal (phi) grid dimension, increasing
NDIM2 accordingly, but they require no special-casing in the MPI or shared
memory layer.  The same shared/private split, broadcast pattern, and
``MPI_Win_allocate_shared`` calls apply unchanged.  The larger NDIM2 simply
means more cells in the shared ``wmain`` and plasma arrays, with
proportionally greater memory savings from the mixed model.
