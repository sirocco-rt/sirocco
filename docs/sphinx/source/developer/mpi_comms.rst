MPI Communication
#################

SIROCCO is parallelised using the Message Passing Interface (MPI). This page contains information on how data is shared
between ranks and should serve as a basic set of instructions for extending or modifying the data communication
routines.

In general, all calls to MPI are isolated from the rest of SIROCCO. Most, if not all, of the MPI code is contained
within five source files, which deal entirely with parallelisation or communication. Currently these files are:

- :code:`communicate_macro.c`
- :code:`communicate_plasma.c`
- :code:`communicate_spectra.c`
- :code:`communicate_wind.c`
- :code:`para_update.c`

Given the names of the files, it should be obvious what sort of code is contained in them. If you need to extend or
implement a new function for MPI, please place it either in one of the above files or create a new file using an
appropriately similar name. Any parallel code should be wrapped by :code:`#ifdef MPI_ON` and :code:`#endif` as shown in
the code example below:

.. code:: c

    void communication_function(void)
    {
    #ifdef MPI_ON
        /* MPI communication could should go between the #ifdef's here */
    #endif
    }

Don't forget to update the Makefile and :code:`templates.h` if you add a new file or function.

Communication pattern: broadcasting data to all ranks
=====================================================

By far the most typical communication pattern in SIROCCO (and, I think, the only pattern) is to broadcast data from one
rank to all other ranks. This is done, for example, to update and synchronise the plasma or macro atom grids in each
rank. As the data structures in SIROCCO are fairly complex and use pointers/dynamic memory allocation, we as forced to
manually pack and unpack a contiguous communication buffer which results in a fairly manual (and error prone?) process
for communicating data.

Calculating the size of the communication buffer
------------------------------------------------

The size of the communication buffer has to be calculated manually, by counting the number of variables being copied
into it and converting this to the appropriate number of bytes. This is done by the :code:`calculate_comm_buffer_size`
function which takes two arguments: 1) the number of :code:`int`'s and 2) the number of :code:`double`'s. We have to
_manually_ count the number of :code:`int` and :code:`double` variables being communicated. Due to the manual nature of
this, greate care has to be taken to ensure the correct number are counted otherwise MPI will cause crash during
communication.

When counting variables, one needs to count the number if _single_ variables of a certain type as well as the number of
elements in an array of that same type. Consider the example below,

.. code:: c

    int my_int;
    int *my_int_arr = malloc(10 * sizeof(int));
    int num_ints = 11;

In this case there are 11 integer variables which will want to be communicated. In practise, calculating the communication
buffer is usually done as in the code example below:

.. code:: c

    /* We need to ensure the buffer is large enough, as soon ranks may be sending a smaller
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

Communication implementation
----------------------------

The general pattern for packing data into a communication buffer and then sharing it between ranks is as follows,

- Loop over all the MPI ranks (in MPI_COMM_WORLD.
- If the loop variable is equal to a rank's ID, that rank will broadcast it's subset of data to the other ranks. This
  rank uses :code:`MPI_Pack` to copy its data into the communication buffer.
- All ranks call :code:`MPI_Bcast`, which sends data from the root rank (this is the rank which has just put its data
  into the communication buffer) and receives it into all non-root ranks.
- Non-root ranks use :code:`MPI_Unpack` to copy data from the communication buffer into the appropriate location.
- This is repeated until all MPI ranks have sent their data root, and have therefore received data from all other ranks.

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

This is likely the most best method to communicate data in SIROCCO, given the complexity of the data structures.
Unfortunately there are not many structures or situations where using a derived data type, to simplify code, is viable
due to none of the structures being contiguous in memory.

Adding a new variable to an existing communication
--------------------------------------------------

- Increment the appropriate variable, or function call to :code:`calculate_comm_buffer_size`, to account for and
  allocate additional space in the communication buffer. For example, if the new variable is an :code:`int` in the
  plasma grid then update :code:`n_cells_max * (20 + 2 * n_phot_total + 1)` to :code:`n_cells_max * (21 + 2 *
  n_phot_total + 1)`
- In the block where :code:`rank == rank_global`, add a new call to :code:`MPI_Pack` using the code which is already
  there as an example.
- In the block where :code:`rank != rank_global`, add a new call to :code:`MPI_Unpack` using the code which is already
  there as an example.

Relationship between sub-structures and communication patterns
==============================================================

The ``plasma_dummy`` and ``macro_dummy`` structures are each divided into three
sub-structures that correspond directly to different MPI communication patterns.
This makes it straightforward to determine which communication function to modify
when adding a new variable:

.. list-table:: Plasma sub-structures and their communication
   :header-rows: 1
   :widths: 20 30 25 25

   * - Sub-structure
     - Contents
     - Communication
     - Key functions
   * - ``state``
     - Thermodynamic state (``ne``, ``t_e``, ``t_r``, ``w``, ``rho``, ``vol``), ion populations (``density``, ``partition``, ``levden``), spectral model parameters, bound-free data
     - Broadcast after wind updates
     - ``broadcast_updated_plasma_properties()``, ``broadcast_plasma_grid()``
   * - ``est``
     - Radiation field estimators (``j``, ``ave_freq``), heating rates (``heat_tot``, ``heat_lines``, etc.), photon counters, flux estimators, cell spectra, ionization estimators
     - Reduced (summed) across ranks after photon transport
     - ``reduce_simple_estimators()``
   * - ``derived``
     - Cooling rates, luminosities, convergence diagnostics, scatter counts, persistent flux averages, ionization parameter (``xi``)
     - Broadcast after wind updates
     - ``broadcast_updated_plasma_properties()``, ``broadcast_wind_luminosity()``, ``broadcast_wind_cooling()``

.. list-table:: Macro-atom sub-structures and their communication
   :header-rows: 1
   :widths: 20 30 25 25

   * - Sub-structure
     - Contents
     - Communication
     - Key functions
   * - ``state``
     - Normalized rate coefficients (``jbar_old``, ``gamma_old``, ``alpha_st_old``, etc.), mode flags
     - Broadcast after wind updates
     - ``broadcast_updated_macro_atom_properties()``
   * - ``est``
     - Raw Sobolev mean intensities (``jbar``), photoionization rates (``gamma``), stimulated recombination rates (``alpha_st``), macro-atom absorption, cooling stores
     - Reduced (summed) across ranks after transport
     - ``reduce_macro_atom_estimators()``
   * - ``derived``
     - Macro-atom emissivities (``matom_emiss``), k-packet rate flags, transition probability matrix
     - Broadcast after computation
     - ``broadcast_macro_atom_emissivities()``

When adding a new variable, place it in the appropriate sub-structure and update
the corresponding communication function.  For ``est`` fields, update the reduction
function.  For ``state`` or ``derived`` fields, update the broadcast function.
In both cases, remember to update the buffer size calculation (the integer and double
counts) to account for the new variable.

MPI-3 shared memory model
-------------------------

When running with more than one MPI rank, SIROCCO uses MPI-3 shared memory windows
to reduce per-node memory consumption.  The key idea is that ranks on the same
physical node share a single copy of data that is read-only during photon transport,
rather than duplicating it across every rank.

During MPI initialisation (in ``sirocco.c``), a *node-local communicator*
is created with ``MPI_Comm_split_type(MPI_COMM_TYPE_SHARED, ...)``.  Three global
variables track the node topology:

- ``node_comm`` — communicator for ranks sharing the same node
- ``node_rank`` — rank index within the node (0 = node leader)
- ``node_size`` — number of ranks on the node

Contiguous block allocation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

All variable-length plasma arrays (density, partition, ioniz, etc.) are allocated
as contiguous blocks in ``calloc_dyn_plasma()`` (in ``gridwind.c``), with each
cell's pointer set to the appropriate offset within the block.  This replaces the
earlier pattern of separate ``calloc`` calls per cell and is a prerequisite for
shared memory, since ``MPI_Win_allocate_shared`` requires contiguous regions.

The allocation is performed by two helper functions, ``alloc_block_double()`` and
``alloc_block_int()``, which accept a ``use_shared`` flag:

- When ``use_shared`` is TRUE and ``np_mpi_global > 1``, only the node leader
  (``node_rank == 0``) allocates memory via ``MPI_Win_allocate_shared``; other
  ranks on the same node obtain a pointer to the same physical memory via
  ``MPI_Win_shared_query``.
- When ``use_shared`` is FALSE, each rank allocates its own private block with
  regular ``calloc``.

The same contiguous block layout is used in non-MPI builds and with a single MPI
rank; the only difference is that ``calloc`` is used unconditionally.

Which arrays are shared
^^^^^^^^^^^^^^^^^^^^^^^

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
     - ``scatters``, ``xscatters``
     - Private
     - Incremented during photon transport (would race in shared memory)

The same shared/private split applies to macro-atom dynamic arrays in
``calloc_estimators()`` (also in ``gridwind.c``).  State and derived arrays
(``jbar_old``, ``gamma_old``, ``matom_emiss``, etc.) are shared, while
estimator arrays (``jbar``, ``gamma``, ``cooling_bf``, etc.) are private.

Block pointer management
^^^^^^^^^^^^^^^^^^^^^^^^

Base pointers for all contiguous blocks are stored in global structs
``plasma_block_ptrs`` (type ``plasma_blocks``) and ``macro_block_ptrs``
(type ``macro_blocks``), declared in ``sirocco.h``.  These structs also hold
the ``MPI_Win`` handles needed to free shared windows and a
``shared_memory_active`` flag that records whether the current allocation
used shared memory.

Synchronisation
^^^^^^^^^^^^^^^

After any broadcast that writes to shared dynamic arrays, an
``MPI_Barrier(node_comm)`` ensures all node-local ranks see the new data
before proceeding.  These barriers appear at the end of:

- ``broadcast_updated_plasma_properties()``
- ``broadcast_plasma_grid()``
- ``broadcast_wind_luminosity()``
- ``broadcast_wind_cooling()``
- ``broadcast_updated_macro_atom_properties()``
- ``broadcast_macro_atom_emissivities()``
- ``reduce_macro_atom_estimators()``

During photon transport, state arrays are read-only so no synchronisation
is required.  The ``sobolev()`` function in ``resonate.c`` previously
modified ``state.density`` temporarily during transport; it now passes a
density override to ``two_level_atom()`` instead, avoiding a race condition
on shared memory.

Cleanup
^^^^^^^

At program exit, ``free_plasma_grid()`` and ``free_macro_grid()`` in
``janitor.c`` free the contiguous blocks.  For shared blocks the memory is
owned by the MPI window, so the pointer is simply NULLed (the MPI runtime
frees it at ``MPI_Finalize``).  Private blocks are freed with ``free()``
as usual.

Memory savings
^^^^^^^^^^^^^^

For a model with *N* plasma cells, *I* ions, and *R* ranks on one node,
the dominant dynamic arrays total roughly ``N * I * 14 * 8`` bytes per rank.
With shared memory the state and derived arrays exist only once per node,
reducing the per-node footprint by approximately ``(R-1)/R`` of the shared
portion.  Estimator arrays remain duplicated across ranks.

In addition to the variable-length dynamic arrays, fixed-size arrays that
were previously embedded in the ``plasma_state`` and ``plasma_derived``
sub-structures (spectral model parameters, persistent flux averages) have
been moved to shared contiguous blocks.  These arrays are declared as
pointers in the struct and point into combined blocks allocated in
``calloc_dyn_plasma()``.  This reduces ``sizeof(plasma_dummy)`` by
approximately 2.3 KB per cell, yielding additional PSS savings of
roughly ``2.3 * N * (R-1)/R`` KB.  The savings scale linearly with
NPLASMA: for a model with 80K cells and 29 ranks, this adds approximately
177 MB of per-rank savings.
