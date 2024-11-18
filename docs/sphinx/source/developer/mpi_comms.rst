MPI Communication
#################

SIROCCO is parallelised using the Message Passing Interface (MPI). This page contains information on how data is shared
between ranks and should serve as a basic set of instructions for extending or modifying the data communication
routines.

In general, all calls to MPI are isolated from the rest of SIROCCO. Most, if not all, of the MPI code is contained
within give source files, which deal entirely with parallelisation or communication. Currently these files are:

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

In this case there are 11 :code:`int`s which will want to be communicated. In practise, calculating the communication
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
