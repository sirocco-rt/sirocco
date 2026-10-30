/***********************************************************/
/** @file  communicate_macro.c
 * @author EJP
 * @date   December 2023
 *
 * @brief Functions for communicating macro atom properties
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

/**********************************************************/
/**
 * @brief  Communicate the macro atom emissivities
 *
 * @param [in] int n_start       The index of the first cell updated by this rank
 * @param [in] int n_stop        The index of the last cell updated by this rank
 * @param [in] int n_cells_rank  The number of cells this rank updated
 *
 * @details
 *
 * The communication pattern and how the size of the communication buffer is
 * determined is documented in
 * $SIROCCO/docs/sphinx/source/developer/mpi_comms.rst.
 * To communicate a new variable, the communication buffer needs to be made
 * bigger and a new `MPI_Pack` and `MPI_Unpack` call need to be added. See the
 * developer documentation for more details.
 *
 **********************************************************/

void
broadcast_macro_atom_emissivities (const int n_start, const int n_stop, const int n_cells_rank)
{
#ifdef MPI_ON
  int i;
  int n_plasma;
  int current_rank;
  int num_comm;
  int position;

  d_xsignal (files.root, "%-20s Begin macro atom emissivity communication\n", "NOK");
  const int n_cells_max = get_max_cells_per_rank (NPLASMA);
  const int comm_buffer_size = calculate_comm_buffer_size (1 + n_cells_max, n_cells_max * (1 + nlevels_macro));

  char *comm_buffer = malloc (comm_buffer_size);
  if (comm_buffer == NULL)
  {
    Error ("broadcast_macro_atom_emissivities: Error in allocating memory for comm_buffer\n");
    Exit (EXIT_FAILURE);
  }

  for (current_rank = 0; current_rank < np_mpi_global; current_rank++)
  {
    position = 0;

    if (rank_global == current_rank)
    {
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
      for (n_plasma = n_start; n_plasma < n_stop; ++n_plasma)
      {
        MPI_Pack (&n_plasma, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.kpkt_emiss, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].derived.matom_emiss, nlevels_macro, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position,
                  MPI_COMM_WORLD);
      }
    }

    MPI_Bcast (comm_buffer, comm_buffer_size, MPI_PACKED, current_rank, MPI_COMM_WORLD);

    position = 0;

    if (rank_global != current_rank)
    {
      MPI_Unpack (comm_buffer, comm_buffer_size, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
      for (i = 0; i < num_comm; i++)
      {
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &n_plasma, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n_plasma].derived.kpkt_emiss, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].derived.matom_emiss, nlevels_macro, MPI_DOUBLE,
                    MPI_COMM_WORLD);
      }
    }
  }

  free (comm_buffer);

  /* Barrier to ensure shared memory writes are visible to all node-local ranks */
  MPI_Barrier (node_comm);

  d_xsignal (files.root, "%-20s Finished macro atom emissivity communication\n", "OK");
#endif
}

/**********************************************************/
/**
 * @brief  Communicate the recombination rates for simple ions between ranks
 *
 * @param [in] int n_start       The index of the first cell updated by this rank
 * @param [in] int n_stop        The index of the last cell updated by this rank
 * @param [in] int n_cells_rank  The number of cells this rank updated
 *
 * @details
 *
 * The communication pattern and how the size of the communication buffer is
 * determined is documented in
 * $SIROCCO/docs/sphinx/source/developer/mpi_comms.rst.
 * To communicate a new variable, the communication buffer needs to be made
 * bigger and a new `MPI_Pack` and `MPI_Unpack` call need to be added. See the
 * developer documentation for more details.
 *
 **********************************************************/

void
broadcast_macro_atom_recomb (const int n_start, const int n_stop, const int n_cells_rank)
{
#ifdef MPI_ON
  int i;
  int n_plasma;
  int current_rank;
  int position;
  int num_comm;

  d_xsignal (files.root, "%-20s Begin macro atom recombination communication\n", "NOK");
  const int n_cells_max = get_max_cells_per_rank (NPLASMA);
  const int comm_buffer_size = calculate_comm_buffer_size (1 + n_cells_max, n_cells_max * (2 * nphot_total));

  char *const comm_buffer = malloc (comm_buffer_size);
  if (comm_buffer == NULL)
  {
    Error ("broadcast_macro_atom_recomb: Error in allocating memory for comm_buffer\n");
    Exit (EXIT_FAILURE);
  }

  for (current_rank = 0; current_rank < np_mpi_global; ++current_rank)
  {
    position = 0;

    if (rank_global == current_rank)
    {
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);   // how many cells to unpack
      for (n_plasma = n_start; n_plasma < n_stop; ++n_plasma)
      {
        MPI_Pack (&n_plasma, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);     // which cell we're working on
        if (nphot_total > 0)
        {
          MPI_Pack (plasmamain[n_plasma].state.recomb_simple, nphot_total, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position,
                    MPI_COMM_WORLD);
          MPI_Pack (plasmamain[n_plasma].state.recomb_simple_upweight, nphot_total, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position,
                    MPI_COMM_WORLD);
        }
      }
    }

    MPI_Bcast (comm_buffer, comm_buffer_size, MPI_PACKED, current_rank, MPI_COMM_WORLD);

    position = 0;

    if (rank_global != current_rank)
    {
      MPI_Unpack (comm_buffer, comm_buffer_size, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
      for (i = 0; i < num_comm; ++i)
      {
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &n_plasma, 1, MPI_INT, MPI_COMM_WORLD);

        if (nphot_total > 0)
        {
          MPI_Unpack (comm_buffer, comm_buffer_size, &position, plasmamain[n_plasma].state.recomb_simple,
                      nphot_total, MPI_DOUBLE, MPI_COMM_WORLD);
          MPI_Unpack (comm_buffer, comm_buffer_size, &position, plasmamain[n_plasma].state.recomb_simple_upweight,
                      nphot_total, MPI_DOUBLE, MPI_COMM_WORLD);
        }
      }
    }
  }

  free (comm_buffer);

  /* Barrier to ensure shared memory writes are visible to all node-local ranks */
  MPI_Barrier (node_comm);

  d_xsignal (files.root, "%-20s Finished macro atom recombination communication\n", "OK");
#endif
}

/**********************************************************/
/**
 * @brief  Communicate the macro atom properties updated in `wind_update`
 *
 * @param [in] int n_start       The index of the first cell updated by this rank
 * @param [in] int n_stop        The index of the last cell updated by this rank
 * @param [in] int n_cells_rank  The number of cells this rank updated
 *
 * @details
 *
 * The communication pattern and how the size of the communication buffer is
 * determined is documented in
 * $SIROCCO/docs/sphinx/source/developer/mpi_comms.rst.
 * To communicate a new variable, the communication buffer needs to be made
 * bigger and a new `MPI_Pack` and `MPI_Unpack` call need to be added. See the
 * developer documentation for more details.
 *
 **********************************************************/

int
broadcast_updated_macro_atom_properties (const int n_start, const int n_stop, const int n_cells_rank)
{
#ifdef MPI_ON
  int i;
  int n_plasma;
  int current_rank;
  int position;
  int num_comm;

  d_xsignal (files.root, "%-20s Begin macro atom updated properties communication\n", "NOK");
  const int n_cells_max = get_max_cells_per_rank (NPLASMA);
  const int comm_buffer_size = calculate_comm_buffer_size (1 + 3 * n_cells_max, n_cells_max * (8 * size_gamma_est + 2 * size_Jbar_est + 1));

  char *const comm_buffer = malloc (comm_buffer_size);
  if (comm_buffer == NULL)
  {
    Error ("broadcast_updated_macro_atom_properties: Error in allocating memory for comm_buffer\n");
    Exit (EXIT_FAILURE);
  }

  for (current_rank = 0; current_rank < np_mpi_global; ++current_rank)
  {
    position = 0;

    if (rank_global == current_rank)
    {
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
      for (n_plasma = n_start; n_plasma < n_stop; ++n_plasma)
      {
        MPI_Pack (&n_plasma, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].est.jbar, size_Jbar_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].state.jbar_old, size_Jbar_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].est.gamma, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].state.gamma_old, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].est.gamma_e, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].state.gamma_e_old, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].est.alpha_st, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].state.alpha_st_old, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].est.alpha_st_e, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].state.alpha_st_e_old, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (&macromain[n_plasma].derived.kpkt_rates_known, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&macromain[n_plasma].derived.matrix_rates_known, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&macromain[n_plasma].est.energy_flow_in, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);

      }
    }

    MPI_Bcast (comm_buffer, comm_buffer_size, MPI_PACKED, current_rank, MPI_COMM_WORLD);

    position = 0;

    if (rank_global != current_rank)
    {
      MPI_Unpack (comm_buffer, comm_buffer_size, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
      for (i = 0; i < num_comm; ++i)
      {
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &n_plasma, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].est.jbar, size_Jbar_est, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].state.jbar_old, size_Jbar_est, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].est.gamma, size_gamma_est, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].state.gamma_old, size_gamma_est, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].est.gamma_e, size_gamma_est, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].state.gamma_e_old, size_gamma_est, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].est.alpha_st, size_gamma_est, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].state.alpha_st_old, size_gamma_est, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].est.alpha_st_e, size_gamma_est, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].state.alpha_st_e_old, size_gamma_est, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &macromain[n_plasma].derived.kpkt_rates_known, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &macromain[n_plasma].derived.matrix_rates_known, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &macromain[n_plasma].est.energy_flow_in, 1, MPI_DOUBLE, MPI_COMM_WORLD);
      }
    }
  }

  free (comm_buffer);

  /* Barrier to ensure shared memory writes are visible to all node-local ranks */
  MPI_Barrier (node_comm);

  d_xsignal (files.root, "%-20s Finished macro atom updated properties communication\n", "OK");
#endif
  return EXIT_SUCCESS;
}

/**********************************************************/
/**
 * @brief communicates the macro-atom B matrices between threads
 *
 *
 * @details communicates the macro-atom B matrices between threads
 * using MPI_Reduce. Only does anything if MPI_ON flag is on,
 * and should only be called if geo.rt_mode == RT_MODE_MACRO
 * and nlevels_macro > 0
 *
 * The communication pattern and how the size of the communication buffer is
 * determined is documented in
 * $SIROCCO/docs/sphinx/source/developer/mpi_comms.rst.
 * To communicate a new variable, the communication buffer needs to be made
 * bigger and a new `MPI_Pack` and `MPI_Unpack` call need to be added. See the
 * developer documentation for more details.
 *
 **********************************************************/

int
broadcast_macro_atom_state_matrix (int n_start, int n_stop, int n_cells_rank)
{
#ifdef MPI_ON
  int n_mpi, n_mpi2, num_comm;
  int n, position;

  d_xsignal (files.root, "%-20s Begin macro atom state matrix communication\n", "NOK");

  /* matom_matrix lives in shared memory on each node, so for a single-node run
   * the writing rank's data is already visible to all node-local ranks.
   * A barrier is sufficient; no cross-rank data transfer is needed. */
  if (num_nodes == 1)
  {
    MPI_Barrier (node_comm);
    d_xsignal (files.root, "%-20s Finished macro atom state matrix communication\n", "OK");
    return (0);
  }

  const int matrix_size = nlevels_macro + 1;
  const int n_cells_max = get_max_cells_per_rank (NPLASMA);
  const int comm_buffer_size = calculate_comm_buffer_size (1 + n_cells_max, n_cells_max * (matrix_size * matrix_size));

  char *comm_buffer = malloc (comm_buffer_size);
  if (comm_buffer == NULL)
  {
    Error ("broadcast_macro_atom_state_matrix: Error in allocating memory for comm_buffer\n");
    Exit (EXIT_FAILURE);
  }

  for (n_mpi = 0; n_mpi < np_mpi_global; n_mpi++)
  {
    position = 0;

    if (rank_global == n_mpi)
    {
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
      for (n = n_start; n < n_stop; n++)
      {
        MPI_Pack (&n, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);

        /* we only communicate the matrix if it is being stored in this cell */
        if (macromain[n].state.store_matom_matrix == TRUE)
        {
          MPI_Pack (macromain[n].derived.matom_matrix[0], matrix_size * matrix_size, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position,
                    MPI_COMM_WORLD);
        }
      }
    }

    MPI_Bcast (comm_buffer, comm_buffer_size, MPI_PACKED, n_mpi, MPI_COMM_WORLD);

    position = 0;

    if (rank_global != n_mpi)
    {
      MPI_Unpack (comm_buffer, comm_buffer_size, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
      for (n_mpi2 = 0; n_mpi2 < num_comm; n_mpi2++)
      {
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &n, 1, MPI_INT, MPI_COMM_WORLD);

        /* we only communicate the matrix if it is being stored in this cell */
        if (macromain[n].state.store_matom_matrix == TRUE)
        {
          MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n].derived.matom_matrix[0], matrix_size * matrix_size, MPI_DOUBLE,
                      MPI_COMM_WORLD);
        }
      }
    }
  }

  free (comm_buffer);

  /* Barrier to ensure shared memory writes are visible to all node-local ranks */
  MPI_Barrier (node_comm);

  d_xsignal (files.root, "%-20s Finished macro atom state matrix communication\n", "OK");
#endif
  return (0);
}

/**********************************************************/
/**
 * @brief Average one contiguous block of estimators over the ranks, in place
 *
 * @param [in,out] double *block  The block, which is replaced by the mean over the ranks
 * @param [in] long n  The number of doubles in the block
 *
 * @details Each value is divided by the number of ranks and the block is then summed
 * over the ranks with MPI_IN_PLACE, so no copy of the block is needed. The sum is done in
 * chunks so that the count passed to MPI fits in an int.
 *
 **********************************************************/

#ifdef MPI_ON
#define REDUCE_CHUNK (1L << 27)

static void
sum_block_in_place (double *block, long n, MPI_Comm comm)
{
  long start, count;

  for (start = 0; start < n; start += REDUCE_CHUNK)
  {
    count = n - start;
    if (count > REDUCE_CHUNK)
      count = REDUCE_CHUNK;
    MPI_Allreduce (MPI_IN_PLACE, block + start, (int) count, MPI_DOUBLE, MPI_SUM, comm);
  }
}

static void
reduce_block_mean_in_place (double *block, long n)
{
  long i;

  if (block == NULL || n <= 0)
    return;

  for (i = 0; i < n; i++)
  {
    block[i] /= np_mpi_global;
  }

  sum_block_in_place (block, n, MPI_COMM_WORLD);
}
#endif

/**********************************************************/
/**
 * @brief
 *
 * @details averages the macro-atom estimators between tasks using MPI_Allreduce.
 *   It should only be called if the MPI_ON flag was present
 *   in compilation, and returns 0 immediately if no macro atom levels.
 *
 *   The estimator arrays are held in contiguous blocks (see calloc_estimators in
 *   gridwind.c), and each block is averaged in place, so no temporary copies are made.
 *
 *   recomb_sp, recomb_sp_e and the k-packet cooling rates are not estimators. They are
 *   calculated from the plasma state at the start of each cycle, are the same on every
 *   rank, and are held in shared memory (see calc_macro_recomb_and_kpkt_rates and #1218).
 *
 **********************************************************/

void
reduce_macro_atom_estimators (void)
{
#ifdef MPI_ON                   // these routines should only be called anyway in parallel but we need these to compile

  int mpi_i;
  double *cell_helper;

  d_xsignal (files.root, "%-20s Begin reduction of macro atom estimators\n", "NOK");

  if (nlevels_macro == 0 && geo.nmacro == 0)
  {
    /* in this case no space would have been allocated for macro-atom estimators */
    Log ("No need to communicate matom estimators as no macro-atoms!\n");
    return;
  }

  /* the two estimators held in the cell structures themselves */
  cell_helper = calloc (sizeof (double), 2 * NPLASMA);
  if (cell_helper == NULL)
  {
    Error ("reduce_macro_atom_estimators: Error allocating cell_helper\n");
    Exit (EXIT_FAILURE);
  }

  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    cell_helper[mpi_i] = plasmamain[mpi_i].est.kpkt_abs;
    cell_helper[mpi_i + NPLASMA] = macromain[mpi_i].est.energy_flow_out;
  }

  reduce_block_mean_in_place (cell_helper, 2L * NPLASMA);

  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    plasmamain[mpi_i].est.kpkt_abs = cell_helper[mpi_i];
    macromain[mpi_i].est.energy_flow_out = cell_helper[mpi_i + NPLASMA];
  }

  free (cell_helper);

  /* the estimator arrays, one contiguous block each */
  reduce_block_mean_in_place (macro_block_ptrs.matom_abs_block, (long) NPLASMA * nlevels_macro);
  reduce_block_mean_in_place (macro_block_ptrs.jbar_block, (long) NPLASMA * size_Jbar_est);
  reduce_block_mean_in_place (macro_block_ptrs.alpha_st_block, (long) NPLASMA * size_gamma_est);
  reduce_block_mean_in_place (macro_block_ptrs.alpha_st_e_block, (long) NPLASMA * size_gamma_est);
  reduce_block_mean_in_place (macro_block_ptrs.gamma_block, (long) NPLASMA * size_gamma_est);
  reduce_block_mean_in_place (macro_block_ptrs.gamma_e_block, (long) NPLASMA * size_gamma_est);

  d_xsignal (files.root, "%-20s Finished reduction of macro atom estimators\n", "OK");
#endif
}

/**********************************************************/
/**
 * @brief Prepare the shared rate arrays before calc_macro_recomb_and_kpkt_rates fills them
 *
 * @details
 * The arrays recomb_sp, recomb_sp_e, cooling_bf, cooling_bf_col and cooling_bb are filled
 * by each rank for its own cells, and then completed by summing over the ranks or nodes
 * (see share_macro_recomb_and_kpkt_rates).  For that the cells calculated elsewhere must
 * be zero.  With shared memory on a single node nothing needs to be done, since every
 * cell is written by a rank on the node.  Otherwise the node leader (or, without shared
 * memory, each rank) sets the arrays to zero.
 *
 **********************************************************/

void
prepare_macro_rate_blocks (void)
{
#ifdef MPI_ON
  int shared = macro_block_ptrs.shared_memory_active;
  long n_alpha = (long) NPLASMA * size_alpha_est;
  long n_phot = (long) NPLASMA * nphot_total;
  long n_lines = (long) NPLASMA * nlines;

  if (np_mpi_global == 1 || macro_block_ptrs.recomb_sp_block == NULL)
    return;

  if (shared)
    MPI_Barrier (node_comm);    /* no rank on the node may still be using the arrays */

  if ((shared == FALSE) || (num_nodes > 1 && node_rank == 0))
  {
    memset (macro_block_ptrs.recomb_sp_block, 0, n_alpha * sizeof (double));
    memset (macro_block_ptrs.recomb_sp_e_block, 0, n_alpha * sizeof (double));
    memset (macro_block_ptrs.cooling_bf_block, 0, n_phot * sizeof (double));
    memset (macro_block_ptrs.cooling_bf_col_block, 0, n_phot * sizeof (double));
    memset (macro_block_ptrs.cooling_bb_block, 0, n_lines * sizeof (double));
  }

  if (shared)
    MPI_Barrier (node_comm);
#endif
}

/**********************************************************/
/**
 * @brief Make the spontaneous recombination and k-packet cooling rates known on all ranks
 *
 * @param [in] int n_start  The first cell calculated by this rank
 * @param [in] int n_stop   One more than the last cell calculated by this rank
 *
 * @details
 * calc_macro_recomb_and_kpkt_rates calculates recomb_sp, recomb_sp_e and the k-packet
 * cooling rates for cells n_start to n_stop - 1 on each rank.
 *
 * The arrays are in shared memory, and each rank has written its own cells into its
 * node's copy.  On a single node a barrier is therefore all that is needed.  With more
 * than one node, the copies were set to zero beforehand (prepare_macro_rate_blocks), so
 * the node leaders complete them by summing over leader_comm, in place.  Without shared
 * memory, every rank sums its private copy over MPI_COMM_WORLD.  As each value comes
 * from exactly one rank, the sums are exact.
 *
 * The cooling totals and kpkt_rates_known are held in each rank's macromain, and are
 * summed over all ranks in the same way, from a temporary array that is zero except for
 * the cells this rank calculated.
 *
 **********************************************************/

#define NKPKT_SCALARS 10

void
share_macro_recomb_and_kpkt_rates (const int n_start, const int n_stop)
{
#ifdef MPI_ON
  int n, shared;
  double *x;
  long n_alpha = (long) NPLASMA * size_alpha_est;
  long n_phot = (long) NPLASMA * nphot_total;
  long n_lines = (long) NPLASMA * nlines;

  if (np_mpi_global == 1 || macro_block_ptrs.recomb_sp_block == NULL)
    return;

  d_xsignal (files.root, "%-20s Begin sharing of macro atom recombination and k-packet rates\n", "NOK");

  shared = macro_block_ptrs.shared_memory_active;

  /* The arrays */
  if (shared)
    MPI_Barrier (node_comm);    /* every rank on the node has written its cells */

  if (shared == FALSE || (num_nodes > 1 && node_rank == 0))
  {
    MPI_Comm comm = shared ? leader_comm : MPI_COMM_WORLD;
    sum_block_in_place (macro_block_ptrs.recomb_sp_block, n_alpha, comm);
    sum_block_in_place (macro_block_ptrs.recomb_sp_e_block, n_alpha, comm);
    sum_block_in_place (macro_block_ptrs.cooling_bf_block, n_phot, comm);
    sum_block_in_place (macro_block_ptrs.cooling_bf_col_block, n_phot, comm);
    sum_block_in_place (macro_block_ptrs.cooling_bb_block, n_lines, comm);
  }

  if (shared)
    MPI_Barrier (node_comm);    /* the leaders have completed the node's copy */

  /* The quantities held in each rank's macromain */
  x = calloc (sizeof (double), (long) NKPKT_SCALARS * NPLASMA);
  if (x == NULL)
  {
    Error ("share_macro_recomb_and_kpkt_rates: Error allocating memory\n");
    Exit (EXIT_FAILURE);
  }

  for (n = n_start; n < n_stop; n++)
  {
    x[n * NKPKT_SCALARS] = macromain[n].derived.cooling_normalisation;
    x[n * NKPKT_SCALARS + 1] = macromain[n].derived.cooling_bbtot;
    x[n * NKPKT_SCALARS + 2] = macromain[n].derived.cooling_bftot;
    x[n * NKPKT_SCALARS + 3] = macromain[n].derived.cooling_bf_coltot;
    x[n * NKPKT_SCALARS + 4] = macromain[n].derived.cooling_ff;
    x[n * NKPKT_SCALARS + 5] = macromain[n].derived.cooling_ff_lofreq;
    x[n * NKPKT_SCALARS + 6] = macromain[n].derived.cooling_adiabatic;
    x[n * NKPKT_SCALARS + 7] = macromain[n].derived.cooling_bb_simple_tot;
    x[n * NKPKT_SCALARS + 8] = macromain[n].derived.kpkt_rates_known;
    x[n * NKPKT_SCALARS + 9] = 1.0;     /* a check that each cell was calculated once */
  }

  sum_block_in_place (x, (long) NKPKT_SCALARS * NPLASMA, MPI_COMM_WORLD);

  for (n = 0; n < NPLASMA; n++)
  {
    if (x[n * NKPKT_SCALARS + 9] != 1.0)
    {
      Error ("share_macro_recomb_and_kpkt_rates: cell %d was calculated %g times\n", n, x[n * NKPKT_SCALARS + 9]);
    }
    macromain[n].derived.cooling_normalisation = x[n * NKPKT_SCALARS];
    macromain[n].derived.cooling_bbtot = x[n * NKPKT_SCALARS + 1];
    macromain[n].derived.cooling_bftot = x[n * NKPKT_SCALARS + 2];
    macromain[n].derived.cooling_bf_coltot = x[n * NKPKT_SCALARS + 3];
    macromain[n].derived.cooling_ff = x[n * NKPKT_SCALARS + 4];
    macromain[n].derived.cooling_ff_lofreq = x[n * NKPKT_SCALARS + 5];
    macromain[n].derived.cooling_adiabatic = x[n * NKPKT_SCALARS + 6];
    macromain[n].derived.cooling_bb_simple_tot = x[n * NKPKT_SCALARS + 7];
    macromain[n].derived.kpkt_rates_known = (x[n * NKPKT_SCALARS + 8] > 0.5) ? TRUE : FALSE;
  }

  free (x);

  d_xsignal (files.root, "%-20s Finished sharing of macro atom recombination and k-packet rates\n", "OK");
#endif
}
