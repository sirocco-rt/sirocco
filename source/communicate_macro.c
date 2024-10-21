/***********************************************************/
/** @file  communicate_macro.c
 * @author EJP
 * @date   December 2023
 *
 * @brief Functions for communicating macro atom properties
 *
 * @TODO: as much as this as possible should use non-blocking communication
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
 * The communication pattern is as outlined in broadcast_updated_plasma_properties.
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

  for (current_rank = 0; current_rank < np_mpi_global; current_rank++)
  {
    position = 0;

    if (rank_global == current_rank)
    {
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
      for (n_plasma = n_start; n_plasma < n_stop; ++n_plasma)
      {
        MPI_Pack (&n_plasma, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].kpkt_emiss, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].matom_emiss, nlevels_macro, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
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
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n_plasma].kpkt_emiss, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].matom_emiss, nlevels_macro, MPI_DOUBLE, MPI_COMM_WORLD);
      }
    }
  }

  free (comm_buffer);
  d_xsignal (files.root, "%-20s Finished macro atom emissivity communication\n", "OK");
#endif
}

/**********************************************************/
/**
 * @brief  Communicate the macro atom recombination properties between ranks
 *
 * @param [in] int n_start       The index of the first cell updated by this rank
 * @param [in] int n_stop        The index of the last cell updated by this rank
 * @param [in] int n_cells_rank  The number of cells this rank updated
 *
 * @details
 *
 * The communication pattern is as outlined in broadcast_updated_plasma_properties.
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
  const int comm_buffer_size = calculate_comm_buffer_size (1 + n_cells_max, n_cells_max * (2 * size_alpha_est + 2 * nphot_total));
  char *const comm_buffer = malloc (comm_buffer_size);

  for (current_rank = 0; current_rank < np_mpi_global; ++current_rank)
  {
    position = 0;

    if (rank_global == current_rank)
    {
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);   // how many cells to unpack
      for (n_plasma = n_start; n_plasma < n_stop; ++n_plasma)
      {
        MPI_Pack (&n_plasma, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);     // which cell we're working on
        if (nlevels_macro > 0)
        {
          MPI_Pack (macromain[n_plasma].recomb_sp, size_alpha_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack (macromain[n_plasma].recomb_sp_e, size_alpha_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        }
        if (nphot_total > 0)
        {
          MPI_Pack (plasmamain[n_plasma].recomb_simple, nphot_total, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack (plasmamain[n_plasma].recomb_simple_upweight, nphot_total, MPI_DOUBLE, comm_buffer, comm_buffer_size,
                    &position, MPI_COMM_WORLD);
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

        if (nlevels_macro > 0)
        {
          MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].recomb_sp, size_alpha_est, MPI_DOUBLE, MPI_COMM_WORLD);
          MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].recomb_sp_e,
                      size_alpha_est, MPI_DOUBLE, MPI_COMM_WORLD);
        }
        if (nphot_total > 0)
        {
          MPI_Unpack (comm_buffer, comm_buffer_size, &position, plasmamain[n_plasma].recomb_simple,
                      nphot_total, MPI_DOUBLE, MPI_COMM_WORLD);
          MPI_Unpack (comm_buffer, comm_buffer_size, &position, plasmamain[n_plasma].recomb_simple_upweight,
                      nphot_total, MPI_DOUBLE, MPI_COMM_WORLD);
        }
      }
    }
  }

  free (comm_buffer);
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
 * The communication pattern is as outlined in broadcast_updated_plasma_properties.
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
  const int comm_buffer_size = calculate_comm_buffer_size (1 + 3 * n_cells_max, n_cells_max * (6 * size_gamma_est + 2 * size_Jbar_est));
  char *const comm_buffer = malloc (comm_buffer_size);

  for (current_rank = 0; current_rank < np_mpi_global; ++current_rank)
  {
    position = 0;

    if (rank_global == current_rank)
    {
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
      for (n_plasma = n_start; n_plasma < n_stop; ++n_plasma)
      {
        MPI_Pack (&n_plasma, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].jbar, size_Jbar_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].jbar_old, size_Jbar_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].gamma, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].gamma_old, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].gamma_e, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].gamma_e_old, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].alpha_st, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].alpha_st_old, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&macromain[n_plasma].kpkt_rates_known, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&macromain[n_plasma].matrix_rates_known, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
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
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].jbar, size_Jbar_est, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].jbar_old, size_Jbar_est, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].gamma, size_gamma_est, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].gamma_old, size_gamma_est, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].gamma_e, size_gamma_est, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].gamma_e_old, size_gamma_est, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].alpha_st, size_gamma_est, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n_plasma].alpha_st_old, size_gamma_est, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &macromain[n_plasma].kpkt_rates_known, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &macromain[n_plasma].matrix_rates_known, 1, MPI_INT, MPI_COMM_WORLD);
      }
    }
  }

  free (comm_buffer);
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
 **********************************************************/

int
broadcast_macro_atom_state_matrix (int n_start, int n_stop, int n_cells_rank)
{
#ifdef MPI_ON
  int n_mpi, n_mpi2, num_comm;
  int n, position;

  d_xsignal (files.root, "%-20s Begin macro atom state matrix communication\n", "NOK");
  const int matrix_size = nlevels_macro + 1;
  const int n_cells_max = get_max_cells_per_rank (NPLASMA);
  const int comm_buffer_size = calculate_comm_buffer_size (1 + n_cells_max, n_cells_max * (matrix_size * matrix_size));
  char *comm_buffer = malloc (comm_buffer_size);

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
        if (macromain[n].store_matom_matrix == TRUE)
        {
          MPI_Pack (macromain[n].matom_matrix[0], matrix_size * matrix_size, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position,
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
        if (macromain[n].store_matom_matrix == TRUE)
        {
          MPI_Unpack (comm_buffer, comm_buffer_size, &position, macromain[n].matom_matrix[0], matrix_size * matrix_size, MPI_DOUBLE,
                      MPI_COMM_WORLD);
        }
      }
    }
  }

  free (comm_buffer);
  d_xsignal (files.root, "%-20s Finished macro atom state matrix communication\n", "OK");
#endif
  return (0);
}

/**********************************************************/
/**
 * @brief
 *
 * @details averages the macro-atom estimators between tasks using MPI_Reduce.
 *   It should only be called if the MPI_ON flag was present
 *   in compilation, and returns 0 immediately if no macro atom levels.
 *   This should probably be improved by working out exactly
 *   what is needed in simple-ion only mode.
 *
 **********************************************************/

void
reduce_macro_atom_estimators (void)
{
#ifdef MPI_ON                   // these routines should only be called anyway in parallel but we need these to compile

  int n, mpi_i;
  double *gamma_helper, *alpha_helper;
  double *level_helper, *cell_helper, *jbar_helper;
  double *gamma_helper2, *alpha_helper2;
  double *level_helper2, *cell_helper2, *jbar_helper2;
  double *cooling_bf_helper, *cooling_bb_helper;
  double *cooling_bf_helper2, *cooling_bb_helper2;

  d_xsignal (files.root, "%-20s Begin reduction of macro atom estimators\n", "NOK");

  if (nlevels_macro == 0 && geo.nmacro == 0)
  {
    /* in this case no space would have been allocated for macro-atom estimators */
    Log ("No need to communicate matom estimators as no macro-atoms!\n");
    return;
  }

  /* allocate helper arrays for the estimators we want to communicate */
  /* the sizes of these arrays should match the allocation in calloc_estimators in gridwind.c */
  /* we need two arrays for each set of variables. Note that we stick all estimators of
     the same size in the same helper array */
  /* could put some error conditions here to check memory allocation worked */
  jbar_helper = calloc (sizeof (double), NPLASMA * size_Jbar_est);
  gamma_helper = calloc (sizeof (double), NPLASMA * 4 * size_gamma_est);
  alpha_helper = calloc (sizeof (double), NPLASMA * 2 * size_alpha_est);
  level_helper = calloc (sizeof (double), NPLASMA * nlevels_macro);
  cell_helper = calloc (sizeof (double), 8 * NPLASMA);
  cooling_bf_helper = calloc (sizeof (double), NPLASMA * 2 * nphot_total);
  cooling_bb_helper = calloc (sizeof (double), NPLASMA * nlines);

  jbar_helper2 = calloc (sizeof (double), NPLASMA * size_Jbar_est);
  gamma_helper2 = calloc (sizeof (double), NPLASMA * 4 * size_gamma_est);
  alpha_helper2 = calloc (sizeof (double), NPLASMA * 2 * size_alpha_est);
  level_helper2 = calloc (sizeof (double), NPLASMA * nlevels_macro);
  cell_helper2 = calloc (sizeof (double), 8 * NPLASMA);
  cooling_bf_helper2 = calloc (sizeof (double), NPLASMA * 2 * nphot_total);
  cooling_bb_helper2 = calloc (sizeof (double), NPLASMA * nlines);

  /* now we loop through each cell and copy the values of our variables
     into our helper arrays */
  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    /* one kpkt_abs quantity per cell */
    cell_helper[mpi_i] = plasmamain[mpi_i].kpkt_abs / np_mpi_global;

    /* each of the cooling sums and normalisations also have one quantity per cell */
    cell_helper[mpi_i + NPLASMA] = macromain[mpi_i].cooling_normalisation / np_mpi_global;
    cell_helper[mpi_i + 2 * NPLASMA] = macromain[mpi_i].cooling_bftot / np_mpi_global;
    cell_helper[mpi_i + 3 * NPLASMA] = macromain[mpi_i].cooling_bf_coltot / np_mpi_global;
    cell_helper[mpi_i + 4 * NPLASMA] = macromain[mpi_i].cooling_bbtot / np_mpi_global;
    cell_helper[mpi_i + 5 * NPLASMA] = macromain[mpi_i].cooling_ff / np_mpi_global;
    cell_helper[mpi_i + 6 * NPLASMA] = macromain[mpi_i].cooling_ff_lofreq / np_mpi_global;
    cell_helper[mpi_i + 7 * NPLASMA] = macromain[mpi_i].cooling_adiabatic / np_mpi_global;



    for (n = 0; n < nlevels_macro; n++)
    {
      level_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].matom_abs[n] / np_mpi_global;
    }

    for (n = 0; n < size_Jbar_est; n++)
    {
      jbar_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].jbar[n] / np_mpi_global;
    }

    for (n = 0; n < size_gamma_est; n++)
    {
      gamma_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].alpha_st[n] / np_mpi_global;
      gamma_helper[mpi_i + ((n + size_gamma_est) * NPLASMA)] = macromain[mpi_i].alpha_st_e[n] / np_mpi_global;
      gamma_helper[mpi_i + ((n + 2 * size_gamma_est) * NPLASMA)] = macromain[mpi_i].gamma[n] / np_mpi_global;
      gamma_helper[mpi_i + ((n + 3 * size_gamma_est) * NPLASMA)] = macromain[mpi_i].gamma_e[n] / np_mpi_global;
    }

    for (n = 0; n < size_alpha_est; n++)
    {
      alpha_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].recomb_sp[n] / np_mpi_global;
      alpha_helper[mpi_i + ((n + size_alpha_est) * NPLASMA)] = macromain[mpi_i].recomb_sp_e[n] / np_mpi_global;
    }

    for (n = 0; n < nphot_total; n++)
    {
      cooling_bf_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].cooling_bf[n] / np_mpi_global;
      cooling_bf_helper[mpi_i + ((n + nphot_total) * NPLASMA)] = macromain[mpi_i].cooling_bf_col[n] / np_mpi_global;
    }

    for (n = 0; n < nlines; n++)
    {
      cooling_bb_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].cooling_bb[n] / np_mpi_global;
    }
  }

  /* because in the above loop we have already divided by number of processes, we can now do a sum
     with MPI_Reduce, passing it MPI_SUM as an argument. This will give us the mean across threads */
  MPI_Allreduce (cell_helper, cell_helper2, NPLASMA * 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (level_helper, level_helper2, NPLASMA * nlevels_macro, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (jbar_helper, jbar_helper2, NPLASMA * size_Jbar_est, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (gamma_helper, gamma_helper2, NPLASMA * 4 * size_gamma_est, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (alpha_helper, alpha_helper2, NPLASMA * 2 * size_alpha_est, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (cooling_bf_helper, cooling_bf_helper2, NPLASMA * 2 * nphot_total, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (cooling_bb_helper, cooling_bb_helper2, NPLASMA * nlines, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /* We now need to copy these reduced variables to the plasma structure in each thread */


  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    /* one kpkt_abs quantity per cell */
    plasmamain[mpi_i].kpkt_abs = cell_helper2[mpi_i];

    /* each of the cooling sums and normalisations also have one quantity per cell */
    macromain[mpi_i].cooling_normalisation = cell_helper2[mpi_i + NPLASMA];
    macromain[mpi_i].cooling_bftot = cell_helper2[mpi_i + 2 * NPLASMA];
    macromain[mpi_i].cooling_bf_coltot = cell_helper2[mpi_i + 3 * NPLASMA];
    macromain[mpi_i].cooling_bbtot = cell_helper2[mpi_i + 4 * NPLASMA];
    macromain[mpi_i].cooling_ff = cell_helper2[mpi_i + 5 * NPLASMA];
    macromain[mpi_i].cooling_ff_lofreq = cell_helper2[mpi_i + 6 * NPLASMA];
    macromain[mpi_i].cooling_adiabatic = cell_helper2[mpi_i + 7 * NPLASMA];


    for (n = 0; n < nlevels_macro; n++)
    {
      macromain[mpi_i].matom_abs[n] = level_helper2[mpi_i + (n * NPLASMA)];
    }

    for (n = 0; n < size_Jbar_est; n++)
    {
      macromain[mpi_i].jbar[n] = jbar_helper2[mpi_i + (n * NPLASMA)];
    }

    for (n = 0; n < size_gamma_est; n++)
    {
      macromain[mpi_i].alpha_st[n] = gamma_helper2[mpi_i + (n * NPLASMA)];
      macromain[mpi_i].alpha_st_e[n] = gamma_helper2[mpi_i + ((n + size_gamma_est) * NPLASMA)] / np_mpi_global;
      macromain[mpi_i].gamma[n] = gamma_helper2[mpi_i + ((n + 2 * size_gamma_est) * NPLASMA)];
      macromain[mpi_i].gamma_e[n] = gamma_helper2[mpi_i + ((n + 3 * size_gamma_est) * NPLASMA)];
    }

    for (n = 0; n < size_alpha_est; n++)
    {
      macromain[mpi_i].recomb_sp[n] = alpha_helper2[mpi_i + (n * NPLASMA)];
      macromain[mpi_i].recomb_sp_e[n] = alpha_helper2[mpi_i + ((n + size_alpha_est) * NPLASMA)];
    }

    for (n = 0; n < nphot_total; n++)
    {
      macromain[mpi_i].cooling_bf[n] = cooling_bf_helper2[mpi_i + (n * NPLASMA)];
      macromain[mpi_i].cooling_bf_col[n] = cooling_bf_helper2[mpi_i + ((n + nphot_total) * NPLASMA)];
    }

    for (n = 0; n < nlines; n++)
    {
      macromain[mpi_i].cooling_bb[n] = cooling_bb_helper2[mpi_i + (n * NPLASMA)];
    }
  }

  free (cell_helper);
  free (level_helper);
  free (jbar_helper);
  free (gamma_helper);
  free (alpha_helper);
  free (cooling_bf_helper);
  free (cooling_bb_helper);

  free (cell_helper2);
  free (level_helper2);
  free (jbar_helper2);
  free (gamma_helper2);
  free (alpha_helper2);
  free (cooling_bf_helper2);
  free (cooling_bb_helper2);

  d_xsignal (files.root, "%-20s Finished reduction of macro atom estimators\n", "OK");
#endif
}
