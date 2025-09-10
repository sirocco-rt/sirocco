/***********************************************************/
/** @file  communicate_wind.c
 * @author EJP
 * @date   December 2023
 *
 * @brief Functions for communicating wind properties
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
 * @brief Broadcast the wind grid, `wmain`, to all ranks
 *
 * @param [in] int n_start       The index of the first cell updated by this rank
 * @param [in] int n_stop        The index of the last cell updated by this rank
 * @param [in] int n_cells_rank  The number of cells this rank updated
 *
 * @details
 *
 * This function should only be called once, after grid initialisation.
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
broadcast_wind_grid (const int n_start, const int n_stop, const int n_cells_rank)
{
#ifdef MPI_ON
  int i;
  int n_wind;
  int current_rank;
  int num_comm;
  int position;
  int bytes_wcone;

  WindPtr cell;

  d_xsignal (files.root, "%-20s Begin communication of wind grid\n", "NOK");

  /* For the wcone field, we should create a derived type to send the struct
   * more efficiently. Although, it may just as easy and quick to communicate
   * each field one by one... but this is the right way to do it */
  MPI_Datatype wcone_derived_type;
  const int count = 2;
  const int block_lengths[] = { 1, 1 };
  const MPI_Datatype block_types[] = { MPI_DOUBLE, MPI_DOUBLE };
  /* We need to find the memory displacements. We'll use the wcone struct in
   * the first cell for this. Each struct should have the same amount of
   * alignment for the fields, so this should be OK */
  MPI_Aint base_address;
  MPI_Aint block_offsets[count];
  MPI_Get_address (&wmain[0].wcone, &base_address);
  MPI_Get_address (&wmain[0].wcone.z, &block_offsets[0]);
  MPI_Get_address (&wmain[0].wcone.dzdr, &block_offsets[1]);
  for (i = 0; i < 2; ++i)
  {
    block_offsets[i] = MPI_Aint_diff (block_offsets[i], base_address);
  }
  MPI_Type_create_struct (count, block_lengths, block_offsets, block_types, &wcone_derived_type);
  MPI_Type_commit (&wcone_derived_type);

  /* Calculate the size of the communication buffer */
  const int n_cells_max = get_max_cells_per_rank (NDIM2);
  // const int num_ints = 5 * n_cells_max + 1;
  // const int num_doubles = n_cells_max * (13 + 3 * 3 + 1 * 9);   // *3 for x, xcen... *9 for v_grad
  MPI_Pack_size (n_cells_max, wcone_derived_type, MPI_COMM_WORLD, &bytes_wcone);
  const int comm_buffer_size = calculate_comm_buffer_size (1 + 5 * n_cells_max, n_cells_max * (13 + 3 * 3 + 1 * 9)) + bytes_wcone;

  char *comm_buffer = malloc (comm_buffer_size);
  if (comm_buffer == NULL)
  {
    Error ("Unable to allocate memory for communication buffer in broadcast_wind_grid\n");
    Exit (EXIT_FAILURE);
  }

  for (current_rank = 0; current_rank < np_mpi_global; current_rank++)
  {
    position = 0;

    if (rank_global == current_rank)
    {
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
      for (n_wind = n_start; n_wind < n_stop; ++n_wind)
      {
        cell = &wmain[n_wind];
        MPI_Pack (&n_wind, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->ndom, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->nwind_dom, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->nplasma, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->x, 3, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->xcen, 3, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->r, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->rcen, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->theta, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->thetacen, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->dtheta, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->dr, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->wcone, 1, wcone_derived_type, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->v, 3, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->v_grad, 9, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->div_v, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->dvds_ave, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->dvds_max, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->vol, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->xgamma, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->xgamma_cen, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->dfudge, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->inwind, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
      }
    }

    MPI_Bcast (comm_buffer, comm_buffer_size, MPI_PACKED, current_rank, MPI_COMM_WORLD);

    position = 0;

    if (rank_global != current_rank)
    {
      MPI_Unpack (comm_buffer, comm_buffer_size, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
      for (i = 0; i < num_comm; i++)
      {
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &n_wind, 1, MPI_INT, MPI_COMM_WORLD);
        cell = &wmain[n_wind];
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->ndom, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->nwind_dom, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->nplasma, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->x, 3, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->xcen, 3, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->rcen, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->theta, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->thetacen, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->dtheta, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->dr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->wcone, 1, wcone_derived_type, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->v, 3, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->v_grad, 9, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->div_v, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->dvds_ave, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->dvds_max, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->vol, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->xgamma, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->xgamma_cen, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->dfudge, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->inwind, 1, MPI_INT, MPI_COMM_WORLD);
      }
    }
  }

  free (comm_buffer);
  MPI_Type_free (&wcone_derived_type);
  d_xsignal (files.root, "%-20s Finished communication of wind grid\n", "NOK");
#endif
}
