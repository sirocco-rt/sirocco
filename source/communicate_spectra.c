/***********************************************************/
/** @file  communicate_wind.c
 * @author EJP
 * @date   December 2023
 *
 * @brief Functions for communicating wind properties
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
 * @brief Normalize the synthetic spectra across ranks.
 *
 * @details
 *
 * Uses an all-to-all reduction to sum up and normalize the synthetic spectra
 * between MPI ranks using an MPI_Allreduce. Works on both the linear and the
 * log spectra.
 *
 **********************************************************/

int
normalize_spectra_across_ranks (void)
{
#ifdef MPI_ON
  int i;
  int j;
  int nspec;
  int size_of_commbuffer;;
  double *spectrum_buffer;

  d_xsignal (files.root, "%-20s Begin spectrum reduction\n", "NOK");

  /* In ionization cycles, we don't need to normalize/reduce an additional
   * geo.nangles, as these are for the observer angles in spectrum cycles */
  if (geo.ioniz_or_extract == CYCLE_EXTRACT)
  {
    nspec = MSPEC + geo.nangles;
  }
  else
  {
    nspec = MSPEC;
  }

  size_of_commbuffer = 4 * nspec * NWAVE_MAX;   // We need space for all 4 separate spectra we are normalizing
  spectrum_buffer = calloc (sizeof (double), size_of_commbuffer);

  for (i = 0; i < NWAVE_MAX; i++)
  {
    for (j = 0; j < nspec; j++)
    {
      spectrum_buffer[i * nspec + j] = xxspec[j].f[i] / np_mpi_global;
      spectrum_buffer[i * nspec + j + (NWAVE_MAX * nspec)] = xxspec[j].lf[i] / np_mpi_global;
      spectrum_buffer[i * nspec + j + (2 * NWAVE_MAX * nspec)] = xxspec[j].f_wind[i] / np_mpi_global;
      spectrum_buffer[i * nspec + j + (3 * NWAVE_MAX * nspec)] = xxspec[j].lf_wind[i] / np_mpi_global;
    }
  }

  MPI_Allreduce (MPI_IN_PLACE, spectrum_buffer, size_of_commbuffer, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (i = 0; i < NWAVE_MAX; i++)
  {
    for (j = 0; j < nspec; j++)
    {
      xxspec[j].f[i] = spectrum_buffer[i * nspec + j];
      xxspec[j].lf[i] = spectrum_buffer[i * nspec + j + (NWAVE_MAX * nspec)];
      xxspec[j].f_wind[i] = spectrum_buffer[i * nspec + j + (2 * NWAVE_MAX * nspec)];
      xxspec[j].lf_wind[i] = spectrum_buffer[i * nspec + j + (3 * NWAVE_MAX * nspec)];
    }
  }

  free (spectrum_buffer);
  d_xsignal (files.root, "%-20s Finished spectrum reduction\n", "OK");
#endif
  return (0);
}
