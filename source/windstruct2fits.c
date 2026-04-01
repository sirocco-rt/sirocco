
/***********************************************************/
/** @file  windstruct2fits.c
 * @author ksl
 * @date   April, 2026
 *
 * @brief  Write wind struct (wmain) for all cells — including
 * lower-hemisphere transport cells — to a FITS binary table.
 *
 * ### Notes ###
 * windsave2fits only covers the plasma cells (upper hemisphere).
 * This program covers all 2*NDIM2 wind cells so that quantities
 * stored per wind cell (position, velocity, velocity gradient,
 * cell boundaries, dvds) can be inspected and compared between
 * the two hemispheres.
 *
 * Usage:
 *   windstruct2fits root
 *
 * Reads root.wind_save and writes root_windstruct.fits
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"
#include "fitsio.h"


/**********************************************************/
/**
 * @brief  Write one binary-table extension with one row per
 *         wind cell (all 2*NDIM2 cells).
 *
 * Columns written:
 *   nwind      J   - cell index in wmain
 *   ndom       J   - domain number
 *   nplasma    J   - plasma cell index (upper+lower share same)
 *   inwind     J   - W_ALL_INWIND / W_NOT_INWIND / ...
 *   x         3E   - inner vertex position (cm)
 *   xcen      3E   - cell centre position (cm)
 *   v         3E   - velocity at inner vertex (cm/s)
 *   v_grad    9E   - velocity gradient tensor (row-major, 1/s)
 *   cell_z_min E   - signed lower z boundary (cm)
 *   cell_z_max E   - signed upper z boundary (cm)
 *   dvds_ave   E   - direction-averaged |dv/ds| (1/s)
 *   dvds_max   E   - maximum |dv/ds| at vertex (1/s)
 *
 **********************************************************/

int
write_wind_table (fitsfile * fptr)
{
  int status = 0;
  int num_rows, num_cols;
  int n;
  WindPtr cell;

  num_rows = 2 * NDIM2;
  num_cols = 12;

  /* Column definitions */
  char *ttype[] = {
    "nwind", "ndom", "nplasma", "inwind",
    "x", "xcen", "v", "v_grad",
    "cell_z_min", "cell_z_max",
    "dvds_ave", "dvds_max"
  };
  char *tform[] = {
    "J", "J", "J", "J",
    "3E", "3E", "3E", "9E",
    "E", "E",
    "E", "E"
  };
  char *tunit[] = {
    "", "", "", "",
    "cm", "cm", "cm/s", "1/s",
    "cm", "cm",
    "1/s", "1/s"
  };

  if (fits_create_tbl (fptr, BINARY_TBL, num_rows, num_cols, ttype, tform, tunit, "windstruct", &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  /* Allocate column buffers */
  int *col_nwind = calloc (num_rows, sizeof (int));
  int *col_ndom = calloc (num_rows, sizeof (int));
  int *col_nplasma = calloc (num_rows, sizeof (int));
  int *col_inwind = calloc (num_rows, sizeof (int));
  float *col_x = calloc (num_rows * 3, sizeof (float));
  float *col_xcen = calloc (num_rows * 3, sizeof (float));
  float *col_v = calloc (num_rows * 3, sizeof (float));
  float *col_v_grad = calloc (num_rows * 9, sizeof (float));
  float *col_cell_z_min = calloc (num_rows, sizeof (float));
  float *col_cell_z_max = calloc (num_rows, sizeof (float));
  float *col_dvds_ave = calloc (num_rows, sizeof (float));
  float *col_dvds_max = calloc (num_rows, sizeof (float));

  if (!col_nwind || !col_ndom || !col_nplasma || !col_inwind ||
      !col_x || !col_xcen || !col_v || !col_v_grad ||
      !col_cell_z_min || !col_cell_z_max || !col_dvds_ave || !col_dvds_max)
  {
    fprintf (stderr, "windstruct2fits: memory allocation failed\n");
    return -1;
  }

  for (n = 0; n < num_rows; n++)
  {
    cell = &wmain[n];

    col_nwind[n] = cell->nwind;
    col_ndom[n] = cell->ndom;
    col_nplasma[n] = cell->nplasma;
    col_inwind[n] = (int) cell->inwind;

    col_x[n * 3 + 0] = (float) cell->x[0];
    col_x[n * 3 + 1] = (float) cell->x[1];
    col_x[n * 3 + 2] = (float) cell->x[2];

    col_xcen[n * 3 + 0] = (float) cell->xcen[0];
    col_xcen[n * 3 + 1] = (float) cell->xcen[1];
    col_xcen[n * 3 + 2] = (float) cell->xcen[2];

    col_v[n * 3 + 0] = (float) cell->v[0];
    col_v[n * 3 + 1] = (float) cell->v[1];
    col_v[n * 3 + 2] = (float) cell->v[2];

    col_v_grad[n * 9 + 0] = (float) cell->v_grad[0][0];
    col_v_grad[n * 9 + 1] = (float) cell->v_grad[0][1];
    col_v_grad[n * 9 + 2] = (float) cell->v_grad[0][2];
    col_v_grad[n * 9 + 3] = (float) cell->v_grad[1][0];
    col_v_grad[n * 9 + 4] = (float) cell->v_grad[1][1];
    col_v_grad[n * 9 + 5] = (float) cell->v_grad[1][2];
    col_v_grad[n * 9 + 6] = (float) cell->v_grad[2][0];
    col_v_grad[n * 9 + 7] = (float) cell->v_grad[2][1];
    col_v_grad[n * 9 + 8] = (float) cell->v_grad[2][2];

    col_cell_z_min[n] = (float) cell->cell_z_min;
    col_cell_z_max[n] = (float) cell->cell_z_max;
    col_dvds_ave[n] = (float) cell->dvds_ave;
    col_dvds_max[n] = (float) cell->dvds_max;
  }

  /* Write columns */
  fits_write_col (fptr, TINT, 1, 1, 1, num_rows, col_nwind, &status);
  fits_write_col (fptr, TINT, 2, 1, 1, num_rows, col_ndom, &status);
  fits_write_col (fptr, TINT, 3, 1, 1, num_rows, col_nplasma, &status);
  fits_write_col (fptr, TINT, 4, 1, 1, num_rows, col_inwind, &status);
  fits_write_col (fptr, TFLOAT, 5, 1, 1, num_rows * 3, col_x, &status);
  fits_write_col (fptr, TFLOAT, 6, 1, 1, num_rows * 3, col_xcen, &status);
  fits_write_col (fptr, TFLOAT, 7, 1, 1, num_rows * 3, col_v, &status);
  fits_write_col (fptr, TFLOAT, 8, 1, 1, num_rows * 9, col_v_grad, &status);
  fits_write_col (fptr, TFLOAT, 9, 1, 1, num_rows, col_cell_z_min, &status);
  fits_write_col (fptr, TFLOAT, 10, 1, 1, num_rows, col_cell_z_max, &status);
  fits_write_col (fptr, TFLOAT, 11, 1, 1, num_rows, col_dvds_ave, &status);
  fits_write_col (fptr, TFLOAT, 12, 1, 1, num_rows, col_dvds_max, &status);

  if (status)
    fits_report_error (stderr, status);

  free (col_nwind);
  free (col_ndom);
  free (col_nplasma);
  free (col_inwind);
  free (col_x);
  free (col_xcen);
  free (col_v);
  free (col_v_grad);
  free (col_cell_z_min);
  free (col_cell_z_max);
  free (col_dvds_ave);
  free (col_dvds_max);

  return status;
}


/**********************************************************/
/**
 * @brief  main
 **********************************************************/

int
main (int argc, char *argv[])
{
  char inroot[LINELENGTH];
  char infile[LINELENGTH];
  char outfile[LINELENGTH];
  fitsfile *fptr;
  int status = 0;

  if (argc < 2)
  {
    fprintf (stderr, "Usage: windstruct2fits root\n");
    fprintf (stderr, "  Reads root.wind_save and writes root_windstruct.fits\n");
    return EXIT_FAILURE;
  }

  get_root (inroot, argv[argc - 1]);
  sprintf (infile, "%.150s.wind_save", inroot);
  sprintf (outfile, "!%.150s_windstruct.fits", inroot);

  printf ("Reading %s\n", infile);
  printf ("Writing %s\n", outfile + 1);   /* skip the leading '!' in the message */

  zdom = calloc (MAX_DOM, sizeof (domain_dummy));
  if (zdom == NULL)
  {
    fprintf (stderr, "windstruct2fits: unable to allocate domain memory\n");
    return EXIT_FAILURE;
  }

  wind_read (infile);

  printf ("NDIM2 = %d  (total wind cells = %d)\n", NDIM2, 2 * NDIM2);

  if (fits_create_file (&fptr, outfile, &status))
  {
    fits_report_error (stderr, status);
    return EXIT_FAILURE;
  }

  /* Empty primary HDU required by FITS standard */
  if (fits_create_img (fptr, FLOAT_IMG, 0, NULL, &status))
  {
    fits_report_error (stderr, status);
  }

  write_wind_table (fptr);

  if (fits_close_file (fptr, &status))
    fits_report_error (stderr, status);

  return EXIT_SUCCESS;
}
