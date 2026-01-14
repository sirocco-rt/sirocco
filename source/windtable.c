/***********************************************************/
/** @file  windtable.c
 * @author Long, Matthews
 * @date   January 2026
 *
 * @brief  Main driver for windtable - exports wind data to FITS binary table
 *
 * windtable reads a Sirocco wind save file and writes all key variables
 * to a single FITS binary table. This provides a simpler alternative to
 * windsave2table which writes multiple ASCII files.
 *
 * Usage: windtable [-e] rootname
 *
 * Options:
 *   -e, --edge    Include edge cells in output (cells at domain boundaries)
 *   --version     Print version information
 *
 * Output: rootname.windtable.fits (or rootname.N.windtable.fits for multi-domain)
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atomic.h"
#include "sirocco.h"
#include "fitsio.h"
#include "windtable.h"


/* External variable for edge cell handling - shared with windsave2table_sub.c */
extern int xedge;


/**********************************************************/
/**
 * @brief  Print program usage information
 *
 ***********************************************************/

void
windtable_usage (void)
{
  printf ("Usage: windtable [-e] rootname\n");
  printf ("\n");
  printf ("Options:\n");
  printf ("  -e, --edge    Include edge cells in output\n");
  printf ("  --version     Print version information\n");
  printf ("\n");
  printf ("Output: rootname.windtable.fits\n");
}


/**********************************************************/
/**
 * @brief  Parse command line arguments
 *
 * @param [in]  argc        Number of command line arguments
 * @param [in]  argv        Command line argument array
 * @param [out] root        Root name of wind save file
 * @param [out] edge_switch Set to 1 if edge cells should be included
 * @return                  0 on success, -1 on error
 *
 ***********************************************************/

int
windtable_parse_args (int argc, char *argv[], char root[], int *edge_switch)
{
  int i;
  int got_root = 0;

  *edge_switch = 0;
  root[0] = '\0';

  for (i = 1; i < argc; i++)
  {
    if (strcmp (argv[i], "--version") == 0)
    {
      printf ("windtable version %s\n", VERSION);
      printf ("Git commit: %s\n", GIT_COMMIT_HASH);
      exit (0);
    }
    else if (strcmp (argv[i], "-e") == 0 || strcmp (argv[i], "--edge") == 0)
    {
      *edge_switch = 1;
    }
    else if (argv[i][0] == '-')
    {
      printf ("Unknown option: %s\n", argv[i]);
      windtable_usage ();
      return -1;
    }
    else
    {
      /* This should be the root name */
      strcpy (root, argv[i]);
      got_root = 1;
    }
  }

  if (!got_root)
  {
    printf ("Error: No root name specified\n");
    windtable_usage ();
    return -1;
  }

  return 0;
}


/**********************************************************/
/**
 * @brief  Write FITS binary table for one domain
 *
 * @param [in] ndom       Domain number
 * @param [in] rootname   Root name for output file
 * @param [in] vars       Array of variable descriptors
 * @param [in] nvars      Number of variables
 * @param [in] ions       Array of ion column descriptors
 * @param [in] nions      Number of ion columns
 * @return                0 on success, non-zero on error
 *
 * @details
 *
 * Creates a FITS file with a single binary table extension containing:
 * - Geometry columns (x, z, r, etc. from wmain[])
 * - Physical variables (from plasmamain[] via get_one())
 * - Ion fractions (from plasmamain[] via get_ion())
 *
 * One row per wind cell in the domain.
 *
 ***********************************************************/

int
windtable_write_fits (int ndom, char *rootname, windtable_var_t * vars, int nvars, windtable_ion_t * ions, int nions)
{
  fitsfile *fptr;
  int status = 0;
  char filename[LINELENGTH];
  int i, j, col;
  int nrows, ncols;
  int ngeom;
  int nstart, ndim2;
  int ii, jj;
  double *data;
  double *dcolumn;
  int *icolumn;
  char **ttype, **tform, **tunit;
  char name_dummy[132];

  /* Number of geometry columns:
   * i, j, nwind, nplasma, inwind (5 integers)
   * x, z, xcen, zcen, r, rcen, theta, thetacen, v_x, v_y, v_z (11 doubles)
   * Total = 16 columns
   */
  ngeom = 16;

  /* Total columns = geometry + variables + ions */
  ncols = ngeom + nvars + nions;

  /* Number of rows = number of cells in domain */
  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;
  nrows = ndim2;

  /* Build output filename - ! prefix overwrites existing file */
  sprintf (filename, "!%s.windtable.fits", rootname);

  Log ("windtable: Writing %s (domain %d, %d rows x %d cols)\n", filename + 1, ndom, nrows, ncols);

  /* Create FITS file */
  if (fits_create_file (&fptr, filename, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  /* Allocate column definition arrays */
  ttype = (char **) calloc (ncols, sizeof (char *));
  tform = (char **) calloc (ncols, sizeof (char *));
  tunit = (char **) calloc (ncols, sizeof (char *));

  for (i = 0; i < ncols; i++)
  {
    ttype[i] = (char *) calloc (WINDTABLE_VAR_NAME_LEN, sizeof (char));
    tform[i] = (char *) calloc (8, sizeof (char));
    tunit[i] = (char *) calloc (WINDTABLE_VAR_UNIT_LEN, sizeof (char));
  }

  /* Define geometry columns - cell identification first */
  col = 0;
  strcpy (ttype[col], "i");
  strcpy (tform[col], "1J");
  strcpy (tunit[col], "");
  col++;
  strcpy (ttype[col], "j");
  strcpy (tform[col], "1J");
  strcpy (tunit[col], "");
  col++;
  strcpy (ttype[col], "nwind");
  strcpy (tform[col], "1J");
  strcpy (tunit[col], "");
  col++;
  strcpy (ttype[col], "nplasma");
  strcpy (tform[col], "1J");
  strcpy (tunit[col], "");
  col++;
  strcpy (ttype[col], "inwind");
  strcpy (tform[col], "1J");
  strcpy (tunit[col], "");
  col++;

  /* Position columns */
  strcpy (ttype[col], "x");
  strcpy (tform[col], "1D");
  strcpy (tunit[col], "cm");
  col++;
  strcpy (ttype[col], "z");
  strcpy (tform[col], "1D");
  strcpy (tunit[col], "cm");
  col++;
  strcpy (ttype[col], "xcen");
  strcpy (tform[col], "1D");
  strcpy (tunit[col], "cm");
  col++;
  strcpy (ttype[col], "zcen");
  strcpy (tform[col], "1D");
  strcpy (tunit[col], "cm");
  col++;
  strcpy (ttype[col], "r");
  strcpy (tform[col], "1D");
  strcpy (tunit[col], "cm");
  col++;
  strcpy (ttype[col], "rcen");
  strcpy (tform[col], "1D");
  strcpy (tunit[col], "cm");
  col++;
  strcpy (ttype[col], "theta");
  strcpy (tform[col], "1D");
  strcpy (tunit[col], "rad");
  col++;
  strcpy (ttype[col], "thetacen");
  strcpy (tform[col], "1D");
  strcpy (tunit[col], "rad");
  col++;

  /* Velocity columns */
  strcpy (ttype[col], "v_x");
  strcpy (tform[col], "1D");
  strcpy (tunit[col], "cm/s");
  col++;
  strcpy (ttype[col], "v_y");
  strcpy (tform[col], "1D");
  strcpy (tunit[col], "cm/s");
  col++;
  strcpy (ttype[col], "v_z");
  strcpy (tform[col], "1D");
  strcpy (tunit[col], "cm/s");
  col++;

  /* Define variable columns */
  for (i = 0; i < nvars; i++)
  {
    strcpy (ttype[col], vars[i].name);
    strcpy (tform[col], vars[i].is_int ? "1J" : "1D");
    strcpy (tunit[col], vars[i].unit);
    col++;
  }

  /* Define ion columns */
  for (i = 0; i < nions; i++)
  {
    strcpy (ttype[col], ions[i].name);
    strcpy (tform[col], "1D");
    strcpy (tunit[col], "");
    col++;
  }

  /* Create binary table extension */
  if (fits_create_tbl (fptr, BINARY_TBL, nrows, ncols, ttype, tform, tunit, "WINDTABLE", &status))
  {
    fits_report_error (stderr, status);
    goto cleanup;
  }

  /* Add header keywords */
  fits_update_key (fptr, TSTRING, "ROOTNAME", rootname, "Sirocco run rootname", &status);
  fits_update_key (fptr, TINT, "NDOM", &ndom, "Domain number", &status);
  fits_update_key (fptr, TSTRING, "VERSION", VERSION, "Sirocco version", &status);
  fits_update_key (fptr, TINT, "NDIM2", &ndim2, "Number of wind cells", &status);
  fits_update_key (fptr, TINT, "COORDTYP", &zdom[ndom].coord_type, "Coordinate type (0=spherical,1=cylind,2=rtheta)", &status);

  /* Allocate data arrays */
  dcolumn = (double *) calloc (nrows, sizeof (double));
  icolumn = (int *) calloc (nrows, sizeof (int));

  /* Write geometry columns */
  col = 1;                      /* FITS columns are 1-indexed */

  /* i - grid index (uses wind_n_to_ij for 2D coords, just n for 1D) */
  for (j = 0; j < nrows; j++)
  {
    if (zdom[ndom].coord_type == SPHERICAL)
    {
      icolumn[j] = j;
    }
    else
    {
      wind_n_to_ij (ndom, nstart + j, &ii, &jj);
      icolumn[j] = ii;
    }
  }
  fits_write_col (fptr, TINT, col++, 1, 1, nrows, icolumn, &status);

  /* j - grid index (0 for 1D models) */
  for (j = 0; j < nrows; j++)
  {
    if (zdom[ndom].coord_type == SPHERICAL)
    {
      icolumn[j] = 0;
    }
    else
    {
      wind_n_to_ij (ndom, nstart + j, &ii, &jj);
      icolumn[j] = jj;
    }
  }
  fits_write_col (fptr, TINT, col++, 1, 1, nrows, icolumn, &status);

  /* nwind */
  for (j = 0; j < nrows; j++)
    icolumn[j] = wmain[nstart + j].nwind;
  fits_write_col (fptr, TINT, col++, 1, 1, nrows, icolumn, &status);

  /* nplasma */
  for (j = 0; j < nrows; j++)
    icolumn[j] = wmain[nstart + j].nplasma;
  fits_write_col (fptr, TINT, col++, 1, 1, nrows, icolumn, &status);

  /* inwind */
  for (j = 0; j < nrows; j++)
    icolumn[j] = wmain[nstart + j].inwind;
  fits_write_col (fptr, TINT, col++, 1, 1, nrows, icolumn, &status);

  /* x */
  for (j = 0; j < nrows; j++)
    dcolumn[j] = wmain[nstart + j].x[0];
  fits_write_col (fptr, TDOUBLE, col++, 1, 1, nrows, dcolumn, &status);

  /* z */
  for (j = 0; j < nrows; j++)
    dcolumn[j] = wmain[nstart + j].x[2];
  fits_write_col (fptr, TDOUBLE, col++, 1, 1, nrows, dcolumn, &status);

  /* xcen */
  for (j = 0; j < nrows; j++)
    dcolumn[j] = wmain[nstart + j].xcen[0];
  fits_write_col (fptr, TDOUBLE, col++, 1, 1, nrows, dcolumn, &status);

  /* zcen */
  for (j = 0; j < nrows; j++)
    dcolumn[j] = wmain[nstart + j].xcen[2];
  fits_write_col (fptr, TDOUBLE, col++, 1, 1, nrows, dcolumn, &status);

  /* r */
  for (j = 0; j < nrows; j++)
    dcolumn[j] = wmain[nstart + j].r;
  fits_write_col (fptr, TDOUBLE, col++, 1, 1, nrows, dcolumn, &status);

  /* rcen */
  for (j = 0; j < nrows; j++)
    dcolumn[j] = wmain[nstart + j].rcen;
  fits_write_col (fptr, TDOUBLE, col++, 1, 1, nrows, dcolumn, &status);

  /* theta */
  for (j = 0; j < nrows; j++)
    dcolumn[j] = wmain[nstart + j].theta;
  fits_write_col (fptr, TDOUBLE, col++, 1, 1, nrows, dcolumn, &status);

  /* thetacen */
  for (j = 0; j < nrows; j++)
    dcolumn[j] = wmain[nstart + j].thetacen;
  fits_write_col (fptr, TDOUBLE, col++, 1, 1, nrows, dcolumn, &status);

  /* v_x */
  for (j = 0; j < nrows; j++)
    dcolumn[j] = wmain[nstart + j].v[0];
  fits_write_col (fptr, TDOUBLE, col++, 1, 1, nrows, dcolumn, &status);

  /* v_y */
  for (j = 0; j < nrows; j++)
    dcolumn[j] = wmain[nstart + j].v[1];
  fits_write_col (fptr, TDOUBLE, col++, 1, 1, nrows, dcolumn, &status);

  /* v_z */
  for (j = 0; j < nrows; j++)
    dcolumn[j] = wmain[nstart + j].v[2];
  fits_write_col (fptr, TDOUBLE, col++, 1, 1, nrows, dcolumn, &status);

  /* Write variable columns using get_one() */
  for (i = 0; i < nvars; i++)
  {
    data = get_one (ndom, vars[i].get_one_name);
    if (data == NULL)
    {
      Error ("windtable: get_one failed for variable '%s'\n", vars[i].get_one_name);
      for (j = 0; j < nrows; j++)
        dcolumn[j] = 0.0;
    }
    else
    {
      for (j = 0; j < nrows; j++)
        dcolumn[j] = data[j];
      free (data);
    }

    if (vars[i].is_int)
    {
      for (j = 0; j < nrows; j++)
        icolumn[j] = (int) dcolumn[j];
      fits_write_col (fptr, TINT, col++, 1, 1, nrows, icolumn, &status);
    }
    else
    {
      fits_write_col (fptr, TDOUBLE, col++, 1, 1, nrows, dcolumn, &status);
    }

    if (status)
    {
      fits_report_error (stderr, status);
      Error ("windtable: Error writing column %s\n", vars[i].name);
    }
  }

  /* Write ion columns using get_ion() */
  for (i = 0; i < nions; i++)
  {
    data = get_ion (ndom, ions[i].z, ions[i].istate, 0, name_dummy);
    if (data == NULL)
    {
      Error ("windtable: get_ion failed for %s (Z=%d, istate=%d)\n", ions[i].name, ions[i].z, ions[i].istate);
      for (j = 0; j < nrows; j++)
        dcolumn[j] = 0.0;
    }
    else
    {
      for (j = 0; j < nrows; j++)
        dcolumn[j] = data[j];
      free (data);
    }

    fits_write_col (fptr, TDOUBLE, col++, 1, 1, nrows, dcolumn, &status);

    if (status)
    {
      fits_report_error (stderr, status);
      Error ("windtable: Error writing ion column %s\n", ions[i].name);
    }
  }

  free (dcolumn);
  free (icolumn);

cleanup:
  /* Free column definition arrays */
  for (i = 0; i < ncols; i++)
  {
    free (ttype[i]);
    free (tform[i]);
    free (tunit[i]);
  }
  free (ttype);
  free (tform);
  free (tunit);

  /* Close FITS file */
  fits_close_file (fptr, &status);

  if (status == 0)
  {
    Log ("windtable: Successfully wrote %s\n", filename + 1);
  }

  return status;
}


/**********************************************************/
/**
 * @brief  Main entry point for windtable
 *
 * @param [in] argc   Number of command line arguments
 * @param [in] argv   Command line argument array
 * @return            0 on success, non-zero on error
 *
 ***********************************************************/

int
main (int argc, char *argv[])
{
  char root[LINELENGTH], rootname[LINELENGTH];
  char windsavefile[LINELENGTH];
  int edge_switch;
  int ndom;
  int nvars, nions;
  windtable_var_t *vars;
  windtable_ion_t *ions;

  /* Set logging verbosity */
  Log_set_verbosity (3);

  /* Parse command line arguments */
  if (windtable_parse_args (argc, argv, root, &edge_switch) < 0)
  {
    return EXIT_FAILURE;
  }

  /* Set global edge flag for get_one/get_ion */
  xedge = edge_switch;

  printf ("windtable: Reading data from %s\n", root);

  /* Allocate domain structure */
  zdom = calloc (MAX_DOM, sizeof (domain_dummy));
  if (zdom == NULL)
  {
    Error ("windtable: Unable to allocate memory for domain\n");
    return EXIT_FAILURE;
  }

  /* Build wind save filename and read */
  sprintf (windsavefile, "%s.wind_save", root);

  if (wind_read (windsavefile) < 0)
  {
    Error ("windtable: Could not open %s\n", windsavefile);
    return EXIT_FAILURE;
  }

  printf ("windtable: Read wind_save file %s\n", windsavefile);
  printf ("windtable: Atomic data from %s\n", geo.atomic_filename);

  /* Get variable registry */
  vars = windtable_get_vars (&nvars);
  printf ("windtable: %d variables in registry\n", nvars);

  /* Generate ion columns from atomic data */
  if (windtable_generate_ion_columns (&ions, &nions) < 0)
  {
    Error ("windtable: Failed to generate ion columns\n");
    return EXIT_FAILURE;
  }
  printf ("windtable: %d ion fraction columns\n", nions);

  /* Write FITS table for each domain */
  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    if (geo.ndomain > 1)
    {
      sprintf (rootname, "%s.%d", root, ndom);
    }
    else
    {
      strcpy (rootname, root);
    }

    if (windtable_write_fits (ndom, rootname, vars, nvars, ions, nions) != 0)
    {
      Error ("windtable: Failed to write FITS file for domain %d\n", ndom);
    }
  }

  /* Cleanup */
  windtable_free_ion_columns (ions);

  return 0;
}
