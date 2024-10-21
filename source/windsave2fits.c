
/***********************************************************/
/** @file  windsave2fits.c 
 * @author ksl
 * @date   August,2024   
 *
 * @brief  Routines to write variious portions of  in a wind structure 
 * to a fits file
 *
 *###Notes###
 * This routine was written primarily to better understand
 * how well the spectra we use for estimating ionization rates
 * model the actual cell spectra; but it can be rather
 * straightforwardly adatpted to create images or table
 * of other portions of a windsave file, that involve
 * large amounts of data.
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"
#include "fitsio.h"


/* The first portions of the program contain some utility routines 
   for working with cfitsio */

char inroot[LINELENGTH], outroot[LINELENGTH], model_file[LINELENGTH], folder[LINELENGTH];
int model_flag, ksl_flag, cmf2obs_flag, obs2cmf_flag;

double line_matom_lum_single (double lum[], PlasmaPtr xplasma, int uplvl);
int line_matom_lum (int uplvl);
int create_matom_level_map ();

// Define a structure to hold spectral data
typedef struct
{
  int num_wavelengths;
  int num_spectra;
  float **data;                 // 2D array to hold spectra data
} Spectra;

/* Next functions are utilities that it should be possitle to call repeatedly*/

/**********************************************************/
/**
 * @brief      Function to perpare 2d image data from a 2d array
 *
 * @details The routine simply returns a flattened version of the 2d array
 * ### Notes ###
 *
 *
 **********************************************************/

float *
prepare_image_data (float **data, int width, int height)
{
  // Allocate memory for a contiguous 1D array to store the 2D image data
  float *image_data = calloc (width * height, sizeof (float));
  if (!image_data)
  {
    fprintf (stderr, "Memory allocation error\n");
    return NULL;
  }

  // Flatten the 2D data into a 1D array for FITS writing
  for (int i = 0; i < height; i++)
  {
    for (int j = 0; j < width; j++)
    {
      image_data[i * width + j] = data[i][j];
    }
  }

  return image_data;
}



/**********************************************************/
/**
 * @brief      Function to write a 2d image data to a fits file with
 * a ginve extension name
 *
 * @details The routine normally returns 0                                
 * ### Notes ###
 * This routine is called after prepare_image data, which will
 * have flattened and converted the original array to a float
 *
 *
 **********************************************************/

int
write_image_extension (fitsfile *fptr, float *image_data, int width, int height, const char *extname)
{
  int status = 0;               // CFITSIO status value MUST be initialized to zero
  long naxes[2];                // Array to store dimension sizes

  // Set the dimensions of the 2D image: [width, height]
  naxes[0] = width;             // Number of columns
  naxes[1] = height;            // Number of rows

  // Create a new image extension of type FLOAT_IMG
  if (fits_create_img (fptr, FLOAT_IMG, 2, naxes, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  // Set the extension name using the EXTNAME keyword
  if (fits_update_key (fptr, TSTRING, "EXTNAME", (void *) extname, NULL, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  // Write the image data to the FITS extension
  if (fits_write_img (fptr, TFLOAT, 1, width * height, image_data, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  return 0;
}




/**********************************************************/
/**
 * @brief      Function to write a 1d image extension with
 * a ginve extension name
 *
 * @details The routine normally returns 0                                   
 * ### Notes ###
 *
 *
 **********************************************************/

int
write_1d_image_extension (fitsfile *fptr, float *data, long length, const char *extname)
{
  int status = 0;               // CFITSIO status value MUST be initialized to zero
  long naxes[1];                // Array to store the dimension size (1D array)

  naxes[0] = length;            // Length of the 1D array

  // Create a new image extension of type FLOAT_IMG
  if (fits_create_img (fptr, FLOAT_IMG, 1, naxes, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  // Set the extension name using the EXTNAME keyword
  if (fits_update_key (fptr, TSTRING, "EXTNAME", (void *) extname, NULL, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  // Write the 1D array data to the FITS image extension
  if (fits_write_img (fptr, TFLOAT, 1, length, data, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  return 0;
}


/* This routine is not longer used but has been kept for reference
// Function to write a binary table with mixed data types to a FITS extension
int
write_mixed_table_extension (fitsfile *fptr, int num_rows, int num_columns, int *int_data, float *float_data, const char *extname)
{
  int status = 0;
  // Define the names, formats, and units for each column
  char *ttype[] = { "INT_COL", "FLOAT_COL" };   // Column names
  char *tform[] = { "J", "E" }; // Formats: 'J' for integer, 'E' for float
  char *tunit[] = { "", "" };   // Units

  // Create a new binary table extension
  if (fits_create_tbl (fptr, BINARY_TBL, num_rows, num_columns, ttype, tform, tunit, NULL, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  // Set the extension name using the EXTNAME keyword
  if (fits_update_key (fptr, TSTRING, "EXTNAME", (void *) extname, NULL, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  // Write the integer data to the first column
  if (fits_write_col (fptr, TINT, 1, 1, 1, num_rows, int_data, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  // Write the float data to the second column
  if (fits_write_col (fptr, TFLOAT, 2, 1, 1, num_rows, float_data, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  return 0;
}
*/

/* This ends the utility routines */

/**********************************************************/
/**
 * @brief      write_spectral_model_table extension to the
 * fits file
 *
 * @details The routine normally returns 0                                   
 * ### Notes ###
 *
 *
 **********************************************************/

int
write_spectra_model_table (fitsfile *fptr)
{

  /* We need one row for each cell and band so
     our table will be NPLASMA*nbands long
   */

  int i, j, k;
  int num_rows, num_cols;
  int nbands;
  int *ichoice, *iplasma, *iband, *nxtot;
  float *exp_w, *exp_temp, *pl_log_w, *pl_alpha;

  nbands = geo.nxfreq;
  num_rows = nbands * NPLASMA;
  printf ("xtest %d %d\n", nbands, NPLASMA);
  num_cols = 8;                 // for now
  char *table_name = "spec_model";

  ichoice = calloc (num_rows, sizeof (int *));
  iplasma = calloc (num_rows, sizeof (int *));
  iband = calloc (num_rows, sizeof (int *));
  exp_w = calloc (num_rows, sizeof (float *));
  exp_temp = calloc (num_rows, sizeof (float *));
  pl_log_w = calloc (num_rows, sizeof (float *));
  pl_alpha = calloc (num_rows, sizeof (float *));
  nxtot = calloc (num_rows, sizeof (int *));

  k = 0;
  for (i = 0; i < NPLASMA; i++)
  {
    for (j = 0; j < nbands; j++)
    {
      iplasma[k] = i;
      iband[k] = j;
      ichoice[k] = plasmamain[i].spec_mod_type[j];
      exp_w[k] = plasmamain[i].exp_w[j];
      exp_temp[k] = plasmamain[i].exp_temp[j];
      pl_log_w[k] = plasmamain[i].pl_log_w[j];
      pl_alpha[k] = plasmamain[i].pl_alpha[j];
      nxtot[k] = plasmamain[i].nxtot[j];
      printf ("nxtot %d\n", plasmamain[i].nxtot[j]);
      k++;
    }
  }


  int status = 0;
  // Define the names, formats, and units for each column
  char *ttype[] = { "nplasma", "band", "spec_mod_type", "exp_w", "exp_temp", "pl_log_w", "pl_alpha", "nxtot" }; // Column names
  char *tform[] = { "J", "J", "J", "E", "E", "E", "E", "J" };   // Formats: 'J' for integer, 'E' for float
  char *tunit[] = { "", "", "", "", "", "", "", "" };   // Units

  // Create a new binary table extension
  if (fits_create_tbl (fptr, BINARY_TBL, num_rows, num_cols, ttype, tform, tunit, NULL, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  // Set the extension name using the EXTNAME keyword
  if (fits_update_key (fptr, TSTRING, "EXTNAME", (void *) table_name, NULL, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  printf ("num_rows: %d\n", num_rows);
  for (i = 0; i < num_rows; i++)
  {
    printf ("%d %d %d\n", i, iplasma[i], ichoice[i]);
  }

  // Write the integer data to the first column
  if (fits_write_col (fptr, TINT, 1, 1, 1, num_rows, iplasma, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }


  // Write the float data to the second column
  if (fits_write_col (fptr, TINT, 2, 1, 1, num_rows, iband, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  // Write the float data to the third  column
  if (fits_write_col (fptr, TINT, 3, 1, 1, num_rows, ichoice, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }


  // Write the float data to the fourth column
  if (fits_write_col (fptr, TFLOAT, 4, 1, 1, num_rows, exp_w, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }



  // Write the float data to the fifth columnmn
  if (fits_write_col (fptr, TFLOAT, 5, 1, 1, num_rows, exp_temp, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }


  // Write the float data to the sixth  column
  if (fits_write_col (fptr, TFLOAT, 6, 1, 1, num_rows, pl_log_w, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }


  // Write the float data to the seventh column
  if (fits_write_col (fptr, TFLOAT, 7, 1, 1, num_rows, pl_alpha, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  // Write the integer data to the eightth column
  if (fits_write_col (fptr, TINT, 8, 1, 1, num_rows, nxtot, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }


  return 0;
}




/**********************************************************/
/**
 * @brief      write_a fits file that contains information
 * releveant to how well the cell spectra are modelled
 * fits file
 *
 * @details The routine normally returns 0                                   
 * ### Notes ###
 *
 * This is the main routine, in the sense that everyting is
 * controlled frrom here
 *
 **********************************************************/



int
make_spec (inroot)
     char *inroot;
{

  fitsfile *fptr;               // Pointer to the FITS file
  int status = 0;               // CFITSIO status value MUST be initialized to zero
  char outfile[LINELENGTH];

  /* The obscure exclamation mark at the beginning of the name means yu want to overwrie 
     and existing file, if it exists */

  sprintf (outfile, "!%s_cellspec.fits", inroot);
  printf ("outfile %s\n", outfile);



  // Create a new FITS file
  //if (fits_create_file (&fptr, "!multi_extension.fits", &status))
  if (fits_create_file (&fptr, outfile, &status))
  {
    fits_report_error (stderr, status); // Print any error message
    return status;
  }

// Create an empty primary HDU
  if (fits_create_img (fptr, FLOAT_IMG, 0, NULL, &status))
  {
    fits_report_error (stderr, status);
  }

  printf ("Hello World %s \n", inroot);
  printf ("Plasma  %d \n", NPLASMA);
  printf ("NBINS in spec %d \n", NBINS_IN_CELL_SPEC);


  Spectra spectra;
  spectra.num_spectra = NPLASMA;
  spectra.num_wavelengths = NBINS_IN_CELL_SPEC;
  spectra.data = calloc (spectra.num_spectra, sizeof (float *));

  for (int i = 0; i < spectra.num_spectra; i++)
  {
    spectra.data[i] = calloc (spectra.num_wavelengths, sizeof (float));

    // Fill with sample data
    for (int j = 0; j < spectra.num_wavelengths; j++)
    {
      // spectra.data[i][j] = (float) (i + j);
      spectra.data[i][j] = (float) plasmamain[i].cell_spec_flux[j];
    }
  }

  // Prepare image data for the first extension
  float *image_data1 = prepare_image_data (spectra.data, spectra.num_wavelengths, spectra.num_spectra);
  if (!image_data1)
  {
    // Free allocated memory before exiting
    for (int i = 0; i < spectra.num_spectra; i++)
    {
      free (spectra.data[i]);
    }
    free (spectra.data);
    return 1;
  }


  status = write_image_extension (fptr, image_data1, spectra.num_wavelengths, spectra.num_spectra, "Jnu");

  Spectra freq;
  freq.num_spectra = 1;
  freq.num_wavelengths = NBINS_IN_CELL_SPEC;
  freq.data = calloc (freq.num_spectra, sizeof (float *));
  freq.data[0] = calloc (freq.num_wavelengths, sizeof (float));

  for (int i = 0; i < NBINS_IN_CELL_SPEC; i++)
  {
    freq.data[0][i] = (float) pow (10., geo.cell_log_freq_min + i * geo.cell_delta_lfreq);
  }


  float *image_data2 = prepare_image_data (freq.data, NBINS_IN_CELL_SPEC, 1);
/*
  double freq[2000][1];
  int i;

  for (i = 0; i < NBINS_IN_CELL_SPEC; i++)
  {
    freq[i][0] = (float) pow (10., geo.cell_log_freq_min + i * geo.cell_delta_lfreq);
  }

  float *image_data2 = prepare_image_data (freq, NBINS_IN_CELL_SPEC, 1);
*/
  // status = write_image_extension (fptr, image_data2, NBINS_IN_CELL_SPEC, 1, "nu");
  status = write_1d_image_extension (fptr, image_data2, NBINS_IN_CELL_SPEC, "nu");

  /* Elimainate this for now

     // Prepare data for the binary table
     int num_rows = 5;
     int int_data[] = { 1, 2, 3, 4, 5 };   // Example integer data
     float float_data[] = { 0.1, 0.2, 0.3, 0.4, 0.5 };     // Example float data

     write_mixed_table_extension (fptr, num_rows, 2, int_data, float_data, "Table_Ext");
   */


/* Now create a 1d extension hold the frequencies of the broad bands */

  int i;
  float xfreq[100];

  for (i = 0; i <= geo.nxfreq; i++)
  {
    xfreq[i] = geo.xfreq[i];
  }

  status = write_1d_image_extension (fptr, xfreq, geo.nxfreq + 1, "nu_model");

/* Now write my table */

  write_spectra_model_table (fptr);

  if (fits_close_file (fptr, &status))
  {
    fits_report_error (stderr, status);
  }

  return 0;
}

/**********************************************************/
/**
 * @brief      parses the command line options
 *
 * @param [in]  int  argc   the number of command line arguments
 * @param [in]  char *  argv[]   The command line arguments
 *
 *
 * ###Notes###
 *
 * The general purpose of each of the command line options
 * should be fairly obvious from reading the code.
 *
 *
 * Although this routine uses the standard Log and Error commands
 * the diag files have not been created yet and so this information
 * is really simply written to the terminal.
 * 
 * WARNING: this has not been updeated. Currently there
 * routine expects only the rootname of a windsave file
 *
 *
 **********************************************************/

int
xparse_command_line (argc, argv)
     int argc;
     char *argv[];
{
  int j = 0;
  int i;
  char dummy[LINELENGTH];
  int mkdir ();
  char *fgets_rc;


  sprintf (outroot, "%s", "new");

  model_flag = ksl_flag = obs2cmf_flag = cmf2obs_flag = 0;

  if (argc == 1)
  {
    printf ("Parameter file name (e.g. my_model.pf, or just my_model):");
    fgets_rc = fgets (dummy, LINELENGTH, stdin);
    if (!fgets_rc)
    {
      printf ("Input rootname is NULL or invalid\n");
      exit (1);
    }
    get_root (inroot, dummy);
  }
  else
  {

    for (i = 1; i < argc; i++)
    {
      if (strcmp (argv[i], "-out_root") == 0)
      {
        if (sscanf (argv[i + 1], "%s", dummy) != 1)
        {
          printf ("sirocco: Expected out_root after -out_root switch\n");
          exit (0);
        }

        get_root (outroot, dummy);
        i++;
        j = i;

      }
      if (strcmp (argv[i], "-model_file") == 0)
      {
        if (sscanf (argv[i + 1], "%s", dummy) != 1)
        {
          printf ("sirocco: Expected a model file containing density, velocity and temperature after -model_file switch\n");
          exit (0);
        }
        get_root (model_file, dummy);
        i++;
        j = i;
        printf ("got a model file %s\n", model_file);
        model_flag = 1;
      }
      else if (strcmp (argv[i], "-ksl") == 0)
      {
        printf ("Carrying out a simple hard wired ion modification\n");
        ksl_flag = 1;
      }
      else if (strcmp (argv[i], "--dry-run") == 0)
      {
        modes.quit_after_inputs = 1;
        j = i;
      }
      else if (strcmp (argv[i], "-cmf") == 0)
      {
        obs2cmf_flag = 1;
      }
      else if (strcmp (argv[i], "-obs") == 0)
      {
        cmf2obs_flag = 1;
      }
      else if (strncmp (argv[i], "-", 1) == 0)
      {
        printf ("sirocco: Unknown switch %s\n", argv[i]);
        exit (0);
      }
    }

    /* The last command line variable is always the windsave file */

    if (j + 1 == argc)
    {
      printf ("All of the command line has been consumed without specifying a parameter file name, so exiting\n");
      exit (1);
    }
    strcpy (dummy, argv[argc - 1]);
    get_root (inroot, dummy);

  }

  return (0);
}



/**********************************************************/
/**
 * @brief      the main routine which carries out the effort
 *
 * @param [in]  int  argc   the number of command line arguments
 * @param [in]  char *  argv[]   The command line arguments
 *
 *
 * ###Notes###
 *
 * This routine oversees the effort.  The basic steps are
 *
 * - parse the command line to get the names of files
 * - call the routine to craetate fits file
 *
 *
 **********************************************************/



int
main (argc, argv)
     int argc;
     char *argv[];
{

  char infile[LINELENGTH], outfile[LINELENGTH];
  FILE *fopen ();
  int mkdir ();


  xparse_command_line (argc, argv);

  sprintf (infile, "%.150s.wind_save", inroot);
  sprintf (outfile, "%.150s.txt", inroot);

  printf ("Reading %s and writing to %s\n", infile, outfile);

  zdom = calloc (MAX_DOM, sizeof (domain_dummy));
  if (zdom == NULL)
  {
    printf ("Unable to allocate memory for domain\n");
    return EXIT_FAILURE;
  }

  wind_read (infile);

  make_spec (inroot);

  return (0);
}
