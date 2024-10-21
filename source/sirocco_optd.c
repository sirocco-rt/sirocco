/** ************************************************************************* */
/**
 * @file  py_optical_depth.c
 * @author  Edward Parkinson
 * @date  February 2021
 *
 * @brief  File containing the main functions defining the program.
 *
 * ************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

#include "atomic.h"
#include "sirocco.h"
#include "sirocco_optd.h"

struct CommandlineArguments
{
  double freq_min;
  double freq_max;
  double n_freq;
  double inclinations[MAX_CUSTOM_ANGLES];
};

/* ************************************************************************* */
/**
 * @brief  Create spectra of tau vs lambda for each observer angle
 *
 * @details
 *
 * This is the main function which will generate the optical depth spectra for
 * each observer angle in xxspec. The algorithm is similar to extract and the
 * tau_diag algorithm which this function is called in.
 *
 * A photon is generated at the central source of the model and is extracted
 * from this location towards the observer where it escapes, where integrate_tau_across_wind
 * returns the integrated optical depth along its path to escape. This is done
 * for a range of photon frequencies to find how optical depth changes with
 * frequency.
 *
 * This processes can take some time compared to tau_evalulate_photo_edges. But,
 * since NUM_FREQUENCY_BINS photons are being generated for each spectrum and the fact
 * that these photons do not interact, the spectra does not actually take that
 * long to complete.
 *
 * ************************************************************************** */

void
create_optical_depth_spectrum (double u_freq_min, double u_freq_max, double *input_inclinations)
{
  int i, j;
  int err;
  double *tau_spectrum;
  double c_optical_depth, c_column_density;
  double c_frequency, freq_min, freq_max, d_freq;
  struct photon photon;

  int n_inclinations;
  SightLines_t *inclinations = initialize_inclination_angles (&n_inclinations, input_inclinations);

  printf ("Creating optical depth spectra:\n");

  tau_spectrum = calloc (n_inclinations * NUM_FREQUENCY_BINS, sizeof *tau_spectrum);
  if (tau_spectrum == NULL)
  {
    errormsg ("cannot allocate %lu bytes for tau_spectrum\n", n_inclinations * NUM_FREQUENCY_BINS * sizeof *tau_spectrum);
    exit (EXIT_FAILURE);
  }

  /*
   * We have a complicated if statement first, though. If a freq_min
   * or a freq_max was provided, then we need to get this first and set
   * the frequency limits appropriately. If neither are defined, then we will
   * use some hardwired limits. The frequency range of the extracted will be
   * used, however if xxpsec is NULL (no observer spectrum exists), then the
   * frequency range will be over a default 100 - 10,000 Angstrom band.
   */

  if (u_freq_min > 0 || u_freq_max > 0)
  {
    if (u_freq_min > 0)
    {
      freq_min = u_freq_min;
    }
    else
    {
      freq_min = VLIGHT / (10000 * ANGSTROM);
    }

    if (u_freq_max > 0)
    {
      freq_max = u_freq_max;
    }
    else
    {
      freq_max = VLIGHT / (100 * ANGSTROM);
    }

    if (freq_max < freq_min)
    {
      errormsg ("frequency range given has set freq_max (%e) < freq_min (%e) \n", freq_max, freq_min);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    if ((geo.nangles == 0 && xxspec == NULL) || (geo.swavemax == 0 && geo.swavemin == 0))
    {
      freq_min = VLIGHT / (10000 * ANGSTROM);
      freq_max = VLIGHT / (100 * ANGSTROM);
    }
    else
    {
      freq_min = VLIGHT / (geo.swavemax * ANGSTROM);
      freq_max = VLIGHT / (geo.swavemin * ANGSTROM);
      if (sane_check (freq_min))
      {
        freq_min = VLIGHT / (10000 * ANGSTROM);
        errormsg ("freq_min has an invalid value setting to %e\n", freq_min);
      }
      if (sane_check (freq_max))
      {
        freq_max = VLIGHT / (100 * ANGSTROM);
        errormsg ("freq_min has an invalid value setting to %e\n", freq_max);
      }
    }
  }

  d_freq = (log10 (freq_max) - log10 (freq_min)) / NUM_FREQUENCY_BINS;
  kbf_need (freq_min, freq_max);

  /*
   * Now create the optical depth spectra for each inclination
   */

  for (i = 0; i < n_inclinations; i++)
  {
    printf ("  - Creating spectrum: %s\n", inclinations[i].name);
    c_frequency = log10 (freq_min);

    for (j = 0; j < NUM_FREQUENCY_BINS; j++)
    {
      c_optical_depth = 0.0;
      c_column_density = 0.0;

      err = create_photon (&photon, pow (10, c_frequency), inclinations[i].lmn);
      if (err == EXIT_FAILURE)
      {
        errormsg ("skipping photon of frequency %e\n", pow (10, c_frequency));
        continue;
      }

      err = integrate_tau_across_wind (&photon, &c_column_density, &c_optical_depth);
      if (err == EXIT_FAILURE)
        continue;

      tau_spectrum[i * NUM_FREQUENCY_BINS + j] = c_optical_depth;
      c_frequency += d_freq;
    }
  }

  write_optical_depth_spectrum (inclinations, n_inclinations, tau_spectrum, freq_min, d_freq);
  free (tau_spectrum);
  free (inclinations);
}

/* ************************************************************************* */
/**
 * @brief Calculate the optical depth for various optical depth edges and
 *        extract the column density.
 *
 * @details
 *
 * This is the main function which will control the procedure for calculating
 * various diagnostic numbers for the optical depth's experienced in the current
 * model. Namely, this function aims to show the total integrated optical depth
 * to each observer angle using (originally) the following optical depths:
 *
 *  - Lymann edge
 *  - Balmer edge
 *  - Helium II edge
 *
 * Once these integrated optical depths have been calculated for each angle, a
 * spectrum of optical depth vs wavelength is created for each angle.
 *
 * The aim of these diagnostic numbers it to provide some sort of quick metric
 * on the optical thickness of the current model.
 *
 * ************************************************************************** */

void
evaluate_photoionization_edges (double *input_inclinations)
{
  int i, j;
  int err;
  double c_frequency, c_optical_depth, c_column_density;
  double *optical_depth_values = NULL, *column_density_values = NULL;
  struct photon photon;
  enum RunModeEnum original_run_mode = RUN_MODE;
  RUN_MODE = RUN_MODE_NO_ES_OPACITY;

  Edges_t edges[] = {
    {"HLymanEdge", 3.387485e+15},
    {"HBalmerEdge", 8.293014e+14},
    {"HeI24eVEdge", 5.9483e+15},
    {"HeII54eVEdge", 1.394384e+16}
  };

  const int n_edges = sizeof edges / sizeof edges[0];

  int n_inclinations;
  SightLines_t *inclinations = initialize_inclination_angles (&n_inclinations, input_inclinations);

  optical_depth_values = calloc (n_inclinations * n_edges, sizeof *optical_depth_values);
  if (optical_depth_values == NULL)
  {
    errormsg ("cannot allocate %lu bytes for optical_depths\n", n_inclinations * n_edges * sizeof *optical_depth_values);
    exit (EXIT_FAILURE);
  }

  column_density_values = calloc (n_inclinations, sizeof *column_density_values);
  if (column_density_values == NULL)
  {
    errormsg ("cannot allocate %lu bytes for column_densities\n", n_inclinations * sizeof *column_density_values);
    exit (EXIT_FAILURE);
  }

  /*
   * Now extract the optical depths and mass column densities. We loop over
   * each PI edge for each inclination angle.
   */

  for (i = 0; i < n_inclinations; i++)
  {
    for (j = 0; j < n_edges; j++)
    {
      c_optical_depth = 0.0;
      c_column_density = 0.0;
      c_frequency = edges[j].freq;

      err = create_photon (&photon, c_frequency, inclinations[i].lmn);
      if (err == EXIT_FAILURE)
      {
        errormsg ("skipping photon of frequency %e\n", c_frequency);
        continue;
      }

      err = integrate_tau_across_wind (&photon, &c_column_density, &c_optical_depth);
      if (err == EXIT_FAILURE)
        continue;               // do not throw extra warning when one is already thrown in integrate_tau_across_wind

      optical_depth_values[i * n_edges + j] = c_optical_depth;
      column_density_values[i] = c_column_density;
    }
  }

  print_optical_depths (inclinations, n_inclinations, edges, n_edges, optical_depth_values, column_density_values);
  free (inclinations);
  free (optical_depth_values);
  free (column_density_values);
  RUN_MODE = original_run_mode;
}

/* ************************************************************************* */
/**
 * @brief   Find the electron scattering photosphere.
 *
 * @details
 *
 * This is the main controlling function for finding the photosphere. The
 * electron scattering optical depth is controlled by the global variable
 * TAU_DEPTH.
 *
 * ************************************************************************** */

void
find_photosphere (void)
{
  int i, err;
  double optical_depth, column_density;
  struct photon photon;
  SightLines_t *inclinations;
  double input_inclinations[MAX_CUSTOM_ANGLES];

  for (i = 0; i < MAX_CUSTOM_ANGLES; ++i)       // required for initialize_inclination_angles, even though unused
    input_inclinations[i] = -1.0;

  int n_inclinations;
  inclinations = initialize_inclination_angles (&n_inclinations, input_inclinations);

  Positions_t *positions = calloc (n_inclinations, sizeof (Positions_t));
  if (positions == NULL)
  {
    errormsg ("unable to allocate memory for the positions array\n");
    exit (EXIT_FAILURE);
  }

  const double test_freq = 8e14;        // todo: this probably need to be a possible input

  printf ("Locating electron scattering photosphere surface for tau_es = %f\n", TAU_DEPTH);

  for (i = 0; i < n_inclinations; i++)
  {
    err = create_photon (&photon, test_freq, inclinations[i].lmn);
    if (err)
    {
      positions[i].x = positions[i].y = positions[i].z = -1.0;
      continue;
    }

    optical_depth = column_density = 0;

    err = integrate_tau_across_wind (&photon, &column_density, &optical_depth);
    if (err)
    {
      positions[i].x = positions[i].y = positions[i].z = -1.0;
      continue;
    }

    positions[i].x = photon.x[0];
    positions[i].y = photon.x[1];
    positions[i].z = photon.x[2];
  }

  write_photosphere_location_to_file (positions, n_inclinations);
  free (inclinations);
  free (positions);
}

/* ************************************************************************* */
/**
 * @brief  Help text for the program.
 *
 * ************************************************************************** */

void
print_help (void)
{
  char *help =
    "A utility program to analyse the optical depth in a Python model.\n\n"
    "usage: py_optical_depth [-h] [-d ndom] [-p tau_stop] [-cion nion] \n"
    "                        [-freq_min min] [-freq_max max] [-i i1 i2 ...]\n"
    "                        [--nonrel] [--smax frac] [--no-es] [--version]\n"
    "                        root\n\n"
    "This program can be used in multiple ways. By default, the integrated continuum\n"
    "optical depth along the defined observer lines of sight for the model are returned.\n"
    "If none of these have been defined, then a set of default lines of sight are used\n"
    "instead or can be specified using -i. This program can also find the surface of\n"
    "constant electron scattering optical using the -p option.\n\n"
    "Please see below for a list of all flags.\n\n"
    "-h               Print this help message and exit\n"
    "-d ndom          Set the domain to launch photons from\n"
    "-p tau_stop      Integrate from outwards to find the electron scattering photosphere\n"
    "-cion nion       Extract the column density for an ion of number nion\n"
    "-freq_min min    The lower frequency boundary for optical depth spectra\n"
    "-freq_max max    The upper frequency boundary for optical depth spectra\n"
    "-i i1 i2 i3 ...  Calculate the optical depth for the given space seperated list of sight lines\n"
    "--nonrel         Use linear frequency transforms, to be used when Python was run\n"
    "                 in non-relativistic mode\n"
    "--smax frac      Set the maximum fraction a photon can move in terms of cell distances\n"
    "--no-es          Do not include opacity contributions from electron scattering\n"
    "--version        Print the version information and exit.\n";
  printf ("%s", help);
}

/* ************************************************************************* */
/**
 * @brief  Parse run time arguments provided at the command line.
 *
 * @param[in]  argc  The number of arguments provided
 * @param[in]  argv  The command line arguments
 *
 * @details
 *
 * Reads the command line arguments. Assumes the final argument is the root name
 * of the model. If no arguments are provided, then the program will run in a
 * default mode and query the user for the root name of the simulation.
 *
 * ************************************************************************** */

struct CommandlineArguments
get_arguments (int argc, char *argv[])
{
  int i;
  int n_args_read;
  char input[LINELENGTH];
  struct CommandlineArguments arguments;

  arguments.freq_min = arguments.freq_max = 0;
  arguments.n_freq = NUM_FREQUENCY_BINS;
  for (i = 0; i < MAX_CUSTOM_ANGLES; ++i)
    arguments.inclinations[i] = -1.0;

  /*
   * If no command line arguments are provided, then use fgets to query the
   * root name from the user and return
   */

  if (argc == 1)
  {
    printf ("Please input the model root name: ");
    if (fgets (input, LINELENGTH, stdin) == NULL)
    {
      printf ("Unable to parse root name\n");
      exit (EXIT_FAILURE);
    }
    get_root (files.root, input);
    return arguments;
  }

  /*
   * Otherwise, loop over the command line arguments and initialize various
   * variables
   */

  n_args_read = 0;
  for (i = 1; i < argc; ++i)
  {
    if (!strcmp (argv[i], "-h"))        //NOTE: print help text
    {
      print_help ();
      exit (EXIT_SUCCESS);
    }
    else if (!strcmp (argv[i], "--nonrel"))    //NOTE: use linear frequency transforms
    {
      printf ("Using linear approximations for Doppler shifts\n");
      rel_mode = REL_MODE_LINEAR;
      n_args_read = i;
    }
    else if (!strcmp (argv[i], "--version"))    //NOTE: print version number
    {
      printf ("Python version %s\n", VERSION);
      printf ("Built from git commit %s\n", GIT_COMMIT_HASH);
      if (GIT_DIFF_STATUS)
        printf ("This version was compiled with %d files with uncommited changes\n", GIT_DIFF_STATUS);
      exit (EXIT_SUCCESS);
    }
    else if (!strcmp (argv[i], "-cion"))        //NOTE: extract column density for specific ion, rather than H
    {
      char *check;
      COLUMN_MODE = COLUMN_MODE_ION;
      COLUMN_MODE_ION_NUMBER = (int) strtol (argv[i + 1], &check, 10);
      if (*check != '\0')
      {
        printf ("Unable to convert argument provided for -cion into an integer\n");
        exit (EXIT_FAILURE);
      }
      if (COLUMN_MODE_ION_NUMBER < 0)
      {
        printf ("Argument for -cion cannot be negative\n");
        exit (EXIT_FAILURE);
      }
      n_args_read = i++;
    }
    else if (!strcmp (argv[i], "-i"))
    {
      int k;
      char *stdtod_check;
      int n_angles_provided = 0;

      for (k = 1; k < MAX_CUSTOM_ANGLES + 2; k++)
      {
        double inclination = strtod (argv[i + k], &stdtod_check);
        if (inclination < 0 || inclination > 90)
        {
          printf ("Invalid inclination angle of %f. Has to be 0 < inclination < 90\n", inclination);
          exit (EXIT_FAILURE);
        }
        if (k > MAX_CUSTOM_ANGLES)
        {
          errormsg ("Maximum number of inclination angles provided: only %d are allowed\n", MAX_CUSTOM_ANGLES);
          exit (EXIT_FAILURE);
        }
        if (strcmp (argv[i + k], stdtod_check) == 0)    // conversion error indicates no more numbers
        {
          break;
        }
        arguments.inclinations[k - 1] = inclination;
        n_angles_provided += 1;
      }
      if (n_angles_provided == 0)
      {
        errormsg ("Invalid input: no inclination were provided for -i option\n");
        exit (EXIT_FAILURE);
      }

      n_args_read = i++;
      i += n_angles_provided - 1;
    }
    else if (!strcmp (argv[i], "-p"))   //NOTE: run in photosphere mode
    {
      RUN_MODE = RUN_MODE_ES_PHOTOSPHERE;
      char *check;
      TAU_DEPTH = strtod (argv[i + 1], &check);
      if (*check != '\0')
      {
        printf ("Unable to convert argument provided for -p into a double\n");
        exit (EXIT_FAILURE);
      }
      n_args_read = i++;
    }
    else if (!strcmp (argv[i], "--smax"))
    {
      char *check;
      SMAX_FRAC = strtod (argv[i + 1], &check);
      if (*check != '\0')
      {
        printf ("Unable to convert argument for -smax into a double");
        exit (EXIT_FAILURE);
      }
      n_args_read = i++;
    }
    else if (!strcmp (argv[i], "-d"))   //NOTE: change the lanching domain for photons
    {
      char *check;
      N_DOMAIN = (int) strtol (argv[i + 1], &check, 10);
      if (*check != '\0')
      {
        printf ("Unable to convert argument provided for -d into an integer\n");
        exit (EXIT_FAILURE);
      }
      n_args_read = i++;
    }
    else if (!strcmp (argv[i], "-freq-min"))    //NOTE: lower frequency boundary for optical depth spectrum
    {
      char *check;
      arguments.freq_min = strtod (argv[i + 1], &check);
      if (*check != '\0')
      {
        printf ("Unable to convert argument provided for -freq_min to a double\n");
        exit (EXIT_FAILURE);
      }
      n_args_read = i++;
    }
    else if (!strcmp (argv[i], "-freq-max"))    //NOTE: upper frequency boundary for optical depth spectrum
    {
      char *check;
      arguments.freq_max = strtod (argv[i + 1], &check);
      if (*check != '\0')
      {
        printf ("Unable to convert argument provided for -freq_max to a double\n");
        exit (EXIT_FAILURE);
      }
      n_args_read = i++;
    }
    else if (!strcmp (argv[i], "--no-es"))
    {
      RUN_MODE = RUN_MODE_NO_ES_OPACITY;
      n_args_read = i++;
    }
    else if (!strncmp (argv[i], "-", 1))
    {
      printf ("Unknown argument %s\n", argv[i]);
      print_help ();
      exit (EXIT_FAILURE);
    }
  }

  if (n_args_read + 1 == argc)
  {
    printf ("All command line arguments have been processed and did not find a root name\n");
    exit (EXIT_FAILURE);
  }

  get_root (files.root, argv[argc - 1]);

  return arguments;
}

/* ************************************************************************* */
/**
 * @brief  The main function of the program.
 *
 * @param  argc  The number of command line arguments
 * @param  argv  The command line arguments
 *
 * @return  EXIT_SUCCESS
 *
 * @details
 *
 * ************************************************************************** */

int
main (int argc, char *argv[])
{
  char windsave_filename[LINELENGTH + 24];
  char specsave_filename[LINELENGTH + 24];

  timer ();

  /*
   * Initialize global variables which are required for photon
   * transport and general program operation
   */

  Log_set_verbosity (2);
  Log_print_max (10);
  Log_quit_after_n_errors ((int) 1e8);
  init_rand ((int) time (NULL));

  rel_mode = REL_MODE_FULL;     // this is updated in get_arguments if required
  SMAX_FRAC = 1.0;
  DENSITY_PHOT_MIN = 1.e-10;
  COLUMN_MODE = COLUMN_MODE_RHO;
  RUN_MODE = RUN_MODE_TAU_INTEGRATE;
  N_DOMAIN = 0;

  struct CommandlineArguments arguments;
  arguments = get_arguments (argc, argv);
  printf ("%-20s Optical depth diagnostics beginning\n", "TAU");
  strcpy (windsave_filename, files.root);
  strcat (windsave_filename, ".wind_save");
  strcpy (specsave_filename, files.root);
  strcat (specsave_filename, ".spec_save");

  /*
   * Read in the wind_save file and initialize the wind cones and DFUDGE which
   * are important for photon transport. The atomic data is also read in at
   * this point (which is also very important)
   */

  zdom = calloc (MAX_DOM, sizeof (domain_dummy));
  if (zdom == NULL)
  {
    printf ("Unable to allocate memory for domain\n");
    return EXIT_FAILURE;
  }

  if (wind_read (windsave_filename) < 0)
  {
    errormsg ("unable to open %s\n", windsave_filename);
    exit (EXIT_FAILURE);
  }

  DFUDGE = setup_dfudge ();
  setup_windcone ();

  /*
   * If a spec_save exists, and there are spectral cycles (possibly a redundant
   * check), then read in the spec_save file.
   */

  if (access (specsave_filename, F_OK) == 0)
  {
    if (geo.pcycle > 0)
    {
      if (spec_read (specsave_filename) < 0)
      {
        errormsg ("Unable to open %s when spectrum cycles have been run\n", specsave_filename);
        exit (EXIT_FAILURE);
      }
    }
  }

  /*
   * Do some further error checking on the COLUMN_MODE_ION to ensure that
   * it is a sensible number and print the ion to be extracted
   */

  printf ("\n");
  if (RUN_MODE != RUN_MODE_ES_PHOTOSPHERE)
  {
    if (COLUMN_MODE == COLUMN_MODE_ION)
    {
      if (COLUMN_MODE_ION_NUMBER < 0 || COLUMN_MODE_ION_NUMBER > nions - 1)
      {
        errormsg ("The ion number %i is an invalid ion number: there are %i ions which have been loaded\n", COLUMN_MODE_ION_NUMBER, nions);
        exit (EXIT_FAILURE);
      }

      printf ("Extracting column density for %s %i\n", ele[ion[COLUMN_MODE_ION_NUMBER].nelem].name, ion[COLUMN_MODE_ION_NUMBER].istate);
    }
    else
    {
      printf ("Extracting mass and Hydrogen column density\n");
    }
  }
  else
  {
    printf ("Finding the photosphere location for electron scattering optical depth %f\n", TAU_DEPTH);
  }

  /*
   * Now we can finally being the optical depth diagnostics... first we close
   * the log and re-open it, with a verbosity of 5. This clears the errors due
   * to atomic data, which we should not need to worry about for this program
   */

  if (RUN_MODE == RUN_MODE_TAU_INTEGRATE || RUN_MODE == RUN_MODE_NO_ES_OPACITY)
  {
    evaluate_photoionization_edges (arguments.inclinations);
    create_optical_depth_spectrum (arguments.freq_min, arguments.freq_max, arguments.inclinations);
  }
  else if (RUN_MODE == RUN_MODE_ES_PHOTOSPHERE)
  {
    find_photosphere ();
  }
  else
  {
    errormsg ("Mode %d is an unknown run mode, not sure how you got here so exiting the program\n", RUN_MODE);
    exit (EXIT_FAILURE);
  }

  printf ("\n%-20s Optical depth diagnostics completed\n", "TAU");
  printf ("Completed optical depth diagnostics. The elapsed TIME was %f\n", timer ());

  error_summary ("end of program");

  return EXIT_SUCCESS;
}
