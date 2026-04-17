/** ********************************************************************************************************************
 *
 *  @file load_model.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date Jan 2024
 *
 *  @brief
 *
 * ****************************************************************************************************************** */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "../atomic.h"
#include "../sirocco.h"

/* rdpar_stat is defined in rdpar.c; reset it to 0 to allow opar to open a new
 * file even if a previous test left it in the "file open" state (rdpar_stat==2). */
extern int rdpar_stat;

#define PATH_SEPARATOR '/'

/** *******************************************************************************************************************
 *
 * @brief
 *
 * @return
 *
 * @details
 *
 *
 * ****************************************************************************************************************** */

const char *
get_sirocco_env_variable (void)
{
  const char *env = getenv ("SIROCCO");
  if (env == NULL)
  {
    fprintf (stderr, "Failed to find SIROCCO environment variable\n");
    return NULL;
  }

  return env;
}

/** *******************************************************************************************************************
 *
 * @brief Free a pointer and set to NULL
 *
 * @param [in]  void** ptr  An address to the pointer to free and set to null
 *
 * @details
 *
 * To use this function, you need to pass the address of a pointer cast as void**, e.g.
 *
 *      free_and_null((void**) &wmain);
 *
 * ****************************************************************************************************************** */

static void
free_and_null (void **ptr)
{
  if (ptr != NULL && *ptr != NULL)
  {
    free (*ptr);
    *ptr = NULL;
  }
}

/** *******************************************************************************************************************
 *
 * @brief Clean up after a model has run as a test case.
 *
 * @param [in]  root_name  The root name of the model in $SIROCCO/source/tests/test_data/define_wind
 *
 * @details
 *
 * ****************************************************************************************************************** */

int
cleanup_model (const char *root_name)
{
  char *SIROCCO_ENV;
  char parameter_filepath[LINELENGTH];

  (void) root_name;

  SIROCCO_ENV = getenv ("SIROCCO");
  if (SIROCCO_ENV == NULL)
  {
    fprintf (stderr, "Failed to find SIROCCO environment variable\n");
    return EXIT_FAILURE;
  }

  snprintf (parameter_filepath, LINELENGTH, "%s/source/tests/test_data/define_wind/%s.pf", SIROCCO_ENV, files.root);
  if (cpar (parameter_filepath) != 1)   /* cpar returns 1 when something is "normal" */
  {
    return EXIT_FAILURE;
  }

  free_domains ();
  free_wind_grid ();
  free_plasma_grid ();
  free_and_null ((void **) &photstoremain);
  free_and_null ((void **) &matomphotstoremain);

  if (nlevels_macro > 0)
  {
    free_macro_grid ();
  }

  NDIM2 = 0;
  NPLASMA = 0;                  /* We probably don't need to do this, but better safe than sorry */

  /* free atomic data elements */
  free_and_null ((void **) &ele);
  free_and_null ((void **) &ion);
  free_and_null ((void **) &xconfig);
  free_and_null ((void **) &line);
  free_and_null ((void **) &auger_macro);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Get the last component in a file path.
 *
 * @details
 *
 * This will return the file name and extension within a file path.
 *
 * ****************************************************************************************************************** */

const char *
get_last_component (const char *path)
{
  const char *last_separator = strrchr (path, PATH_SEPARATOR);
  return last_separator ? last_separator + 1 : path;
}

/** *******************************************************************************************************************
 *
 * @brief Find the location of the atomic data, relative to a model.
 *
 * @details
 *
 * This function exists to modify the atomic data in the parameter file to make it an absolute path to where
 * it has been set at suite initialisation.
 *
 * ****************************************************************************************************************** */

int
set_atomic_data_filename (const char *atomic_data_location)
{
  char temp_filepath[LINELENGTH];
  char atomic_data_filepath[LINELENGTH];

  if (strlen (atomic_data_location) + strlen (get_last_component (geo.atomic_filename)) >= LINELENGTH)
  {
    perror ("New atomic data filepath will be too long for internal buffer");
    return EXIT_FAILURE;
  }

  strncpy (temp_filepath, atomic_data_location, LINELENGTH - 1);
  temp_filepath[LINELENGTH - 1] = '\0';

  size_t len1 = strlen (temp_filepath);
  if (len1 > 0 && temp_filepath[len1 - 1] != PATH_SEPARATOR)
  {
    strncat (temp_filepath, "/", LINELENGTH - len1 - 1);
  }
  snprintf (atomic_data_filepath, LINELENGTH, "%s%s", temp_filepath, get_last_component (geo.atomic_filename));

  if (strlen (atomic_data_filepath) >= LINELENGTH)
  {
    perror ("Buffer overflow when creating new atomic data filepath");
    return EXIT_FAILURE;
  }

  strcpy (geo.atomic_filename, atomic_data_filepath);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Initialise all the necessary parameters to define the wind
 *
 * @param [in]  root_name  The root name of the model in $SIROCCO/source/tests/test_data/define_wind
 * @param [in]  atomic_data_location  The location of atomic data
 *
 * @details
 *
 * All the required parameters to define the wind, plasma and macro grid are initialised. There are a number of things
 * which aren't initialised (e.g. binary_basics, get_spec_type, init_observers, etc.) as these are not required to
 * initialise the grids.
 *
 * ****************************************************************************************************************** */

int
setup_model_grid (const char *root_name, const char *atomic_data_location)
{
  int n_dom;
  char *SIROCCO_ENV;
  char rdchoice_answer[LINELENGTH];
  char rdchoice_choices[LINELENGTH];
  char parameter_filepath[LINELENGTH];

  SIROCCO_ENV = getenv ("SIROCCO");
  if (SIROCCO_ENV == NULL)
  {
    fprintf (stderr, "Failed to find SIROCCO environment variable\n");
    return EXIT_FAILURE;
  }

  geo.run_type = RUN_TYPE_NEW;

  /* Reset rdpar input state in case a previous test left opar open (e.g. via a
   * fatal CUnit assertion that skipped cleanup_model). */
  rdpar_stat = 0;

  /* Set up parameter file, that way we can get all the parameters from that
   * instead of defining them manually */
  strcpy (files.root, root_name);
  snprintf (parameter_filepath, LINELENGTH, "%s/source/tests/test_data/define_wind/%s.pf", SIROCCO_ENV, files.root);
  if (opar (parameter_filepath) != 2)   /* opar returns 2 when reading for the parameter file */
  {
    fprintf (stderr, "Unable to read from parameter file %s.pf", files.root);
    return EXIT_FAILURE;
  }

  zdom = calloc (MAX_DOM, sizeof (domain_dummy));       /* We'll allocate MAX_DOM to follow sirocco */
  if (zdom == NULL)
  {
    fprintf (stderr, "Unable to allocate space for domain structure\n");
    return EXIT_FAILURE;
  }

  strncpy (rdchoice_answer, "star", LINELENGTH);
  snprintf (rdchoice_choices, LINELENGTH, "%d,%d,%d,%d,%d", SYSTEM_TYPE_STAR, SYSTEM_TYPE_CV, SYSTEM_TYPE_BH, SYSTEM_TYPE_AGN,
            SYSTEM_TYPE_PREVIOUS);
  geo.system_type = rdchoice ("System_type(star,cv,bh,agn,previous)", rdchoice_choices, rdchoice_answer);
  init_geo ();

  /* These routines don't seem to depend on the atomic data or anything which
   * depends on the atomic data */
  const double star_lum = get_stellar_params ();
  get_disk_params ();
  get_bl_and_agn_params (star_lum);

  rdint ("Wind.number_of_components", &geo.ndomain);
  if (geo.ndomain > MAX_DOM)
  {
    fprintf (stderr, "Using more domains (%d) in model than MAX_DOM (%d)\n", geo.ndomain, MAX_DOM);
    return EXIT_FAILURE;
  }

  for (n_dom = 0; n_dom < geo.ndomain; ++n_dom)
  {
    get_domain_params (n_dom);
  }

  /* We have to be a bit creative with the atomic data, to munge the correct
   * filepath with what's in the parameter file */
  init_ionization ();
  rdstr ("Atomic_data", geo.atomic_filename);
  if (set_atomic_data_filename (atomic_data_location))
  {
    return EXIT_FAILURE;
  }
  setup_atomic_data (geo.atomic_filename);

  /* We should now be able to initialise everything else which seems to have
   * some dependence on the ionisation settings or atomic data */
  for (n_dom = 0; n_dom < geo.ndomain; ++n_dom)
  {
    get_wind_params (n_dom);
  }

  if (geo.system_type == SYSTEM_TYPE_CV)
  {
    binary_basics ();
  }

  setup_windcone ();
  DFUDGE = setup_dfudge ();

  return EXIT_SUCCESS;
}
