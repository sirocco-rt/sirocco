/** ********************************************************************************************************************
 *
 *  @file test_run_mode.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date November 2023
 *
 *  @brief Unit tests for `define_wind`
 *
 * ****************************************************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <CUnit/CUnit.h>

#include "../assert.h"
#include "../unit_test.h"
#include "../../atomic.h"
#include "../../sirocco.h"

/** *******************************************************************************************************************
 *
 * @brief Test a SV wind model for a CV system.
 *
 * @details
 *
 * ****************************************************************************************************************** */

static void
test_previous_wind (void)
{
  geo.run_type = RUN_TYPE_PREVIOUS;
}

/** *******************************************************************************************************************
 *
 * @brief Test a SV wind model for a CV system.
 *
 * @details
 *
 * ****************************************************************************************************************** */

static void
test_restart_model (void)
{
  geo.run_type = RUN_TYPE_RESTART;
}

/** *******************************************************************************************************************
 *
 * @brief Test a SV wind model for a CV system.
 *
 * @details
 *
 * ****************************************************************************************************************** */

static int
suite_init (void)
{
  char windsave_filepath[LINELENGTH];

  /* Check that the wind_save was generated in an earlier test */
  const char *sirocco_loc = get_sirocco_env_variable ();
  if (sirocco_loc == NULL)
  {
    return EXIT_FAILURE;
  }
  snprintf (windsave_filepath, LINELENGTH, "%s/source/tests/test_data/define_wind/restart_cv.wind_save", sirocco_loc);
  if (access (windsave_filepath, F_OK) == -1)
  {
    fprintf (stderr, "Failed to open %s for test\n", windsave_filepath);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Test a SV wind model for a CV system.
 *
 * @details
 *
 * ****************************************************************************************************************** */

static int
suite_teardown (void)
{
  char windsave_filepath[LINELENGTH];

  /* Delete the test file */
  const char *sirocco_loc = get_sirocco_env_variable ();
  if (sirocco_loc == NULL)
  {
    return EXIT_FAILURE;
  }
  snprintf (windsave_filepath, LINELENGTH, "%s/source/tests/test_data/define_wind/restart_cv.wind_save", sirocco_loc);
  int err = remove (windsave_filepath);
  if (err != 0)
  {
    fprintf (stderr, "Failed to remove %s after test\n", windsave_filepath);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Test a SV wind model for a CV system.
 *
 * @details
 *
 * ****************************************************************************************************************** */

void
create_run_mode_test_suite (void)
{
  CU_pSuite suite = CU_add_suite ("Run Mode", suite_init, suite_teardown);

  if (suite == NULL)
  {
    fprintf (stderr, "Failed to create `Run Mode` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }

  if ((CU_add_test (suite, "Restart Model", test_restart_model) == NULL)
      || (CU_add_test (suite, "Previous Model", test_previous_wind) == NULL))
  {
    fprintf (stderr, "Failed to add tests to `Run Mode` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
