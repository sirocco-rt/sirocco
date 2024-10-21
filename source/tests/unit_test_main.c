/** ********************************************************************************************************************
 *
 *  @file main.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief Main function for unit tests
 *
 *  https://gitlab.com/cunity/cunit
 *
 * ****************************************************************************************************************** */

#ifdef MPI_ON
#include <mpi.h>
#endif

#include <time.h>
#include <stdlib.h>
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include <gsl/gsl_errno.h>

/* Test suite prototypes */
#include "tests/suites.h"

/* Python logging prototypes */
#include "../atomic.h"
#include "../sirocco.h"

/** *******************************************************************************************************************
 *
 *  @brief Entry point for unit tests
 *
 * ****************************************************************************************************************** */

int
main (int argc, char **argv)
{
#ifdef MPI_ON
  int mpi_err = MPI_Init (&argc, &argv);
  if (mpi_err != EXIT_SUCCESS)
  {
    Error ("Failed to initialise MPI\n");
    Exit (mpi_err);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &rank_global);
  MPI_Comm_size (MPI_COMM_WORLD, &np_mpi_global);
#else
  rank_global = 0;
  np_mpi_global = 1;
#endif

  if (np_mpi_global != 1)
  {
    fprintf (stderr, "Using multiple MPI ranks is not supported\n");
    return EXIT_FAILURE;
  }

  /* We don't want this to get in the way, like in Python */
  gsl_set_error_handler_off ();

  /* I don't want any output from Python, thank you very much */
  Log_init ("logfile");
  Log_set_verbosity (0);
  rdpar_set_verbose (0);
  Log_set_mpi_rank (rank_global, np_mpi_global);

  /* Initialize the CUnit test registry */
  if (CU_initialize_registry () != CUE_SUCCESS)
  {
    return CU_get_error ();
  }

  /* Initialise some stuff for Python -- e.g. RNG */
  init_rand ((int) time (NULL));

  /* Create test suites */
  create_matrix_test_suite ();
  create_compton_test_suite ();
  create_define_wind_test_suite ();
//  create_run_mode_test_suite ();
  create_translate_test_suite ();

  /* Run the test suites */
  CU_basic_set_mode (CU_BRM_VERBOSE);
  CU_basic_run_tests ();

  const int num_tests_failed = (int) CU_get_number_of_tests_failed ();

  if (num_tests_failed > 0)
  {
    printf ("\033[1;31m%d test(s) failed\n\033[1;0m", num_tests_failed);
  }
  else
  {
    printf ("\033[1;32mAll tests ran successfully\n\033[1;0m");
  }

  /* Clean up the CUnit registry */
  CU_cleanup_registry ();
  Log_close ();

#ifdef MPI_ON
  MPI_Finalize ();
#endif

  return num_tests_failed;
}
