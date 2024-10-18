/** ********************************************************************************************************************
 *
 *  @file test_matrix.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief Unit tests for matrix functions in Python
 *
 * ****************************************************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <CUnit/CUnit.h>

#include "../../atomic.h"
#include "../../sirocco.h"
#include "../assert.h"

#define BUFFER_LENGTH 512


/** *******************************************************************************************************************
 *
 * @brief Load in test data for the matrix inversion tests
 *
 * @param [in] matrix_path filepath to the input matrix for the tests
 * @param [in] inverse_path filepath to the input matrix for verification
 * @param [out] matrix pointer to pointer for matrix input
 * @param [out] inverse point to pointer for inverse used for verification
 * @param [out] size the size of the matrices
 *
 * @return an error status
 *
 * @details
 *
 * The matrix and inverse variables are pointers to pointers, as the memory is allocated for those variables in this
 * function.
 *
 * ****************************************************************************************************************** */

static int
get_invert_matrix_test_data (const char *matrix_path, const char *inverse_path, double **matrix, double **inverse, int *size)
{
  FILE *fp_matrix = fopen (matrix_path, "r");
  if (!fp_matrix)
  {
    return EXIT_FAILURE;
  }

  int size_from_matrix = 0;
  fscanf (fp_matrix, "%d", &size_from_matrix);
  *size = size_from_matrix;

  *matrix = malloc (size_from_matrix * size_from_matrix * sizeof (double));
  *inverse = malloc (size_from_matrix * size_from_matrix * sizeof (double));

  int i;
  for (i = 0; i < size_from_matrix * size_from_matrix; ++i)
  {
    fscanf (fp_matrix, "%le", &(*matrix)[i]);
  }
  fclose (fp_matrix);

  FILE *fp_inverse = fopen (inverse_path, "r");
  if (!fp_inverse)
  {
    return EXIT_FAILURE;
  }

  fscanf (fp_inverse, "%*d");
  for (i = 0; i < size_from_matrix * size_from_matrix; ++i)
  {
    fscanf (fp_inverse, "%le", &(*inverse)[i]);
  }

  fclose (fp_inverse);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Load in test data for the solve matrix tests
 *
 * @param [in] a_path the filepath to the A matrix for input
 * @param [in] b_path the filepath to the b vector for input
 * @param [in] x_path the filepath to the x vector for verification
 * @param [out] a_matrix pointer to pointer for A matrix for input
 * @param [out] b_vector pointer to pointer for b vector for input
 * @param [out] x_vector pointer to pointer for x vector for verification
 * @param [out] size the size of the vectors
 *
 * @return an error status
 *
 * @details
 *
 * Pointer to pointers are used for a_matrix and etc. because the memory for those pointers/arrays are allocated here.
 *
 * ****************************************************************************************************************** */

static int
get_solve_matrix_test_data (const char *a_path, const char *b_path, const char *x_path, double **a_matrix, double **b_vector,
                            double **x_vector, int *size)
{
  FILE *fp_a = fopen (a_path, "r");
  if (!fp_a)
  {
    return EXIT_FAILURE;
  }

  int size_from_a = 0;
  fscanf (fp_a, "%d", &size_from_a);
  *size = size_from_a;

  *a_matrix = malloc (size_from_a * size_from_a * sizeof (double));
  *b_vector = malloc (size_from_a * sizeof (double));
  *x_vector = malloc (size_from_a * sizeof (double));

  int i;
  for (i = 0; i < size_from_a * size_from_a; ++i)
  {
    fscanf (fp_a, "%le", &(*a_matrix)[i]);
  }
  fclose (fp_a);

  FILE *fp_b = fopen (b_path, "r");
  if (!fp_b)
  {
    return EXIT_FAILURE;
  }

  FILE *fp_x = fopen (x_path, "r");
  if (!fp_x)
  {
    return EXIT_FAILURE;
  }

  fscanf (fp_b, "%*d");
  fscanf (fp_x, "%*d");
  for (i = 0; i < size_from_a; ++i)
  {
    fscanf (fp_b, "%le", &(*b_vector)[i]);
    fscanf (fp_x, "%le", &(*x_vector)[i]);
  }

  fclose (fp_b);
  fclose (fp_x);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Test case for invert_matrix
 *
 * @param [in] test_name the test name
 *
 * @return an error status
 *
 * @details
 *
 * Calls invert matrix with given input and checks the absolute tolerance.
 *
 * ****************************************************************************************************************** */

static int
call_invert_matrix (const char *test_name)
{
  const char *sirocco_path = getenv ((const char *) "SIROCCO");
  if (sirocco_path == NULL)
  {
    CU_FAIL_FATAL ("$SIROCCO has not been set");
  }

  double *matrix;
  double *inverse;
  char matrix_filepath[BUFFER_LENGTH];
  char inverse_filepath[BUFFER_LENGTH];

  sprintf (matrix_filepath, "%s/source/tests/test_data/matrix/%s/matrix.txt", sirocco_path, test_name);
  sprintf (inverse_filepath, "%s/source/tests/test_data/matrix/%s/inverse.txt", sirocco_path, test_name);

  int matrix_size;
  const int get_err = get_invert_matrix_test_data (matrix_filepath, inverse_filepath, &matrix, &inverse, &matrix_size);
  if (get_err)
  {
    CU_FAIL_MSG_FATAL ("Unable to load test data");
  }

  double *test_inverse = malloc (matrix_size * matrix_size * sizeof (double));
  const int matrix_err = invert_matrix (matrix, test_inverse, matrix_size);
  if (matrix_err)
  {
    CU_FAIL_MSG_FATAL ("`invert_matrix` failed with error");
  }

  CU_ASSERT_DOUBLE_ARRAY_EQUAL_FATAL (test_inverse, inverse, matrix_size, EPSILON);

  free (matrix);
  free (inverse);
  free (test_inverse);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Test case for solve_matrix
 *
 * @param [in] test_name the test name
 *
 * @return an error status
 *
 * @details
 *
 * Calls solve_matrix with given input and checks the absolute tolerance.
 *
 * ****************************************************************************************************************** */

int
call_solve_matrix (const char *test_name)
{
  const char *sirocco_path = getenv ((const char *) "SIROCCO");
  if (sirocco_path == NULL)
  {
    CU_FAIL_FATAL ("$SIROCCO has not been set");
  }

  double *matrix_a;
  double *vector_b;
  double *vector_x;
  char matrix_a_filepath[BUFFER_LENGTH];
  char vector_b_filepath[BUFFER_LENGTH];
  char vector_x_filepath[BUFFER_LENGTH];

  sprintf (matrix_a_filepath, "%s/source/tests/test_data/matrix/%s/A.txt", sirocco_path, test_name);
  sprintf (vector_b_filepath, "%s/source/tests/test_data/matrix/%s/b.txt", sirocco_path, test_name);
  sprintf (vector_x_filepath, "%s/source/tests/test_data/matrix/%s/x.txt", sirocco_path, test_name);

  int vector_size;
  const int get_err =
    get_solve_matrix_test_data (matrix_a_filepath, vector_b_filepath, vector_x_filepath, &matrix_a, &vector_b, &vector_x, &vector_size);
  if (get_err)
  {
    CU_FAIL ("Unable to load test data");
  }

  double *test_vector_x = malloc (vector_size * sizeof (double));
  const int matrix_err = solve_matrix (matrix_a, vector_b, vector_size, test_vector_x, -1);
  if (matrix_err)
  {
    CU_FAIL ("`solve_matrix` failed with error");
  }

  CU_ASSERT_DOUBLE_ARRAY_EQUAL_FATAL (test_vector_x, vector_x, vector_size, EPSILON);

  free (matrix_a);
  free (vector_b);
  free (vector_x);
  free (test_vector_x);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Tests for `solve_matrix`
 *
 * @details
 *
 * This function is called by the "Matrix Functions" tests suite. Each test for `solve_matrix` should be put into here.
 *
 * ****************************************************************************************************************** */

void
test_solve_matrix (void)
{
  call_solve_matrix ("small_matrix");
  call_solve_matrix ("matrix_ion");
}

/** *******************************************************************************************************************
 *
 * @brief Tests for `invert_matrix`
 *
 * @details
 *
 * This function is called by the "Matrix Functions" tests suite. Each test for `invert_matrix` should be put into here.
 *
 * ****************************************************************************************************************** */

void
test_invert_matrix (void)
{
  call_invert_matrix ("inverse_small");
  call_invert_matrix ("inverse_macro");
}

/** *******************************************************************************************************************
 *
 * @brief Initialise the program for the matrix test suite
 *
 * @return An error code status
 *
 * @details
 *
 * Either initialises cuSolver or turns off the GSL error handler.
 *
 * ****************************************************************************************************************** */

int
matrix_suite_init (void)
{
  int error = EXIT_SUCCESS;

#ifdef CUDA_ON
  error = cusolver_create ();
#endif

  return error;
}

/** *******************************************************************************************************************
 *
 * @brief Clean up after the matrix test suite
 *
 * @return An error code status
 *
 * @details
 *
 * Either destroys the cuSolver handle, or turns GSL error handling back on.
 *
 * ****************************************************************************************************************** */

int
matrix_suite_teardown (void)
{
  int error = EXIT_SUCCESS;

#ifdef CUDA_ON
  error = cusolver_destroy ();
#endif

  return error;
}

/** *******************************************************************************************************************
 *
 * @brief Create a CUnit test suite for matrix functions
 *
 * @details
 *
 * This function will create a test suite for the matrix functions.
 *
 * The design philosophy of the tests is to only test what Python is compiled for, e.g. does it use GSL or CUDA? This is
 * done to keep the code and build system as simple as possible.
 *
 * ****************************************************************************************************************** */

void
create_matrix_test_suite (void)
{
#ifdef CUDA_ON
  char *suite_name = "Matrix Functions: cuSolver";
#else
  char *suite_name = "Matrix Functions: GNU Science Library";
#endif

  CU_pSuite suite = CU_add_suite (suite_name, matrix_suite_init, matrix_suite_teardown);

  if (suite == NULL)
  {
    fprintf (stderr, "Failed to create `Matrix Functions Suite` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }

  /* Add CPU tests to suite */
  if ((CU_add_test (suite, "Solve Matrix", test_solve_matrix) == NULL) ||
      (CU_add_test (suite, "Invert Matrix", test_invert_matrix) == NULL))
  {
    fprintf (stderr, "Failed to add tests to `Matrix Functions Suite: CPU`\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
