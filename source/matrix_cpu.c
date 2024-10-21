/* ****************************************************************************************************************** */
/**
 *  @file matrix_cpu.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>

#include "atomic.h"
#include "sirocco.h"

#ifndef CUDA_ON                 /* removes a compiler warning about unused functions */

/* ****************************************************************************************************************** */
/**
 *  @brief Check the return status of a GSL function
 *
 *  @param [in] status  the status to check
 *
 *  ***************************************************************************************************************** */

#define GSL_CHECK(status)                                                                                              \
  do {                                                                                                                 \
    if (status != GSL_SUCCESS) {                                                                                       \
      Error("GSL Error: %s (%d)\n", gsl_strerror(status), status);                                                     \
      return status;                                                                                                   \
    }                                                                                                                  \
  } while (0)

/* ****************************************************************************************************************** */
/**
 * @brief Get the error string for a error code from the GSL matrix solver
 *
 * @param [in] error_code the error code returned from the GSL matrix solver function
 *
 * @return The error string
 *
 * @details
 *
 * @todo
 *
 * These error codes should be enumerators.
 *
 *  ***************************************************************************************************************** */


static const char *
cpu_matrix_error_string (int error_code)
{
  switch (error_code)
  {
  case 2:
    return "some matrix rows failing relative error check";
  case 3:
    return "some matrix rows failing absolute error check";
  case 4:
    return "Unsolvable matrix! Determinant is zero. Defaulting to no change";
  default:
    return "bad return from solve_matrix";
  }
}

/* ****************************************************************************************************************** */
/**
 * @brief  Solve the linear system A x = b, for the vector x
 *
 * @param  [in]  a_matrix a square matrix on the LHS
 * @param  [in]  b_vector the B resultant vector
 * @param  [in]  matrix-size the number of rows (and columns) in the square matrix matrix and vectors
 * @param  [out] x_vector the x vector on the RHS
 * @param  [in]  nplasma the index of the plasma cell we are working on
 *
 * @return an integer representing the error state
 *
 * @details
 * Performs LU decomposition to solve for x in the linear system A x = b. The calculation is perform, in serial, on
 * the CPU using GSL.
 *
 * ### Notes ###
 *
 * https://www.gnu.org/software/gsl/doc/html/linalg.html
 *
 *  ***************************************************************************************************************** */

static int
cpu_solve_matrix (double *a_matrix, double *b_vector, int matrix_size, double *x_vector, int nplasma)
{
  int i, gsl_err, signnum;
  double test_val;
  double m_determinant;
  int n_error = 0;
  gsl_permutation *p;
  gsl_matrix_view m;
  gsl_vector_view b;
  gsl_vector *test_vector;
  gsl_vector *solution_vector;
  gsl_matrix *test_matrix;

  /* create gsl matrix/vector views of the arrays of rates.
   * This is the structure that gsl uses do define an array.
   * It contains not only the data but the dimensions, etc.*/
  m = gsl_matrix_view_array (a_matrix, matrix_size, matrix_size);

  /* gsl_vector_view_array creates the structure that gsl uses to define a vector
   * It contains the data and the dimension, and other information about where
   * the vector is stored in memory etc.
   */
  b = gsl_vector_view_array (b_vector, matrix_size);
  solution_vector = gsl_vector_alloc (matrix_size);

  /* permutations are special structures that contain integers 0 to size-1, which can
   * be manipulated -- in the GPU example, this is the pivot array */
  p = gsl_permutation_alloc (matrix_size);

  /* these are used for testing the solution below */
  test_matrix = gsl_matrix_alloc (matrix_size, matrix_size);
  test_vector = gsl_vector_alloc (matrix_size);
  GSL_CHECK (gsl_matrix_memcpy (test_matrix, &m.matrix));       // create copy for testing

  /* Decomposes a_matrix (now m as a GSL matrix) into its LU Components. The L and U components are stored in m.
   */
  gsl_err = gsl_linalg_LU_decomp (&m.matrix, p, &signnum);      /* this modifies a_matrix for some reason */
  if (gsl_err)
  {
    Error ("Solve_matrix: gsl_linalg_LU_decomp failure %d for cell %i \n", gsl_err, nplasma);
    Exit (EXIT_FAILURE);
  }

  gsl_err = gsl_linalg_LU_solve (&m.matrix, p, &b.vector, solution_vector);
  if (gsl_err)
  {
    m_determinant = gsl_linalg_LU_lndet (&m.matrix);    // get the determinant to report to user
    Error ("Solve_matrix: gsl_linalg_LU_solve failure (error %d determinant %.3e) for plasma cell %i \n", gsl_err, m_determinant, nplasma);
    return (4);
  }

  gsl_permutation_free (p);

  /* Check that the populations vector we have just created really is a solution to
     the matrix equation. The following line does the matrix multiplication test_vector = 1.0 * test_matrix * populations.
     The CblasNoTrans statement just says we do not do anything to test_matrix, and the 0.0 means we do not add a
     second matrix to the result. If the solution has worked, then test_vector should be equal to b_temp */

  gsl_err = gsl_blas_dgemv (CblasNoTrans, 1.0, test_matrix, solution_vector, 0.0, test_vector); /* calculates matrix vector product */
  if (gsl_err != 0)
  {
    Error ("Solve_matrix: bad return (%d) when testing matrix solution to rate equations in plasma cell %d\n", gsl_err, nplasma);
  }

  /* now cycle through and check the solution to y = m * populations really is (1, 0, 0 ... 0) */
  for (i = 0; i < matrix_size; i++)
  {
    /* get the element of the vector we want to check */
    test_val = gsl_vector_get (test_vector, i);

    /* b_vector is (1, 0, 0, 0 ...) when we do matom rates. test_val is normally something like
       1e-16 if it's supposed to be 0. We have a different error check if b_vector[i] is 0 */
    if (b_vector[i] > 0.0)
    {
      if (n_error == 0 && fabs ((test_val - b_vector[i]) / test_val) > EPSILON)
      {
        Error ("Solve_matrix: test solution fails relative error for row %i %e != %e frac_error %3 in plasma cell %d\n", i, test_val,
               b_vector[i], fabs ((test_val - b_vector[i]) / test_val), nplasma);
        gsl_err = 2;
        n_error += 1;
      }
    }
    else if ((n_error) == 0 && fabs (test_val - b_vector[i]) > EPSILON) // if b_vector is 0, check absolute error
    {
      Error ("Solve_matrix: test solution fails absolute error for row %i %e != %e in plasma cell %d\n", i, test_val, b_vector[i], nplasma);
      gsl_err = 3;
      n_error += 1;
    }
  }

  if (n_error > 1)
  {
    Error ("Solve_matrix: There were %d row errors in all for plasma cell %d\n", n_error, nplasma);
  }

  /* copy the solution vector to the output array */
  for (i = 0; i < matrix_size; i++)
  {
    x_vector[i] = gsl_vector_get (solution_vector, i);
  }

  /* free memory */
  gsl_vector_free (test_vector);
  gsl_matrix_free (test_matrix);
  gsl_vector_free (solution_vector);

  return gsl_err;
}

/* ****************************************************************************************************************** */
/**
 * @brief Calculate the inverse of the square matrix `matrix`
 *
 * @param  [in]  matrix the matrix to compute the inverse for
 * @param  [out] inverse_inverse the inverse of `matrix`
 * @param  [in]  matrix_size  the size of the matrix
 *
 * @return an integer representing the error state
 *
 * @details
 *
 * Performs LU decomposition to compute the inverse of a matrix. Uses the GSL matrix library, which is (I think) a
 * wrapper around lapack.
 *
 *  ***************************************************************************************************************** */

static int
cpu_invert_matrix (double *matrix, double *inverse_out, int matrix_size)
{
  int s;
  int i, j;
  gsl_matrix_view N;
  gsl_matrix *inverse;
  gsl_permutation *p;

  N = gsl_matrix_view_array (matrix, matrix_size, matrix_size);
  p = gsl_permutation_alloc (matrix_size);
  inverse = gsl_matrix_alloc (matrix_size, matrix_size);

  if (!p || !inverse)
  {
    Error ("unable to allocate memory for matrix inversion\n");
    return GSL_ENOMEM;
  }

  GSL_CHECK (gsl_linalg_LU_decomp (&N.matrix, p, &s));
  GSL_CHECK (gsl_linalg_LU_invert (&N.matrix, p, inverse));
  gsl_permutation_free (p);

  for (i = 0; i < matrix_size; ++i)     /* i is mm in macro_accelerate.c */
  {
    for (j = 0; j < matrix_size; ++j)   /* j is nn in macro_accelerate.c */
    {
      inverse_out[i * matrix_size + j] = gsl_matrix_get (inverse, i, j);
    }
  }

  gsl_matrix_free (inverse);

  return EXIT_SUCCESS;
}

#endif

/* ****************************************************************************************************************** */
/**
 * @brief Get the error string for a error code for `solve_matrix`
 *
 * @param [in] error_code the error code
 *
 * @return The error string
 *
 * @details
 *
 *  ***************************************************************************************************************** */

const char *
get_matrix_error_string (int error_code)
{
#ifdef CUDA_ON
  return cusolver_get_error_string (error_code);
#else
  return cpu_matrix_error_string (error_code);
#endif
}

/* ****************************************************************************************************************** */
/**
 * @brief Solve the linear system A x = b, for the vector x
 *
 * @param  [in]  a_matrix - a square matrix on the LHS
 * @param  [in]  b_vector - the B resultant vector
 * @param  [in]  size - the number of rows (and columns) in the square matrix matrix and vectors
 * @param  [out] x_vector - the x vector on the RHS
 * @param  [in]  nplasma - the index of the plasma cell we are working on
 *
 * @return an integer representing the error state
 *
 * @details
 * Performs LU decomposition to solve for x in the linear system A x = b.
 *
 * This function is a wrapper around two others: `cpu_solve_linear_system` and `gpu_solve_linear_system`.
 *
 * ### Notes ###
 *
 * `nplasma` is only used to indicate a cell number, if the calculation fails.
 *
 *  ***************************************************************************************************************** */

#include <time.h>
#include <sys/time.h>

int
solve_matrix (double *a_matrix, double *b_matrix, int size, double *x_matrix, int nplasma)
{
  int error;

#ifdef CUDA_ON
  error = gpu_solve_matrix (a_matrix, b_matrix, size, x_matrix);
#else
  error = cpu_solve_matrix (a_matrix, b_matrix, size, x_matrix, nplasma);
#endif

  return error;
}

/* ****************************************************************************************************************** */
/**
 * @brief
 *
 * @param  [in]  a_matrix
 * @param  [out] a_inverse
 * @param  [in]  num_rows
 *
 * @return an integer representing the error state
 *
 * @details
 *
 *  ***************************************************************************************************************** */

int
invert_matrix (double *matrix, double *inverted_matrix, int num_rows)
{
  int error;

#ifdef CUDA_ON
  error = gpu_invert_matrix (matrix, inverted_matrix, num_rows);
#else
  error = cpu_invert_matrix (matrix, inverted_matrix, num_rows);
#endif

  return error;
}
