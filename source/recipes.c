
/***********************************************************/
/** @file  recipes.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  These are collection of math utility routines
 * used in Python
 *
 *
 * The routines have been modified to make them ansi compatable
 * Some of these routines were originally Numerical Recipes 
 * Routines, but these have been replaced with gsl
 *
 * These routines should be kept SEPARATE from routines that require the Python specific
 * structures in sirocco.h
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>

#include "constants.h"
#include "math_struc.h"
#include "math_proto.h"
#include "log.h"





/**********************************************************/
/**
 * @brief      A wrapper function that carries out numerical integration on a supplied function between limits
 *
 *
 * @param [in] func - a function to integrate, needs to be of the form (double (double,void*))
 * @param [in] a - lower bound
 * @param [in] b - upper bound
 * @param [in] eps - the relative accuracy desired.
 * @return   The integral                          .
 *
 * @details
 * This routine carries out numerical integration of the function func from a to b. It currently
 * uses the gsl implementation of the ROMBERG integration scheme. Initial tests using QAG showed that
 * the type of integrals we do - often with exponential tails - is not well treated. There is an 
 * internal test to ensure that the user is not asking for an integral where the function is zero everywhere.
 * The gsl routine does not handle this well, and takes ages to run and (sometimes) return zero.
 *
 * ### Notes ###
 *
 **********************************************************/


double
num_int (func, a, b, eps)
     double (*func) (double, void *);
     double a, b;
     double eps;
{
  double result, error, result2;
  double alpha = 0.0;
  void *test = NULL;
  double delta;
  int zflag, i;
  int status = 0;
  int status2 = 0;

  size_t neval;
  gsl_function F;
  F.function = func;
  F.params = &alpha;
  zflag = 1;
  if (func (a, test) == 0.0 && func (b, test) == 0.0)
  {
    zflag = 0;
    delta = (b - a) / 101;
    for (i = 0; i < 100; i++)
    {
      if (func (a + delta * i, test) != 0)
        zflag = 1;
    }
  }
  if (zflag == 1)
  {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
    status = gsl_integration_qags (&F, a, b, 0, eps, 1000, w, &result, &error);
    if (status)
    {
      if (status == GSL_EROUND) //The rounding error has been triggered - try a different integrator
      {
        gsl_integration_workspace_free (w);
        gsl_integration_romberg_workspace *w = gsl_integration_romberg_alloc (30);
        status2 = gsl_integration_romberg (&F, a, b, 0, eps, &result2, &neval, w);
        gsl_integration_romberg_free (w);
        if (status2)
        {
          Error ("num_init: some kind of error in romberg and qags integration\n");
        }
      }
    }
    else
    {
      gsl_integration_workspace_free (w);
    }
  }
  else
  {
    result = 0.0;
  }

  return (result);
}

/**********************************************************/
/**
* @brief      A wrapper function that finds the root of a function between two limits
*
*
* @param [in] *func - the function we want to find the root of
* @param [in] double x_lo - lower bound
* @param [in] double x_hi - upper bound
* @param [in] double bol - the relative accuracy desired.
* @param [out] int * ierr An error return, TRUE if an error
* @return   The location between a and b of the zero point                         .
*
* @details
* This routine finds the root  of the function func from x_lo to x_hi. It currently
* uses the gsl implementation of the BRENT root finding scheme. This replaced the 
* zbrentnumerical recipes, and the call is intended to be identical
*
* ### Notes ###
*
*  The function needs to be of the form 
*
* double (* function) (double x, void * params)
*
* See the documentation for gsl_function, but historically passed
* the information needed via external variables and so typically
* nothing is passed in the parameters.  This is something that
* ultimately should be changed.
*
**********************************************************/


double
zero_find (func, x1, x2, tol, ierr)
     double (*func) (double, void *);
     double x1, x2;
     double tol;
     int *ierr;
{
  double result;
  double x_below, x_above;
  double alpha = 0.0;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  int iter = 0, max_iter = 100;
  int status;
  double x_lo, x_hi;

// If necessary, reorder the inputs to the gsl function

  if (x1 < x2)
  {
    x_lo = x1;
    x_hi = x2;
  }
  else
  {
    x_lo = x2;
    x_hi = x1;
  }



  *ierr = FALSE;
  // Check that the interval is bracketed

  x_below = func (x_lo, (void *) &alpha);
  x_above = func (x_hi, (void *) &alpha);

  if (x_below * x_above > 0.0)
  {
    Error ("zero_find: function not bracketed x_lo %e -> %e, x_hi %e -> %e\n", x_lo, x_below, x_hi, x_above);
    *ierr = TRUE;
  }


  gsl_function F;
  F.function = func;
  F.params = &alpha;



  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);



  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi, tol, 0);



  }
  while (status == GSL_CONTINUE && iter < max_iter);

  if (status != GSL_SUCCESS)
  {
    x_below = func (x_lo, (void *) &alpha);
    x_above = func (x_hi, (void *) &alpha);
    Error ("zero_find failed: %d of %d brackets x_lo %e -> %e, x_hi %e -> %e\n", iter, max_iter, x_lo, x_below, x_hi, x_above);
    *ierr = TRUE;
  }

  result = (x_lo + x_hi) / 2.0;

  gsl_root_fsolver_free (s);

  return (result);
}


/**********************************************************/
/**
* @brief      Find the mimimum value of a function in an interval
*
*
* @param [in] a - lower bound
* @param [in] b - guess for the minumum
* @param [in] c - upper bound
* @param [in] f  the function to minimize
* @param [in] tol A tolerance
* @param [out] * xmin - The place where the minimum occors
* @return   The minimum of the function f                           .
*
* @details
* This is a wrapper routine which uses gsl routines to 
* find the minimum of a function
*
* ### Notes ###
*
* At of  212222, this routine is only use by the roche lobe 
* realated routines
* 
**********************************************************/



double
find_function_minimum (a, m, b, func, tol, xmin)
     double (*func) (double, void *);
     double a, m, b;
     double tol, *xmin;

{
  int status = 0;
  void *test = NULL;

  int iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;

  gsl_function F;
  F.function = func;
  F.params = 0;

  T = gsl_min_fminimizer_brent;

  s = gsl_min_fminimizer_alloc (T);
  status = gsl_min_fminimizer_set (s, &F, m, a, b);
  if (status)
  {
    if (status == GSL_EINVAL)   //THere is no minimum 
    {
      return fmin (func (a, test), func (b, test));     //keep old behaviour
    }
  }

  do
  {
    iter++;
    status = gsl_min_fminimizer_iterate (s);
    m = gsl_min_fminimizer_x_minimum (s);
    a = gsl_min_fminimizer_x_lower (s);
    b = gsl_min_fminimizer_x_upper (s);
    status = gsl_min_test_interval (a, b, 0.0, tol);
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  *xmin = m;

  gsl_min_fminimizer_free (s);

  return func (m, test);
}





/**********************************************************/
/**
 * @brief Perform linear/logarithmic interpolation of an array
 *
 * @param [in] double  value   The value used for interpoaltion
 * @param [in] double  array[]   An array containing a set of asceding values
 * @param [in] int  npts   The size of the array
 * @param [out] int *  ival  The lower index to use in the interpolation
 * of the array which need to be used to interpolate on
 * @param [out] double *  f   The fraction of the upper point to use in the interpolation
 * @param [in] int  mode  A switch to choose linear(0)  or lograrithmic(1) interpolation
 * @return     Usually returns 0, but returns -1 if the input value is less
 * than the first elememnt in the array, and 1 if it is greater than the
 * last element int he array.  In either of these cases, the fractions are
 * set up only to access one array element
 *
 * @details
 *
 *
 * Typically one has two parallel arrays, one containing a set of values,
 * in the original case frequencies, and another containing some function
 * of those values.  This routine finds the fraction of the function at
 * two point in the data array to interpolate.
 *
 * The values for the input array can be interpolated linear or logarithmically.
 * that is the fraction that is returned are based on the values in the
 * array or the logarithm of them.
 *
 *
 * The routine uses bisection
 *
 *
 *
 * ### Notes ###
 *
 * fraction is a utility written to speed up the search for the bracketing
 * topbase photionization x-sections, but in principle it should work in
 * other situations within sirocco (which at one time at least contributed
 * significantly to the runtime for Python).
 *
 * Today, fraction is called directly in the routines that are used to set up coordinate
 * systems, and indirectly through linterp (below) which should be inspected
 * to see how the logarithmic interpolation is supposed to work.
 *
 * The routine is similar to the numerical * recipes routine locate.
 * There may be a gsl routine as well.  The routine
 * should probably be replaced.
 *
 *
 **********************************************************/

int
fraction (value, array, npts, ival, f, mode)
     double array[];            // The array in we want to search
     int npts, *ival;           // ival is the lower point
     double value;              // The value we want to index
     double *f;                 // The fractional "distance" to the next point in the array
     int mode;                  // 0 = compute in linear space, 1=compute in log space
{
  int imin, imax, ihalf;

  if (value < array[0])
  {
    *ival = 0;
    *f = 0.0;
    return (-1);
  }

  imax = npts - 1;
  if (value > array[imax])
  {
    *ival = npts - 2;
    *f = 1.0;
    return (1);
  }



  imin = 0;

/* In what follows, there is a specific issue as to how to treat
the situation where the value is exactly on an array element. Right
now this is set to try to identify the array element below the one
on which the value sits, and to set the fraction to 1.  This was
to reflect the behavior of the search routine in where_in_grid. */

  while (imax - imin > 1)
  {
    ihalf = (imin + imax) >> 1; // Compute a midpoint >> is a bitwise right shift
    if (value > array[ihalf])
    {
      imin = ihalf;
    }
    else
      imax = ihalf;
  }

// So array[imin] just <= value

  if (mode == 0)
    *f = (value - array[imin]) / (array[imax] - array[imin]);   //linear interpolation
  else if (mode == 1)
    *f = (log (value) - log (array[imin])) / (log (array[imax]) - log (array[imin]));   //log interpolation
  else
  {
    Error ("Fraction - unknown mode %i\n", mode);
    exit (0);
    return (0);
  }

  *ival = imin;

  return (0);
}







/**********************************************************/
/**
 * @brief      Perform a linear interpolation on two parallel arrays, the first
 * of which contains a set of values to be interpolated and the second of
 * which has the function at those values
 *
 * @param [in] double  x   A value
 * @param [in] double  xarray[]   The array that is interpolated
 * @param [in] double  yarray[]   The array that contains a function of the values in xarray
 * @param [in] int  xdim   The length of the two arrays
 * @param [out] double *  y   The resulting interpolated value
 * @param [in] int  mode   A switch to choose linear(0) or "logarithmic" (1)
 * interpolation
 *
 * @return     The number of the array element that is used for the lower of the two
 * elements that are interpolated on.
 *
 * @details
 * Given a number x, and an array of x's in xarray, and functional
 * values y = f(x) in yarray and the dimension of xarray and yarray,
 * linterp calculates y, the linearly interpolated value of y(x). It
 * also returns the element in the xarray so that one can consider
 * not doing a search to find the right array element in certain
 * circumstances

 *
 * ### Notes ###
 *
 * For mode 0, the value that is returned is
 *
 * (1-f)*y[nelem]+f*[nelem+1)
 *
 * For mode 1, the value returned is
 * exp ((1. - f) * log (y[nelem]) + f * log (y[nelem + 1]))
 *
 *
 **********************************************************/

int
linterp (x, xarray, yarray, xdim, y, mode)
     double x;                  // The value that we wish to index i
     double xarray[], yarray[];
     int xdim;
     double *y;
     int mode;                  //0 = linear, 1 = log
{
  int nelem = 0;
  double frac;


  fraction (x, xarray, xdim, &nelem, &frac, mode);

  if (yarray[nelem] == yarray[nelem + 1])
  {
    // Prevent round-off errors when the numbers are identical`
    *y = yarray[nelem];
  }
  else if (mode == 0)
    *y = (1. - frac) * yarray[nelem] + frac * yarray[nelem + 1];
  else if (mode == 1)
    *y = exp ((1. - frac) * log (yarray[nelem]) + frac * log (yarray[nelem + 1]));
  else
  {
    Error ("linterp - unknown mode %i\n", mode);
    exit (0);
  }

  return (nelem);

}
