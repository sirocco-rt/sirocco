/***********************************************************/
/** @file  cdf.c
 * @author ksl,nsh
 * @date   March, 2018
 *
 * @brief  Routines to create a cumulatimve distribution (CDF) from a function or array, and then to
 * sample that cdf in a fleximble manner.
 *
 * @details
 *
 * CDFs can be generated either from a function using cdf_gen_from_func or from an array using
 * cdf gen_from_array.  The function or the array that is input are used as the proability densities
 * for the cdfs.  The cdfs are stored in a Cdf structure (See sirocco.h)
 *
 * The procedure is either
 * to call cdf_gen_from_func if you have a function whose probability distribution you
 * want (as might be the case if you were generating bremsstrahlung photons and you knew
 * the formula)  or cdf_gen_from_array (as might be the case if you were trying to
 * generate Monte Carlo spectra from precalculated Kurucz models).
 *
 * Once the cdfs are generated one can sample the full distribution distritution function
 * with cdf_get_rand, or one can sample a part of the distribution by setting the
 * part that one wants with cdf_limit and then sampling the distribution with cdf_get_rand_limit
 *
 * There are a number of helper functions that are internal to the generation of the cdfs,
 * and verification that the cdfs are readonable.
 *
 * ###Notes###
 *
 * In generating the CDFs, one must be careful of places where the pdf is discontinuous, or more
 * generally changing in a way that will not be well represented via linear interpolations.  When
 * generating a CDF from a function it is possible to specify positions where one wishes to have
 * points, e.g on either side of a discontinuity.  When generating a CDF from a function, the
 * discontinuities should be well sampled in the array that is provided.
 *
 * @bug For reasons, which are currently unclear there are differences in the number of points
 * maintained in the cdfs for different generation methods.
 *
 * These routines should be kept SEPARATE from routines that require the Python specific
 * structures in sirocco.h 
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <math.h>
#include "constants.h"
#include "math_struc.h"
#include "math_proto.h"
#include "log.h"
//#include "sirocco.h"
//#include "models.h"

/// This is the initial value of PDFSTEPS
#define PDFSTEPS 100000

/// This is the value of pdfsteps at this point in time         
int pdf_steps_current;

/// Flag to say whether structures to hold a cdf have been initialised
int init_pdf = 0;

/// Pointer to pdf_array - made external because it is used in varous routines
double *pdf_array;


/**********************************************************/
/**
 * @brief      Generate a cumulative distrubution function (cdf) from a function

 * @param [in, out] CdfPtr  cdf    	A ptr to a cdf structure
 * @param [in] double *func  		The probablility density function to be integrated to create the pdf
 * @param [in] double xmin,xmax		The range over which the function is to be integrated
 * @param [in] int njumps			The number of points at which we want to force a break in the cdf
 * @param [in] double jump			an array of points in the x-corrdinates at which we want to force points
 * @return     0 on successful completion,
 *
 * @details
 * This routine allows one to generate a CDF from a function. It is used
 * to make photons from blackbody and bremstrahlung sources and
 * also to sample things like cosine functions for directions of generated photons.
 *
 *
 * ### Notes ###
 *
 * The jumps are used to force points in the CDF at certain values of X.
 *
 *
 **********************************************************/

int
cdf_gen_from_func (cdf, func, xmin, xmax, njumps, jump)
     CdfPtr cdf;
     double (*func) (double, void *);
     double xmin, xmax;
     double jump[];
     int njumps;
{
  double xstep;
  double y;
  int j, m, mm, n;
  int njump_min, njump_max;
  int icheck, pdfsteps;
  double delta;

  njump_min = njump_max = 0;

  /* Check the input data before proceeding

     First check to see that xmin is actually less than xmax

     Then  check that, if there are jumps, they are in order

     Next find the indices into the jump array that contains jumps in between xmin and xmax

     At this point njump_min will point to the first jump which is betewen xmin and
     xmax or it will equal to njumps in which case there were no jumps which were greater
     than xmin.

     Similarly njump_max will be the point to the first jump above xmax or if there are
     no jumps above xmax, then it will be njumps.
   */

  if (xmax <= xmin)
  {
    Error ("pdf_gen_from_func: xmin %g <= xmax %g\n", xmin, xmax);
    Exit (1);
  }
  if (njumps > 0)
  {
    for (j = 1; j < njumps; j++)
    {
      if (jump[j] <= jump[j - 1])
      {
        Error ("pdf_gen_from_func: jump[%d]=%g <=jump[%d]=%g out of order\n", j, jump[j], j - 1, jump[j - 1]);
        Exit (1);
      }
    }
    njump_min = 0;
    while (njump_min < njumps && jump[njump_min] <= xmin)
      njump_min++;
    njump_max = 0;
    while (njump_max < njumps && jump[njump_max] < xmax)
      njump_max++;


    njumps = njump_max - njump_min;
  }

  /* OK, all the input data seems OK */

  /* Construct what is effectively is the definite integral from xmin to x. Note
     however that currently cdf_array[0] corresponds awkwardly to the integral
     from xmin to xstep.

     There code increases the size of the array if the steps appear so large that they
     will affect the digitization.

     Problems might occur if there are very large jumps in the function.

     10oct - ksl - Incorprated into the code in attempting to create a power law.
   */

  delta = 1.;
  n = 0;
  pdfsteps = PDFSTEPS;

  while (n < 3)
  {
    delta = gen_array_from_func (func, xmin, xmax, pdfsteps);
    if (delta < 0.1 / FUNC_CDF) //FUNC_CDF is set in sirocco.h - it is the number of steps we want in the final CDF
      break;
    pdfsteps *= 10;
    n = n + 1;
  }


  xstep = (xmax - xmin) / pdfsteps;

/* At this point cdf_array (an external variable) contains a finely discretised normalized version of
  the CDF for the function it contains pdfsteps points, with a seperation in the x variable of xstep.
  we now 'downsample' the pdf to an array with FUNC_CDF points. We ensure that we have points at locations given
  by the jumps. The resaon for doing this is to limit the size of the final CDF - it will be used a great
  many times, and the original reason for a small array is to speed up searching in the cdf_get_rand bit.
  */


  /* Set the first points in the cdf */

  cdf->x[0] = xmin;
  cdf->y[0] = 0;

  n = 0;                        //This is the position in pdf_array - we start at 0
  mm = 1;                       //This is the index to a desired value of y at the end of the CDF we want it to equal 1!
  j = njump_min;                //This refers to the jumps
  for (m = 1; m < FUNC_CDF; m++)        //Loop over the number of required points in the downsampled CDF
  {
    y = (float) mm / (FUNC_CDF - njumps);       // Desired value of y for the next point in the down sampled CDF

    while (pdf_array[n] < y && n < pdfsteps)    // Work one's way through pdf_array, until we hit the next required value of y
    {
      if (j < njump_max && jump[j] <= xmin + (n + 1) * xstep)   //We have passed a jump before hitting our next required value of y - we want to make a point here
      {
        cdf->x[m] = xmin + (n + 1) * xstep;     //Not exactly jump but close- - we set the next x point in the CDF to the X point we just passed
        cdf->y[m] = pdf_array[n];       //And we set the y point in the CDF to the corresponding Y point
        j++;                    //increment the jump number
        m++;                    //increment the pdf structure number
      }
      n++;                      //This is the index increment for the pdf_array we are going thruogh
    }
    // The while has triggered - we have passed out next required value of y - so we make a point in the new down sampled CDF
    /* So at this point pdf_array[n-1] < x and pdf_array[n]>x */
    cdf->x[m] = xmin + (n + 1) * xstep;
    cdf->y[m] = pdf_array[n];
    mm++;                       // increment the number associated with the desired y ignoring jumps
    /* So pdf->y will contain numbers from 0 to 1 */

  }

  cdf->x[FUNC_CDF] = xmax;      //Set the last x point in the downsampled CDF
  cdf->y[FUNC_CDF] = 1.0;       //And the last y point must equal 1.
  cdf->norm = 1.;               /* pdf_gen_from array produces a properly nomalized cdf and so the
                                   normalization is 1.  110629 ksl */
  cdf->ncdf = FUNC_CDF;         //The number of points is

/* Set the limits for cdf_get_rand_limit to the full array */
  cdf->limit1 = 0;
  cdf->limit2 = 1.0;
  cdf->x1 = cdf->x[0];
  cdf->x2 = xmax;

/* Calculate the gradients */
  if (calc_cdf_gradient (cdf))
  {
    Error ("cdf_gen_from_func: Error returned from calc_cdf_gradient\n");
  }


  /* Check the cdf */
  if ((icheck = cdf_check (cdf)) != 0)
  {
    Error ("cdf_gen_from_function: error %d on cdf_check\n", icheck);
  }
  return (icheck);

}


/**********************************************************/
/**
 * @brief      generates a cdf from a function
 *
 * @param [in] double *func The probablility density function to be integrated  to create the pdf
 * @param [in] double xmin,xmax;	The range over which the function is to be integrated
 * @param [in] int pdfsteps;	number of steps that are calculated
 * @return     delta - the largest change in the cdf between any two points in the grid.
 *
 * @details
 * This is a routine which is called by cdf_gen_from_func which simply calculates the cumulative
 * distribution of the function in equally spaced steps between xmin and xmax.  The CDF is properly
 * normalized.
 *
 * The function to represent the PDF is required to have arguments double and void *, due to the GSL integrator used
 * in num_int. As the void * parameter is unused, we pass NULL instead of an array of extra parameters.
 *
 * ### Notes ###
 *
 * 10oct - ksl -This is rather brute force.  An alternative would have been to have increased the density of points
 * only in the region where it was needed.  This could be done by looking at the regions where the
 * sampling was poor in the initial calculation and addeding new points in in a new calculations.  But
 * this would require one to pass both the values of the CDF and the positions where the values were
 * taken.  This would be  not be hard but would require more extensive modifications to pdf_gen_from func than
 * currently.
 *
 *
 **********************************************************/

double
gen_array_from_func (func, xmin, xmax, pdfsteps)
     double (*func) (double, void *);
     double xmin, xmax;
     int pdfsteps;
{

  double x, z, xstep;
  double sum;
  int m, n;
  int idiag;
  double delta;

  xstep = (xmax - xmin) / pdfsteps;     //The size of the step in the x variable
  idiag = 0;

  /* First, allocate an array for internal use, if that
     has not already been done */


  if (init_pdf == 0 || pdf_steps_current < pdfsteps)
  {

    if (pdf_array != NULL)
      free (pdf_array);

    if ((pdf_array = calloc (sizeof (x), pdfsteps)) == NULL)
    {
      Error ("pdf: Could not allocate space for pdf_array\n");
      Exit (1);
    }
    init_pdf = 1;
    pdf_steps_current = pdfsteps;
  }


  for (n = 0; n < pdfsteps; n++)        //Loop over the required number of steps in the cdf
  {
    x = xmin + (n + 0.5) * xstep;       /*The next value of x - it is midway between points on the required cdf because we
                                           will be adding the new z value onto the cumlative total, we assume the function is best approximated over the
                                           whole rage from x[n] to x[n+1] by using the value of the function between the two */
    if ((z = (*func) (x, NULL)) < 0 || z > VERY_BIG || sane_check (z))  //check the function return is sensible
    {
      Error ("pdf_gen_from_func: probability density %g < 0 at %g\n", z, x);
    }
    if (n == 0)                 //We are at the start of the cdf
      pdf_array[0] = z;
    else
      pdf_array[n] = pdf_array[n - 1] + z;      //increment the cdf by the value of the function at the modpoint between points
    /* Check to see if the integral seems too large */
    if (pdf_array[n] > 1.0e100 && idiag == 0)
    {
      idiag = 1;
      Log ("pdf_gen_from_func: ZOWIE  n %d z %g pdf_array[n] %g x %g\n", n, z, pdf_array[n], x);
      for (m = 0; m < n; m += 100)
      {
        z = (*func) (x, NULL);
        Log ("pdf_gen_from_func: zowie n %d x %g z %g pdf_array %g\n", m, x = xmin + (0.5 * m) * xstep, z, pdf_array[m]);
      }
    }
  }
  /* Thus, pdf_array is proportional to  the definite integral from xmin to x
     (where x=xmin+(n+1)*xstep) */

  sum = pdf_array[pdfsteps - 1];        //The total integral is just the last point

  if (sane_check (sum))
  {
    Error ("pdf_gen_from_func:sane_check Sum %f is NaN\n", sum);
  }

  /* Renormalize the array so that this really is a cumulative distribution function */
  delta = 0;
  for (n = 0; n < pdfsteps; n++)
  {
    pdf_array[n] /= sum;
    if (n > 0)
    {
      x = pdf_array[n] - pdf_array[n - 1];
      if (x > delta)
      {
        delta = x;              /*This is the difference in cumlative probability between
                                   adjacent points - it is a measure of the discretisation. If it is too
                                   large then we may want to use a finer grid to get a smoother cdf. */
      }
    }

  }


  return (delta);
}


#define PDF_ARRAY  28000

double pdf_x[PDF_ARRAY], pdf_y[PDF_ARRAY], pdf_z[PDF_ARRAY];
int pdf_n;



/**********************************************************/
/**
 * @brief      Generate a cumulative distribution function (cdf) from an array of values
 *
 * @param [in] cdf	Pointer to the cdf that will be populated here
 * @param [in] x	Locations of points in supplied data - usually wavelength or freuency
 * @param [in] y	Data - e.g. cross section or emissivity
 * @param [in] n_xy	Number of data points
 * @param [in] xmin	Lower limit to required cdf
 * @param [in] xmin	Upper limit to required cdf
 * @return     echeck	the return value from cdf_check
 *
 * @details
 * Takes two arrays of some values of a function y at locations
 * x and turns it into a cumlative density function to allow
 * random values to be obtained from the function. There is no
 * resampling.
 *
 * ###Notes###
 *
 * cdf_gen_from_array does not have the concept of jumps. All of
 * this needs to be taken care of by the routine that generates
 * the array.  This means the x's should be monotonic, without
 * duplicates, and all the values of y should be 0 or positive.
 * The program exits if these conditions are not fulfilled.
 *
 * In constructing the cdf the routine collapses regions of the
 * pdf that have zero probability of occuring.
 *
 **********************************************************/

int
cdf_gen_from_array (cdf, x, y, n_xy, xmin, xmax)
     CdfPtr cdf;
     double x[], y[];
     int n_xy;
     double xmin, xmax;
{
  int allzero;
  int nmin, nmax, cdf_n;
  int n, nn;
  double sum;
  int echeck, zcheck = 0;


/* Perform various checks on the inputs */


  if (n_xy > NCDF)              //If the supplied data array is larger than the output array - fixed as NCDF
  {
    Error ("cdf_gen_from_array: supplied data %i is larger than default aray size %i - increase NCDF\n", n_xy, NCDF);
    Exit (1);
  }


  if (xmax < xmin)              //This must be a mistake, the limits are reversed
  {
    Error ("cdf_gen_from_array: xmin %g <= xmax %g\n", xmin, xmax);
    Exit (1);
  }

  /*We are now going to crawl down the input array, checking for mistakes.
     Firstly we want to see if all the array values are equal to zero.
   */


  allzero = 0;

  for (n = 1; n < n_xy; n++)    //Go over all elements
  {
    if (x[n] <= x[n - 1])       //One x-point is less than the one above - this would cause problems later
    {
      for (nn = 0; nn < n_xy; nn++)
      {
        Log ("cdf_gen_from_array: %5d %11.6e %11.6e\n", nn, x[nn], y[nn]);
      }
      Error ("cdf_gen_from_array: input x not in ascending order at element %5d/%5d  %11.6e %11.6e\n", n, n_xy, x[n - 1], x[n]);

      Exit (1);
    }
    if (y[n] < 0)               //A negative point!
    {
      for (nn = 0; nn < n_xy; nn++)
      {
        Log ("cdf_gen_from_array: %5d %11.6e %11.6e\n", nn, x[nn], y[nn]);
      }
      Error ("cdf_gen_from_array: probability density %g < 0 at element %d\n", y[n], n);
      Exit (1);
    }
    else if (y[n] > 0)          //At least one point is positive
    {
      allzero = 1;              //So they are not all zero.
    }
  }



  /* ksl - Our requirements are as follows.  We want the cdf to have boundaries given by xmin and xmax, unless
   * the array does not encompass this, in which case we will only include what the array contains.  We also
   * want to suppress regions at the end of the cdf where the pdf was 0
   */


  /*First, find the point in the input array that matches the required start point
     we step up through the input data until we reach a point where the x-point is greater
     than xmin, or we reach the end of the array */



  nmin = 0;
  while (x[nmin] <= xmin && nmin < n_xy)
  {
    nmin++;
  }

  /*After this loop, nmin will point to an x array element larger than xmin. So, to
     interpolate we interpolate between elements xmin-1 and xmin */

  /*This next deals with what happens if we ran off the end of the array in the last loop */

  if (nmin == n_xy)
  {
    Error ("cdf_gen_from_array: xmin greater than all array values\n");
    Exit (1);
  }

  /* if xmin is equal to one of the values in the x array then nmin will point to that value, but otherwise it
     will be slightly larger and we will need to fix this */


  /* The next loop finds the location of the uppermost required x-value in the supplied array */

  nmax = nmin;
  while (x[nmax] < xmax && nmax < n_xy)
  {
    nmax++;
  }


  /* In general, nmax should be one past the last required element so we will interpolate
     between nmax-1 and nmax */

  /*now deal with a few pathological cases */

  if (nmax == nmin)
  {
    cdf_inputs_to_file (x, y, n_xy, xmin, xmax, "diag_cdf_inputs.txt");
//OLD    cdf_to_file (cdf, "cdf_gen_from_array");
    Error ("cdf_gen_from_array: nmin and nmax are identical which is not desirable\n");
//OLD    Exit (1);
  }

  if (nmin == 0 && xmin < x[0]) /*We are requesting a CDF that starts below where we have data
                                   essentially we can assume that the CDF is zero between our originally requested xmin, 
                                   so we change xmin to the lowest point for which we have information */
  {
    xmin = x[0];                //nmin is already set to zero
  }

  if (nmax > n_xy - 1)          /*This is a problem - we have run off the end of the array into places 
                                   where there is no data. If the last array element is within the range we want, then we
                                   can recover by resetting xmax to the last array with data, and setting nmax to n_xy-1 */
  {
    if (x[nmax - 1] <= xmax)    /* hopefully the last X value is less than or equal to xmax */
    {
      nmax = n_xy - 1;          //point to the last point in the array of daya
      xmax = x[n_xy - 1];       //reset the maximum value of x to that          
    }
  }



  /* We now need to deal with the possibility that xmin and xmax are between elements.
     we do this by slightly moving the relvant points in the x/y arrays to make points
     that are exactly at the required values of x. We use linear interpolation. */


  if (nmin > 0 && x[nmin] > xmin)       /*test to check we can sensibly interpolate */
  {
    nmin--;                     //We will be resetting the lower bracketing point 
    linterp (xmin, x, y, n_xy, &y[nmin], 0);
    x[nmin] = xmin;
  }

  if (x[nmax] > xmax)           /*If we are off the end of the data, we cannot interpolate - 
                                   this should really have already been caught tho... */
  {
    linterp (xmax, x, y, n_xy, &y[nmax], 0);
    x[nmax] = xmax;
  }




  /* So at this point we want the array running from nmin to nmax */

// We have changed nmax and nmin, so just check once more to ensure we havent messed up
  if (nmax == nmin)
  {
    Error ("cdf_gen_from_array - modified nmax=nmin followin interpolation at ends\n");
//XXX    Exit (1);
  }


  if (xmax < x[0] || xmin > x[n_xy - 1] || allzero == 0)
  {                             // These are special (probably nonsensical) cases
    cdf->x[0] = xmin;
    cdf->y[0] = 0.;
    cdf->x[1] = xmax;
    cdf->y[1] = 1.;
    cdf->norm = 1.0;
    cdf->ncdf = 1;
    Error ("cdf_gen_from_array: all y's were zero or xmin xmax out of range of array x-- returning uniform distribution %d\n", allzero);
    zcheck = -1;

  }
  else
  {
/* So at this point, have unscaled probability density in pdf_x, pdf_y for the points
 * specified by the input array but we want the cumulative distribution
 * We also have assured that there is one value of pdf_x that corresponds to all of the jumps
 */

    // The following lines perform an integration via the trapezoid rule - each point contains
    // the integral up to that poont, so it starts at 0 and ends at the total

    cdf_n = (nmax - nmin);      //The number of points in the integration
    cdf->x[0] = x[nmin];        //The initial x - point in the cdf
    cdf->y[0] = 0.0;            //The initial y value of the cdf - must be zero at the start
    for (n = 1; n < cdf_n + 1; n++)     //Loop over all the CDF
    {
      cdf->x[n] = x[nmin + n];  //The next point in the cdf is just the next point in the input array - no regridding
      cdf->y[n] = cdf->y[n - 1] + 0.5 * (y[nmin + n - 1] + y[nmin + n]) * (x[nmin + n] - x[nmin + n - 1]);      //Just the average value of y over the interval.
    }
    sum = cdf->y[cdf_n];        //the total integrated pdf
    for (n = 1; n <= cdf_n; n++)        //Loop over the CDF and divide each point by the max (final) point
    {
      cdf->y[n] /= sum;         //this is now a cdf - we go from 0 to 1.        
    }
    /*These next few lines are belt and braces - to ensure the last point is in the right place - it should already be */
    cdf->ncdf = cdf_n;
    cdf->x[cdf->ncdf] = x[nmax];
    cdf->y[cdf->ncdf] = 1.0;
    cdf->norm = sum;

  }


/* Calculate the gradients */
  if (calc_cdf_gradient (cdf))
  {
    Error ("cdf_gen_from_array: Error returned from calc_cdf_gradient\n");
  }

  if ((echeck = cdf_check (cdf)) != 0)
  {
    Error ("cdf_gen_from_array: error %d on cdf_check - check CDF_err.diag\n", echeck);
    cdf_to_file (cdf, "CDF_err.diag");  //output the CDF to a file
//XXX    Exit (1);
  }

  if (nmax == nmin)
  {
    cdf_to_file (cdf, "cdf_gen_from_array with only two elements");
//OLD    Exit (1);
  }

/* Set the limits for cdf_get_rand_limit to the full array */
  cdf->limit1 = 0;
  cdf->limit2 = 1.0;
  cdf->x1 = cdf->x[0];
  cdf->x2 = x[nmax];



  if (zcheck)
  {
    return (zcheck);            // Trap the case where values were initially all zeros
  }

  return (echeck);

}

/**********************************************************/
/**
 * @brief      Generate a single sample from a cdf
 *
 * @param [in] CdfPtr  cdf   a structure which contains a cumulative distribution function.
 *
 * @return   x  a random value fronm the cdf between xmin and xmax
 *
 * @details
 *
 * ### Notes ###
 *
 * The quadratic equation that is solved is intended to assure
 * that the recreated pdf is continuous at the boundaries
 * between the points that ae contained in the cdf structure
 *
 **********************************************************/

double
cdf_get_rand (cdf)
     CdfPtr cdf;
{
  double x, r;
  int i, j;
  double q;
  double a, b, c, s[2];
  int quadratic ();

/* Gnerate a random number and then find the interval n the cdf in which x lies */
  r = random_number (0.0, 1.0); //This *excludes* 0.0 and 1.0.
  i = gsl_interp_bsearch (cdf->y, r, 0, cdf->ncdf);

/* Now calculate a place within that interval - we use the gradient of the 
 * CDF to get a more accurate value between the CDF points
 * 
 * To avoid round-off errors and wasted CPU cycles we check that there
 * is that the gradient change is significant
 */
  q = random_number (0.0, 1.0);

  if (fabs (a = (cdf->d[i + 1] - cdf->d[i])) > 1.e-6)
  {
    a *= 0.5;
    b = cdf->d[i];
    c = (-0.5) * (cdf->d[i + 1] + cdf->d[i]) * q;

    if ((j = quadratic (a, b, c, s)) < 0)
    {
      Error ("cdf_get_rand: No positive roots found %d\n", j);
    }
    else
    {
      q = s[j];
      if (q < 0 || q > 1)
      {
        Error ("cdf_get_rand: q out of range %d  %f\n", j, q);
      }
    }
  }

  x = cdf->x[i] * (1. - q) + cdf->x[i + 1] * q;

  if (!(cdf->x[0] <= x && x <= cdf->x[cdf->ncdf]))
  {
    Error ("cdf_get_rand: %g %d %g %g\n", r, i, q, x);
  }

  return (x);
}





/**********************************************************/
/**
 * @brief      sets limit1 and limit2 so that one can generate distributions
 * from a limited range within a previously created cdf.
 *
 * @param [in] CdfPtr  cdf   A ptr to a cdf structure
 * @param [in] double  xmin   The minimum value to return
 * @param [in] double  xmax   The maximum value to return
 * @return     Always returns 0
 *
 * @details
 *
 * We have defined the cdf in such a way that one can generate a random number 0 and
 * 1, find the place in the y array that corresponds to this and map back onto the
 * x values (which is what is returned. We now want to limit the range of the returned
 * x values, but we need to know in terms of the original array where these values
 * are and the generate a random number between some minimum (>0) and maximum (<1) whichi
 * will allow one to find the x that corresponds. The point of this routine is to
 * find that minimum and maximum value. Having done this one can call cdf_get_rand_limit
 *
 * ### Notes ###
 * Sometimes it is expensive computationally to recalculate the distribution function
 * every time you change the limits.   An example of this is the distribution function
 * for a bb where the distribution function is quite smooth (and can be represented in a
 * dimensionless manner) but calculating the function and the distribution function
 * can be expensive.  This routine together with cdf_gen_rand_limit allow one to sample a part of a
 * cdf.  Note that cdf_gen-from_func or cdf_gen_from_array should have been
 * called previously!
 *
 * xmin and xmax will cause cdf[].limit1 and cdf[].limit2 to be set
 * in such a way that when you call pdf_get_rand_limit the distribution
 * will be sampled only between xmin and xmax.
 *
 * If these routine is called improperly, e.g if
 * xmin and xmax lie outside of the original value of xmin
 * and xmax then the distribution will not extend this far,
 * and the progrem will exit.
 *
 **********************************************************/

int
cdf_limit (cdf, xmin, xmax)
     CdfPtr cdf;
     double xmin, xmax;
{
  int i;
  double q;
  if (cdf->y[cdf->ncdf] != 1.0)
  {
    Error ("cdf_limit: cdf not defined!)");
    Exit (1);
  }
  if (xmin >= cdf->x[cdf->ncdf])
  {
    Error ("cdf_limit: xmin %g > cdf->x[cdf->ncdf] %g\n", xmin, cdf->x[cdf->ncdf]);
  }
  if (xmax <= cdf->x[0])
  {
    Error ("cdf_limit: xmax %g < cdf->x[0] %g\n", xmax, cdf->x[0]);
    Exit (1);
  }

/* Set the limits for the minimum */

  if (xmin <= cdf->x[0])
  {
    cdf->limit1 = 0;
    cdf->x1 = cdf->x[0];
  }
  else
  {
    cdf->x1 = xmin;
    i = 0;
    while (xmin > cdf->x[i])
    {
      i++;
    }
    q = (xmin - cdf->x[i - 1]) / (cdf->x[i] - cdf->x[i - 1]);
    cdf->limit1 = cdf->y[i - 1] + q * (cdf->y[i] - cdf->y[i - 1]);
  }

/* Now set the limits for the maximum */

  if (xmax >= cdf->x[cdf->ncdf])
  {
    cdf->limit2 = 1.0;
    cdf->x2 = cdf->x[cdf->ncdf];
  }
  else
  {
    cdf->x2 = xmax;
    i = cdf->ncdf;
    while (xmax <= cdf->x[i])
    {
      i--;
    }
    q = (xmax - cdf->x[i]) / (cdf->x[i + 1] - cdf->x[i]);
    cdf->limit2 = cdf->y[i] + q * (cdf->y[i + 1] - cdf->y[i]);
  }

  return (0);
}


/**********************************************************/
/**
 * @brief      get a random number for a cdf inside limits set by cdf_limit.
 *
 * @param [in, out] CdfPtr  cdf   A ptr to a cdf structure
 * @return     A random number drawn from the Cdf, but within inside limits which have been set previously
 *
 * @details
 *
 * The basic idea here is that we have created a cdf that covers a broad range of, for example, frequency
 * space.  However, at this point we want photons that only cover some portion of the broad range.
 * Instead of creating a new cdf that just covers a desired frequency space, we only sample a portion of
 * the broad range as defined by cdf_limit (for this particular cdf)
 *
 * ### Notes ###
 *
 **********************************************************/

double
cdf_get_rand_limit (cdf)
     CdfPtr cdf;
{
  double x, r;
  int i, j;
  double q;
  double a, b, c, s[2];
  int quadratic ();
  r = random_number (0.0, 1.0);

  r = r * cdf->limit2 + (1. - r) * cdf->limit1;
  i = r * cdf->ncdf;
  while (cdf->y[i + 1] < r && i < cdf->ncdf - 1)
    i++;
  while (cdf->y[i] > r && i > 0)
    i--;
  while (TRUE)
  {
    q = random_number (0.0, 1.0);
    a = 0.5 * (cdf->d[i + 1] - cdf->d[i]);
    b = cdf->d[i];
    c = (-0.5) * (cdf->d[i + 1] + cdf->d[i]) * q;
    if ((j = quadratic (a, b, c, s)) < 0)
    {
      Error ("pdf_get_rand_limit: No positive roots %d\n", j);
    }
    else
    {
      q = s[j];
      if (q < 0 || q > 1)
      {
        Error ("cdf_get_rand_limit: q out of range %d  %f\n", j, q);
      }
    }

    x = cdf->x[i] * (1. - q) + cdf->x[i + 1] * q;
    if (cdf->x1 < x && x < cdf->x2)
      break;
  }

  return (x);
}

int cdf_write_init = 0;
/**********************************************************/
/**
 * @brief      Write the full structure of the cumulative distribution function to a file
 *
 * @param [in] CdfPtr  cdf   A ptr to a cdf structure
 * @param [in] char  comment[]   A string with a comment
 * @return     Always returns 0
 *
 * @details
 * This is a diagnostic routine that writes a cdf to a file for analysis
 *
 * This routine should really be parallelized so that it only
 * writes to a file from thread 0, but cdf.c has
 * been written so it can be removed from sirocco as a whole, and
 * tested separately
 *
 * ### Notes ###
 **********************************************************/

int
cdf_to_file (cdf, comment)
     CdfPtr cdf;
     char comment[];
{
  FILE *fopen (), *fptr;
  int n;
  if (cdf_write_init == 0)
  {
    fptr = fopen ("cdf_diag.txt", "w");
    cdf_write_init++;
  }
  else
  {
    fptr = fopen ("cdf_diag.txt", "a");
  }
  fprintf (fptr, "# %s\n", comment);
  fprintf (fptr, "# limits (portion.to.sample)   %10.4g %10.4g\n", cdf->limit1, cdf->limit2);
  fprintf (fptr, "# x1 x2  Range(to.be.returned) %10.4g %10.4g\n", cdf->x1, cdf->x2);
  fprintf (fptr, "# norm   Scale.factor          %10.4g \n", cdf->norm);
  fprintf (fptr, "# n x y  1-y d\n");
  for (n = 0; n <= cdf->ncdf; n++)
    fprintf (fptr, "%6d %14.8e %14.8e %14.8e %14.8e\n", n, cdf->x[n], cdf->y[n], 1. - cdf->y[n], cdf->d[n]);
  fclose (fptr);
  return (0);
}

/**********************************************************/
/**
 * @brief      Write the array inputs used for generating a cdf to a file
 *
 * @param [in] double   x[]     An array containing the x variable
 * @param [in] double   y[]     An array containing the y variable 
 * @param [in] double   xmin    The min value in x to be used to generate the cdf
 * @param [in] double   xmax    The max value in x to be used to generate the cdf
 * @param [in] char  filename[]   The name of the file to which the cdf should be written
 * @return     Always returns 0
 *
 * @details
 * This is a diagnostic routine
 *
 * ### Notes ###
 **********************************************************/

int
cdf_inputs_to_file (x, y, n_xy, xmin, xmax, filename)
     double x[], y[];
     int n_xy;
     double xmin, xmax;
     char filename[];
{
  FILE *fopen (), *fptr;
  int n;
  fptr = fopen (filename, "w");
  fprintf (fptr, "# Number of samples in array %d\n", n_xy);
  fprintf (fptr, "# xmin xmax  Range(to.be.returned) %10.4g %10.4g\n", xmin, xmax);
  for (n = 0; n <= n_xy; n++)
    fprintf (fptr, "%5d %14.8e %14.8e\n", n, x[n], y[n]);
  fclose (fptr);
  return (0);
}

/**********************************************************/
/**
 * @brief      a simple routine to check some basics of the cumulative distribution function.
 *
 *
 * @param [in, out] CdfPtr  cdf   A ptr to a cdf structure
 * @return     0 if eveything is OK, a positive number otherwise
 *
 * @details
 * This is a diagnostic designed to see whether a cdf is reasonable
 *
 * ### Notes ###
 *
 **********************************************************/

int
cdf_check (cdf)
     CdfPtr cdf;
{
  int n;
  double x, y;
  int hcheck, icheck, jcheck, kcheck, fcheck;
  hcheck = icheck = jcheck = kcheck = fcheck = 0;
  x = cdf->x[0];
  y = cdf->y[0];
  if (y != 0.0)
  {
    Error ("cdf_check: cumulative distribution function should start at 0 not %e\n", y);
    hcheck = 1;
    Log ("%e %e %e %e\n", cdf->x[0], cdf->y[0], cdf->x[1], cdf->y[1]);
  }
  if (cdf->y[cdf->ncdf] != 1.0)
  {
    Error ("cdf_check: cumulative distribution function should end at 1 not %e\n", cdf->y[cdf->ncdf - 1]);
    icheck = 1;
  }

  for (n = 1; n < cdf->ncdf + 1; n++)
  {
    // Note the equal sign here
    if (x <= cdf->x[n])
      x = cdf->x[n];
    else
    {
      jcheck = 1;
      Error ("cdf_check: x problem n %d x %f pdf->x[n] %f\n", n, x, cdf->x[n]);
    }
    // Note the equal sign here
    if (y <= cdf->y[n])
      y = cdf->y[n];
    else
    {
      kcheck = 1;
      Error ("cdf_check: y problem n %d y %f pdf->y[n] %f\n", n, y, cdf->y[n]);
    }
  }
  if (jcheck == 1)
  {
    Error ("cdf_check: x values in pdf should be monotonic\n");
  }
  if (kcheck == 1)
  {
    Error ("cdf_check: cumulative distribution function should be monotonic\n");
  }

  if (hcheck != 0 || icheck != 0 || jcheck != 0 || kcheck != 0)
  {
    cdf_to_file (cdf, "cdf.diag");
    fcheck = hcheck * 1000 + icheck * 100 + jcheck * 10 + kcheck;
    Error ("cdf_check %d %d %d %d %d\n", fcheck, hcheck, icheck, jcheck, kcheck);
  }

/* Next section is a dangerous attempt to repair the cdf  */

  if (hcheck != 0)
    cdf->y[0] = 0.0;
  if (icheck != 0)
    cdf->y[cdf->ncdf] = 1.0;
  if (jcheck != 0)
  {
    for (n = 0; n < cdf->ncdf; n++)
    {
      if (cdf->x[n] >= cdf->x[n + 1])
        cdf->x[n + 1] = cdf->x[n] + 1.e-20;
    }
  }
  if (kcheck != 0)
  {
    for (n = 0; n < cdf->ncdf; n++)
    {
      if (cdf->y[n] >= cdf->y[n] + 1)
        cdf->y[n + 1] = cdf->y[n] + 1.e-20;
    }
  }

  return (fcheck);
}




/**********************************************************/
/**
 * @brief      Calculate gradients for a cdf to be used to better approximate
 * a cdf calculated at n points
 *
 * @param [in] CdfPtr  cdf   A ptr to a cdf structure
 * @return     Normally returns 0; a positive return indicates a problem
 *
 * @details
 * Calculate the gradient of a cdf at each point in the cdf array in order
 * to allow one to calculate a first order correction
 * to a uniform distibution between the points
 * where the gradient has been calculated.
 *
 * ### Notes ###
 * The gradients are calculated by differencing the cdf.  Thus
 * 	this routine should work both for cdf's calculated
 * 	from functions and from arrays
 * 	the way that the ends are dealt with is not
 * 	great - we should really try to come up with an
 * 	extrapolation rather than just fill in the same gradients
 * 	for the second and penultimate cells into the
 * 	first and last cells.
 *
 *
 **********************************************************/

int
calc_cdf_gradient (cdf)
     CdfPtr cdf;
{
  int n, istat;
  double dx1, dx2, dy1, dy2;

  istat = 0;
  for (n = 1; n < cdf->ncdf; n++)
  {
    dy1 = cdf->y[n] - cdf->y[n - 1];
    dx1 = cdf->x[n] - cdf->x[n - 1];
    dy2 = cdf->y[n + 1] - cdf->y[n];
    dx2 = cdf->x[n + 1] - cdf->x[n];
    if (dx1 != 0.0 && dx2 != 0.0)
    {
      cdf->d[n] = 0.5 * (dy1 / dx1 + dy2 / dx2);
    }
    else if (dx1 != 0.0)
    {
      cdf->d[n] = dy1 / dx1;
    }
    else if (dx2 != 0.0)
    {
      cdf->d[n] = dy2 / dx2;
    }
    else
    {
      cdf->d[n] = 0.0;
      Error ("calc_cdf_gradient: dx1 and dx2 both 0 at  %3d %11.6e\n", n, cdf->x[n]);
      istat = 1;
    }

  }
  /* Fill in the ends */
  //Simple versuion

  cdf->d[0] = cdf->d[1];

  cdf->d[cdf->ncdf] = cdf->d[cdf->ncdf - 1];


  //more complex inear interpolation - might be risky so currently commented out.
//  cdf->d[0] = cdf->d[1] - (cdf->d[2] - cdf->d[1]) / (cdf->x[2] - cdf->x[1]) * (cdf->x[1] - cdf->x[0]); //

//  cdf->d[cdf->ncdf] =
//    cdf->d[cdf->ncdf - 1] + (cdf->d[cdf->ncdf - 2] - cdf->d[cdf->ncdf - 1]) / (cdf->x[cdf->ncdf - 2] -
//                                                                               cdf->x[cdf->ncdf - 1]) * (cdf->x[cdf->ncdf] -
//                                                                                                         cdf->x[cdf->ncdf - 1]);

  return (istat);
}



/**********************************************************/
/**
 * @brief      Given two paralel arrays, x and y, this routine reorders the array
 * in ascending order of x
 *
 * @param [in, out] double *  x   The array that is sorted
 * @param [in, out] double *  y   The parallel array that is sorted in order of x
 * @param [in, out] int  n_xy   The length of the two arrays
 * @return     The length of the input array
 *
 * @details
 *
 * This routine is a routine that is intended to fix a potential
 * problem with an input pdf, which one is intending to convert into
 * a cdf. The cdf genearion routines assume that an input routine is in ascending
 * order of the x variable, often a freqency.  The y variable is usually an unnormalized
 * pdf. This routine checks that the
 * x variables is indeed in ascending order and if not reorders the two arrays.
 *
 * ### Notes ###
 *
 * The routine makes an assumption  about what one wants if the x array contains
 * two x values that are identical.  It assumes one wants to the value with the lower
 * corresponding pdf value to be first in the output arrays (as is genearlly
 * the case whn one has a sharp jump for a photoionization x-section.  This will
 * generate a at the upper edge of a square wave.
 *
 **********************************************************/

int
cdf_array_fixup (x, y, n_xy)
     double *x, *y;
     int n_xy;
{
  int n, m;
  size_t *order;
  double *xx, *yy;

  order = calloc (sizeof (size_t), n_xy);
  xx = calloc (sizeof (double), n_xy);
  yy = calloc (sizeof (double), n_xy);

  gsl_sort_index (order, x, 1, n_xy);

  for (n = 0; n < n_xy; n++)
  {
    xx[n] = x[order[n]];
    yy[n] = y[order[n]];
  }

  /* So now I know the order, but there could be two x values that are the same */

  m = 0;
  x[0] = xx[0];
  y[0] = yy[0];
  for (n = 0; n < n_xy; n++)
  {
    if (m == 0 && xx[n] > 0)
    {
      x[0] = x[n];
      y[0] = y[n];
      m++;
    }
    else if (xx[n] > xx[n - 1])
    {
      x[m] = xx[n];
      y[m] = yy[n];
      m++;
    }
  }
  free (order);
  free (xx);
  free (yy);



  return (m);
}
