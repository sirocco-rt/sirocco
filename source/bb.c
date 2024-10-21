
/***********************************************************/
/** @file  bb.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  Routines concerning the generation of photons from blackbodies
 *
 * ### Notes ##
 *
 * The two main routines are emittance_bb or planck.  The rest of the routines
 * are helper routines that should not normally be called directly.
 *
 * The first call to either of these routines (surely emittance_bb) results in a call to integ_planck_d,
 * which in turn calls integ_planck_init.  This initialization routine in turn populates the
 * array integ_planck , which contains the integral of the bb function in an array.
 * bb_emittance continues to access the array integ_plank through integ_planck_d every
 * future time is is called.

 * planck does the same thing albeit more indirectly. It sets up a cdf each time new frequency
 * limits are placed on it.  planck therefore really uses the
 * cdf.
 *
 * Both emittence_bb and planck accept inputs in physical units, frequencies and temperatures. Internally
 * howeve the routines convert physical to dimensionless units.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "atomic.h"
#include "sirocco.h"

#define ALPHAMIN 0.4            // Region below which we will use a low frequency approximation
#define ALPHAMAX 30.            // Region above which we will use a high frequency approximation
#define ALPHABIG 100.           // Region over which can maximmally integrate the Planck function
#define NJUMPS 30               //The number of 'jumps' - places in the CDF where we want to force points
#define NMAX 		1000    //The number of points at which the planck function is integrated between ALPHAMIN and ALPHAMAX for storage/

int ninit_planck = 0;           //A flag to say wether we have computed our stored blackbody integral

double bb_old_t = 0;
double bb_old_freqmin = 0;
double bb_old_freqmax = 0;
double alphamin, alphamax;
double cdf_bb_lo, cdf_bb_hi, cdf_bb_tot;        // The precise boundaries in the the bb cdf
double cdf_bb_ylo, cdf_bb_yhi;  // The places in the CDF defined by freqmin & freqmax
double lo_freq_alphamin, lo_freq_alphamax, hi_freq_alphamin, hi_freq_alphamax;  //  the limits to use for the low and high frequency values


/**********************************************************/
/**
 * Array that cdf_gen_from_func uses to establish the
 * specific points in the cdf of the dimensionless bb function.
 *
 * These are what we call 'jumps' and are used by cdf_gen_from_func to
 * ensure important parts of the CDF have points
 **********************************************************/
double bb_set[] = {

  10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
  19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 29.99999
};


int error_bb_hi = 0;
int error_bb_lo = 0;


/**********************************************************/
/**
 * @brief      returns the frequency for a photon which follows a Planck distribution
 * within defined frequency limits
 *
 * @param [in] double  t   The temperature of the bb
 * @param [in] double  freqmin   The minimum frequency for the photon
 * @param [in] double  freqmax   The maximum frequency for the photon
 * @return     The frequency drawn randomly from a BB function with
 *      temperature T
 *
 * @details
 *
 * The frequency is pseudo random in the following limited sense.  The photons are
 * returned are weighted so that the energy distribution of a  function is approximately
 * reproduced. Within an energy bin the photon frequency is uniform.
 *
 * ### Notes ###
 *
 * The first time one enters the program, a cdf for a diminensionless
 * BB function is createda.  The cdf is created for values between ALPHAMIN and
 * ALPHAMAX where ALPHA=h nu/kT. The number of points in the cdf is determined by NMAX
 *
 * On subseqent entries, when the temperature
 * or frequency limits are changed, we use standard routines to limit
 * what portion of the dimensionless cdf to use.
 *
 * If the frequency range and temperature for a photon falls outside of ALPHAMIN
 * and ALPHAMAX special routines are used to sample the distribution there.
 * so that the spectrum falls off as a power law at low frequencies and an
 * exponential at high frequencies.
 *
 *
 **********************************************************/

double
planck (t, freqmin, freqmax)
     double t, freqmin, freqmax;
{
  double freq, alpha, y;
  int echeck;


  if (t <= 0)
  {
    Error ("planck: A value of %e for t is unphysical\n", t);
    return (freqmin);
  }
  /*First time through create the array containing the proper boundaries for the integral
   * of the BB function, Note calling cdf_gen_from func also defines ylo and yhi */

  if (ninit_planck == 0)
  {                             /* First time through p_alpha must be initialized */
    if ((echeck = cdf_gen_from_func (&cdf_bb, &planck_d_2, ALPHAMIN, ALPHAMAX, 21, bb_set)) != 0)
    {
      Error ("Planck: on return from cdf_gen_from_func %d\n", echeck);
    }
    /* We need the integral of the bb function outside of the regions of interest as well
     *
     * cdf_bb_lo is the position in the full cdf of the low frequcny boundary
     * cdf_bb_hi is position in the full cdf of the hi frequcny boundary
     * */


    cdf_bb_tot = num_int (planck_d, 0, ALPHABIG, 1e-8);
    cdf_bb_lo = num_int (planck_d, 0, ALPHAMIN, 1e-8) / cdf_bb_tot;
    cdf_bb_hi = 1. - num_int (planck_d, ALPHAMAX, ALPHABIG, 1e-8) / cdf_bb_tot;


    ninit_planck++;

  }


/* If temperatures or frequencies have changed since the last call to planck
 * redefine various limits, including the portion of the cdf to be used
*/

  if (t != bb_old_t || freqmin != bb_old_freqmin || freqmax != bb_old_freqmax)
  {

    alphamin = PLANCK * freqmin / (BOLTZMANN * t);
    alphamax = PLANCK * freqmax / (BOLTZMANN * t);

    bb_old_t = t;
    bb_old_freqmin = freqmin;
    bb_old_freqmax = freqmax;

    cdf_bb_ylo = cdf_bb_yhi = 1.0;

    if (alphamin < ALPHABIG)    //check to make sure we get a sensible number - planck_d(ALPHAMAX is too small to sensibly integrate)
    {
      cdf_bb_ylo = num_int (planck_d, 0, alphamin, 1e-8) / cdf_bb_tot;  //position in the full cdf of current low frequency boundary

      if (cdf_bb_ylo > 1.0)
        cdf_bb_ylo = 1.0;
    }
    if (alphamax < ALPHABIG)    //again, check to see that the integral will be sensible
    {
      cdf_bb_yhi = num_int (planck_d, 0, alphamax, 1e-8) / cdf_bb_tot;  //position in the full cdf of currnet hi frequency boundary

      if (cdf_bb_yhi > 1.0)
        cdf_bb_yhi = 1.0;
    }

/* These variables are not always used */

    lo_freq_alphamin = alphamin;        //Set the minimum frequency to use the low frequency approximation to the lower band limit
    lo_freq_alphamax = alphamax;        //Set to a default value

    if (lo_freq_alphamax > ALPHAMIN)    //If the upper alpha for this band is above the loew frequency approximation lower limit
      lo_freq_alphamax = ALPHAMIN;      //Set the maximum alpha we will use the low frequency approximation to the default value

    hi_freq_alphamax = alphamax;        //Set the maximum frequency to use the high frequency approximation to to the upper band limit
    hi_freq_alphamin = alphamin;        //Set to a default value
    if (hi_freq_alphamin < ALPHAMAX)    //If the lower band limit is less than the high frequency limit
      hi_freq_alphamin = ALPHAMAX;      //Se the minimum alpha value to use the high frequency limit to the default value


/* Check whether the limits for alpha min and max are within the 'normal' bb range and if so
* set the portion of the full cdf to use.
*
* Note that alphamin is always less than alphamax.
*/

    if (alphamin < ALPHAMAX && alphamax > ALPHAMIN)
    {
      cdf_limit (&cdf_bb, alphamin, alphamax);
    }

  }
  /* End of section redefining limits */


  y = random_number (0.0, 1.0); //We get a random number between 0 and 1 (excl)

  y = cdf_bb_ylo * (1. - y) + cdf_bb_yhi * y;   // y is now in an allowed place in the cdf


/* There are 3 cases to worry about
 *	The case where everything is in the low frequency limit
 *	The case where everything is in the normal limit
 *	The case where some photons are in the low regime and some are
 *	in the normal regime
*/

  if (y <= cdf_bb_lo || alphamax < ALPHAMIN)    //we are in the low frequency limit
  {
    alpha = get_rand_pow (lo_freq_alphamin, lo_freq_alphamax, 2.);
  }
  else if (y >= cdf_bb_hi || alphamin > ALPHAMAX)       //We are in the high frequency limit
  {
    alpha = get_rand_exp (hi_freq_alphamin, hi_freq_alphamax);
  }
  else
  {
    alpha = cdf_get_rand_limit (&cdf_bb);       //We are in the region where we use the BB function
  }

  freq = BOLTZMANN * t / PLANCK * alpha;
  if (freq < freqmin || freqmax < freq)
  {
    Error ("planck: freq %g out of range %g %g\n", freq, freqmin, freqmax);
  }
  return (freq);
}



/**********************************************************/
/**
 * @brief      obtains a random number between x1 and x2
 *  	for a power law densiity distribution with index alpha
 *
 * @param [in, out] double  x1   the minimum allowed value to return
 * @param [in, out] double  x2   the maximum allowed value to return
 * @param [in, out] double  alpha   the index of the power law
 * @return     A single value taken from a power law distribution
 *
 * @details
 *
 *
 * ### Notes ###
 *
 * It is straight fowrward to integrate a power law
 * in closed form and given this it one can easily
 * find the value of the upper limit to the definite
 * integral that is a given fraction of the total interval
 *
 * The case with alpha==1 is a special case, but this
 * is also easy to integrate.
 *
 **********************************************************/

double
get_rand_pow (x1, x2, alpha)
     double x1, x2, alpha;
{
  double r;
  double a;

  r = random_number (0.0, 1.0); //This produces a random number between 0 and 1 excl

  if (alpha == -1)
  {
    x1 = log (x1);
    x2 = log (x2);
    a = (1. - r) * x1 + r * x2;
    a = exp (a);
  }
  else
  {
    x1 = pow (x1, alpha + 1.);
    x2 = pow (x2, alpha + 1.);

    a = (1. - r) * x1 + r * x2;

    a = pow (a, 1. / (alpha + 1.));
  }
  return (a);
}



/**********************************************************/
/**
 * @brief      obtain a random number between alpha_min
 * and alpha_max from an distribution of the form e**-alpha
 *
 * @param [in] double  alpha_min   The lower limit for the random number
 * @param [in] double  alpha_max   The upper limit for the random number
 * @return     A double precision random number
 *
 * @details
 *
 * ### Notes ###
 *
 * The cdf for an exponential distribution can be easily
 * shown to be given by solving this equation for alpha
 *
 * r*(exp(-alpha_min)-exp(-alpha_max))=(exp(-alpha_min)-exp(alpha))
 *
 * but you can recast this and solve for delta_alpha
 *
 * exp(-delta_alpha)=(1-R)+R*exp(-(alpha_max-alpha_min))
 *
 * This has the advantage that it depends only on the
 * difference of alpha_max and alpha_min and not their
 * actual values, and as long as the exp of a very
 * large number turns out to be zero and not
 * not a number, it shuld not generate NANs
 *
 **********************************************************/

double
get_rand_exp (alpha_min, alpha_max)
     double alpha_min, alpha_max;
{
  double r;
  double x;
  double a, aa;
  double delta_alpha;

  r = random_number (0.0, 1.0);

  x = exp (alpha_min - alpha_max);


  aa = (1. - r) + r * x;
  delta_alpha = -(log (aa));

  a = alpha_min + delta_alpha;

  if (sane_check (a))
  {
    Error ("get_rand_exp:sane_check: alpha_min %e alpha_max %e --> %e %e %e %e %e\n", alpha_min, alpha_max, a, aa, delta_alpha, x, r);
    return (alpha_min);
  }
  return (a);
}

// The dimensionless planck function integrated from ALPHAMIN to a range of values of alpha
double integ_planck[NMAX + 1];

// A flag to say whether we have initialised integ_planck.
int i_integ_planck_d = 0;

/**********************************************************/
/**
 * @brief      Obtains the integral of the dimensionless blackbody function
 *    between alphamin and alphamax
 *
 * @param [in] double  alphamin   The minimum value of alpha (h nu/ kT) considered
 * @param [in] double  alphamax   The maximum value of alpha to be considered
 * @return     The value of the integral
 *
 * @details
 *
 * To save computing time, the routine actually accesses an array that contains the
 * the integral of the bb function from 0 to alpha for a number of values of alpha.
 * By differencing (and intepolating) this array for the element corrsponding
 * to alphamax and alphamin, the routine returns the integral between alphamax
 * and alphamin
 *
 *
 *
 * ### Notes ###
 *
 * The first time the function is called the integ_planck[] is filled.
 *
 *
 * integ_plank[n]  contains the integral of the
 * dimensionless planck function from  ALPHAMIN to ALPHAMAX.  Therefore if one
 * wants to obtain the integral of the dimensionless bb function, one simply interpolates
 * this array.
 *
 **********************************************************/

double
integ_planck_d (alphamin, alphamax)
     double alphamin, alphamax;
{

  double x, z1, z2;
  int n;
  int init_integ_planck_d ();

  /* If this is the first time, integ_plank_d is called, initilaize the integ_planck array */
  if (i_integ_planck_d == 0)
  {
    init_integ_planck_d ();
    i_integ_planck_d++;
  }


  /* Find the array elements that are associated with alphamin and interpolate to find
   * the value of the integrated planck functiion at alphamin
   *
   * If alphamin is off the bottom of the tabulated values we set the integral for x1 to 0
   * and continue on to x2
   *
   * If alphamin is off the top end of the array we assume result is 0, and return what is
   * essentially an error condition, since alphamax will also be too high to use as well.
   */

  x = (alphamin - ALPHAMIN) / (ALPHAMAX - ALPHAMIN) * NMAX;
  if (x <= 0.0)
    z1 = 0.0;
  else if (x >= (NMAX))
  {
    return (0.0);
  }
  else
  {
    n = x;                      //n is the array element
    x -= n;                     //x is now the fractional distance between array elements
    z1 = integ_planck[n] * (1. - x) + integ_planck[n + 1] * x;  //Interpolate
  }

  /* Now find the array elements associated with alphamax
   *
   * If the array element associated with alphamax is negative, the highest frequencey
   * is below the bottom of the tabulated values so return 0, which is essentially
   * an error.
   *
   * if the array element is above the maximum alpha, we use the maximum alpha.  Otherwise
   * (the normal case) we interpolate
   *
   */

  x = (alphamax - ALPHAMIN) / (ALPHAMAX - ALPHAMIN) * NMAX;

  if (x < 0.0)
  {
    return (0.0);
  }

  else if (x >= (NMAX))
  {
    z2 = integ_planck[NMAX];
  }
  else
  {
    n = x;
    x -= n;
    z2 = integ_planck[n] * (1. - x) + integ_planck[n + 1] * x;  //Interpolate
  }

  /* z1 and z2 are the integrals of the dimensionless BB at alphamin and alphamax
     Return the difference between the integral at alphamax and alphamin */

  return (z2 - z1);
}




/**********************************************************/
/**
 * @brief      calculates integrals of the dimensionless bb function from
 * 0 to alpha and the stores the reusults in an array used by integ_plank_d
 *
 * @return     Always returns 0
 *
 * @details
 *
 * This is a helper routine which is only called once from integ_planck in order to initialize
 * the array integ_planck.
 *
 * It calculates the integral of the dimensionless bb function between
 * 0 and alpha for a series of alpha values values between ALPHAMIN and
 * ALPHAMAX.
 *
 * ### Notes ###
 *
 * ALPHAMIN and ALPHAMAX are  hardcoded
 * which are hardcorded.
 *
 *
 **********************************************************/

int
init_integ_planck_d ()
{
  double x;
  double planck_d (), qromb ();
  int n;
  for (n = 0; n < NMAX + 1; n++)
  {
    x = ALPHAMIN + n * (ALPHAMAX - ALPHAMIN) / NMAX;
    integ_planck[n] = num_int (planck_d, 0.0, x, 1e-7);
  }

  return (0);
}



#define EPSILON	1.e-6


/**********************************************************/
/**
 * @brief      The value of the dimensioless BB function at alpha
 *
 * @param [in] double  alpha
 * @param [in] void  params   An extra (unused) variable to make it paletable for the gsl integrator

 * @return   The value of the dimensionless BB function at x
 *
 * @details
 *
 * The BB emittence is given by
 *
 * F_nu= 2*PI* (kT/h**3)*planck_d(h*freq/ k T)
 *
 * alpha is simply h nu/kT
 *
 * ### Notes ###
 *
 **********************************************************/

double
planck_d (double alpha, void *params)
{
  return (planck_d_2 (alpha, params));
}


double
planck_d_2 (double alpha, void *params)
{
  double x;
  (void) params;
  if (alpha < EPSILON || alpha > ALPHABIG)
    return (0);
  x = (alpha * alpha * alpha) / (exp (alpha) - 1);
  return (x);
}


/**********************************************************/
/**
 * @brief      Calculate the emittance of a bb between freqmin and freqmax
 *
 * @param [in] double  freqmin   The minimum frequency
 * @param [in] double  freqmax   The maximum frequency
 * @param [in] double  t   The temperature at which the calculation is made
 * @return     the emmittance between freqmin and freqmax
 *
 * @details
 *
 * ### Notes ###
 * This should integrate to sigma  if freqmin and freqmax go effeictively
 * from 0 to infinity
 *
 **********************************************************/

double
emittance_bb (freqmin, freqmax, t)
     double freqmin, freqmax, t;
{
  double alphamin, alphamax, q1;
  double integ_planck_d ();
  q1 = 2. * PI * (BOLTZMANN * BOLTZMANN * BOLTZMANN * BOLTZMANN) / (PLANCK * PLANCK * PLANCK * VLIGHT * VLIGHT);

  alphamin = PLANCK * freqmin / (BOLTZMANN * t);
  alphamax = PLANCK * freqmax / (BOLTZMANN * t);


  if (alphamin > ALPHAMIN && alphamax < ALPHAMAX)
  {
    return (q1 * t * t * t * t * integ_planck_d (alphamin, alphamax));
  }
  else if (alphamax > ALPHABIG)
  {
    if (alphamin > ALPHABIG)
      return (0);
    else
      return (q1 * t * t * t * t * num_int (planck_d, alphamin, ALPHABIG, 1e-7));
  }
  else
  {
    return (q1 * t * t * t * t * num_int (planck_d, alphamin, alphamax, 1e-7));

  }
}


/**********************************************************/
/**
 * @brief      decides whether a maximum frequency requested for an integral is sensible.
 * If it is too far off the end of the planck function, the numerical integration routine  will malfunction. We
 * just have to set it to a frequency where the BB function is tiny, say where hnu/kT =100.
 *
 * @param [in] double  freq_max   The maximum frequency
 * @param [in] double  temp   The temperature
 * @return   A frequency which is the maximum value for which one should try to evaluate the
 * BB function
 *
 * @details
 * We use ALPHABIG to define the place in the BB spectrum where we want to give up
 *
 * ### Notes ###
 *
 * This helper routine was written by NSH in August 2012. We
 * were having lots of problems with trying to integrate cross sections
 * of a planck function at temperatures where there is no flux at fmax. 
 *
 * The routine just checks if fmax is going to be more than 100xt/(h/k)
 * which will return tiny numbers for planck, and may cause the numerical
 * integration routine to return return
 * nonsense.
 * 
 * The routine is called from function calc_pi_rate
 *
 *
 **********************************************************/

double
check_freq_max (freq_max, temp)
     double freq_max, temp;
{
  double bblim;


  bblim = ALPHABIG * (temp / H_OVER_K);

  if (bblim < freq_max)
  {
    freq_max = bblim;
  }

  return (freq_max);

}


#undef NMAX
#undef ALPHAMIN
#undef ALPHAMAX
