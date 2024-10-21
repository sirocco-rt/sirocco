
/***********************************************************/
/** @file  brem.c
 * @author nsh
 * @date   October, 2015
 *
 * @brief  Functions to allow a Bremsstrahlung type input spectrum
 *
 * This was added as part of the Rad-hydro effort and is intended
 * to allow a spectrum to match the assumed spectrum in the Blondin
 * heating and cooling rates. Note that this spectrum is very
 * similar in some ways to a blackbody spectrum - essentially
 * a power law turning into a high frequency exponential drop off
 * and so there is a great deal of similarity with routines in
 * bb.c There is scope here for simplification.
 ***********************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief      The integrand used in num_int to compute the luminosity of a
 *             Bremsstrahlung source
 *
 * @param [in] double  freq		 The frequency at which to compute the
 *                             Bremsstrahlung luminosity
 * @param [in] void *  params  UNUSED parameters required for num_int
 *
 * @return     				 The luminosity at frequency freq
 *
 * @details
 * Since geo.const_agn is a luminosity, this function returns a luminosity
 * with units ergs/s
 *
 * params is unused, but the function pointer num_int requires it. It is re-cast
 * to void  to avoid compiler warnings.
 *
 * ### Notes ###
 * 10/15 - Written by NSH
 *
 **********************************************************/

double
integ_brem (double freq, void *params)
{
  double answer;
  (void) params;
  answer = geo.const_agn * pow (freq, geo.brem_alpha) * exp ((-1.0 * PLANCK * freq) / (BOLTZMANN * geo.brem_temp));
  return (answer);
}



/**********************************************************/
/**
 * @brief      The integrand for integrating a dimensionless Bremsstrahlung
 *             spectrum
 *
 * @param [in] double  alpha	  h*freq/k_b/T -
 * @param [in] void *  params   UNUSED parameters required for num_int
 *
 * @return     					luminosity of Bremsstrahlung function at alpha
 *
 * @details
 * This produces a dimensionless Bremsstrahlung spectrum so that one
 * can break fast changing parts of the spectrum (i.e. the exponential
 * drop off) into chunks to ensure we get a properly defnined spectrum
 * of randomly generated photons.
 *
 * params is unused, but the function pointer num_int requires it. It is re-cast
 * to void  to avoid compiler warnings.
 *
 * ### Notes ###
 * 10/15 - Written by NSH
 *
 **********************************************************/

double
brem_d (double alpha, void *params)
{
  double answer;
  (void) params;
  answer = pow (alpha, geo.brem_alpha) * exp (-1.0 * alpha);
  return (answer);
}



/*These parameters are to set up breaks in the spectrum. The low frequency approximation
 is just a power law, and the high freuqnecy aspproximation is just an exponential. One
 can invert simple functions for these ends to ontasin a random photon frequency. In between
 one needs to explicitly integrate the function */

#define BREM_ALPHAMIN 0.01      // Region below which we will use a low frequency approximation
#define BREM_ALPHAMAX 2.        // Region above which we will use a high frequency approximation
#define BREM_ALPHABIG 100.      //  Region over which can maximmally integrate the Bremsstrahlung function


/* These variables are used to store details of a previously made cdf. If we are getting
 random frequency photons fron a Bremsstrahlung spectrum, we will not want to re-create the cdf
 every time we need a new photon - which could be millions of times - so we only remake the cdf
 if the temperature of the spectrum has changed (unlikely) or the frequency bands have
 changed (this will happen as we move thruogh the photon generation bands) */

int ninit_brem = 0;             // This is a flag to say wether a cdf has already been made
double old_brem_alpha_tiny = 0;
double old_brem_t = 0;          // This is the temperature last used to make a cdf
double old_brem_freqmin = 0;    // This is the lower frequency last used to make a cdf
double old_brem_freqmax = 0;    // This is the lower frequency last used to make a cdf


/* These variables are used in the code, but are made global so they persist and can be re-used
 they are only refedined if the frequency limits, or the temperature has changed */

double brem_alphamin, brem_alphamax;    // The input frequency range in dimensionless values
double cdf_brem_lo, cdf_brem_hi, cdf_brem_tot;  // The precise boundaries in the the bb cdf
double cdf_brem_ylo, cdf_brem_yhi;      // The places in the CDF defined by freqmin & freqmax
double brem_lo_freq_alphamin, brem_lo_freq_alphamax, brem_hi_freq_alphamin, brem_hi_freq_alphamax;      //  the limits to use for the low and high frequency values

/**********************************************************/
/** brem_set is the array that cdf_gen_from_func uses to esablish the
 * specific points in the cdf of the dimensionless Bremsstrahlung function.
 * The intention is get a smoooth spectrum. These used to be called 'jumps'
 **********************************************************/
double brem_set[] = {
  0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
  0.9, 1.9
};





/**********************************************************/
/**
 * @brief      Obtain a random frequency photon from a Bremsstrahlung spectrum
 *
 * @param [in] double  freqmin   Minimum frequency to be generated
 * @param [in] double  freqmax   Maximum frequency to be generated
 * @return     					A random photon frequency between freqmin andf freqmax
 *
 * @details
 * This subroutine is used solely to generate random photon frequencies
 * from a Bremsstrahlung type spectrum. The temperature of the spectrum (defining
 * the high frequency exponential dropoff exp(-hnu/kT)) and the spectrasl index
 * of the low frequency part of the spectrum are defined in the geo structure.
 * This code is heavily based upon planck and shares a great deal of code
 * with it. There is a clear possibility of generalising the code to make some
 * kind of generalised power law plus exponential dropoff spectrum.
 *
 *
 * ### Notes ###
 * 10/15 - Written by NSH
 *
 **********************************************************/

double
get_rand_brem (freqmin, freqmax)
     double freqmin, freqmax;
{
  double freq, alpha, y, brem_alpha_tiny;
  int echeck;

  brem_alpha_tiny = PLANCK * xband.f1[0] / BOLTZMANN / geo.brem_temp;

  if (brem_alpha_tiny > BREM_ALPHAMIN)
    brem_alpha_tiny = BREM_ALPHAMIN / 10.;

  /*The first time of calling this function we produce a CDF which runs from BREM_ALPHAMIN to BREM_ALPHAMAX
     the dmiensinless frequency range over which we cannot use the power law or exponential approximations.
     This only needs to be done one, ans so a flag is set once the CDF has been made */

  if (ninit_brem == 0 || brem_alpha_tiny != old_brem_alpha_tiny)
  {                             /* First time through p_alpha must be initialized */
    if ((echeck = cdf_gen_from_func (&cdf_brem, &brem_d, BREM_ALPHAMIN, BREM_ALPHAMAX, 10, brem_set)) != 0)
    {
      Error ("get_rand_brem: on return from cdf_gen_from_func %d\n", echeck);
    }

    /* We need the integral of the brem function outside of the regions of interest as well */

//    cdf_brem_tot = qromb (brem_d, brem_alpha_tiny, BREM_ALPHABIG, 1e-8);
//    cdf_brem_lo = qromb (brem_d, brem_alpha_tiny, BREM_ALPHAMIN, 1e-8) / cdf_brem_tot;  //position in the full cdf of low frequcny boundary
//    cdf_brem_hi = 1. - qromb (brem_d, BREM_ALPHAMAX, BREM_ALPHABIG, 1e-8) / cdf_brem_tot;       //postion in fhe full hi frequcny boundary

    cdf_brem_tot = num_int (brem_d, brem_alpha_tiny, BREM_ALPHABIG, 1e-8);
    cdf_brem_lo = num_int (brem_d, brem_alpha_tiny, BREM_ALPHAMIN, 1e-8) / cdf_brem_tot;        //position in the full cdf of low frequcny boundary
    cdf_brem_hi = 1. - num_int (brem_d, BREM_ALPHAMAX, BREM_ALPHABIG, 1e-8) / cdf_brem_tot;     //postion in fhe full hi frequcny boundary



    ninit_brem++;               //Set the flag to tell the code we have made the CDF.
    old_brem_alpha_tiny = brem_alpha_tiny;

  }

/* We need to define limits - this is only needed to be done once per frequency band. There is some
  integrations needed, so we check to see if they have changed before doing a load of work.
*/

  if (geo.brem_temp != old_brem_t || freqmin != old_brem_freqmin || freqmax != old_brem_freqmax)
  {

/* set alphamin and alphamax - the dimensionless versions of the frequency range	*/
    brem_alphamin = PLANCK * freqmin / (BOLTZMANN * geo.brem_temp);
    brem_alphamax = PLANCK * freqmax / (BOLTZMANN * geo.brem_temp);

/* set the parameters for which these calculations have been done, so we dont redo them */
    old_brem_t = geo.brem_temp;
    old_brem_freqmin = freqmin;
    old_brem_freqmax = freqmax;


/* we now compute the location in the cdf where these limits occur. We will be using a random number from 0 to 1 to
	select a frequency, so we need to know where 0 and 1 occur! */
    cdf_brem_ylo = cdf_brem_yhi = 1.0;
    if (brem_alphamin < BREM_ALPHABIG)  //There is *some* emission
    {
      //     cdf_brem_ylo = qromb (brem_d, brem_alpha_tiny, brem_alphamin, 1e-8) / cdf_brem_tot;       //The position in full CDF of the upper frequency bound
      cdf_brem_ylo = num_int (brem_d, brem_alpha_tiny, brem_alphamin, 1e-8) / cdf_brem_tot;     //The position in full CDF of the upper frequency bound

      if (cdf_brem_ylo > 1.0)
        cdf_brem_ylo = 1.0;
    }
    if (brem_alphamax < BREM_ALPHABIG)
    {
//      cdf_brem_yhi = qromb (brem_d, brem_alpha_tiny, brem_alphamax, 1e-8) / cdf_brem_tot;       //position in the full cdf of currnt hi frequcny boundary
      cdf_brem_yhi = num_int (brem_d, brem_alpha_tiny, brem_alphamax, 1e-8) / cdf_brem_tot;     //position in the full cdf of currnt hi frequcny boundary

      if (cdf_brem_yhi > 1.0)
        cdf_brem_yhi = 1.0;
    }




    brem_lo_freq_alphamin = brem_alphamin;
    brem_lo_freq_alphamax = brem_alphamax;
    if (brem_lo_freq_alphamax > BREM_ALPHAMIN)  //If the upper frequency bound is in or past the region of the CDF where we need to use the full brem sprecrum
      brem_lo_freq_alphamax = BREM_ALPHAMIN;    //Set an upper bound to the range where we can use the power law approximation

    brem_hi_freq_alphamax = brem_alphamax;
    brem_hi_freq_alphamin = brem_alphamin;
    if (brem_hi_freq_alphamin < BREM_ALPHAMAX)  //If the lower frequency bound is in or below the region of the CDF where we need to use the full brem sprecrum
      brem_hi_freq_alphamin = BREM_ALPHAMAX;    //Set a lower bound to the range where we can use the exponential approximation

/* This test is baciscally asking if there is any part of the frequency range that falls in the
	range whwere we need to use the proper Bremsstrahlung spectrum */

    if (brem_alphamin < BREM_ALPHAMAX && brem_alphamax > BREM_ALPHAMIN)
    {
      cdf_limit (&cdf_brem, brem_alphamin, brem_alphamax);      //We limit the cdf because we might not need the full extent
    }

  }

/* End of section (re)defining limits */

  y = random_number (0.0, 1.0); //Get a random number beween 0 and 1.

  y = cdf_brem_ylo * (1. - y) + cdf_brem_yhi * y;       // y is now in an allowed place in the cdf

/* There are 3 cases to worry about
	The case where everything is in the low frequency limit
	The case where everything is in the normal limit
	The case where some photons are in the low regime and some are
	in the normal regime
*/

  if (y <= cdf_brem_lo || brem_alphamax < BREM_ALPHAMIN)        //We will be selecting a frequency using the lower frequency (PL) approximation
  {
    alpha = get_rand_pow (brem_lo_freq_alphamin, brem_lo_freq_alphamax, geo.brem_alpha);
  }
  else if (y >= cdf_brem_hi || brem_alphamin > BREM_ALPHAMAX)
  {
    alpha = get_rand_exp (brem_hi_freq_alphamin, brem_hi_freq_alphamax);        //We will be selecting a frequency using the high frequency (exp) approximation
  }
  else
  {
    alpha = cdf_get_rand_limit (&cdf_brem);     //We will be using the full CDF approach because we are in the regime where PL and exp are inappropriate
  }

  freq = BOLTZMANN * geo.brem_temp / PLANCK * alpha;    //Get a frequency back from the dimenionless alpha parameter
  if (freq < freqmin || freqmax < freq)
  {
    Error ("get_rand_brem: freq %g out of range %g %g\n", freq, freqmin, freqmax);
  }
  return (freq);
}
