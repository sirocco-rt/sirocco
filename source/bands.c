
/***********************************************************/
/** @file  bands.c
 * @author ksl
 * @date   March  , 2018
 *
 * @brief  Setup the frequency bands used for stratified sampling and
 * for recording spectra in wind cells.
 *
 * The subroutines here deal with two related issues, how to choose how
 * many photons or photon bundles to create in various wavelength ranges,
 * and how to record the spectrum that pass through individual cells.
 *
 * Photons, or photon bundles, in Python are created with certain
 * frequency and carry a certain amount of energy.  Without
 * any kind of stratified sampling all of the photons would be created
 * with the same energy content and frequencies which sample the spectrum
 * in equal energy intervals.  The sum of energies of all of the
 * photons bundles originating in, for example, a stellar photosphere
 * would equal that of the luminosity of the star.
 *
 * Normally, however, Python uses a form of stratified sampling to assure
 * we have photons in various energy bands.  This is important, especially
 * for calculating ionization rates since one needs enough photons
 * above ionization edges to estimate the rate accurately.  In banding,
 * one chooses a number of photons in a given frequency range and changes
 * the energy of the photon bundle so that the proper normalization is
 * obtained.  
 *
 * Once one has decided how many photons in what wave bands to generate
 * one also needs to capture a crude spectrum in each cell of the grid
 * at least for many of the ionization modes we use.  The resolution
 * needs to be good enough to allow for various important ionization 
 * edges, but not so fine that one cannot characterize the slope of the
 * spectrum in each interval.  
 *
 * The routines init_bands in this file sets up the bands to for creating
 * photons based on various user inputs; it then calls freqs_init that 
 * sets up the frequency boundaries for recording the spectra.
 *
 * The two tasks are only loosely related, and in fact while there
 * is a lot of flexibility in for stratified sampling, the spectra
 * in the cells are recorded largely on the basis of the system type.
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

/* Actual structures are in sirocco.h.  Here for reference only.

#define NBANDS 10
struct xbands
{
	double f1[NBANDS],f2[NBANDS];
	double min_fraction[NBANDS];
	double nat_fraction[NBANDS];          // The fraction of the accepted luminosity in this band
	double used_fraction[NBANDS];
	double flux[NBANDS];                  //The "luminosity" within a band
	double weight[NBANDS];
	int nphot[NBANDS];
	int nbands;           // Actual number of bands in use
}
xband;
*/


/**********************************************************/
/** 
 * @brief      This is the routine that initializes the frequency bands    
 * used to ensure that photons are created in sufficient numbers across
 * a large range of frequency during ionization cycles.
 *
 * Python uses a form of stratified sampling in an attempt to assure
 * that there are photon bundles at (high, generally) frequencies
 * in sufficient numbers for accurate ionization equilibria to be
 * calculated.
 * 
 * This is the routine that initializes the bands, and the fraction of photons
 * to be generated in each band.    There are
 * a number of possiblilities for setting up the bands
 *
 * @param [in] int  imode   A switch used for determining how the bands are to be populated
 * @param [out] struct xbands *  band   The structure that is populated with the band information 
 * @return     The routine itself simply returns 0 on success
 *
 * The outputs are passed to other routines through the pointer
 * to xbands.  
 *
 * @details
 *
 * The currently allow modes are
 * *0	Use temperature to define a single band
 * *1	Use f1 and f2 to define a single band
 * *2	Use t,f1 and f2, and hardwired bands to define multiple bands which have
 *    been tuned to be relevant to CVs and other systems where a hot disk
 *    (T_eff of about 100,000 K is assumed)
 * *3 Bands tuned for yso
 * *4 Query the user to specify an arbitrary set of bands
 * *5 Bands set up by nsh to test cloudy
 * *6 Hardwired very wide bands for testing
 * *7 Bands hardwired for AGN paper1 by nsh
 * *8  Define bands in logarithmic intervals
 *
 *
 * ### Notes ###
 *
 * Various choices for banding have been developed over time, often in
 * the process of working on different types of simulations.  In general,
 * it is up to the user (by pure thought or by experimentation) to find
 * a set of frequency ranges and photon fractions in each band so that a 
 * calculation can be carried out efficiently.  Choosing the wrong banding
 * will make the program slow to converge/noisy because one does not have
 * enough photons in the necessary intervals, and in extreme situations 
 * to get the ionization structure wrong.  Evidently, option 4 allows 
 * the user the most flexibility.
 *
 * This routine also calls the routine freqs_init which initializes the
 * boundaries for recording coarse versions of the spectra in each
 * cell. This process is also sometimes referred to (confusingly) as
 * banding.  At present this is more or less hardwired.
 *
 **********************************************************/

int
bands_init (imode, band)
     int imode;                 // A switch used for determining how the bands are to be populated
     struct xbands *band;

{
  int mode;
  int nband;
  double xx;
  double tmax, freqmin, freqmax;
  double t;                     // A temperature which can be used to set absolute limits on the bands
  double f1, f2;                // frequency limits that can overide any other limits
  double fmax;
  double f1_log, f2_log, df;
  int ii;
  char answer[LINELENGTH];
  double band_min_frac;
  double total_min_frac = 0;


  freqmin = VLIGHT / 12000e-8;  /*20000 A */

  tmax = 30000.;                /* This sets a floor on freqmax */

  for (ii = 0; ii < geo.ndomain; ii++)
  {
    if (zdom[ii].twind > tmax)
    {
      tmax = zdom[ii].twind;
    }
  }

  if (geo.tstar > tmax)
    tmax = geo.tstar;
  if (geo.t_bl > tmax && geo.lum_bl > 0.0)
    tmax = geo.t_bl;
  if ((0.488 * teff (1.)) > tmax)
    tmax = 0.488 * teff (1.);
  freqmax = BOLTZMANN * tmax / PLANCK * 10.;
  if (freqmax < 2.0 * 54.418 / HEV)
  {
    Log ("Increasing maximum frequency to twice the Helium edge\n");
    freqmax = 2.0 * 54.418 / HEV;
  }
  else
    Log ("Maximum frequency %8.2e determined by T %8.2e\n", freqmax, tmax);
  geo.tmax = tmax;              /*this a global variable so it is available to the code to make informed guesses as to the possible 
                                   location of any BB driven exponential dropoff in the spectrum */
  t = tmax;
  f1 = freqmin;
  f2 = freqmax;

  /* end of import */

  if (imode == -1)
  {
    // TODO: ideally, the args for rdchoice should be constructed automatically using the enums somehow, to avoid mistakes
    strcpy (answer, "cv");
    mode =
      rdchoice ("Photon_sampling.approach(T_star,cv,yso,AGN,tde_bb,min_max_freq,user_bands,cloudy_test,wide,logarithmic)",
                "0,2,3,7,9,1,4,5,6,8", answer);
  }
  else
  {
    mode = imode;
  }

  if (mode == T_STAR_BAND)
  {
    /* Mode 0 is sets a single band based on the temperature given */
    band->nbands = 1;
    band->f1[0] = BOLTZMANN * t / PLANCK * 0.05;
    band->f2[0] = BOLTZMANN * t / PLANCK * 20.;
    band->min_fraction[0] = 1.0;
  }
  else if (mode == MIN_MAX_FREQ_BAND)
  {
    /* Mode 1 sets a single wide band defined by f1 and f2 */
    band->nbands = 1;
    band->f1[0] = f1;
    band->f2[0] = f2;
    band->min_fraction[0] = 1.0;
  }
  else if (mode == CV_BAND)     /* Traditional cv setup */
  {
    band->nbands = 4;
    band->f1[0] = f1;
    band->f2[0] = band->f1[1] = 13.599 / HEV;
    band->f2[1] = band->f1[2] = 24.588 / HEV;
    band->f2[2] = band->f1[3] = 54.418 / HEV;
    band->f2[3] = f2;
    band->min_fraction[0] = 0;
    band->min_fraction[1] = 0.1;
    band->min_fraction[2] = 0.1;
    band->min_fraction[3] = 0.1;
    if (f1 > band->f2[0])
    {
      Error ("bands_init: f1 (%e) > 13.599/HEV)\n", f1);
      Exit (0);
    }
    if (f2 < band->f2[2])
    {
      Error ("bands_init: f2 (%e) < 54.418/HEV)\n", f2);
      Exit (0);
    }

  }
  else if (mode == YSO_BAND)    /* YSO setup */
  {
    band->nbands = 4;
    band->f1[0] = f1;
    band->f2[0] = band->f1[1] = 1.511 / HEV;
    band->f2[1] = band->f1[2] = 3.3998 / HEV;
    band->f2[2] = band->f1[3] = 6.0000 / HEV;
    band->f2[3] = f2;
    band->min_fraction[0] = 0;
    band->min_fraction[1] = 0.1;
    band->min_fraction[2] = 0.2;
    band->min_fraction[3] = 0.2;
    if (f1 > band->f2[0])
    {
      Error ("bands_init: f1 (%e) > 13.599/HEV)\n", f1);
      Exit (0);
    }
    if (f2 < band->f2[2])
    {
      Error ("bands_init: f2 (%e) < 54.418/HEV)\n", f2);
      Exit (0);
    }

  }
  else if (mode == USER_DEF_BAND)
  {
    rdint ("Photon_sampling.nbands", &band->nbands);

    if (band->nbands > NBANDS)
    {
      Error ("bands: Asking for more bands than allowed (%d). Reducing to maximum value.\n", NBANDS);
      band->nbands = NBANDS;
    }

    Log ("Lowest photon energy is ev (freq) is %f (%.2e)\n", f1 * HEV, f1);
    Log ("Highest photon energy is ev (freq) is %f (%.2e)\n", f2 * HEV, f2);
    Log ("Enter band boundaries in increasing eV, and assure they are between lowest and highest energy\n");

    rddoub ("Photon_sampling.low_energy_limit(eV)", &xx);
    f1 = xx / HEV;

    rddoub ("Photon_sampling.high_energy_limit(eV)", &xx);
    f2 = xx / HEV;

    if (f2 <= f1)
    {
      Error ("bands_int: high energy limit must be greater than low energy limit\n");
      Exit (0);
    }

    Log ("Lowest photon energy is ev (freq) is %f (%.2e)\n", f1 * HEV, f1);
    Log ("Highest photon energy is ev (freq) is %f (%.2e)\n", f2 * HEV, f2);

    band->f1[0] = f1;

    for (nband = 0; nband < band->nbands - 1; nband++)
    {
      rddoub ("Photon_sampling.band_boundary(eV)", &xx);
      band->f2[nband] = band->f1[nband + 1] = xx / HEV;

    }
    band->f2[nband] = f2;

    Log ("Enter minimum fraction of photons in each band. The total must be <= 1\n");

    for (nband = 0; nband < band->nbands; nband++)
    {
      rddoub ("Photon_sampling.band_min_frac", &band_min_frac);
      total_min_frac += band_min_frac;
      if (total_min_frac > 1.0)
      {
        Log ("bands_init: total minimum fraction for all bands > 1: total_min_frac = %f\n", total_min_frac);
        Exit (1);
      }
      band->min_fraction[nband] = band_min_frac;
    }
    for (nband = 0; nband < band->nbands; nband++)
    {
      Log ("For band %i, f1=%10.3e, f2=%10.3e, frac=%.2f\n", nband, band->f1[nband], band->f2[nband], band->min_fraction[nband]);
    }



  }
  else if (mode == CLOUDY_TEST_BAND)    /* Set up to compare with cloudy power law table command note
                                           that this also sets up the weight and photon index for the PL, to ensure a continuous distribution */
  {
    if (geo.agn_ion_spectype != SPECTYPE_CL_TAB)
    {
      Error ("Trying to use a broken power law banding without setting spectype to broken power law - must set spectype to 4\n");
      Exit (0);
    }
    rddoub ("Photon_sampling.low_energy_limit(eV)", &xx);

    if (xx > geo.agn_cltab_low)
    {
      xx = geo.agn_cltab_low / 10.0;
      Log ("Lowest  frequency reset to 1/10 of low frequency break\n");
    }
    f1 = xx / HEV;
    rddoub ("Photon_sampling.high_energy_limit(eV)", &xx);

    if (xx < geo.agn_cltab_hi)
    {
      xx = geo.agn_cltab_hi * 10.0;
      Log ("highest  frequency reset to 10x high frequency break\n");
    }
    f2 = xx / HEV;

    if (f2 <= f1)
    {
      Error ("bands_int: high energy limit must be greater than low energy limit\n");
      Exit (0);
    }
    Log ("Lowest photon energy is ev (freq) is %f (%.2e)\n", f1 * HEV, f1);
    Log ("Highest photon energy is ev (freq) is %f (%.2e)\n", f2 * HEV, f2);


    band->nbands = 12;

    band->f1[0] = (geo.agn_cltab_low / HEV) / 1000.0;
    band->f2[0] = band->f1[1] = (geo.agn_cltab_low / HEV) / 100.0;
    band->f2[1] = band->f1[2] = (geo.agn_cltab_low / HEV) / 10.0;
    band->f2[2] = (geo.agn_cltab_low / HEV);

/* Now set up a set of log spaced bands in the range over the central range */
    f1_log = log10 (geo.agn_cltab_low / HEV);
    f2_log = log10 (geo.agn_cltab_hi / HEV);
    df = (f2_log - f1_log) / (6);

    ii = 3;
    while (ii < 9)
    {
      band->f1[ii] = pow (10., f1_log + (ii - 3) * df);
      band->f2[ii] = pow (10., f1_log + ((ii - 3) + 1) * df);
      ii++;
    }

    f1_log = log10 (geo.agn_cltab_hi / HEV);
    f2_log = log10 (f2);
    df = (f2_log - f1_log) / (3);

    ii = 9;
    while (ii < 12)
    {
      band->f1[ii] = pow (10., f1_log + (ii - 9) * df);
      band->f2[ii] = pow (10., f1_log + ((ii - 9) + 1) * df);
      ii++;
    }

    //Set number of photons in each band

    band->min_fraction[0] = 0.0666;
    band->min_fraction[1] = 0.0666;
    band->min_fraction[2] = 0.0666;

    band->min_fraction[3] = 0.1;
    band->min_fraction[4] = 0.1;
    band->min_fraction[5] = 0.1;
    band->min_fraction[6] = 0.1;
    band->min_fraction[7] = 0.1;
    band->min_fraction[8] = 0.1;

    band->min_fraction[9] = 0.0666;
    band->min_fraction[10] = 0.0666;
    band->min_fraction[11] = 0.0666;

    //Set alpha for each band

    band->alpha[0] = geo.agn_cltab_low_alpha;
    band->alpha[1] = geo.agn_cltab_low_alpha;
    band->alpha[2] = geo.agn_cltab_low_alpha;
    band->alpha[3] = geo.alpha_agn;
    band->alpha[4] = geo.alpha_agn;
    band->alpha[5] = geo.alpha_agn;
    band->alpha[6] = geo.alpha_agn;
    band->alpha[7] = geo.alpha_agn;
    band->alpha[8] = geo.alpha_agn;
    band->alpha[9] = geo.agn_cltab_hi_alpha;
    band->alpha[10] = geo.agn_cltab_hi_alpha;
    band->alpha[11] = geo.agn_cltab_hi_alpha;

    //Set the constant for each band to ensure continuous distribution

    band->pl_const[0] = geo.const_agn * pow ((band->f2[2]), geo.alpha_agn) / pow ((band->f2[2]), band->alpha[0]);
    band->pl_const[1] = geo.const_agn * pow ((band->f2[2]), geo.alpha_agn) / pow ((band->f2[2]), band->alpha[0]);
    band->pl_const[2] = geo.const_agn * pow ((band->f2[2]), geo.alpha_agn) / pow ((band->f2[2]), band->alpha[0]);
    band->pl_const[3] = geo.const_agn;
    band->pl_const[4] = geo.const_agn;
    band->pl_const[5] = geo.const_agn;
    band->pl_const[6] = geo.const_agn;
    band->pl_const[7] = geo.const_agn;
    band->pl_const[8] = geo.const_agn;

    band->pl_const[9] = geo.const_agn * pow ((band->f2[8]), geo.alpha_agn) / pow ((band->f2[8]), band->alpha[9]);

    band->pl_const[10] = geo.const_agn * pow ((band->f2[8]), geo.alpha_agn) / pow ((band->f2[8]), band->alpha[9]);

    band->pl_const[11] = geo.const_agn * pow ((band->f2[8]), geo.alpha_agn) / pow ((band->f2[8]), band->alpha[9]);


    for (nband = 0; nband < band->nbands; nband++)
      Log ("f1=%e,f2=%e,alpha=%e,const=%e,lum1=%e,lum2=%e\n",
           band->f1[nband], band->f2[nband], band->alpha[nband],
           band->pl_const[nband],
           band->pl_const[nband] * pow (band->f1[nband],
                                        band->alpha[nband]), band->pl_const[nband] * pow (band->f2[nband], band->alpha[nband]));
  }

  else if (mode == WIDE_BAND)   //Test for balance to have a really wide frequency range
  {
    tmax = geo.tstar;
    fmax = tmax * WIEN;         //Use wiens law to get peak frequency
    band->nbands = 17;


    band->f1[0] = 1e10;
    band->f2[0] = band->f1[1] = fmax * 0.01;
    band->f2[1] = band->f1[2] = fmax * 0.1;
    band->f2[2] = band->f1[3] = fmax;
    band->f2[3] = band->f1[4] = fmax * 1.5;
    band->f2[4] = band->f1[5] = fmax * 2;
    band->f2[5] = band->f1[6] = fmax * 2.5;
    band->f2[6] = band->f1[7] = fmax * 3;
    band->f2[7] = band->f1[8] = fmax * 4;
    band->f2[8] = band->f1[9] = fmax * 6;
    band->f2[9] = band->f1[10] = fmax * 8;
    band->f2[10] = band->f1[11] = fmax * 10;
    band->f2[11] = band->f1[12] = fmax * 12;
    band->f2[12] = band->f1[13] = fmax * 14;
    band->f2[13] = band->f1[14] = fmax * 16;
    band->f2[14] = band->f1[15] = fmax * 18;
    band->f2[15] = band->f1[16] = fmax * 20;
    band->f2[16] = 1e20;

    band->min_fraction[0] = 0.1;
    band->min_fraction[1] = 0.1;
    band->min_fraction[2] = 0.1;
    band->min_fraction[3] = 0.05;
    band->min_fraction[4] = 0.05;
    band->min_fraction[5] = 0.05;
    band->min_fraction[6] = 0.05;
    band->min_fraction[7] = 0.05;
    band->min_fraction[8] = 0.05;
    band->min_fraction[9] = 0.05;
    band->min_fraction[10] = 0.05;
    band->min_fraction[11] = 0.05;
    band->min_fraction[12] = 0.05;
    band->min_fraction[13] = 0.05;
    band->min_fraction[14] = 0.05;
    band->min_fraction[15] = 0.05;
    band->min_fraction[16] = 0.05;
  }

  else if (mode == AGN_BAND)    //Test for balance matching the bands we have been using for AGN runs
  {

    band->nbands = 10;
    band->f1[0] = 1e14;
    band->f2[0] = band->f1[1] = 1e15;
    band->f2[1] = band->f1[2] = 3.162e15;
    band->f2[2] = band->f1[3] = 1e16;
    band->f2[3] = band->f1[4] = 3.162e16;
    band->f2[4] = band->f1[5] = 1e17;
    band->f2[5] = band->f1[6] = 3.162e17;
    band->f2[6] = band->f1[7] = 1e18;
    band->f2[7] = band->f1[8] = 3.162e18;
    band->f2[8] = band->f1[9] = 1e19;
    band->f2[9] = 1e20;

    band->min_fraction[0] = 0.1;
    band->min_fraction[1] = 0.1;
    band->min_fraction[2] = 0.1;
    band->min_fraction[3] = 0.1;
    band->min_fraction[4] = 0.1;
    band->min_fraction[5] = 0.1;
    band->min_fraction[6] = 0.1;
    band->min_fraction[7] = 0.1;
    band->min_fraction[8] = 0.1;
    band->min_fraction[9] = 0.1;

  }
  else if (mode == TDE_BB_BAND)
  {
    band->nbands = 7;
    band->f1[0] = 1e14;
    band->f2[0] = band->f1[1] = 1e15;
    band->f2[1] = band->f1[2] = 3.162e15;
    band->f2[2] = band->f1[3] = 1e16;
    band->f2[3] = band->f1[4] = 3.162e16;
    band->f2[4] = band->f1[5] = 1e17;
    band->f2[5] = band->f1[6] = 3.162e17;
    band->f2[6] = band->f1[7] = 1e18;

    for (ii = 0; ii < band->nbands; ++ii)
      band->min_fraction[ii] = (int) (1.0 / band->nbands);
  }
  else if (mode == LOG_USER_DEF_BAND)   /* Generalized method to set up logarithmic bands */
  {
    Log ("Lowest photon energy is ev (freq) is %f (%.2e)\n", f1 * HEV, f1);
    Log ("Highest photon energy is ev (freq) is %f (%.2e)\n", f2 * HEV, f2);

    band->nbands = 5;
    rdint ("Photon_sampling.nbands", &band->nbands);

    if (band->nbands > NBANDS)
    {
      Error ("bands: Asking for more bands than allowed (%d). Reducing to maximum value\n", NBANDS);
      band->nbands = NBANDS;
    }

    xx = f1 * HEV;
    rddoub ("Photon_sampling.low_energy_limit(eV)", &xx);
    f1 = xx / HEV;

    xx = f1 * HEV;
    rddoub ("Photon_sampling.high_energy_limit(eV)", &xx);
    f2 = xx / HEV;

    if (f2 <= f1)
    {
      Error ("bands_int: high energy limit must be greater than low energy limit\n");
      Exit (0);
    }

    f1_log = log10 (f1);
    f2_log = log10 (f2);
    df = (f2_log - f1_log) / (band->nbands);
    ii = 0;
    while (ii < band->nbands)
    {
      band->f1[ii] = pow (10., f1_log + ii * df);
      band->f2[ii] = pow (10., f1_log + (ii + 1) * df);
      band->min_fraction[ii] = 1. / band->nbands;
      ii++;
    }

  }
  else
  {
    Error ("bands_init: Unknown mode %d\n", mode);
    Exit (0);
  }

  if (geo.agn_ion_spectype == SPECTYPE_CL_TAB && mode != CLOUDY_TEST_BAND)
  {
    /* There's quite a lot of hardwired behaviour that means you need to use cloudy_test banding 
       if you want a Cloudy broken power-law SED */
    Error ("Using Cloudy broken-law: this only works with Photon_sampling.approach set to cloudy_test.\n");
    Exit (0);
  }


  Log ("bands_init: There are %d bands\n", band->nbands);
  for (nband = 0; nband < band->nbands; nband++)
  {
    Log ("bands_init: band %i,  f1=%10.3e,  f2=%10.3e, frac=%.2f\n", nband, band->f1[nband], band->f2[nband], band->min_fraction[nband]);
    Log ("bands_init: band %i, eV1=%10.3e, eV2=%10.3e, frac=%.2f\n", nband,
         band->f1[nband] * HEV, band->f2[nband] * HEV, band->min_fraction[nband]);
    Log ("bands_init: band %i, alpha1=%f, alpha2=%f, frac=%.2f\n", nband,
         band->f1[nband] * PLANCK / (BOLTZMANN * tmax), band->f2[nband] * PLANCK / (BOLTZMANN * tmax), band->min_fraction[nband]);
  }

  check_appropriate_banding (band, mode);

  /* we used to call freqs init here, but now the photon generation bands are 
     also tied to the ionization bands, see e.g. gh issue #1084 */
  ion_bands_init (mode, band->f1[0], band->f2[band->nbands - 1], band);

  /* Now define the freqquency boundaries for the cell spectra */
  geo.cell_log_freq_min = log10 (band->f1[0]);
  geo.cell_log_freq_max = log10 (band->f2[band->nbands - 1]);
  geo.cell_delta_lfreq = (geo.cell_log_freq_max - geo.cell_log_freq_min) / NBINS_IN_CELL_SPEC;

  return (0);
}




/**********************************************************/
/** 
 * @brief      This is the routine where the frequency
 * 	binning for coarse spectra in each plasma cell is established
 * 
 * @param [in] int  mode   The banding mode
 * @param [in] double  freqmin   The minimum frequency
 * @param [in] double  freqmax   The maximum frequency
 * @param [in] struct xbands *band bands structure to check
 * @return     Always returns 0.  
 *
 * The frequency intervals are stored
 * in geo.xfreq.  The total number of frequencies is geo.nxfreq
 *
 * @details
 * 
 * In order to approximate the radiation field in each plasma cell, one
 * needs to record a coarse spectrum.  There is no point in creaating 
 * a detailed spectrum since this will (usually) not be justified by the
 * photon statistics and in any event wouuld require to much memory. Memory
 * restrictions also mean we cannot simple record the effects of individual
 * photons.  
 *
 * So that one can sensibly record a coarse spectum, one needs to define
 * the frequency boundaries for the course spectra.  That is the purpose
 * of this routine.
 *
 * ### Notes ###
 *
 * At present everything is hardwired.  In the case of SPECTYPE_CL_TAB,
 * a special broken power law used for comparisons, a different set
 * of bands are used.
 *
 * freqmin and freqmax are used in order to limit the total range of
 * the spectral bands.
 *
 * Although some the description of this routine refer to banding
 * this routine does have anything to do with the stratified sampling
 * used to create photons.  
 *
 **********************************************************/
#define MIN_N_IONBANDS 7

int
ion_bands_init (mode, freqmin, freqmax, band)
     int mode;
     double freqmin, freqmax;
     struct xbands *band;
{
  int i, n, ngood, good[NXBANDS];
  double xfreq[NXBANDS];
  int nxfreq, nband;


/* At present set up a single energy band for 2 - 10 keV */
/*NSH 70g - bands set up to match the bands we are currently using in the.pf files. This should probably end up tied together in the long run! */
/* nxfreq = 7;
 xfreq[0] = 1.0 / HEV;
 xfreq[1] = 13.6 / HEV;
 xfreq[2] = 54.42 / HEV;
 xfreq[3] = 392. / HEV;
 xfreq[4] = 739. / HEV;
 xfreq[5] = 2000 / HEV;
 xfreq[6] = 10000 / HEV;
 xfreq[7] = 50000 / HEV;*/

  /* bands to match the cloudy table spectrum - needed to cover all frequencies to let induced compton work OK */
  if (mode == CLOUDY_TEST_BAND)
  {
    nxfreq = 3;
    xfreq[0] = 0.0001 / HEV;
    xfreq[1] = geo.agn_cltab_low / HEV;
    xfreq[2] = geo.agn_cltab_hi / HEV;
    xfreq[3] = 100000000 / HEV;
  }
  /* if we have a sufficient number of user-selected bands, then just use those band boundaries */
  else if (band->nbands >= MIN_N_IONBANDS)
  {
    for (nband = 0; nband < band->nbands; nband++)
    {
      xfreq[nband] = band->f1[nband];
    }
    nxfreq = band->nbands;
    xfreq[band->nbands] = band->f2[band->nbands - 1];
  }

  /* if we don't have enough bands, then use the old hardwired bands */
  else
  {
    nxfreq = 10;
    xfreq[0] = freqmin;         //We need the whole range to be modelled for induced compton heating to work
    xfreq[1] = 1e15;            //This should be below the lowest threshold frequency of any element in our model
    xfreq[2] = 3.162e15;
    xfreq[3] = 1e16;
    xfreq[4] = 3.162e16;
    xfreq[5] = 1e17;
    xfreq[6] = 3.162e17;
    xfreq[7] = 1e18;
    xfreq[8] = 3.162e18;
    xfreq[9] = 1.2e19;          //This is the highest frequency defined in our ionization data
    xfreq[10] = freqmax;

    Log ("ion_bands_init: only %d generation bands, so adopting default 10 ionization bands from %8.4e Hz to %8.4e Hz\n", band->nbands,
         freqmin, freqmax);
  }

  Log ("ion_bands_init: %d ionization bands from %8.2feV (%8.2eHz) to %8.2feV (%8.2eHz)\n", nxfreq, freqmin * HEV, freqmin, freqmax * HEV,
       freqmax);

  ngood = 0;
  for (i = 0; i < nxfreq; i++)
  {
    if (freqmin < xfreq[i] && xfreq[i] < freqmax)
    {
      good[i] = 1;
      ngood++;
    }
    else if (freqmin < xfreq[i + 1] && xfreq[i + 1] < freqmax)
    {
      good[i] = 1;
      ngood++;
    }
    else
    {
      good[i] = 0;
    }
  }

  Log_silent ("ion_bands_init: Of %d starting intervals, %d will have photons\n", nxfreq, ngood);

  n = 0;
  for (i = 0; i < nxfreq; i++)
  {
    if (good[i] == 1)
    {
      geo.xfreq[n] = xfreq[i];
      geo.xfreq[n + 1] = xfreq[i + 1];
      n++;
    }
  }
  geo.nxfreq = n;

  /* OK at this point we know at least some photons will be generated in each interval, but we still don't know
   * that the we are going to have a possibilty of photons throughout the first and last intervals.
   */

  if (freqmin > geo.xfreq[0])
  {
    geo.xfreq[0] = freqmin;
  }

  if (freqmax < geo.xfreq[geo.nxfreq])
  {
    geo.xfreq[geo.nxfreq] = freqmax;
  }


  Log_silent ("ion_bands_init: There were %d final intervals\n", geo.nxfreq);
  for (n = 0; n < geo.nxfreq; n++)
  {
    Log_silent ("ion_bands_init: %8.2f (%8.2e)    %8.2f (%8.2e)  \n", geo.xfreq[n] * HEV, geo.xfreq[n], geo.xfreq[n + 1] * HEV,
                geo.xfreq[n + 1]);
  }



  return (0);

}

/**********************************************************/
/** 
 * @brief Helper routine for checking the photon generation
 * banding is reasonable. At this stage fairly rudimentary
 *
 * @param [in] struct xbands *band bands structure to check
 * @param [in] int  mode   The banding mode
 *
 * @details
 * 
 **********************************************************/

void
check_appropriate_banding (band, mode)
     struct xbands *band;
     int mode;
{
  if (geo.system_type == SYSTEM_TYPE_AGN)
  {
    if (mode == CV_BAND)
      Error ("Using CV banding for AGN system. Not recommended!\n");
    else if (mode == YSO_BAND)
      Error ("Using YSO banding for AGN system. Not recommended!\n");
    else if (mode == T_STAR_BAND)
      Error ("Using Tstar banding for AGN system. Not recommended!\n");
    else if (band->nbands < 4)
      Error ("You only have %d photon generation bands for AGN system. Not recommended!\n");
  }
  else if (geo.system_type == SYSTEM_TYPE_CV)
  {
    if (mode == YSO_BAND)
      Error ("Using YSO banding for CV system. Not recommended!\n");
    else if (band->nbands < 4)
      Error ("You only have %d photon generation bands for CV system. Not recommended!\n");
  }
  else if (geo.system_type == SYSTEM_TYPE_BH)
  {
    if (mode == CV_BAND)
      Error ("Using CV banding for BH system. Not recommended!\n");
    else if (mode == YSO_BAND)
      Error ("Using YSO banding for BH system. Not recommended!\n");
    else if (mode == T_STAR_BAND)
      Error ("Using Tstar banding for BH system. Not recommended!\n");
    else if (band->nbands < 4)
      Error ("You only have %d photon generation bands for BH system. Not recommended!\n");
  }
}
