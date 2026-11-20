
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
 * Photons, or photon bundles, in Sirocco are created with certain
 * frequency and carry a certain amount of energy.  Without
 * any kind of stratified sampling all of the photons would be created
 * with the same energy content and frequencies which sample the spectrum
 * in equal energy intervals.  The sum of energies of all of the
 * photons bundles originating in, for example, a stellar photosphere
 * would equal that of the luminosity of the star.
 *
 * Normally, however, Sirocco uses a form of stratified sampling to assure
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
 * The routine bands_init in this file sets up the bands for creating
 * photons based on various user inputs; it then calls ion_bands_init,
 * which sets up the frequency boundaries for recording the spectra.
 *
 * The two tasks are only loosely related.  The bands in which the
 * spectra in the cells are recorded are chosen separately, with
 * Wind.ionization.band_model, and are the photon generation bands
 * only if that is set to photon_bands.
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
 * Sirocco uses a form of stratified sampling in an attempt to assure
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
 * This routine also calls the routine ion_bands_init which initializes the
 * boundaries for recording coarse versions of the spectra in each
 * cell. This process is also sometimes referred to (confusingly) as
 * banding.  How those boundaries are chosen is set by
 * Wind.ionization.band_model.
 *
 **********************************************************/

int
bands_init (int imode, struct xbands *band)
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

  if (cell_spec_on_thresholds ())
  {
    cell_spec_threshold_grid ();
    return (0);
  }

  geo.cell_delta_lfreq = (geo.cell_log_freq_max - geo.cell_log_freq_min) / geo.nbins_in_cell_spec;

  if (geo.cell_freq != NULL)
    free (geo.cell_freq);
  geo.cell_freq = calloc (geo.nbins_in_cell_spec + 1, sizeof (double));

  for (ii = 0; ii <= geo.nbins_in_cell_spec; ii++)
  {
    geo.cell_freq[ii] = pow (10., (geo.cell_log_freq_min + ii * geo.cell_delta_lfreq));
  }

  return (0);
}


/**********************************************************/
/**
 * @brief      Whether the cell spectra have bin boundaries at every
 * photoionization threshold
 *
 * @return     TRUE if the boundaries include the thresholds
 *
 * @details
 * The cell spectra (cell_spec_flux) are binned on a uniform log grid,
 * except when they are used to calculate photoionization rates: with
 * band model cell_spectrum, and for the development table written with
 * -pi_rate_table.  Then every ground-state and inner-shell threshold is
 * also a bin boundary, so that no ion's own edge falls inside a bin.
 **********************************************************/

int
cell_spec_on_thresholds ()
{
  if (geo.band_model == BAND_MODEL_CELL_SPECTRUM || modes.pi_rate_table == TRUE)
  {
    return (TRUE);
  }
  return (FALSE);
}


/**********************************************************/
/**
 * @brief      Compare two doubles, for qsort
 *
 * @param [in] const void *a   The first double
 * @param [in] const void *b   The second double
 * @return     -1, 0 or 1 as the first is less than, equal to or greater than the second
 **********************************************************/

static int
cmp_double (const void *a, const void *b)
{
  double x;
  double y;

  x = *(const double *) a;
  y = *(const double *) b;
  if (x < y)
  {
    return (-1);
  }
  if (x > y)
  {
    return (1);
  }
  return (0);
}


/**********************************************************/
/**
 * @brief      The photoionization thresholds inside the range of
 * the cell spectra
 *
 * @param [out] double *thresh   The thresholds, sorted, or NULL just to count them
 * @return     The number of distinct thresholds
 *
 * @details
 * The thresholds are those of the ground-state cross section of each ion
 * (as used by calc_pi_rate) and of every inner-shell cross section, that
 * lie strictly inside the range of the cell spectra.  Thresholds within
 * 1e-6 dex of each other count as one.  The number depends only on the
 * atomic data and the frequency range.
 **********************************************************/

static int
cell_spec_thresholds (double *thresh)
{
  double *t;
  double freq;
  double fmin;
  double fmax;
  int nt;
  int nion;
  int n;
  int nunique;

  fmin = pow (10., geo.cell_log_freq_min);
  fmax = pow (10., geo.cell_log_freq_max);

  t = calloc (nions + n_inner_tot + 1, sizeof (double));
  nt = 0;

  for (nion = 0; nion < nions; nion++)
  {
    freq = -1.0;
    if (ion[nion].phot_info > 0 && ion[nion].ntop_ground >= 0)
    {
      freq = phot_top[ion[nion].ntop_ground].freq[0];
    }
    else if (ion[nion].phot_info == 0 && ion[nion].nxphot >= 0)
    {
      freq = phot_top[ion[nion].nxphot].freq[0];
    }
    if (freq > fmin && freq < fmax)
    {
      t[nt] = freq;
      nt++;
    }
  }

  for (n = 0; n < n_inner_tot; n++)
  {
    freq = inner_cross[n].freq[0];
    if (freq > fmin && freq < fmax)
    {
      t[nt] = freq;
      nt++;
    }
  }

  qsort (t, nt, sizeof (double), cmp_double);

  nunique = 0;
  for (n = 0; n < nt; n++)
  {
    if (nunique == 0 || log10 (t[n] / t[nunique - 1]) > 1e-6)
    {
      t[nunique] = t[n];
      nunique++;
    }
  }

  if (thresh != NULL)
  {
    for (n = 0; n < nunique; n++)
    {
      thresh[n] = t[n];
    }
  }

  free (t);
  return (nunique);
}


/**********************************************************/
/**
 * @brief      Set the cell spectrum bins to a log grid with every
 * photoionization threshold added as a boundary
 *
 * @return     Always 0
 *
 * @details
 * For a new run, geo.nbins_in_cell_spec is the number of log bins
 * (-cell_spec_dim, or 1000 with band model cell_spectrum); the thresholds
 * are added to these.  On a restart, geo.nbins_in_cell_spec is the total
 * number of bins (read from the wind save, or as reset by -cell_spec_dim or
 * band model cell_spectrum), and the number of log bins is that total less
 * the number of thresholds.  The program stops if that leaves no log bins,
 * or if the total exceeds NBINS_IN_CELL_SPEC.
 **********************************************************/

int
cell_spec_threshold_grid ()
{
  double *thresh;
  double *logbins;
  int nthresh;
  int nlog;
  int ntot;
  int i;
  int j;
  int k;

  nthresh = cell_spec_thresholds (NULL);

  if (geo.run_type == RUN_TYPE_NEW)
  {
    nlog = geo.nbins_in_cell_spec;
  }
  else
  {
    nlog = geo.nbins_in_cell_spec - nthresh;
  }

  ntot = nlog + nthresh;
  if (nlog < 1 || ntot > NBINS_IN_CELL_SPEC)
  {
    Error
      ("cell_spec_threshold_grid: %d log bins and %d thresholds do not fit (maximum %d bins in all). A restart with band model cell_spectrum or -pi_rate_table needs a saved run that used one of them\n",
       nlog, nthresh, NBINS_IN_CELL_SPEC);
    Exit (1);
  }

  thresh = calloc (nthresh + 1, sizeof (double));
  cell_spec_thresholds (thresh);

  geo.cell_delta_lfreq = (geo.cell_log_freq_max - geo.cell_log_freq_min) / nlog;
  logbins = calloc (nlog + 1, sizeof (double));
  for (i = 0; i <= nlog; i++)
  {
    logbins[i] = pow (10., (geo.cell_log_freq_min + i * geo.cell_delta_lfreq));
  }

  if (geo.cell_freq != NULL)
  {
    free (geo.cell_freq);
  }
  geo.cell_freq = calloc (ntot + 1, sizeof (double));

  /* merge the log boundaries and the thresholds, which are both sorted */
  i = 0;
  j = 0;
  k = 0;
  while (k <= ntot)
  {
    if (j >= nthresh || (i <= nlog && logbins[i] <= thresh[j]))
    {
      geo.cell_freq[k] = logbins[i];
      i++;
    }
    else
    {
      geo.cell_freq[k] = thresh[j];
      j++;
    }
    k++;
  }

  geo.nbins_in_cell_spec = ntot;

  Log ("cell_spec_threshold_grid: cell spectra have %d bins: %d log bins plus %d photoionization thresholds\n", ntot, nlog, nthresh);

  free (logbins);
  free (thresh);
  return (0);
}


#define MIN_N_IONBANDS 7

/* Edges used to place ionization band boundaries when there are fewer
   than MIN_N_IONBANDS photon generation bands.  Each entry is (z, istate)
   of the ion whose ground-state threshold is the edge; (-1, 0) is the
   Balmer edge.  The order is the order of priority: an edge is used only
   if it lies inside the photon range, and edges are added until there are
   MIN_N_IONBANDS bands. */

static int ion_band_edge_list[][2] = {
  {1, 1}, {2, 1}, {2, 2}, {-1, 0},      /* H I, He I, He II, Balmer */
  {6, 3}, {6, 4}, {8, 6}, {7, 5}, {14, 4},      /* C III, C IV, O VI, N V, Si IV */
  {8, 7}, {8, 8}, {6, 5}, {6, 6}, {7, 6}, {7, 7},       /* O VII, O VIII, C V, C VI, N VI, N VII */
  {10, 9}, {10, 10}, {12, 11}, {14, 13}, {14, 14}, {16, 15},    /* Ne IX, Ne X, Mg XI, Si XIII, Si XIV, S XV */
  {26, 25}, {26, 26}            /* Fe XXV, Fe XXVI */
};

#define BALMER_EDGE_EV 3.40
#define ION_BAND_MIN_DEX 0.05   /* no edge is added closer than this (in log nu) to an existing boundary */


/**********************************************************/
/**
 * @brief      The ground-state photoionization threshold of an ion
 *
 * @param [in] int  z   The element, or -1 for the Balmer edge
 * @param [in] int  istate   The ionization state
 * @return     The threshold frequency in Hz, or -1 if the ion, or its
 * photoionization cross section, is not in the atomic data
 **********************************************************/

static double
ion_band_edge_freq (int z, int istate)
{
  int nion;
  int n;

  if (z < 0)
  {
    return (BALMER_EDGE_EV / HEV);
  }

  for (nion = 0; nion < nions; nion++)
  {
    if (ion[nion].z == z && ion[nion].istate == istate)
    {
      /* For macro-atom ions, find the cross section whose lower level is the ground
         level directly from the levels, so that the result does not depend on
         ion[].ntop_ground (which was wrong for macro-atom ions before the fix in
         get_atomic_data) */
      if (ion[nion].macro_info == 1)
      {
        for (n = 0; n < nphot_total; n++)
        {
          if (phot_top[n].nion == nion && phot_top[n].nlev == ion[nion].first_nlte_level)
          {
            return (phot_top[n].freq[0]);
          }
        }
        return (-1.0);
      }
      if (ion[nion].phot_info > 0 && ion[nion].ntop_ground >= 0)
      {
        return (phot_top[ion[nion].ntop_ground].freq[0]);
      }
      if (ion[nion].phot_info == 0 && ion[nion].nxphot >= 0)
      {
        return (phot_top[ion[nion].nxphot].freq[0]);
      }
      return (-1.0);
    }
  }
  return (-1.0);
}


/**********************************************************/
/**
 * @brief      Add a band boundary, keeping the boundaries in order
 *
 * @param [in, out] double *edges   The boundaries, in increasing order
 * @param [in] int  nedges   The number of boundaries
 * @param [in] double  freq   The boundary to add
 * @return     The new number of boundaries
 *
 * @details
 * The boundary is not added if it lies outside the first and last
 * boundaries, or within ION_BAND_MIN_DEX of an existing boundary.
 **********************************************************/

static int
ion_band_add_edge (double *edges, int nedges, double freq)
{
  int i;
  int k;

  if (freq <= edges[0] || freq >= edges[nedges - 1])
  {
    return (nedges);
  }

  for (i = 0; i < nedges; i++)
  {
    if (fabs (log10 (edges[i] / freq)) < ION_BAND_MIN_DEX)
    {
      return (nedges);
    }
  }

  k = nedges;
  while (k > 0 && edges[k - 1] > freq)
  {
    edges[k] = edges[k - 1];
    k--;
  }
  edges[k] = freq;

  return (nedges + 1);
}


/**********************************************************/
/**
 * @brief      Halve the widest band until there are enough bands
 *
 * @param [in, out] double *xfreq   The band boundaries, in increasing order
 * @param [in] int  nedges   The number of boundaries
 * @param [in] int  nbands   The number of bands wanted
 * @return     The new number of boundaries, nbands + 1
 *
 * @details
 * The widest band, measured in log nu, is split at its geometric mean
 * frequency, and this is repeated until there are nbands bands.
 **********************************************************/

static int
ion_band_split_to (double *xfreq, int nedges, int nbands)
{
  int k;
  int kmax;
  double freq;
  double width;
  double wmax;

  while (nedges - 1 < nbands)
  {
    wmax = -1.0;
    kmax = 0;
    for (k = 0; k < nedges - 1; k++)
    {
      width = log (xfreq[k + 1] / xfreq[k]);
      if (width > wmax)
      {
        wmax = width;
        kmax = k;
      }
    }
    freq = sqrt (xfreq[kmax] * xfreq[kmax + 1]);
    for (k = nedges; k > kmax + 1; k--)
    {
      xfreq[k] = xfreq[k - 1];
    }
    xfreq[kmax + 1] = freq;
    nedges++;
  }

  return (nedges);
}


/**********************************************************/
/**
 * @brief      Ionization bands from the generation bands and the
 * important photoionization edges
 *
 * @param [in] struct xbands *band   The photon generation bands
 * @param [out] double *xfreq   The band boundaries
 * @return     The number of bands
 *
 * @details
 * Used when there are fewer than MIN_N_IONBANDS generation bands.
 * The boundaries of the generation bands are kept, since photons are
 * generated in each of them.  Boundaries are then added at the edges in
 * ion_band_edge_list, in order of priority, if they lie inside the photon
 * range, until there are MIN_N_IONBANDS bands.  If there are still too few
 * bands, the widest band (in log nu) is halved until there are enough.
 *
 * Placing the boundaries on edges matters because the spectrum in
 * each band is described by a smooth model, which cannot follow a sharp
 * edge inside a band (#1084, #1169).
 **********************************************************/

static int
ion_bands_from_edges (struct xbands *band, double *xfreq)
{
  int nedges;
  int nband;
  int i;
  int nlist;
  double freq;

  nedges = 0;
  xfreq[nedges] = band->f1[0];
  nedges++;
  for (nband = 0; nband < band->nbands; nband++)
  {
    xfreq[nedges] = band->f2[nband];
    nedges++;
  }

  nlist = sizeof (ion_band_edge_list) / sizeof (ion_band_edge_list[0]);
  for (i = 0; i < nlist && nedges - 1 < MIN_N_IONBANDS; i++)
  {
    freq = ion_band_edge_freq (ion_band_edge_list[i][0], ion_band_edge_list[i][1]);
    if (freq > 0)
    {
      nedges = ion_band_add_edge (xfreq, nedges, freq);
    }
  }

  nedges = ion_band_split_to (xfreq, nedges, MIN_N_IONBANDS);

  return (nedges - 1);
}


/* The elements whose ground-state thresholds become band boundaries in
   ion_bands_from_elements, in order of priority (roughly the order of
   their abundance): H, He, O, C, N.  Each element's ions are taken from
   the neutral up. */

static int ion_band_elements[] = { 1, 2, 8, 6, 7 };

#define ION_BANDS_MAX_THRESH (NXBANDS / 2)      /* at most this many boundaries are placed on thresholds */


/**********************************************************/
/**
 * @brief      Ionization bands from the boundaries given by the user
 * (Wind.ionization.band_model user_bands)
 *
 * @param [in] struct xbands *band   The photon generation bands
 * @param [out] double *xfreq   The band boundaries
 * @return     The number of bands
 *
 * @details
 * The bands run from the bottom of the first generation band to the top
 * of the last, with the user's inner boundaries (geo.band_user_freq)
 * between.  A boundary outside that range is dropped with a warning.
 * There is no splitting: the user gets the bands asked for.
 **********************************************************/

static int
ion_bands_from_user (struct xbands *band, double *xfreq)
{
  int n;
  int nedges;
  double fmin, fmax;

  fmin = band->f1[0];
  fmax = band->f2[band->nbands - 1];
  nedges = 0;
  xfreq[nedges] = fmin;
  nedges++;
  for (n = 0; n < geo.band_user_n - 1; n++)
  {
    if (geo.band_user_freq[n] > fmin && geo.band_user_freq[n] < fmax)
    {
      xfreq[nedges] = geo.band_user_freq[n];
      nedges++;
    }
    else
    {
      Error ("ion_bands_from_user: band boundary %.3f eV is outside the photon range %.3f to %.3f eV, so it is dropped\n",
             geo.band_user_freq[n] * HEV, fmin * HEV, fmax * HEV);
    }
  }
  xfreq[nedges] = fmax;
  Log ("ion_bands_from_user: %d bands from the boundaries given\n", nedges);
  return (nedges);
}


/**********************************************************/
/**
 * @brief      Ionization bands from the thresholds of the abundant
 * elements, with the remaining bands made by splitting the widest bands
 *
 * @param [in] struct xbands *band   The photon generation bands
 * @param [out] double *xfreq   The band boundaries
 * @return     The number of bands, NXBANDS
 *
 * @details
 * The bands cover the whole photon generation range, from the bottom of
 * the first generation band to the top of the last; the generation band
 * boundaries themselves are not kept.
 *
 * Boundaries are first placed on the ground-state photoionization
 * thresholds of H, He, O, C and N, element by element and within each
 * element from the neutral ion up, skipping any threshold outside the
 * range or within ION_BAND_MIN_DEX of a boundary already placed, until
 * there are ION_BANDS_MAX_THRESH of them or the list is used up.  The
 * widest band (in log nu) is then halved until there are NXBANDS bands.
 *
 * This is the rule "elem+split" of the banding study (Sirocco-2611-Banding):
 * offline, against the fine cell histograms of an AGN shell, a star and a
 * CV, it gave photoionization rates of the H, He, C, N and O ions within
 * 0.1 dex (90th percentile), against up to 3 dex for the generation bands.
 * The bands are global and fixed for the whole run, and need no knowledge
 * of the spectrum in the cells.
 **********************************************************/

static int
ion_bands_from_elements (struct xbands *band, double *xfreq)
{
  int nedges;
  int nelem;
  int i;
  int istate;
  double freq;

  nedges = 0;
  xfreq[nedges] = band->f1[0];
  nedges++;
  xfreq[nedges] = band->f2[band->nbands - 1];
  nedges++;

  nelem = sizeof (ion_band_elements) / sizeof (ion_band_elements[0]);
  for (i = 0; i < nelem; i++)
  {
    for (istate = 1; istate <= ion_band_elements[i]; istate++)
    {
      if (nedges - 2 >= ION_BANDS_MAX_THRESH)
      {
        break;
      }
      freq = ion_band_edge_freq (ion_band_elements[i], istate);
      if (freq > 0)
      {
        nedges = ion_band_add_edge (xfreq, nedges, freq);
      }
    }
  }

  Log ("ion_bands_from_elements: %d boundaries on photoionization thresholds of H, He, O, C and N\n", nedges - 2);

  nedges = ion_band_split_to (xfreq, nedges, NXBANDS);

  return (nedges - 1);
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
 * How the bands are chosen is set by geo.band_model
 * (Wind.ionization.band_model): on the thresholds of H, He, O, C and N
 * with the rest made by splitting (ion_bands_from_elements; the default,
 * and also used with cell_spectrum), the photon generation bands (with
 * edges added if there are fewer than MIN_N_IONBANDS,
 * ion_bands_from_edges), or the boundaries given by the user
 * (ion_bands_from_user).  In the case of SPECTYPE_CL_TAB, a special
 * broken power law used for comparisons, a different set of bands is
 * used.
 *
 * freqmin and freqmax are used in order to limit the total range of
 * the spectral bands.
 *
 * Although some the description of this routine refer to banding
 * this routine does have anything to do with the stratified sampling
 * used to create photons.  
 *
 **********************************************************/
int
ion_bands_init (int mode, double freqmin, double freqmax, struct xbands *band)
{
  int i, n, ngood, good[NXBANDS];
  double xfreq[NXBANDS + 1];
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
  /* band model user_bands: the boundaries given by the user, within the range of the generation bands */
  else if (geo.band_model == BAND_MODEL_USER_BANDS)
  {
    nxfreq = ion_bands_from_user (band, xfreq);
  }
  /* the default: boundaries on the thresholds of the abundant elements, the rest by splitting */
  else if (geo.band_model != BAND_MODEL_PHOTON_BANDS)
  {
    nxfreq = ion_bands_from_elements (band, xfreq);
  }
  /* band model photon_bands: if there are enough photon generation bands, use their boundaries */
  else if (band->nbands >= MIN_N_IONBANDS)
  {
    for (nband = 0; nband < band->nbands; nband++)
    {
      xfreq[nband] = band->f1[nband];
    }
    nxfreq = band->nbands;
    xfreq[band->nbands] = band->f2[band->nbands - 1];
  }

  /* if we don't have enough generation bands, place the ionization bands on the important edges */
  else
  {
    nxfreq = ion_bands_from_edges (band, xfreq);

    Log ("ion_bands_init: only %d generation bands, so adding boundaries at photoionization edges\n", band->nbands);
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


  Log ("ion_bands_init: There were %d final intervals\n", geo.nxfreq);
  for (n = 0; n < geo.nxfreq; n++)
  {
    Log ("ion_bands_init: %8.2f (%8.2e)    %8.2f (%8.2e)  \n", geo.xfreq[n] * HEV, geo.xfreq[n], geo.xfreq[n + 1] * HEV, geo.xfreq[n + 1]);
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
check_appropriate_banding (struct xbands *band, int mode)
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
