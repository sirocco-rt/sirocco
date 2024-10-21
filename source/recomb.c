
/***********************************************************/
/** @file  recomb.c
 * @author ksl,nsh
 * @date   January, 2018
 *
 * @brief  The routines in this file all have to do with
 * photoionization and/or
 * recombination rates and emissivities.
 *
 * ###Notes###
 *
 * These routines are quite complex representing a lot of work by
 * various of us over a long period of time, and it is not entirely
 * clear that they are quite what we want at present.  It would be
 * worthwhile to think about
 * their structure the next time we make consider adding a physical
 * process involving recombination.
 *
 * Some of the complication associated with these routines arises
 * from run time considerations which may no longer be very relevant
 * today with parallel processing etc.
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "sirocco.h"

/***************************FBSTRUC ***********************************/
/* The next section contains the freebound structures that can be used for both the
 * specific emissivity of a free-bound transition, and for the recombination coefficient
 * assuming the array has been initialized, which can take a few minutes
*/

#define NTEMPS  60              // The number of temperatures which are stored in each fbstruct
                                /* NSH this was increased from 30 to 60 to take account of 3 extra OOM
                                   intemperature we wanted to have in fb */
#define NFB 20                  // The maximum number of frequency intervals for which the fb emission is calculated


struct fbstruc
{
  double f1, f2;
  double cool[NIONS][NTEMPS];   //cooling rate due to radiative recombination
  double lum[NIONS][NTEMPS];    //emissivity due to radiative recombinaion
  double cool_inner[NIONS][NTEMPS];     //cooling rate due to recombinations to inner shells
}
freebound[NFB];



double xnrecomb[NIONS][NTEMPS]; // There is only one set of recombination coefficients
double xninnerrecomb[NIONS][NTEMPS];    // There is only one set of recombination coefficient

double fb_t[NTEMPS];
int nfb = 0;


/** FBEMISS was calculated as follows:
 * x= 2. * PI * MELEC * BOLTZMANN / (H*H);
 * x=pow(x,-1.5);
 * x*=8. * PI / (C*C);
 * x*= H;
 */
#define FBEMISS   7.67413e-62   // Calculated with constants.c
#define LOG_FBEMISS -140.7224208336316


/* These are external structures used primarily because we need to call
Numerical Recipes routines from fb_verner and fb_topbase */

///Topbase description of a photoionization x-section
struct topbase_phot *fb_xtop;

/// Temperature (and log) at which the emissivity is calculated 
double fbt, log_fbt;

/// fb_choice (see above)
int fbfr;






/**********************************************************/
/**
 * @brief      returns the partial (for a specific ion) emissivity or
 * recombination rate for ions described in terms of Topbase photoionization x-sections.
 *
 * @param [in] double  freq   The freqeuncy of interest
 * @return     An emissivity or a recombination rate
 *
 * What the routine returns depends on the external variable fbfr. The
 * choices are:
 *
 * * FB_FULL         Calculate fb emissivity including energy associated with the threshold
 * * FB_REDUCED      Calculate the fb emissivity without the threshold energy
 * * FB_RATE         Calulate the fb recombinarion rate
 *
 * @details
 *
 *
 * ### Notes ###
 *
 * This routine used to be used for integrations, the wrapper routine fb_topbase_partial2 is
 * now used for that purpose - this is only used directly now.
 * Recast in log space for speedup purposes
 *
 *
 *  231128 - working with log_sigma_phot removed. ksl
 *
 **********************************************************/

double
fb_topbase_partial (freq)
     double freq;
{
  int nion;
  double partial, log_freq;
  double logx;
  double log_gn, log_gion;
  double fthresh;
  double test;
  double x;

  log_freq = log (freq);        //Go into log space

  fthresh = fb_xtop->freq[0];
  if (freq < fthresh)
    return (0.0);               // No recombination at frequencies lower than the threshold freq occur

  nion = fb_xtop->nion;

  /* JM -- below lines to address bug #195 */
//  gn = 1;
  log_gn = 0;
  if (ion[nion].phot_info > 0)  // it's a topbase record
  {
    log_gn = xconfig[fb_xtop->nlev].log_g;
  }

  else if (ion[nion].phot_info == 0)    // it's a VFKY record, so shouldn't really use levels
  {
    log_gn = ion[nion].log_g;
  }
  else
  {
    Error
      ("fb_topbase_partial: Did not understand cross-section type %i for ion %i (z=%i, istate %i). Setting multiplicity to zero!\n",
       ion[nion].phot_info, nion, ion[nion].z, ion[nion].istate);
    log_gn = -999.;             //Not really sure what to do here - probably return a zero
  }

  log_gion = ion[nion + 1].log_g;       // Want the g factor of the next ion up in log space

  x = sigma_phot (fb_xtop, freq);
  if (x < VERY_SMALL)
  {
    return (0.0);
  }

  logx = log (x);

  if (sane_check (logx))
  {
    Error ("fb_topbase_partial: x %e logx %e\n", x, logx);
    return (0.0);
  }



  // Now calculate emission using Ferland's expression - recast in log space 2022 for speed.

//  partial = FBEMISS * gn / (2. * gion) * pow (freq * freq / fbt, 1.5) * exp (H_OVER_K * (fthresh - freq) / fbt) * x;

  test =
    LOG_FBEMISS + log_gn - 0.6931471805599453 - log_gion + 1.5 * (2.0 * log_freq - log_fbt) + (H_OVER_K * (fthresh - freq) / fbt) + logx;


  partial = exp (test);         //Go back to linear space


  // 0=emissivity, 1=heat loss from electrons, 2=photons emissivity

  if (fbfr == FB_REDUCED)
    partial *= (freq - fthresh) / freq;
  else if (fbfr == FB_RATE)
    partial /= (PLANCK * freq);


  if (sane_check (partial))
  {
    Error ("fb_topbase_partial: Failed test %e partial %e\n", test, partial);
    Error ("log_gn %e log_gion %e log_freq %e log_fbt %e fthresh %e  freq %e  fbt %e logx %e\n",
           log_gn, log_gion, log_freq, log_fbt, fthresh, freq, fbt, logx);
    partial = 0.0;

  }

  return (partial);
}





/**********************************************************/
/**
 * @brief      This is a wrapper for fb_topbase_partial to allow it to be used for integrations
 *
 * @param [in] double  freq   The freqeuncy of interest
 * @param [in] void  params   An extra (unused) variable to make it paletable for the gsl integrator

 * @return     An emissivity or a recombination rate
 *
 * What the routine returns depends on the external variable fbfr. The
 * choices are:
 *
 * * FB_FULL         Calculate fb emissivity including energy associated with the threshold
 * * FB_REDUCED      Calculate the fb emissivity without the threshold energy
 * * FB_RATE         Calulate the fb recombinarion rate
 *
 * @details
 *
 *
 * ### Notes ###
 *
 * This routine is integrated over frequency using a gsl routine wrapped in the num_int
 * routine for that purpsoe.  Much of the information is passed externally for historical reasons
 * In Princible the extranl information could be contrained in the parameters.
 *
 *
 *
 **********************************************************/


double
fb_topbase_partial2 (double freq, void *params)
{
  double partial;

  partial = fb_topbase_partial (freq);

  return (partial);
}




/**********************************************************/
/**
 * @brief      calculates the integrated emissivity of the plasma, or the number of
 * recombinations per second of a particular ion within defined frequency bands.
 *
 * @param [in] double  t   The temperature at which the emissivity
 * @param [in] double  f1   The minimum frequency
 * @param [in] double  f2   The maximum frequency
 * @param [in, out] int  nion   The ion for which the emissivity is returned
 * @param [in, out] int  fb_choice   A switch which determines exactly what is to be returned
 * @param [in, out] int  mode   A switch which indicates whether one is interested in normal
 * radiative recombination (OUTER_SHELL) or inner shell recombinaiton (INNER_SHELL)
 * @return     The routine returns the specific emissivity, e.g. the emissivity and
 * 	or recombination rate per electron and per ion.
 *
 * 	The options are as follows:
 *
 * * FB_FULL         Calculate fb emissivity including energy associated with the threshold
 * * FB_REDUCED      Calculate the fb emissivity without the threshold energy
 * * FB_RATE         Calulate the fb recombinarion rate
 *
 *
 * @details
 *
 * The routine returns the integral of the emissivity (or recombination rate) over the
 * frequency interval.  The routine uses pre-calculated values if that is possible, or
 * calculate a new values if that is not possible. (In that case xinteg_fb updates the
 * set of pre-calculated values so that they can be used in a future call.
 *
 * ### Notes ###
 *
 *
 **********************************************************/

double
integ_fb (t, f1, f2, nion, fb_choice, mode)
     double t;                  // The temperature at which to calculate the emissivity
     double f1, f2;             // The frequencies over which to integrate the emissivity
     int nion;                  // The ion for which the "specific emissivity" is calculateed
     int fb_choice;             // 0=full, 1=reduced, 2= rate
     int mode;                  // 1- outer shell 2-inner shell
{
  double fnu;
  int n;

  if (mode == OUTER_SHELL)
  {

    if (fb_choice == FB_FULL)
    {
      for (n = 0; n < nfb; n++)
      {
        /* See if the frequencies correspond to one previously calculated */
        if (f1 == freebound[n].f1 && f2 == freebound[n].f2)
        {
          fnu = get_fb (t, nion, n, fb_choice, mode);
          return (fnu);
        }
      }
      /* If not calculate it here */
      fnu = xinteg_fb (t, f1, f2, nion, fb_choice);
      return (fnu);
    }
    else if (fb_choice == FB_REDUCED)
    {
      for (n = 0; n < nfb; n++)
      {
        /* See if the frequencies correspond to one previously calculated */
        if (f1 == freebound[n].f1 && f2 == freebound[n].f2)
        {
          fnu = get_fb (t, nion, n, fb_choice, mode);
          return (fnu);
        }
      }
      /* If not calculate it here */
      fnu = xinteg_fb (t, f1, f2, nion, fb_choice);
      return (fnu);
    }
    else if (fb_choice == FB_RATE)
    {
      /* See if the frequencies correspond to one previously calculated */
      if (nfb > 0)
      {
        fnu = get_nrecomb (t, nion, mode);
        return (fnu);
      }
      /* If not calculate it here */
      fnu = xinteg_fb (t, f1, f2, nion, fb_choice);
      return (fnu);
    }
    Error ("integ_fb: Unknown fb_choice(%d)\n", fb_choice);
    Exit (0);
  }

  else if (mode == INNER_SHELL) // inner shell
  {
    if (fb_choice == FB_FULL)
    {
      for (n = 0; n < nfb; n++)
      {
        /* See if the frequencies correspond to one previously calculated */
        if (f1 == freebound[n].f1 && f2 == freebound[n].f2)
        {
          fnu = get_fb (t, nion, n, fb_choice, mode);
          return (fnu);
        }
      }
      fnu = xinteg_inner_fb (t, f1, f2, nion, fb_choice);
      return (fnu);
    }
    else if (fb_choice == FB_REDUCED)
    {
      for (n = 0; n < nfb; n++)
      {
        /* See if the frequencies correspond to one previously calculated */
        if (f1 == freebound[n].f1 && f2 == freebound[n].f2)
        {
          fnu = get_fb (t, nion, n, fb_choice, mode);
          return (fnu);
        }
      }
      fnu = xinteg_inner_fb (t, f1, f2, nion, fb_choice);
      return (fnu);
    }
    else if (fb_choice == FB_RATE)
    {
      if (nfb > 0)
      {
        fnu = get_nrecomb (t, nion, mode);
        return (fnu);
      }
      fnu = xinteg_inner_fb (t, f1, f2, nion, fb_choice);
      return (fnu);
    }
    Error ("integ_fb: Unknown fb_choice(%d)\n", fb_choice);
    Exit (0);
  }

  Error ("integ_fb: Unknown mode(%d)\n", mode);
  Exit (0);
  return (0);
}







/**********************************************************/
/**
 * @brief      returns the energy lost from the plasma due to fb emission in a
 * 	single wind cell at a temperature t between the frequncy limits f1 and f2.
 * 	The energy lost is just the kinetic energy lost from the plasma
 * 	because it does not include the ionization potential
 * 	associated with each recombination.  Python tracks effectively the kinetic energy
 * 	of the plasma (not the potential energy available if everything recombined).
 *
 * @param [in] PlasmaPtr  xplasma   The plasma cell of interest
 * @param [in] double  t   The temperature of the cell
 * @param [in] double  f1   The minimum frequency
 * @param [in] double  f2   The maximum frequency
 * @param [in] int  fb_choice   A switch controlling whether emissivites or rates are
 * returned
 * @param [in] int  mode   A switch denoting whether normal recombination (OUTER_SHELL) or inner shell recombination (INNER_SHELL)
 * @return     The routine calculates an emissivity for normal recombination or a cooling rate for normal or dielectronic recombination
 *
 * If OUTER_SHELL is chosen then the options are
 *
 * * FB_FULL         Calculate fb emissivity including energy associated with the threshold
 * * FB_REDUCED      Calculate the fb emissivity without the threshold energy
 *
 * If INNER_SHELL is chosen the recombination rate is returned.
 *
 * @details
 *
 *
 * ### Notes ###
 * The results in this routine are all based on the Milne relation
 *
 * The routine is unusual in that one passes a the wind cell and
 * not the plasma cell.  This probably should be changed.
 *
 * Note that FB_REDUCED is not an option for this routine
 *
 * Question: What is preventing us from calculating a dielectronic emission rate?
 *
 * Total free bound emission in this function is calculated in the CMF. 
 * 
 **********************************************************/

double
total_fb (xplasma, t, f1, f2, fb_choice, mode)
     PlasmaPtr xplasma;
     double t, f1, f2;
     int fb_choice;
     int mode;
{
  double total;
  int nion;

  if (t < 100. || f2 < f1)
    t = 100.;                   /* Set the temperature to 100 K so that if there are free electrons emission by this process continues */

  if (mode == OUTER_SHELL)
    init_freebound (100., 1.e9, f1, f2);


// Calculate the number of recombinations whenever calculating the fb_luminosities
  num_recomb (xplasma, t, mode);

  total = 0;
  xplasma->cool_rr_metals = 0.0;
  xplasma->lum_rr_metals = 0.0;

  /*
   * This loop is now over nions - 1, to avoid out of bounds access and above
   * the final ion there is nothing for it to recombine from.
   */

  for (nion = 0; nion < nions - 1; nion++)
  {
    if (xplasma->density[nion] > DENSITY_PHOT_MIN)
    {
      if (mode == OUTER_SHELL)
      {
        if (fb_choice == FB_FULL)       // we are calculating a luminosity
        {
          total += xplasma->lum_rr_ion[nion] =
            xplasma->vol * xplasma->ne * xplasma->density[nion + 1] * integ_fb (t, f1, f2, nion, fb_choice, mode);
          if (ion[nion].z > 2)
            xplasma->lum_rr_metals += xplasma->lum_rr_ion[nion];
        }
        else                    // we are calculating a cooling rate
        {
          total += xplasma->cool_rr_ion[nion] =
            xplasma->vol * xplasma->ne * xplasma->density[nion + 1] * integ_fb (t, f1, f2, nion, fb_choice, mode);
          if (ion[nion].z > 2)
            xplasma->cool_rr_metals += xplasma->cool_rr_ion[nion];
        }
      }
      else if (mode == INNER_SHELL)     // at present we do not compute a luminosity from DR
        total += xplasma->cool_dr_ion[nion] =
          xplasma->vol * xplasma->ne * xplasma->density[nion + 1] * integ_fb (t, f1, f2, nion, fb_choice, mode);
    }

  }
  return (total);
}


double fb_x[NCDF], fb_y[NCDF];
/// There is at most one jump per level
double fb_jumps[NLEVELS];
/// This is just a dummy array that parallels fb_jumpts
double xfb_jumps[NLEVELS];
int fb_njumps = (-1);

// WindPtr ww_fb;
double one_fb_f1, one_fb_f2, one_fb_te; /* Old values */


/**********************************************************/
/**
 * @brief      generates one free bound photon with specific frequency limits
 *
 * @param [in] PlasmaPtr  xplasma   The wind cell in which the photon is being
 * @param [in] double  f1   The minimum frequency
 * @param [in] double  f2   The frequency limits
 * @return     The frequency of the fb photon that was generated.
 *
 * @details
 *
 * 	Although the routine returns only a single photon, it stores more
 * 	that one for the same conditions. If the routine is called mulitiple
 * 	times from the same cell with the same frequency limits it will use the pre-
 * 	stored values.  This was intended to avoid the process of having to
 * 	generate cdfs multiple times
 *
 * ### Notes ###
 *
 *
 * 	@bug This routine still assumes the possibility of jumps
 *      even though this possibility has been removed from the cdf generation
 *      routines.
 *
 **********************************************************/

double
one_fb (xplasma, f1, f2)
     PlasmaPtr xplasma;         /* a single cell */
     double f1, f2;             /* freqmin and freqmax */
{
  double freq, tt, delta;
  int n, nn, nnn;
  double fthresh, dfreq;
  int nplasma;
  PhotStorePtr xphot;
  nplasma = xplasma->nplasma;
  xphot = &photstoremain[nplasma];

  if (f2 < f1)
  {
    Error ("one_fb: f2 %g < f1 %g Something is rotten  t %g\n", f2, f1, xplasma->t_e);
    Exit (0);
  }

/* Check if an appropriate photon frequency has already been generated, 
   and use that instead if possible 
 */

  tt = xplasma->t_e;
  if (xphot->n < NSTORE && xphot->f1 == f1 && xphot->f2 == f2 && xphot->t == tt)
  {
    freq = xphot->freq[xphot->n];
    (xphot->n)++;
    return (freq);
  }


  delta = tt / 100;             // Fudge factor to prevent generation of a CDF if t has changed only slightly

  /* Check to see if we have already generated a cdf */
  if (tt > (one_fb_te + delta) || tt < (one_fb_te - delta) || f1 != one_fb_f1 || f2 != one_fb_f2)
  {
/* Then need to generate a new cdf */

//OLD    ww_fb = one;

    /* Create the fb_array */

    /* Determine how many intervals are between f1 and f2.  These need to be
       put in increasing frequency order */

    if (f1 != one_fb_f1 || f2 != one_fb_f2)
    {                           // Regenerate the jumps
      fb_njumps = 0;
      for (n = 0; n < nphot_total; n++)
      {
        fthresh = phot_top_ptr[n]->freq[0];
        if (f1 < fthresh && fthresh < f2)
        {
          fb_jumps[fb_njumps] = fthresh;
          fb_njumps++;
        }
      }



      /* The next line sorts the fb_jumps by frequency and eliminates
       * duplicate frequencies which is what was causing the error in
       * cdf.c when more than one jump was intended
       */

      if (fb_njumps > 1)        //We only need to sort and compress if we have more than one jump
      {
        fb_njumps = sort_and_compress (fb_jumps, xfb_jumps, fb_njumps);
        for (n = 0; n < fb_njumps; n++)
        {
          fb_jumps[n] = xfb_jumps[n];
        }
      }


    }



    /*NSH 1707 - modified the loop below to ensure we have points just below and above any jumps */

    nnn = 0;                    //Zero the index for elements in the flux array
    nn = 0;                     //Zero the index for elements in the jump array
    n = 0;                      //Zero the counting element for equally spaced frequencies
    dfreq = (f2 - f1) / (ARRAY_PDF - 1);        //This is the frequency spacing for the equally spaced elements
    while (n < (ARRAY_PDF) && nnn < NCDF)       //We keep going until n=ARRAY_PDF-1, which will give the maximum required frequency
    {
      freq = f1 + dfreq * n;    //The frequency of the array element we would make in the normal run of things
      if (freq > fb_jumps[nn] && nn < fb_njumps)        //The element we were going to make has a frequency abouve the jump
      {
        fb_x[nnn] = fb_jumps[nn] * (1. - DELTA_V / (2. * VLIGHT));      //We make one frequency point DELTA_V cm/s below the jump
        fb_y[nnn] = fb (xplasma, xplasma->t_e, fb_x[nnn], nions, FB_FULL);      //And the flux for that point
        nnn = nnn + 1;          //increase the index of the created array
        fb_x[nnn] = fb_jumps[nn] * (1. + DELTA_V / (2 * VLIGHT));       //And one frequency point just above the jump
        fb_y[nnn] = fb (xplasma, xplasma->t_e, fb_x[nnn], nions, FB_FULL);      //And the flux for that point
        nn = nn + 1;            //We heave dealt with this jump - on to the next one
        nnn = nnn + 1;          //And we will be filling the next array element next time
      }
      else                      //We haven't hit a jump
      {
        if (freq > fb_x[nnn - 1])       //Deal with the unusual case where the upper point in our 'jump' pair is above the next regular point
        {
          fb_x[nnn] = freq;     //Set the next array element frequency
          fb_y[nnn] = fb (xplasma, xplasma->t_e, fb_x[nnn], nions, FB_FULL);    //And the flux
          n = n + 1;            //Increment the regular grid counter
          nnn = nnn + 1;        //Increment the generated array counter
        }
        else                    //We dont need to make a new point, the upper frequency pair of the last jump did the trick
        {
          n = n + 1;            //We only need to increment our regualr grid counter
        }
      }
    }

    //Ensure the last point lines up exatly with f2

    fb_x[nnn - 1] = f2;
    fb_y[nnn - 1] = fb (xplasma, xplasma->t_e, f2, nions, FB_FULL);


    if (nnn > NCDF)
    {
      Error ("one _fb: Overflow of working array\n");
      Exit (0);
    }


    /* At this point, the variable nnn stores the number of points */


    if (cdf_gen_from_array (&cdf_fb, fb_x, fb_y, nnn, f1, f2) != 0)
    {
      Error ("one_fb after cdf_gen_from_array error: f1 %g f2 %g te %g ne %g nh %g vol %g\n",
             f1, f2, xplasma->t_e, xplasma->ne, xplasma->density[1], xplasma->vol);
      Error ("Giving up\n");
      Exit (0);
    }
    one_fb_te = xplasma->t_e;
    one_fb_f1 = f1;
    one_fb_f2 = f2;             /* Note that this may not be the best way to check for a previous cdf */
  }

/* OK, generate photons */

/* First generate the photon we need */
  freq = cdf_get_rand (&cdf_fb);
  if (freq < f1 || freq > f2)
  {
    Error ("one_fb:  freq %e  freqmin %e freqmax %e out of range\n", freq, f1, f2);
  }

/* Now create and store for future use a set of additonal photons */

  for (n = 0; n < NSTORE; n++)
  {
    xphot->freq[n] = cdf_get_rand (&cdf_fb);
    if (xphot->freq[n] < f1 || xphot->freq[n] > f2)
    {
      Error ("one_fb:  freq %e  freqmin %e freqmax %e out of range\n", xphot->freq[n], f1, f2);
    }

  }
  xphot->n = 0;
  xphot->t = tt;
  xphot->f1 = f1;
  xphot->f2 = f2;

  return (freq);
}



/**********************************************************/
/**
 * @brief      calculates the total number of recombinations
 * for all of the ions in a cell
 *
 * @param [in,out] PlasmaPtr  xplasma   The plasma cell of interest
 * @param [in] double  t_e   The temperarture of interest
 * @param [in] int  mode   A switch indicating whether one wants normal radiative recombination (OUTER_SHELL) or
 * dielectronic recombination (INNER_SHELL) rates to be calculated.
 * @return   Always returns 0; the results are stored in xplasma->recomb or xplasma->inner_recomb, depending
 * on the mode
 *
 * @details
 * The routine calculates recombination rates for all the ions in a single cell.  It stores the
 * recombination rates in the arrays recomb (for normal recombination) or inner_recomb (for dielectronic
 * recombination)
 *
 * ### Notes ###
 * The calculation is made purely direct recombination using the photoionization
 * x-sections.
 *
 * The units are #/cm**3/s
 *
 *
 **********************************************************/

int
num_recomb (xplasma, t_e, mode)
     PlasmaPtr xplasma;
     double t_e;
     int mode;
{
  int nelem;
  int i, imin, imax;
  for (nelem = 0; nelem < nelements; nelem++)
  {
    imin = ele[nelem].firstion;
    imax = ele[nelem].lastion;
    for (i = imin; i < imax; i++)
    {
      if (xplasma->density[i] > DENSITY_PHOT_MIN)
      {
        if (mode == OUTER_SHELL)        //outer shell
          xplasma->recomb[i] = xplasma->ne * xplasma->density[i + 1] * integ_fb (t_e, 0.0, VERY_BIG, i, FB_RATE, mode);
        else if (mode == INNER_SHELL)   //innershell
          xplasma->inner_recomb[i] = xplasma->ne * xplasma->density[i + 1] * integ_fb (t_e, 0.0, VERY_BIG, i, FB_RATE, mode);

      }
    }
    xplasma->recomb[imax] = 0.0;        // Can't recombine to highest i-state
    xplasma->inner_recomb[imax] = 0.0;  // Can't recombine to highest i-state

  }

  return (0);
}




/**********************************************************/
/**
 * @brief      calculates either the free_bound emissivity at a specific frequency or 
 * depending on inputs, the recombination rate 
 *
 * @param [in] PlasmaPtr  xplasma   A plasma cell
 * @param [in] double  t   The temperature at which to calculate the emisivity
 * @param [in] double  freq   The frequency at which to calculate the emissivity
 * @param [in] int  ion_choice   Either the total or the emissivity for a specific ion
 * @param [in] int  fb_choice   determines whether what is returned is the emissivity a specific frecuency 0
 * @return     The returns depend on fb_choice
 *
 * If ion_choice is a number less than the number of ions then the value returned for that specific
 * ion.  However if ion_choice is set to the number of ions or greater, the the value returned is 
 * for the sum of all the ions.  
 *
 *
 * The choices are:
 * * FB_FULL         Calculate fb emissivity including energy associated with the threshold
 * * FB_REDUCED      Calculate the fb emissivity without the threshold energy
 * * FB_RATE         Calulate the fb recombination rate
 *
 * @details
 *
 * ### Notes ###
 *
 * This routine calls subroutines that have variablew which are are transmitted
 * externally, e.g fbt,fbfr, rather than though calls.  The reason for this is some
 * of the routines are also used as integrands where only one variable can be transimitted
 * directly.  
 *
 * This routine does not have choices associated with inner shell recombination.
 *
 **********************************************************/

double
fb (xplasma, t, freq, ion_choice, fb_choice)
     PlasmaPtr xplasma;         // A cell with all its associated density data
     double t;                  // The temperature at which to calculate the emissivity
     double freq;               // The frequency at which to calculate the emissivity
     int ion_choice;            // Selects which ions the emissivity is to be calculated for (see above)
     int fb_choice;             // 0=emissivity in the standard sense, 1 heat loss from electons, 2 number of photons
{
  int n;
  double fnu, x;
  int nmin, nmax;               // These are the photo-ionization xsections that are used
  int nion, nion_min, nion_max;

  nion_min = nion_max = 0;
  if (ion_choice < nions)       //Get emissivity for this specific ion_number
  {
    nion_min = ion_choice;
    nion_max = ion_choice + 1;
  }
  else if (ion_choice == nions) // Get the total emissivity
  {
    nion_min = 0;
    nion_max = nions;
  }
  else
  {
    Error ("fb: This choice %d for ion_choice is not supported\n", ion_choice);
    Exit (0);
  }


  fbt = t;                      /* Externally transmitted variable */
  log_fbt = log (t);
  fbfr = fb_choice;             /* Externally transmitted variable */

  fnu = 0.0;                    /* Initially set the emissivity to zero */


  for (nion = nion_min; nion < nion_max; nion++)
  {
    if (ion[nion].phot_info > 0)        // topbase or VFKY+topbase
    {
      nmin = ion[nion].ntop_first;
      nmax = nmin + ion[nion].ntop;
    }
    else if (ion[nion].phot_info == 0)  // VFKY
    {
      nmin = ion[nion].nxphot;
      nmax = nmin + 1;
    }
    else
      nmin = nmax = 0;          // no XS / ionized - don't do anything

    //Debug("in fb for ion %i info %i, nmin nmax %i, %i\n", nion, ion[nion].phot_info, nmin, nmax);

    x = 0.0;

    /* Loop over relevent Topbase photoionization x-sections.  If
       an ion does not have Topbase photoionization x-sections then
       ntmin and ntmax are the same and the loop will be skipped. */

    for (n = nmin; n < nmax; n++)
    {
      fb_xtop = &phot_top[n];   /*Externally transmited to fb_topbase_partial */
      /* We don't want to include fb transitions associated with macro atoms here
         - they are separated out for now. (SS, Apr 04). "If" statement added. */
      if (fb_xtop->macro_info == FALSE || geo.macro_simple == TRUE || geo.rt_mode == RT_MODE_2LEVEL)
      {
        x += fb_topbase_partial (freq);
      }



    }


    /* x is the emissivity from this ion. Add it to the total */

    fnu += xplasma->density[nion + 1] * x;      // nion+1, the ion doing the recombining
  }

  fnu *= xplasma->ne;           // Correct from specific emissivity to the total fb emissivity

  return (fnu);

}



/** Indicates the total number of freebound sets that could be used */
int init_freebound_nfb;

/**********************************************************/
/**
 * @brief      initializes the structure fb_struc as well as some
 * associated arrays and variables (found in sirocco.h) that describe
 * recombination rates and band-limited luminosities.
 *
 * @param [in] double  t1   A lower limit for the temperature
 * @param [in] double  t2   An upper limit for the temperature
 * @param [in] double  f1   The lower limit for a frequency interval
 * @param [in] double  f2   The upper limit for the frequency interval
 * @return     The routine generally returns 0
 *
 * @details
 * Python typically calculates photons in frequency ranges (in order
 * to enable stratified sampling).  For this to work, one needs
 * freebound emissivities and cooling rates corresponding to these
 * freqency ranges.  Since we retrun to the same frequency ranges every cycle,
 * Python stores the necessary information in structures.
 *
 * This routine is responsible for populating these structures, so
 * that they can be accessed later via the routine get_fb.
 *
 *
 *
 * ### Notes ###
 *
 * The first time the routine is called, both recombination
 * rates and band-limited luminosities are calculated.  On
 * subsequent calls the routine checks to see whether it has
 * already calculated the band-limited freebound emissivities,
 * and if so returns without redoing the calculation.  However,
 * if a new frequency interval is provided, the new luminosities
 * are added to the free-bound structure.  To force a
 * re-initialization nfb must be set to 0.
 *
 * The routine allows for the possibility that there are more
 * frequency intervals than place to store data and will recylce
 * the structure if this occurs (indicating this with several error
 * messages).  This allows the program to proceed, but if this
 * happens often then the variable NFB in sirocco.h should be increased.
 *
 **********************************************************/

int
init_freebound (t1, t2, f1, f2)
     double t1, t2, f1, f2;
{
  double t;
  int i, j, nion;
  double ltmin, ltmax, dlt;
  double xinteg_fb ();
  int nput;



  if (nfb == 0)
  {
    if (t2 < t1)
    {
      Error ("init_freebound: t2(%g)<t1(%g)\n", t2, t1);
      Exit (0);
    }

    ltmin = log10 (t1);
    ltmax = log10 (t2);
    dlt = (ltmax - ltmin) / (NTEMPS - 1);

    for (j = 0; j < NTEMPS; j++)
    {
      fb_t[j] = pow (10., ltmin + dlt * j);
    }

    Log ("init_freebound: Creating recombination coefficients\n");
    for (nion = 0; nion < nions; nion++)
    {

      for (j = 0; j < NTEMPS; j++)
      {
        t = fb_t[j];
        xnrecomb[nion][j] = xinteg_fb (t, 0.0, VERY_BIG, nion, FB_RATE);
        xninnerrecomb[nion][j] = xinteg_inner_fb (t, 0.0, VERY_BIG, nion, FB_RATE);
      }
    }
  }
  else if (fabs (fb_t[0] - t1) > 10. || fabs (fb_t[NTEMPS - 1] - t2) > 1000.)
  {
    Error ("init_freebound: Cannot initialize to new temps without resetting nfb");
    Exit (0);

  }

/* Now check to see whether the freebound information has already
been calculated for these conditions, and if so simply return.
*/
  i = 0;
  while ((freebound[i].f1 != f1 || freebound[i].f2 != f2) && i < nfb)
  {

    i++;
  }

  if (i < nfb)
  {
    return (0);
  }

/* We have to calculate a new set of freebound data */
  if (i == NFB)
  {
    /* We've filled all the available space in freebound so we start recycling elements, assuming that the latest
     * ones are still likelyt to be needed
     */
    nput = init_freebound_nfb % NFB;
    init_freebound_nfb++;

    Error ("init_freebound: Recycling freebound, storage for NFB (%d), need %d to avoid \n", NFB, init_freebound_nfb);

  }
  else
  {
    nput = init_freebound_nfb = nfb;
    nfb++;
  }


/* Having reached this point, a new set of fb emissivities
must be calculated.  Note that old information is not destroyed
unless nfb had been set to 0.  The new set is added to the old
on the assumption that the fb information will be reused.
*/


  Log ("init_freebound: Creating recombination emissivities between %e and %e in structure element %d\n", f1, f2, nput);


  freebound[nput].f1 = f1;
  freebound[nput].f2 = f2;

  for (nion = 0; nion < nions; nion++)
  {
    for (j = 0; j < NTEMPS; j++)
    {                           //j covers the temps
      t = fb_t[j];
      freebound[nput].lum[nion][j] = xinteg_fb (t, f1, f2, nion, FB_FULL);
      freebound[nput].cool[nion][j] = xinteg_fb (t, f1, f2, nion, FB_REDUCED);
      freebound[nput].cool_inner[nion][j] = xinteg_inner_fb (t, f1, f2, nion, FB_REDUCED);

    }
  }

  return (0);
}



/**********************************************************/
/**
 * @brief      Return the recombination coefficient
 *
 * @param [in] double  t   The temperature
 * @param [in] int  nion   The ion of interest
 * @param [in] int  mode   A switch to choose normal (OUTER_SHELL) or inner shell (INNER_SHELL) recombination
 * @return     The recombination coefficient
 *
 * @details
 * Uses data from Badnell or another source to get a recombination rate
 *
 * ### Notes ###
 *
 **********************************************************/

double
get_nrecomb (t, nion, mode)
     double t;
     int nion;
     int mode;
{
  int linterp ();
  double x = -99.;
  if (mode == OUTER_SHELL)
    linterp (t, fb_t, xnrecomb[nion], NTEMPS, &x, 0);   //Interpolate in linear space
  else if (mode == INNER_SHELL)
    linterp (t, fb_t, xninnerrecomb[nion], NTEMPS, &x, 0);      //Interpolate in linear space
  else
  {
    Error ("Get_nrecomb: Unkonwn mode/type %i of recombination coefficient", mode);
    Exit (0);
  }
  return (x);
}




/**********************************************************/
/**
 * @brief      Interpolate from a set of stored emissivities
 *
 * @param [in] double  t   The temperure of interest
 * @param [in] int  nion   The ion of interest
 * @param [in] int  narray   The number of the array calculated for a particular frequency range
 * @param [in] int  fb_choice   A switch used only in the case of normal recombination
 * @param [in] int  mode   A switch which indicates whether one is interested 
 * in normal (OUTER_SHELL) or innershell (INNER_SHELL)  recombination
 * @return     The program generally returns an emissivity or a cooling rate, 
 * and the choices are determined by the
 * combination of fb_choice and mode
 *
 * if the mode is set for normal recombiantion, then the possiblities are:
 * * FB_FULL         Calculate fb emissivity including energy associated with the threshold
 * * FB_REDUCED      Calculate the fb emissivity without the threshold energy
 *
 * if the mode is set for inner shell then the cooling rate (effectively FB_FULL) is returned.
 *
 *
 * @details
 * In an effort to save time, information needed to calculate free-bound 
 * emissivities and cooling are calculated when
 * a freqency interval is defined, and this routine is used to retrieve 
 * these values.  The fb structures are populated by init_freebound
 *
 *
 * ### Notes ###
 *
 **********************************************************/

double
get_fb (t, nion, narray, fb_choice, mode)
     double t;
     int nion;
     int narray;
     int fb_choice;
     int mode;
{
  int linterp ();
  double x = -99.;
  if (mode == OUTER_SHELL)
  {
    if (fb_choice == FB_REDUCED)
      linterp (t, fb_t, &freebound[narray].cool[nion][0], NTEMPS, &x, 0);       //Interpolate in linear space
    else if (fb_choice == FB_FULL)
      linterp (t, fb_t, &freebound[narray].lum[nion][0], NTEMPS, &x, 0);        //Interpolate in linear space
    else
    {
      Error ("Get_fb - unexpected mode %i", mode);
      Exit (0);
    }
  }
  else if (mode == INNER_SHELL)
    linterp (t, fb_t, &freebound[narray].cool_inner[nion][0], NTEMPS, &x, 0);   //Interpolate in linear space

  else
  {
    Error ("Get_fb - unkonwn mode %i", mode);
    Exit (0);
  }
  return (x);
}





/**********************************************************/
/**
 * @brief      calculates the integrated emissivity of
 *   an ion in the plasma.
 *
 * @param [in, out] double  t   The temperature of interest
 * @param [in, out] double  f1   The minimum frequency
 * @param [in, out] double  f2   The maximum frequency
 * @param [in, out] int  nion   The ion of interest
 * @param [in, out] int  fb_choice   A switch which determined exactly what is returned
 * @return     The integrated emissivity, depending on the fb_choice.  The possibilites are
 *
 * * FB_FULL         Calculate fb emissivity including energy associated with the threshold
 * * FB_REDUCED      Calculate the fb emissivity without the threshold energy
 * * FB_RATE         Calulate the fb recombinarion rate
 *
 * @details
 *
 * The routine performs a numerical integration over the partial emissivities
 *
 *
 * ### Notes ###
 * This routine is called by integ_fb.  It is not intended to be called
 * directly, which uses a Numerical Recipes routine to integrate over frequncy.
 *
 *
 **********************************************************/

double
xinteg_fb (t, f1, f2, nion, fb_choice)
     double t;                  // The temperature at which to calculate the emissivity
     double f1, f2;             // The frequencies overwhich to integrate the emissivity
     int nion;                  // The ion for which the "specific emissivity is calculateed
     int fb_choice;             // 0=full, otherwise reduced
{
  int n;
  double fnu;
  double dnu;                   //NSH 140120 - a parameter to allow one to restrict the integration limits.
  double fthresh, fmax;
  double den_config ();
  int nmin, nmax;               // These are the limits over which number xsections we will use
  double qromb ();


  dnu = 0.0;                    //Avoid compilation errors.

  nmin = nmax = 0;
  if (-1 < nion && nion < nions)        //Get emissivity for this specific ion_number
  {
    if (ion[nion].phot_info > 0)        // topbase or hybrid
    {
      nmin = ion[nion].ntop_first;
      nmax = nmin + ion[nion].ntop;
    }
    else if (ion[nion].phot_info == 0)  // VFKY
    {
      nmin = ion[nion].nxphot;
      nmax = nmin + 1;
    }
    else
      // the ion is a fully ionized ion  and doesn't have a cross-section, so return 0
      return (0.0);
  }
  else                          // Get the total emissivity
  {
    Error ("integ_fb: %d is unacceptable value of nion\n", nion);
    Exit (0);
  }

  // Put information where it can be used by the integrating function
  fbt = t;
  log_fbt = log (t);

  fbfr = fb_choice;

  /* Limit the frequency range to one that is reasonable before integrating */

  if (f1 < 3e12)
    f1 = 3e12;                  // 10000 Angstroms
  if (f2 > 3e18)                // Set a maximum value for the maximum frequncy
    f2 = 3e18;                  // This is 1 Angstroms
  if (f2 < f1)
    return (0.0);               /* Because there is nothing to integrate */

  fnu = 0.0;


  for (n = nmin; n < nmax; n++)
  {
    // loop over relevent Topbase or VFKY photoionzation x-sections
    fb_xtop = &phot_top[n];

    /* Adding an if statement here so that photoionization that's part of a macro atom is
       not included here (these will be dealt with elsewhere). (SS, Apr04) */
    if (fb_xtop->macro_info == FALSE || geo.macro_simple == TRUE || geo.rt_mode == RT_MODE_2LEVEL)      //Macro atom check. (SS)
    {
      fthresh = fb_xtop->freq[0];
      fmax = fb_xtop->freq[fb_xtop->np - 1];    // Argues that this should be part of structure
      if (f1 > fthresh)
        fthresh = f1;
      if (f2 < fmax)
        fmax = f2;
      // Now calculate the emissivity as long as fmax exceeds xthreshold and there are ions to recombine
      if (fmax > fthresh)
      {
        //NSH 140120 - this is a test to ensure that the exponential will not go to zero in the integrations
        dnu = 100.0 * (fbt / H_OVER_K);
        if (fthresh + dnu < fmax)
        {
          fmax = fthresh + dnu;
        }
//        fnu += qromb (fb_topbase_partial, fthresh, fmax, 1.e-4);
        fnu += num_int (fb_topbase_partial2, fthresh, fmax, 1.e-4);

      }
    }
  }


  return (fnu);
}





/**********************************************************/
/**
 * @brief      calculates the integrated fb emissivity of inner
 *   shell transitions in an ion at a given temperature
 *
 * @param [in] double  t   The temperarure of interest
 * @param [in] double  f1   The minimum frequency
 * @param [in] double  f2   The maximum frequency
 * @param [in] int  nion   The ion of interest
 * @param [in] int  fb_choice   ???
 * @return     Exactly what the the routine returns is determined by fb_choice,
 * which is passed as an an external variable to fb_topbase_partial
 *
 * The possibilites are
 *
 * * FB_FULL         Calculate fb emissivity including energy associated with the threshold
 * * FB_REDUCED      Calculate the fb emissivity without the threshold energy
 * * FB_RATE         Calulate the fb recombinarion rate
 *
 * @details
 *
 *
 *
 * ### Notes ###
 * This routine is  virtual copy of xinteg_fb but considers only inner shells .
 * The integration that is requied to integrate over frequencies is carried out using
 * the Numerical Recipes routine qromb.  This is why a number of variables are
 * passed as external variables, include fb_choice, fmin and fmax.
 *
 *
 *
 **********************************************************/

double
xinteg_inner_fb (t, f1, f2, nion, fb_choice)
     double t;                  // The temperature at which to calculate the emissivity
     double f1, f2;             // The frequencies overwhich to integrate the emissivity
     int nion;                  // The ion for which the "specific emissivity is calculateed
     int fb_choice;             // 0=full, otherwise reduced
{
  int n, nn;
  double fnu;
  double dnu;                   // a parameter to allow one to restrict the integration limits.
  double fthresh, fmax;
  double den_config ();


  dnu = 0.0;                    //Avoid compilation errors.
  fnu = 0.0;
  nn = -1;


  if (f1 < 3e12)
    f1 = 3e12;                  // 10000 Angstroms
  if (f2 > 3e18)                // increase upper limits to include  highly ionised ions that we are now seeing in x-ray illuminated nebulas.
    f2 = 3e18;                  // This is 1 Angstroms
  if (f2 < f1)
    return (0.0);               /* Because there is nothing to integrate */

  for (n = 0; n < n_inner_tot; n++)
  {
    if (inner_cross[n].nion == nion)
    {
      nn = n;
      fbt = t;
      log_fbt = log (t);

      fbfr = fb_choice;

      /* Limit the frequency range to one that is reasonable before integrating */



      // loop over relevent Topbase or VFKY photoionization x-sections
      fb_xtop = &inner_cross[nn];

      /* Adding an if statement here so that photoionization that's part of a macro atom is
         not included here (these will be dealt with elsewhere). (SS, Apr04) */
      if (fb_xtop->macro_info == FALSE || geo.macro_simple == TRUE || geo.rt_mode == RT_MODE_2LEVEL)    //Macro atom check. (SS)
      {
        fthresh = fb_xtop->freq[0];
        fmax = fb_xtop->freq[fb_xtop->np - 1];  // Argues that this should be part of structure
        if (f1 > fthresh)
          fthresh = f1;
        if (f2 < fmax)
          fmax = f2;

        // Now calculate the emissivity as long as fmax exceeds xthreshold and there are ions to recombine
        if (fmax > fthresh)
        {
          //NSH 140120 - this is a test to ensure that the exponential will not go to zero in the integrations
          dnu = 100.0 * (fbt / H_OVER_K);
          if (fthresh + dnu < fmax)
          {
            fmax = fthresh + dnu;
          }
          fnu += num_int (fb_topbase_partial2, fthresh, fmax, 1.e-4);
        }

      }
    }
  }



  return (fnu);
}





/**********************************************************/
/**
 * @brief      get the total rediative recombination rate
 *
 * @param [in] int  nion   The ion number of interest
 * @param [in] double  T   The temperature for which the rate is calculated
 * @return     The total recombination rate for this ion
 * 		for this temperature
 *
 * @details
 * Generates a total recombination rate for
 * a given ion at a given temperature using
 * data, obtained from sources such as Badnell
 * or Shull.
 *
 * If these are not present, an error is generated, and a
 * value using the Milne relation is returned.
 *
 * ### Notes ###
 * Recombination rates can be calculated from the Milne relation
 * but in most cases this does not give one the most accurate
 * recombination
 * rate because one does not have all of the relevant levels.
 *
 **********************************************************/

double
total_rrate (nion, T)
     int nion;
     double T;
{


  double rate;                  //The returned rate
  double rrA, rrB, rrT0, rrT1, rrC, rrT2;       //The parameters
  double term1, term2, term3;   //Some temporary parameters to make calculation simpler


  rate = 0.0;                   /* NSH 130605 to remove o3 compile error */


  if (ion[nion].total_rrflag == 1)      /*We have some kind of total radiative rate data */
  {
    if (total_rr[ion[nion].nxtotalrr].type == RRTYPE_BADNELL)
    {
      rrA = total_rr[ion[nion].nxtotalrr].params[0];
      rrB = total_rr[ion[nion].nxtotalrr].params[1];
      rrT0 = total_rr[ion[nion].nxtotalrr].params[2];
      rrT1 = total_rr[ion[nion].nxtotalrr].params[3];
      rrC = total_rr[ion[nion].nxtotalrr].params[4];
      rrT2 = total_rr[ion[nion].nxtotalrr].params[5];


      rrB = rrB + rrC * exp ((-1.0 * rrT2) / T);        //If C=0, this does nothing


      term1 = sqrt (T / rrT0);
      term2 = 1.0 + sqrt (T / rrT0);
      term2 = pow (term2, (1 - rrB));
      term3 = 1.0 + sqrt (T / rrT1);
      term3 = pow (term3, (1 + rrB));


      rate = pow ((term1 * term2 * term3), -1.0);
      rate *= rrA;
    }
    else if (total_rr[ion[nion].nxtotalrr].type == RRTYPE_SHULL)
    {
      rate = total_rr[ion[nion].nxtotalrr].params[0] * pow ((T / 1.0e4), -1.0 * total_rr[ion[nion].nxtotalrr].params[1]);
    }
    else
    {
      Error ("total_rrate: unknown parameter type for ion %i\n", nion);
      Exit (0);
    }
  }
  else                          /* We dont have coefficients, so use the Milne relation. Note
                                   that the Milne relation uses the lower ion of a pair, and so nion-1
                                   is correct
                                 */
  {
    Error ("total_rrate: No T_RR parameters for ion %i - using Milne relation\n", nion);
    rate = xinteg_fb (T, 3e12, 3e18, nion - 1, FB_RATE);
  }


  return (rate);

}



/**********************************************************/
/**
 * @brief      get the recombination rate to the ground state for an
 * ion at a particular temperature
 *
 * @param [in] int  nion   The ion of interest
 * @param [in] double  T   A temperature
 * @return     the rate
 *
 * @details
 * This routine generates a recombination rate to the ground state for
 * a given ion at a given temperature using data from Badnell or other
 * sources (with similar data formats).
 *
 * If these recombination rates a are not
 * available for the given ion the Milne relation is
 * used.
 *
 * ### Notes ###
 * nion is the ion that is recombining to a less ionizaed
 * state.  There is thus no rate for HI or He I, etc.
 *
 * The routine is used to produce recombination rate coefficients
 * for the matrix ionization
 *
 **********************************************************/

double
gs_rrate (nion, T)
     int nion;
     double T;
{
  double rate, drdt, dt;
  int i, imin, imax;
  double rates[BAD_GS_RR_PARAMS], temps[BAD_GS_RR_PARAMS];
  int ntmin;
  double fthresh, fmax, dnu;


  imin = imax = 0;              /* NSH 130605 to remove o3 compile error */



  if (ion[nion].bad_gs_rr_t_flag == 1 && ion[nion].bad_gs_rr_r_flag == 1)       //We have tabulated gs data

  {
    for (i = 0; i < BAD_GS_RR_PARAMS; i++)
    {
      rates[i] = bad_gs_rr[ion[nion].nxbadgsrr].rates[i];
      temps[i] = bad_gs_rr[ion[nion].nxbadgsrr].temps[i];
    }

    if (T < temps[0])           //we are below the range of GS data
    {
      Log_silent ("bad_gs_rr: Requested temp %e is below limit of data for ion %i(Tmin= %e)\n", T, nion, temps[0]);
      //      rate = rates[0];
      imax = 1;
      imin = 0;
    }

    else if (T >= temps[BAD_GS_RR_PARAMS - 1])  //we are above the range of GS data
    {
      Log_silent
        ("bad_gs_rr: Requested temp %e is above limit (%e) of data for ion %i\n",
         T, nion, bad_gs_rr[ion[nion].nxbadgsrr].temps[BAD_GS_RR_PARAMS - 1]);
      imax = BAD_GS_RR_PARAMS - 1;
      imin = BAD_GS_RR_PARAMS - 2;
      //We will try to extrapolate.

    }
    else                        //We must be within the range of tabulated data
    {
      for (i = 0; i < BAD_GS_RR_PARAMS - 1; i++)
      {
        if (temps[i] <= T && T < temps[i + 1])  //We have bracketed the correct temperature
        {
          imin = i;
          imax = i + 1;
        }
      }
    }
    /* interpolate in log space */
    drdt = (log10 (rates[imax]) - log10 (rates[imin])) / (log10 (temps[imax]) - log10 (temps[imin]));
    dt = (log10 (T) - log10 (temps[imin]));
    rate = pow (10, (log10 (rates[imin]) + drdt * dt));
  }

  /* Use the Milne relation -
     NB - this is different from using xinteg_fb because
     that routine does recombination to all excited levels (at least for topbase ions).
   */
  else
  {
    rate = 0.0;                 /* NSH 130605 to remove o3 compile error */

    fbt = T;
    log_fbt = log (T);
    fbfr = FB_RATE;

    if (ion[nion - 1].phot_info > 0)    //topbase or hybrid
    {
      ntmin = ion[nion - 1].ntop_ground;
      fb_xtop = &phot_top[ntmin];
    }
    else if (ion[nion - 1].phot_info == 0)      //vfky
    {
      fb_xtop = &phot_top[ion[nion - 1].nxphot];
    }

    fthresh = fb_xtop->freq[0];
    fmax = fb_xtop->freq[fb_xtop->np - 1];
    dnu = 100.0 * (fbt / H_OVER_K);

    if (fthresh + dnu < fmax)
    {
      fmax = fthresh + dnu;
    }

    rate = num_int (fb_topbase_partial2, fthresh, fmax, 1e-5);

  }

  return (rate);
}





/**********************************************************/
/**
 * @brief      sort an array into numerical order elimination duplicates
 *
 * @param [in] double *  array_in   The input array
 * @param [out] double *  array_out   The output array
 * @param [in] int  npts   The number of points in the input array
 * @return     The number of valid element in the output array
 *
 * @details
 * The routine uses the GSL routine qsort to sort the array in place,
 * and then copies unique elements of the sorted array into the output
 * array.
 *
 * ### Notes ###
 *
 * The routine is used in the creation of cdfs, which need arrays
 * which are sorted into numerical order, and for which one really
 * does not wish duplicated values.
 *
 **********************************************************/

int
sort_and_compress (double *array_in, double *array_out, int npts)
{
  double *values;
  int n, nfinal;

  values = calloc (sizeof (double), npts);
  for (n = 0; n < npts; n++)
  {
    values[n] = array_in[n];
  }

  /* Sort the array in place */
  qsort (values, npts, sizeof (double), compare_doubles);


  array_out[0] = values[0];     //Copy the first jump into the output array

  nfinal = 1;
  for (n = 1; n < npts; n++)    //Loop over the remaining jumps in the array
  {
    if (values[n] > array_out[nfinal - 1])      //In the next point in the array is larger than the last one (i.e. not equal)
    {
      array_out[nfinal] = values[n];    //Put the next point into the array
      nfinal += 1;              //Increment the size of the array
    }
  }

  free (values);

  return (nfinal);
}



/**********************************************************/
/**
 * @brief      A routine used by qsort in sort_and_compress
 *
 * @param [in] const void *  a   A double precision number
 * @param [in] const void *  b   A second double precision nubmer
 * @return     1 if a is greater than b, 0 otherwise
 *
 * @details
 * This routine just compares two double precision numbers and
 * returns 1 if a is greate than b, and 0 otherwise.  It is
 * needed by qsort which sorts a double precision array into
 * numberical order.
 *
 * ### Notes ###
 *
 **********************************************************/

int
compare_doubles (const void *a, const void *b)
{
  if (*(double *) a > *(double *) b)
    return 1;
  else if (*(double *) a < *(double *) b)
    return -1;
  else
    return 0;
}



/**********************************************************/
/** 
 * @brief selects the frequency of bf macro atom emission
 * 
 * 
 * @param [in]     WindPtr w   the ptr to the structure defining the wind
 * @param [in]     int nconf   the index into phot_top that identifies the continuum we wish to sample
 * @return freq    double freq the frequency of the packet to be emitted
 *
 * a scaled down version of one_fb for use with macro atom implementation. The objective is to select the 
 * emission frequency of a bound free photon that is to be generated by a specific continuum process
 * 
 *
 * ###Notes###
 * 
 * To improve the overall speed of Python, the routine generates
 * Mulitiple bf photons for a transition, and stores them in the 
 * matomxphot structure.  (It is not entirely clear that this
 * represents a signficant savings.)
 *
 * Prior to the creation of this routine, the macro atom routines always did 
 * this using an analytic hydrogenic approximation.
 * (SS/JM 1Aug2018)
***********************************************************/
double
matom_select_bf_freq (WindPtr one, int nconf)
{
  double f1, f2;
  double dfreq, freq;
  double te;
  PlasmaPtr xplasma;
  MatomPhotStorePtr matomxphot;

  int n;

  fbfr = FB_FULL;               //set external variable to sample the full emissivity of this process
  fb_xtop = &phot_top[nconf];   //set external pointer to the right bf process

  xplasma = &plasmamain[one->nplasma];
  te = xplasma->t_e;            //electron temperature in cell
  fbt = te;                     //set external temperature to the right value
  log_fbt = log (te);           //set external temperature to the right value

  //If hydrogenic ion use analytic expression
  if (ion[phot_top[nconf].nion].istate == ion[phot_top[nconf].nion].z)
  {
    return (phot_top[nconf].freq[0] - (log (1. - random_number (0.0, 1.0)) * te / H_OVER_K));
  }



//Check to see if we have some stored ones to use from previous pass
  matomxphot = &matomphotstoremain[one->nplasma];
  if (matomxphot->n < NSTORE && matomxphot->nconf == nconf && matomxphot->t == te)
  {
    freq = matomxphot->freq[matomxphot->n];
    (matomxphot->n)++;
    return (freq);
  }


  //make a new cdf and sample
  f1 = phot_top[nconf].freq[0]; //threshold frequency = minimum frequency for emission
  f2 = phot_top[nconf].freq[phot_top[nconf].np - 1];    //last frequency in list

  if ((H_OVER_K * (f2 - f1) / fbt) > ALPHA_MATOM_NUMAX_LIMIT)
  {
    //flast is currently very far into the exponential tail: so reduce flast to limit value of h nu / k T.
    f2 = f1 + fbt * ALPHA_MATOM_NUMAX_LIMIT / H_OVER_K;
  }



  dfreq = (f2 - f1) / (ARRAY_PDF - 1);  //This is the frequency spacing for the equally spaced elements

  for (n = 0; n < ARRAY_PDF; n++)       //We keep going until n=ARRAY_PDF-1, which will give the maximum required frequency
  {
    freq = f1 + dfreq * n;      //The frequency of the array element we would make in the normal run of things
    fb_x[n] = freq;             //Set the next array element frequency
    fb_y[n] = fb_topbase_partial (freq);        /* should return proportional to the total emissivity from this SINGLE bf. 
                                                   Note we don't need to multiply by n_e or n_ion since we only want a CDF 
                                                   for one process: so these factors will scale out
                                                 */
  }


//OLD  if (ARRAY_PDF > NCDF)
//OLD  {
//OLD    Error ("matom_select_bf_freq: Overflow of working array\n");
//OLD    Exit (0);
//OLD  }


  /* At this point, the variable nnn stores the number of points */


  if (cdf_gen_from_array (&cdf_fb, fb_x, fb_y, ARRAY_PDF, f1, f2) != 0)
  {
    Error ("matom_select_bf_freq after cdf_gen_from_array: f1 %g f2 %g te %g \n", f1, f2, xplasma->t_e);
    Error ("matomc_selct_fb_freg: Printing inputs to macro_recomb.txt\n");
    FILE *fptr;
    fptr = fopen ("macro_recomb.txt", "w");
    fprintf (fptr, "# fmin %e fmax %e\n", f1, f2);
    for (n = 0; n < ARRAY_PDF; n++)
    {
      fprintf (fptr, "%12.6e  %12.6e \n", fb_x[n], fb_y[n]);
    }
    fclose (fptr);

    Error ("Giving up\n");
    Exit (0);
  }


/* OK, generate photons */

/* First generate the photon we need */
  freq = cdf_get_rand (&cdf_fb);
  if (freq < f1 || freq > f2)
  {
    Error ("matom_select_bf_freq:  freq %e  freqmin %e freqmax %e out of range\n", freq, f1, f2);
  }


/* Now create and store for future use a set of additonal photons */

  for (n = 0; n < NSTORE; n++)
  {
    matomxphot->freq[n] = cdf_get_rand (&cdf_fb);
    if (matomxphot->freq[n] < f1 || matomxphot->freq[n] > f2)
    {
      Error ("matom_select_bf_freq:  freq %e  freqmin %e freqmax %e out of range\n", matomxphot->freq[n], f1, f2);
    }

  }
  matomxphot->n = 0;
  matomxphot->t = te;
  matomxphot->nconf = nconf;

  return (freq);


}
