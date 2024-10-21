

/***********************************************************/
/** @file  spectra.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief
 *
 * The subroutines in this file handle allocating, incrementing, and writing out the
 * spectrum arrays for Python
 *
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

int spec_initialized = FALSE;

/**********************************************************/
/**
 * @brief  Allocate memory for the arrays in xxspec.
 * @param  nspec  The number of spectra to allocate.
 *
 * @details
 *
 * NWAVE_MAX bins is always allocated. 
 * The number bins actually used may differ from this, and
 * moreover the number of bins during the ionization and spectral
 * cycles can differ.                               
 *
 ********************************************************/

void
spectrum_allocate (int nspec)
{
  int i;

  for (i = 0; i < nspec; ++i)
  {
    if ((xxspec[i].f = calloc (NWAVE_MAX, sizeof (*xxspec[i].f))) == NULL)
    {
      Error ("Unable to allocate memory for xxspec[%d].f array with %d bins\n", i, NWAVE_MAX);
      Exit (EXIT_FAILURE);
    }
    if ((xxspec[i].lf = calloc (NWAVE_MAX, sizeof (*xxspec[i].lf))) == NULL)
    {
      Error ("Unable to allocate memory for xxspec[%d].lf array with %d bins\n", i, NWAVE_MAX);
      Exit (EXIT_FAILURE);
    }
    if ((xxspec[i].f_wind = calloc (NWAVE_MAX, sizeof (*xxspec[i].f_wind))) == NULL)
    {
      Error ("Unable to allocate memory for xxspec[%d].f_wind array with %d bins\n", i, NWAVE_MAX);
      Exit (EXIT_FAILURE);
    }
    if ((xxspec[i].lf_wind = calloc (NWAVE_MAX, sizeof (*xxspec[i].lf_wind))) == NULL)
    {
      Error ("Unable to allocate memory for xxspec[%d].lf_wind array with %d bins\n", i, NWAVE_MAX);
      Exit (EXIT_FAILURE);
    }
  }
}




/**********************************************************/
/**
 * @brief      allocates space for and initializes the spectrum arrays, storing
 * the criteria that are later used to create each spectrum.
 *
 * @param [in] double  f1   The minimum frequency for the spectra to be created
 * @param [in] double  f2   The maximum frequency in the spectra
 * @param [in] int  nangle   The number of different inclination angles (or more properly the number of spectra) to be created
 * @param [in] double  angle[]   The inclination angles associated with each spectrum
 * @param [in] double  phase[]   The orbital phase associated with each spectrum
 * @param [in] int  scat_select[]   A parameter for each spectrum which allows one to construct spectra with specifid numbers of scatters
 * @param [in] int  top_bot_select[]   A code which allows one to select photons only from below or above the disk
 * @param [in] int  select_extract   FALSE for Live or Die option, TRUE  for a normal extraction
 * @param [in] double  rho_select[]   Rho coordinates for extracting only photons in a particular region
 * @param [in] double  z_select[]   Z cooordiante for selecting only photons that scattered or were created in a praticular region
 * @param [in] double  az_select[]  Aximuthal angle for selecting only photons in a particular region
 * @param [in] double  r_select[]   Radius of the region from which to select photons
 * @return     Always returns 0
 *
 * @details
 * The first time spectrum_init  is called (i.e. if ispec_start=0), it allocates memory
 * for the various spectrum arrays.  (This is done one time, so one needs to allocate
 * space for all the arrays even though they are not all used in the ionization step).
 * The total number of spectra created is nangle+MSPEC.)
 *
 * On subsequent calls to spectrum_init, it rezeros all the spectrum information and
 * calculates the other information associated with each spectrum, such as the
 * angle cosines and the names of each spectrum.
 *
 * ### Notes ###
 * angle[],phase[] and scat_select[] only apply to the spectra extracted at
 * specific angles.  angle[0], phase[0], and scat_select[0] all affect spec[3]
 *
 * scat_select allows one to select spectra with a specific number of scattere or range og
 * scatters.  If nscat > 999 select all.  This is the normal case. The rest are used
 * largely for diagnostic purposes.  If 0 <nscat < MAXScat selected only photons which
 * have scattered nscat times.  If nscat is negattive, then all photons with more than |nscat| re
 * included.
 *
 * If top_bot_select is 0, then all photons are counted in the spectrum; if > 0, then only
 * photons above the disk plane are selected, if <0 only those from below the disk plane are selected.
 *
 * rho_select,z_select,az_select,r_select define a region of space and can be used
 * to create spectra only from photons that scattered or were created in these particular
 * regions.  It is normally only used for diagnostic reasons
 *
 *
 * Warning - Do not put anything in this routine that does anything but initialize
 * or reinitialize the spectrum structures. This is important because this routine
 * is not accessed if one is continuing an old calculation of the detailed spectrum.
 * It is still used on a restart where the detailed spectral cycles have not begun
 * because in that case the spectra are not saved.
 *
 **********************************************************/
int
spectrum_init (f1, f2, nangle, angle, phase, scat_select, top_bot_select, select_extract, rho_select, z_select, az_select, r_select)
     double f1, f2;
     int nangle;
     double angle[], phase[];
     int scat_select[], top_bot_select[];
     int select_extract;
     double rho_select[], z_select[], az_select[], r_select[];
{
  int i, n;
  int nspec;
  double freqmin, freqmax, dfreq;
  double lfreqmin, lfreqmax, ldfreq;
  double x1, x2;
  char dummy[20];

  nspec = nangle + MSPEC;

  /* Allocate memory for the spectrum arrays the first time routine is called */

  if (spec_initialized == FALSE)
  {
    nspectra = nspec;           /* Note that nspectra is a global variable */
    xxspec = calloc (sizeof (spectrum_dummy), nspec);
    if (xxspec == NULL)
    {
      Error ("spectrum_init: Could not allocate memory for %d spectra with %d wavelengths\n", nspec, NWAVE_EXTRACT);
      Exit (EXIT_FAILURE);
    }
    spectrum_allocate (nspec);
    spec_initialized = TRUE;
  }

  /* Zero or rezero all spectral matrices */

  for (i = 0; i <= MAXSCAT; i++)
    nscat[i] = nres[i] = 0;

  for (i = 0; i < NSTAT; i++)
    nstat[i] = 0;

  /* Setup bins for linear and logarithmic spectrum. During the ionization
   * and spectral cycles, we use a different number of wavelength bins. */

  if (geo.ioniz_or_extract == CYCLE_IONIZ)
  {
    NWAVE_NOW = NWAVE_IONIZ;
  }
  else
  {
    NWAVE_NOW = NWAVE_EXTRACT;
  }

  freqmin = f1;
  freqmax = f2;
  dfreq = (freqmax - freqmin) / NWAVE_NOW;

  lfreqmin = log10 (freqmin);
  lfreqmax = log10 (freqmax);
  ldfreq = (lfreqmax - lfreqmin) / NWAVE_NOW;

  for (n = 0; n < nspec; n++)
  {
    xxspec[n].lmn[0] = xxspec[n].lmn[1] = xxspec[n].lmn[2] = 0.0;
    xxspec[n].renorm = 1.0;
    xxspec[n].freqmin = freqmin;
    xxspec[n].freqmax = freqmax;
    xxspec[n].dfreq = dfreq;
    xxspec[n].lfreqmin = lfreqmin;
    xxspec[n].lfreqmax = lfreqmax;
    xxspec[n].ldfreq = ldfreq;

    for (i = 0; i < NSTAT; i++)
    {
      xxspec[n].nphot[i] = 0;
    }

    for (i = 0; i < NWAVE_MAX; i++)
    {
      xxspec[n].f[i] = 0;
      xxspec[n].lf[i] = 0;
      xxspec[n].f_wind[i] = 0;
      xxspec[n].lf_wind[i] = 0;
    }
  }

  strcpy (xxspec[SPEC_CREATED].name, "Created");
  strcpy (xxspec[SPEC_CWIND].name, "WCreated");
  strcpy (xxspec[SPEC_EMITTED].name, "Emitted");
  strcpy (xxspec[SPEC_CENSRC].name, "CenSrc");
  strcpy (xxspec[SPEC_DISK].name, "Disk");
  strcpy (xxspec[SPEC_WIND].name, "Wind");
  strcpy (xxspec[SPEC_HITSURF].name, "HitSurf");
  strcpy (xxspec[SPEC_SCATTERED].name, "Scattered");

  for (n = MSPEC; n < nspec; n++)
  {

    /*
       We want to set up the direction cosines for extractions.  We have to be careful
       about the sense of the orbital phase within the program.  Viewed from the "north"
       pole in a system with the secondary along the x axis the observer moves clockwise
       and phases just before 0 should be in the +x + y quadrant (since we have the disk and
       wind rotating counter clockwise as viewed from the top.  Another way of saying this
       is at phases just before 0, e.g. 0.98, the observer sees the receding side of the
       disk. The minus sign in the terms associated with phase are to make this happen.
       02feb ksl
     */

    sprintf (xxspec[n].name, "A%02.0f", angle[n - MSPEC]);
    xxspec[n].lmn[0] = sin (angle[n - MSPEC] / RADIAN) * cos (-phase[n - MSPEC] * 360. / RADIAN);
    xxspec[n].lmn[1] = sin (angle[n - MSPEC] / RADIAN) * sin (-phase[n - MSPEC] * 360. / RADIAN);
    xxspec[n].lmn[2] = cos (angle[n - MSPEC] / RADIAN);
    Log_silent ("Angle %e Angle cosines:%e %e %e\n", angle[n - MSPEC], xxspec[n].lmn[0], xxspec[n].lmn[1], xxspec[n].lmn[2]);

    /* Initialize variables needed for live or die option.

       There are  various issues associated with extracting a spectrum
       at inclinations near 0 or 90 deg in the live or die mode

       At 90 degrees, the isseus arise, because we normally
       extract on both sides of the disk, that is to say if we want to
       get the flux at 45d, we actually use bands at 45 and 135 degrees,
       explicitly assuming that the program only deals with winds which are
       biconical.

       But if we choose 90 degrees for extraction we are extracting
       basically from 88-92 degrees, not as in the case of 45, from 43-47, and
       133-137.

       Similar issues occur at very low inclination angles near the poles.

     */


    x1 = angle[n - MSPEC] - DANG_LIVE_OR_DIE;
    x2 = angle[n - MSPEC] + DANG_LIVE_OR_DIE;

    /* Adjust when the direction one wishes to extract is close to the poles */
    if (x1 < 0.)
      x1 = 0;
    if (x2 > 180.)
      x2 = 180.;

    /* Adjust when the direction one wishes to extract is close to the disk or xy plane */

    if (x1 < 90 && x2 > 90)
    {
      if (90 - x1 < x2 - 90)
      {
        x1 = 90;
      }
      else
      {
        x2 = 90;
      }
    }

    x1 = fabs (cos (x1 / RADIAN));
    x2 = fabs (cos (x2 / RADIAN));
    if (x1 > x2)
    {
      xxspec[n].mmax = x1;
      xxspec[n].mmin = x2;
    }
    else
    {
      xxspec[n].mmax = x2;
      xxspec[n].mmin = x1;
    }

    if (select_extract == FALSE)        // We are in Live or Die mode
    {
      xxspec[n].renorm = 1. / (xxspec[n].mmax - xxspec[n].mmin);

    }
    else                        // No renormalization in extract mode
      xxspec[n].renorm = 1.;

    /* Completed initialization of variables for live or die */

    strcpy (dummy, "");
    sprintf (dummy, "P%04.2f", phase[n - MSPEC]);
    strcat (xxspec[n].name, dummy);

    xxspec[n].nscat = scat_select[n - MSPEC];

    if (xxspec[n].nscat < MAXSCAT)
    {                           /* Then conditions have been placed on the
                                   number of scatters to be included so update the names */
      strcpy (dummy, "");
      if (xxspec[n].nscat > MAXSCAT)
        sprintf (dummy, "_sc:all");
      if (xxspec[n].nscat >= 0)
        sprintf (dummy, "_sc:%d", xxspec[n].nscat);
      else
        sprintf (dummy, "_sc:>%d", -xxspec[n].nscat);
      strcat (xxspec[n].name, dummy);
    }
    xxspec[n].top_bot = top_bot_select[n - MSPEC];
    if (xxspec[n].top_bot != 0)
    {                           /* Then conditions have been placed on the last
                                   location of the photon so update the names */
      strcpy (dummy, "");
      if (xxspec[n].top_bot == 1)
        sprintf (dummy, "_top");
      else if (xxspec[n].top_bot == -1)
        sprintf (dummy, "_bot");
      else if (xxspec[n].top_bot == 2)
      {
        sprintf (dummy, "_pos");
        xxspec[n].x[0] = rho_select[n - MSPEC] * cos (az_select[n - MSPEC] / RADIAN);
        xxspec[n].x[1] = rho_select[n - MSPEC] * sin (az_select[n - MSPEC] / RADIAN);
        xxspec[n].x[2] = z_select[n - MSPEC];
        xxspec[n].r = r_select[n - MSPEC];

      }
      else
      {
        Error ("spectrum_init: Unknown option %d\n", xxspec[n].top_bot);
        Exit (0);

      }
      strcat (xxspec[n].name, dummy);
    }
  }

  return (0);
}




/**********************************************************/
/**
 * @brief      Increments the spectrum arrays
 *  	after each flight of photons is processed (during ionization
 *  	cycles and for detailed spectra in the Live or Die option).
 *
 * @param [in] PhotPtr  p   A flight of photons
 * @param [in] int  nangle  The number of different angles and phases for which to create detailed spectra
 * @param [in] int  select_extract   Parameter to select whether to use the Live or Die (0) or extract option
 * @return     Always returns 0
 *
 * @details
 * This routine increments the total spectrum arrays based on what has happened to each
 * photon during ionization cycles.  In the Live or Die option, the spectra at specific angles
 * are also created here when detailed spectra are created.
 *
 * The routine is called after each batch of photons has been transported through the wind and
 * prints some intermediate results to assure the user that the program is still running.
 *
 * ### Notes ###
 * Summing up of the spectra in the "extract" option is done in extract.c
 *
 * spectrum_init has to be called prior to spectrum_create to set
 * up the frequency limits, etc. of the spectra.
 *
 * Two versions of the spectra are created, one with linear binning and one with
 * with logarithmic binning.
 *
 *
 **********************************************************/

int
spectrum_create (p, nangle, select_extract)
     PhotPtr p;
     int nangle;
     int select_extract;

{
  int nphot, istat, j, k, k1, n;
  int nspec, nwave, spectype;
  double freqmin, freqmax, dfreq, ldfreq;
  double x1;
  int mscat, mtopbot;
  double delta;
  double nlow, nhigh;
  int k_orig, k1_orig;
  int iwind;
  int max_scat, max_res;

  if (geo.ioniz_or_extract == CYCLE_IONIZ)
  {
    nwave = NWAVE_IONIZ;
  }
  else
  {
    nwave = NWAVE_EXTRACT;
  }

  /* Setup frequency boundaries, etc using values set in spectrum_int */

  freqmin = xxspec[SPEC_CREATED].freqmin;
  freqmax = xxspec[SPEC_CREATED].freqmax;
  dfreq = xxspec[SPEC_CREATED].dfreq;
  ldfreq = xxspec[SPEC_CREATED].ldfreq;

  nspec = nangle + MSPEC;
  nlow = 0.0;                   // variable to store the number of photons that have frequencies which are too low
  nhigh = 0.0;                  // variable to store the number of photons that have frequencies which are too high
  delta = 0.0;                  // fractional frequency error allowod

  for (nphot = 0; nphot < NPHOT; nphot++)
  {
    if ((j = p[nphot].nscat) < 0 || j > MAXSCAT)
      nscat[MAXSCAT]++;
    else
      nscat[j]++;

    if ((j = p[nphot].nrscat) < 0 || j > MAXSCAT)
      nres[MAXSCAT]++;
    else
      nres[j]++;

    /*
     * Determine whether this is a wind photon, that is was it created in the
     * wind or scattered by the wind
     */

    iwind = FALSE;
    if (p[nphot].origin == PTYPE_WIND || p[nphot].origin == PTYPE_WIND_MATOM || p[nphot].nscat > 0)
    {
      iwind = TRUE;
    }

    /* Find the bins that need to be incremented in the logarithimic and linear grid. This
     * has to be done twice, because we want this for the original and final frequencey of
     * the photoon
     */

    k1 = (int) ((log10 (p[nphot].freq) - log10 (freqmin)) / ldfreq);
    if (k1 < 0)
    {
      k1 = 0;
    }
    if (k1 > nwave - 1)
    {
      k1 = nwave - 1;
    }

    k1_orig = (int) ((log10 (p[nphot].freq_orig) - log10 (freqmin)) / ldfreq);
    if (k1_orig < 0)
    {
      k1_orig = 0;
    }
    if (k1_orig > nwave - 1)
    {
      k1_orig = nwave - 1;
    }


    k = (int) ((p[nphot].freq - freqmin) / dfreq);
    if (k < 0)
    {
      if (((1. - p[nphot].freq / freqmin) > delta) && (geo.rt_mode != RT_MODE_MACRO))
        nlow = nlow + 1;
      k = 0;
    }
    else if (k > nwave - 1)
    {
      if (((1. - freqmax / p[nphot].freq) > delta) && (geo.rt_mode != RT_MODE_MACRO))
        nhigh = nhigh + 1;
      k = nwave - 1;
    }

    k_orig = (int) ((p[nphot].freq_orig - freqmin) / dfreq);
    if (k_orig < 0)
    {
      if (((1. - p[nphot].freq_orig / freqmin) > delta) && (geo.rt_mode != RT_MODE_MACRO))
        nlow = nlow + 1;
      k_orig = 0;
    }
    else if (k_orig > nwave - 1)
    {
      if (((1. - freqmax / p[nphot].freq_orig) > delta) && (geo.rt_mode != RT_MODE_MACRO))
        nhigh = nhigh + 1;
      k_orig = nwave - 1;
    }

    /* Having worked out what spectral bins to increment, we now actually increment the various spectra */

    istat = p[nphot].istat;

    if (p[nphot].origin == PTYPE_WIND || p[nphot].origin == PTYPE_WIND_MATOM)
    {
      xxspec[SPEC_CWIND].f[k_orig] += p[nphot].w_orig;
      xxspec[SPEC_CWIND].lf[k1_orig] += p[nphot].w_orig;
      xxspec[SPEC_CWIND].nphot[istat]++;
    }
    else
    {
      xxspec[SPEC_CREATED].f[k_orig] += p[nphot].w_orig;
      xxspec[SPEC_CREATED].lf[k1_orig] += p[nphot].w_orig;
      xxspec[SPEC_CREATED].nphot[istat]++;
    }

    if (iwind)
    {
      xxspec[SPEC_CREATED].f_wind[k_orig] += p[nphot].w_orig;
      xxspec[SPEC_CREATED].lf_wind[k1_orig] += p[nphot].w_orig;
    }

    if (istat == P_ESCAPE)
    {
      xxspec[SPEC_EMITTED].f[k] += p[nphot].w;
      xxspec[SPEC_EMITTED].lf[k1] += p[nphot].w;
      if (iwind)
      {
        xxspec[SPEC_EMITTED].f_wind[k] += p[nphot].w;
        xxspec[SPEC_EMITTED].lf_wind[k1] += p[nphot].w;
      }
      xxspec[SPEC_EMITTED].nphot[istat]++;
      spectype = p[nphot].origin;

      /* When a photon that originated for example in the BL which has a type of PTYPE_BL is scattered in the wind by 
       * a macro atom it's type is increased by 10.  When we want to construct a spectrum for photons originating
       * from the boundary layer we need to subtract 10 from the type.    See sirocco.h 
       */
      if (spectype >= 10)
        spectype -= 10;

      if (p[nphot].nmacro == 0 && (spectype == PTYPE_STAR || spectype == PTYPE_BL || spectype == PTYPE_AGN))
      {
        xxspec[SPEC_CENSRC].f[k] += p[nphot].w;
        xxspec[SPEC_CENSRC].lf[k1] += p[nphot].w;
        if (iwind)
        {
          xxspec[SPEC_CENSRC].f_wind[k] += p[nphot].w;
          xxspec[SPEC_CENSRC].lf_wind[k1] += p[nphot].w;
        }
        xxspec[SPEC_CENSRC].nphot[istat]++;
      }
      else if (p[nphot].nmacro == 0 && spectype == PTYPE_DISK)
      {
        xxspec[SPEC_DISK].f[k] += p[nphot].w;
        xxspec[SPEC_DISK].lf[k1] += p[nphot].w;
        if (iwind)
        {
          xxspec[SPEC_DISK].f_wind[k] += p[nphot].w;
          xxspec[SPEC_DISK].lf_wind[k1] += p[nphot].w;
        }
        xxspec[SPEC_DISK].nphot[istat]++;
      }
      else if (spectype == PTYPE_WIND || p[nphot].nmacro > 0)
      {
        /* In macro atom mode a photon is regarded as being in the wind if it has had a macro atom interaction */
        xxspec[SPEC_WIND].f[k] += p[nphot].w;   /* wind spectrum */
        xxspec[SPEC_WIND].lf[k1] += p[nphot].w; /* logarithmic wind spectrum */
        if (iwind)
        {
          xxspec[SPEC_WIND].f_wind[k] += p[nphot].w;    /* emitted spectrum */
          xxspec[SPEC_WIND].lf_wind[k1] += p[nphot].w;  /* logarithmic emitted spectrum */
        }
        xxspec[SPEC_WIND].nphot[istat]++;
      }
      else
      {
        Error ("spectrum_create: Unknown photon type %d\n", spectype);
      }

      /* For Live or Die option, increment the spectra here */
      if (select_extract == FALSE)
      {
        x1 = fabs (p[nphot].lmn[2]);
        for (n = MSPEC; n < nspec; n++)
        {

          /* Complicated if statement to allow one to choose whether to construct the spectrum
             from all photons or just from photons which have scattered a specific number
             of times.  */

          if (((mscat = xxspec[n].nscat) >= MAXSCAT ||
               p[nphot].nscat == mscat ||
               (mscat < 0 && p[nphot].nscat >= (-mscat))) && ((mtopbot = xxspec[n].top_bot) == 0 || (mtopbot * p[nphot].x[2]) > 0))
          {
            if (xxspec[n].mmin < x1 && x1 < xxspec[n].mmax)
            {
              xxspec[n].f[k] += p[nphot].w;
              xxspec[n].lf[k1] += p[nphot].w;
              if (iwind)
              {
                xxspec[n].f_wind[k] += p[nphot].w;
                xxspec[n].lf_wind[k1] += p[nphot].w;
              }
            }
          }
        }
      }
    }
    else if (istat == P_HIT_STAR || istat == P_HIT_DISK)
    {
      xxspec[SPEC_HITSURF].nphot[istat]++;
    }

    if (istat == P_ESCAPE && (p[nphot].nscat > 0 || p[nphot].nrscat > 0))
    {
      xxspec[SPEC_SCATTERED].f[k] += p[nphot].w;
      xxspec[SPEC_SCATTERED].lf[k1] += p[nphot].w;
      if (p[nphot].w > 0 && p[nphot].w < 1e-100)
      {
        Log ("spectrum_create: very small weight %e (%e) for phot %d\n", p[nphot].w, p[nphot].w_orig, nphot);
      }
      if (iwind)
      {
        xxspec[SPEC_SCATTERED].f_wind[k] += p[nphot].w;
        xxspec[SPEC_SCATTERED].lf_wind[k1] += p[nphot].w;
      }
      if (istat < 0 || istat > NSTAT - 1)
        xxspec[SPEC_SCATTERED].nphot[NSTAT - 1]++;
      else
        xxspec[SPEC_SCATTERED].nphot[istat]++;
    }

    if (istat < 0 || istat > NSTAT - 1)
      nstat[NSTAT - 1]++;
    else
      nstat[istat]++;
  }

  /* At this point all of the spectra have been incremented and so we performe a simple check on the number of
     photons were lost and then we print out some statistics having to do with the number of 
     scatters each photon has undergone.
   */

  if ((nlow / nphot > 0.05) || (nhigh / nphot > 0.05))
  {
    Error ("spectrum_create: Fraction of photons lost: %4.2f wi/ freq. low, %4.2f w/freq hi\n", nlow / nphot, nhigh / nphot);
  }
  else
  {
    Log ("spectrum_create: Fraction of photons lost:  %4.2f wi/ freq. low, %4.2f w/freq hi\n", nlow / nphot, nhigh / nphot);
  }

  max_scat = max_res = 0;

  for (j = 1; j < MAXSCAT; j++)
  {
    if (nscat[j] > 0)
    {
      max_scat = j;
    }
    if (nres[j] > max_res)
    {
      max_res = j;
    }
  }

  Log ("\nNo. of photons which have scattered n times.     The max number of scatters seen was %d\n", max_scat);

  for (j = 0; j <= max_scat; j++)
  {
    Log ("%-9.3g", (double) nscat[j]);
    if ((j % 10) == 9)
      Log ("\n");
  }

  Log ("\nNumber of photons resonantly scattering n times.  The max number of scatters seen was %d\n", max_res);
  for (j = 0; j <= max_res; j++)
  {
    Log ("%-9.3g", (double) nres[j]);
    if ((j % 10) == 9)
      Log ("\n");
  }
  Log ("\nNo of photons and their fates\n!!PhotFate: ");
  for (j = 0; j < NSTAT; j++)
  {
    Log ("%-9.3g", (double) nstat[j]);
    if ((j % 10) == 9)
      Log ("\n");
  }
  Log ("\n");



  Log ("Photons contributing to the various spectra\n");
  Log ("                      Inwind    Scat      Esc       Star      >nscat    err       Absorb    Disk     Sec        Adiab(matom)\n");
  for (n = 0; n < nspectra; n++)
  {
    Log ("%20s ", xxspec[n].name);
    for (j = 0; j < NSTAT; j++)
      Log (" %-9.3g", (double) xxspec[n].nphot[j]);
    Log ("\n");
  }


  return (0);

}





/**********************************************************/
/**
 * @brief      add a single photon to a specific spectrum
 * on the fly  
 *
 * @param [in] PhoPtr p  A single photon
 * @param [in] int spec_id  The spectrum to be incremented   
 * @return     Always returns 0
 *
 * @details
 *
 * This routine allows a spectrum to be incremented during 
 * during its flight through the wind.  It is necessary
 * for example if one wishes to record when a photon
 * hits the disk (and is reflected).
 *
 * ### Notes ###
 * The routine assumes one wants to use the current 
 * frequency/weight of the photon.  If for certain spectra
 * one wants to record the original frequency or weight
 * some additional logic is required.  
 *
 * If spectra are accumulated during the flight of photons
 * through the plasma one needs to avoid double counting
 * in spectrum_create
 *
 **********************************************************/


int
spec_add_one (p, spec_type)
     PhotPtr p;
     int spec_type;
{
  int k;
  int iwind;
  double freq;

  freq = p->freq;


  iwind = FALSE;
  if (p->origin == PTYPE_WIND || p->origin == PTYPE_WIND_MATOM || p->nscat > 0)
  {
    iwind = TRUE;
  }


  k = (freq - xxspec[spec_type].freqmin) / xxspec[spec_type].dfreq;

  if (k > NWAVE_NOW - 1)
    k = NWAVE_NOW - 1;
  else if (k < 0)
    k = 0;

  xxspec[spec_type].f[k] += p->w;

  if (iwind)
  {
    xxspec[SPEC_SCATTERED].f_wind[k] += p->w;
  }

  k = (log10 (freq) - xxspec[spec_type].lfreqmin) / xxspec[spec_type].ldfreq;

  if (k > NWAVE_NOW - 1)
    k = NWAVE_NOW - 1;
  else if (k < 0)
    k = 0;

  xxspec[spec_type].lf[k] += p->w;
  if (iwind)
  {
    xxspec[SPEC_SCATTERED].lf_wind[k] += p->w;
  }

  return (0);

}




/**********************************************************/
/**
 * @brief      writes out the spectra to a file
 *
 * @param [in] char  filename[]   The name of the file to write
 * @param [in] int  nspecmin   The number of the first spectrum to write
 * @param [in] int  nspecmax   The number of the last spectrum to write              .
 * @param [in] int  select_spectype   The type of spectral file you want to create
 * @param [in, out] double  renorm   This is renormalization which incrementally decreases to
 * one as the detailed spectral calculation goes forward.  It
 * was added to allow one to print out the spectrum at the
 * end of each cycle, rather than the end of the entire
 * calculation.
 * @param [in] int loglin FALSE to print the spectrum out in linear units, TRUE in log units
 * @param [in] int iwind If FALSE, print out the normal spectrum; if TRUE, print
 * out only photons that were scattered or created in the wind.
 *
 * @return     Always returns 0, unless there is a major problem in which case the program
 * exits
 *
 * @details
 * This simple routine simply writes the spectra to a file in an easily interpretable
 * ascii format.  The spectra will have already been created (using spectrum_create). 
 *
 * Normally one would write all of the spectra in one go, but  one can use
 * spectrum summary to write various spectra to various files by using the variables
 * nspecmin and nspecmax.
 *
 * Normally s[0],s[1],and s[2] will be the escaping, scattered, and absorbed spectrum.
 * The rest will be those which have been "extracted".
 *
 * It is called multiple times. In Python, it is currently called at two different
 * locations in the code, once at the
 * end of each ionization cycle  and at the end of each spectrum cycle.               
 * 
 * For the the spectrum cycles, the spectrum is 
 * and renormalized so that the overall normalization of the source
 * does not change.  This enables one to plot  the spectrum as the routine
 * is continuing to calculate the detatiled spectra, improving the statistics
 *
 * The rooutine is called separately to write the spectra out with a linear and
 * logarithmic binning.
 *
 * The file includes a copy of the .pf file inputs for the program so that one can
 * track exactly what inputs wre used to create the spectrum.
 *
 * ### Notes ###
 *
 **********************************************************/

int
spectrum_summary (filename, nspecmin, nspecmax, select_spectype, renorm, loglin, iwind)
     char filename[];
     int loglin;
     int nspecmin, nspecmax;
     int select_spectype;
     double renorm;
     int iwind;

{
  FILE *fopen (), *fptr;
  int i, n;
  int nwave;
  char string[LINELENGTH];
  double freq, freqmin, dfreq, freq1;
  double lfreqmin, ldfreq;
  double x, dd;

  if (geo.ioniz_or_extract == CYCLE_IONIZ)
  {
    nwave = NWAVE_IONIZ;
  }
  else
  {
    nwave = NWAVE_EXTRACT;
  }

  if ((fptr = fopen (filename, "w")) == NULL)
  {
    Error ("spectrum_summary: Unable to open %s for writing\n", filename);
    Exit (0);
  }

  if (nspecmin < 0 || nspecmax < 0 || nspecmin > nspecmax)
  {
    Error ("spectrum_summary: nspecmin %d or nspecmax %d not reasonable \n", nspecmin, nspecmax);
    Exit (0);
  }

  /* Construct and write a header string  for the output file */
  fprintf (fptr, "# Python Version %s\n", VERSION);
  fprintf (fptr, "# Git commit hash %s\n", GIT_COMMIT_HASH);

  get_time (string);
  fprintf (fptr, "# Date	%s\n#  \n", string);

  if (select_spectype == SPECTYPE_RAW)
  {
    fprintf (fptr, "\n# Units: L_nu spectrum (erg/s/Hz)\n\n");
  }
  else if (select_spectype == SPECTYPE_FLAMBDA)
  {
    fprintf (fptr, "\n# Units: flambda spectrum (erg/s/cm^2/A) at %.1f parsecs\n\n", D_SOURCE);
  }
  else if (select_spectype == SPECTYPE_FNU)
  {
    fprintf (fptr, "\n# Units: Fnu spectrum (erg/s/cm^2/Hz) at %.1f parsecs\n\n", D_SOURCE);
  }
  else
  {
    Error ("spectrum_summary: Unknown select_spectype %d\n", select_spectype);
    Exit (0);
  }


  /* Save all of the parameter file information to the spectrum file */

  rdpar_save (fptr);


  /* Write the rest of the header for the spectrum file */

  fprintf (fptr, "# \nFreq.             Lambda    ");

  for (n = nspecmin; n <= nspecmax; n++)
  {
    fprintf (fptr, " %-10s", xxspec[n].name);
  }


  fprintf (fptr, "\n");


  /* Ignore the end bins because they include all photons outside the frequency range and there may be some
     as a result of the fact that the bb function generate some IR photons */

  dd = 4. * PI * (D_SOURCE * PC) * (D_SOURCE * PC);

  if (loglin == FALSE)          /* Then write the linear version of the spectra */
  {
    freqmin = xxspec[nspecmin].freqmin;
    dfreq = xxspec[nspecmin].dfreq;
    for (i = 1; i < nwave - 1; i++)
    {
      freq = freqmin + i * dfreq;
      fprintf (fptr, "%-8e %10.5e ", freq, VLIGHT * 1e8 / freq);
      for (n = nspecmin; n <= nspecmax; n++)
      {
        x = xxspec[n].f[i] * xxspec[n].renorm;
        if (iwind)
        {
          x = xxspec[n].f_wind[i] * xxspec[n].renorm;
        }


        if (select_spectype == SPECTYPE_FLAMBDA)
        {
          x *= (freq * freq * 1e-8) / (dfreq * dd * VLIGHT);
        }
        else if (select_spectype == SPECTYPE_FNU)
        {
          x /= (dfreq * dd);
        }
        else if (select_spectype == SPECTYPE_RAW)
        {                       /*generated spectrum */
          x /= (dfreq);
        }
        fprintf (fptr, " %10.5g", x * renorm);
      }

      fprintf (fptr, "\n");
    }
  }
  else if (loglin == TRUE)
  {
    freq1 = lfreqmin = xxspec[nspecmin].lfreqmin;
    ldfreq = xxspec[nspecmin].ldfreq;

    for (i = 1; i < nwave - 1; i++)
    {
      freq = pow (10., (lfreqmin + i * ldfreq));
      dfreq = freq - freq1;
      fprintf (fptr, "%-8e %-8.4g ", freq, VLIGHT * 1e8 / freq);
      for (n = nspecmin; n <= nspecmax; n++)
      {
        x = xxspec[n].lf[i] * xxspec[n].renorm;
        if (iwind)
        {
          x = xxspec[n].lf_wind[i] * xxspec[n].renorm;
        }

        if (select_spectype == SPECTYPE_FLAMBDA)
        {
          x *= (freq * freq * 1e-8) / (dfreq * dd * VLIGHT);
        }
        else if (select_spectype == SPECTYPE_FNU)
        {
          x /= (dfreq * dd);
        }
        else if (select_spectype == SPECTYPE_RAW)
        {
          x /= (dfreq);
        }
        fprintf (fptr, " %10.5g", x * renorm);
      }

      fprintf (fptr, "\n");
      freq1 = freq;
    }
  }
  fclose (fptr);

  return (0);

}









/**********************************************************/
/**
 * @brief      renormalizes the detailed spectra in case
 * of a restart
 *
 * @param [in] int  nangle   The number of discrete angles
 * @return     Always returns 0
 *
 * @details
 *
 * This routine deals with a very special case when one is
 * restarting a previous run, and adding spectral cycles to
 * obtain a detailed spectrum with higher statistics.
 *
 * In the previous run, the spectral will have been normalized
 * so that one gets the correct flux for the simulation.  We
 * are now computing additional spectral cycles, so we have
 * to reduce the fluxes for the spectral cycles that are read
 * in, so we can continue.
 * 	renorm_factor = (cycles in old pf file) / (cycles in new pf file)
 * 
 * Restarts can also occur when one is just completing a previous
 * run, without increasing the number of cycles, and in this case
 * no renormalization is required.
 * ### Notes ###
 * see Issue #134 and more importantly #503
 *
 **********************************************************/

int
spectrum_restart_renormalise (nangle)
     int nangle;
{
  double renorm_factor;
  int n, m, nspec;

  if (geo.pcycles == geo.pcycles_renorm)
  {
    return (0);
  }

  nspec = nangle + MSPEC;

  /* If we have gooten to this point, then the number of pcycles to which the
   * spectrum has been renormalized has changed, and so we must renormalize
   * the spectrum.  The original spectrum current is normalized to g
   * geo.pcycle/geo.pcycles_renorm of the final value.  We want the psectrum
   * to be geo.pcycles/geo.pcycles of the final value, so the renormlization
   * factor is geo.pcycles_renorm/geo.pcycles.
   */


  renorm_factor = ((double) geo.pcycles_renorm) / ((double) geo.pcycles);


  Log ("spectrum_restart_renormalise:  %d %d %d %f\n", geo.pcycle, geo.pcycles_renorm, geo.pcycles, renorm_factor);

  /* loop over each spectrum column and each wavelength bin */
  for (n = 0; n < nspec; n++)
  {
    for (m = 0; m < NWAVE_EXTRACT; m++)
    {
      xxspec[n].f[m] *= renorm_factor;
      xxspec[n].lf[m] *= renorm_factor;
    }
  }

  return (0);
}
