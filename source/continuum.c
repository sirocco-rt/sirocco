
/***********************************************************/
/** @file  continuum.c
 * @author ksl
 * @date May, 2018
 *
 * @brief  Rotines related speciically to generating MC spectra for
 * stellar atmospheres like grids of spectra
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atomic.h"
#include "sirocco.h"
#include "models.h"

double old_t, old_g, old_freqmin, old_freqmax;
//OLD double jump[] = { 913.8 };


/**********************************************************/
/**
 * @brief      get a photon frequency from a grid of stored spectra
 * from a continuum grid
 *
 * @param [in] int  spectype   An index to the continuum grid from which one wants to obtaing a photon
 * @param [in] double  t   The associated temperature (or first variable) defining the grid
 * @param [in] double  g   The associted gravity (or second variable)
 * @param [in] double  freqmin   The minimum frequency desiered for a photon
 * @param [in] double  freqmax   The maximum frequecny for a photon
 * @return   A frequency for a photon
 *
 * @details
 * This is a routine that is used to select a function from a grid of spectra
 * stored in the model structure.  For the purpose of this routine, the
 * model grid is assumed to be based on two variables t (temperature) and
 * gravity (g).  This reflects the fact that most of the grids used in Python
 * to date (e.g Kurucz models) are these kinds of grids.  The routine linearly
 * interpolates  on these two variables to produce a spectrum at the desired
 * t and g.  It then creates a comulative distribution of the portion of
 * the spectrum bewtween fmin and fmax and then uses a random number to
 * obtain a frequncy for a single photon. 
 *
 *
 * ### Notes ###
 *
 * The model structure is general in the sense that it was intended
 * for situation with any number of variables.
 * However what is written here
 * is specific to the case of two variables. This is a problem, as we have
 * found in trying to use this in other situations.  A more general routine is needed.
 * The first step in doing this is to replace much of the code here with
 * a call to the routine model (spectype, par) in models.c
 * This is issue #539, which as of 221109 remains open. The
 * note there states correctly that we should not attack this problem
 * until we have a specific example.
 *
 **********************************************************/


double
one_continuum (spectype, t, g, freqmin, freqmax)
     int spectype;
     double t, g, freqmin, freqmax;
{
  double lambdamin, lambdamax;
  double w_local[NCDF], f_local[NCDF];
  double f, y;
  int n, nwave;
  double par[2];
  int model ();

  /* Check if the parameters are the same as the stored ones, otherwise initialise */
  if (old_t != t || old_g != g || old_freqmin != freqmin || old_freqmax != freqmax)
  {                             /* Then we must initialize */
    if (comp[spectype].nmods == 1)      //If we only have one model (as is the case of an AGN model SED then dont interpolate
    {
      comp[spectype].xmod = mods[comp[spectype].modstart];
      old_t = t;                //These are dummies, but should prevent unwanted regeneration of the array
      old_g = g;
    }
    else if (t != comp[spectype].xmod.par[0] || g != comp[spectype].xmod.par[1])
    {
      par[0] = t;
      par[1] = g;
      model (spectype, par);
      old_t = t;
      old_g = g;
    }

    lambdamin = VLIGHT * 1e8 / freqmax;
    lambdamax = VLIGHT * 1e8 / freqmin;
    nwave = 0;

    /* Create the first element of the array, precisely at lambdamin if that is possible
       Specifically, if the first wavelength in the model is below the wavelength range in the simulation,
       interpolate on the model flux to get the flux at lambdamin. copy relevant wavelengths and
       fluxes to w_local and f_local  */
    if (comp[spectype].xmod.w[0] < lambdamin && lambdamin < comp[spectype].xmod.w[comp[spectype].nwaves - 1])
    {
      w_local[nwave] = lambdamin;
      linterp (lambdamin, comp[spectype].xmod.w, comp[spectype].xmod.f, comp[spectype].nwaves, &y, 0);
      f_local[nwave] = y;
      nwave++;
    }

    /* loop over rest of model wavelengths and fluxes and copy to w_local and f_local.
       This does not include the end points
     */

    for (n = 0; n < comp[spectype].nwaves; n++)
    {
      if (comp[spectype].xmod.w[n] > lambdamin && comp[spectype].xmod.w[n] < lambdamax)
      {
        w_local[nwave] = comp[spectype].xmod.w[n];
        f_local[nwave] = comp[spectype].xmod.f[n];
        nwave++;
      }
    }

    /* No add a point at lambdamax.  Specificly,  check if upper bound is beyond lambdamax, and if so, 
       interpolate to get appropriate flux
       at lambda max. copy to w_local and f_local */

    if (comp[spectype].xmod.w[0] < lambdamax && lambdamax < comp[spectype].xmod.w[comp[spectype].nwaves - 1])
    {
      w_local[nwave] = lambdamax;
      linterp (lambdamax, comp[spectype].xmod.w, comp[spectype].xmod.f, comp[spectype].nwaves, &y, 0);
      f_local[nwave] = y;
      nwave++;
    }

    /* There are two pathological cases to deal with, when we only have one non zero point,
       we need to make an extra point just up/down from the penultimate/second point so we
       can make a sensible CDF. */

    if (f_local[nwave - 2] == 0.0)      //We have a zero just inside the end
    {
      nwave++;
      w_local[nwave - 1] = w_local[nwave - 2];
      f_local[nwave - 1] = f_local[nwave - 2];
      w_local[nwave - 2] = w_local[nwave - 3] / (1. - DELTA_V / (2. * VLIGHT));
      linterp (w_local[nwave - 2], comp[spectype].xmod.w, comp[spectype].xmod.f, comp[spectype].nwaves, &y, 0);
      f_local[nwave - 2] = y;
    }

    /* we should now have our arrays w_local and f_local which can be used to generate a cdf */

    /*  Get_model returns wavelengths in Ang and flux in ergs/cm**2/Ang */



    if (cdf_gen_from_array (&comp[spectype].xcdf, w_local, f_local, nwave, lambdamin, lambdamax) != 0)
    {
      Error ("One_continuum: after return from cdf_gen_from_array\n");
    }

    old_freqmin = freqmin;
    old_freqmax = freqmax;

  }

  /* generate the frequency from the CDF that has been built up from the model fluxes */

  f = (VLIGHT * 1.e8 / cdf_get_rand (&comp[spectype].xcdf));

  /* check if the frequency is too small or too large, and default to simulation limits */
  if (f > freqmax)
  {
    Error ("one_continuum: f too large %e\n");
    f = freqmax;
  }
  if (f < freqmin)
  {
    Error ("one_continuum: f too small %e\n");
    f = freqmin;
  }
  return (f);
}


/**********************************************************/
/** 
 * @brief      get the surface flux for Hubeny/Kurucz like
 * stellar models
 *
 * @param [in] int  spectype   an integer identifying the set of models to use
 * @param [in] double  freqmin   The minimum frequency
 * @param [int] double  freqmax   The maximum frequency
 * @param [in] double  t   A temperature
 * @param [in] double  g   A gravity
 * @return     The surface flux for a star-like spectrum read from a grid
 *
 * @details
 * The routine gets the band-limited flux per unit area for Hubeny or Kurucz models
 * which are both in units of the Eddington flux (H).
 * To allow the same code to work with a much simpler model SED, where there is only
 * one model - there is a switch to disable interpolation. 
 *
 * ### Notes ###
 *
 * This is not a general purpose routine.  Not all spectral
 * models for stars are calculated in terms of the Eddington flux.  The
 * factor of 4 pi is the conversion from Eddington flux to physical flux
 *
 * For an explanation see Hubeny & Mihalas - Theory of Stellar Atmospheres,
 * equation 3.70 (or Mihalas, Stellar Atmospheres 2nd ed, eq (1-27).
 *
 * Note that the observed flux at a distance is given by
 *
 * f/H = 4pi * R**2/d**2
 *
 * and
 *
 * L=f*4pi*R**2 = 16 pi**2 R**2 H
 *
 *
 **********************************************************/

int integ_spectype;             //External variable pointing to the model for our Romburg interpolation.

double
emittance_continuum (spectype, freqmin, freqmax, t, g)
     int spectype;
     double freqmin, freqmax, t, g;
{
  int nwav, n;
  double x, lambdamin, lambdamax, w1, w2, f_interp;

  double par[2];
  int model ();

  lambdamin = VLIGHT / (freqmax * ANGSTROM);
  lambdamax = VLIGHT / (freqmin * ANGSTROM);

  /* Get the spectrum for which we want to calculate the emittance, interpolating if there is
     more than one model.
   */
  if (comp[spectype].nmods == 1)
  {
    comp[spectype].xmod = mods[comp[spectype].modstart];
  }
  else
  {
    par[0] = t;
    par[1] = g;
    model (spectype, par);
  }
  nwav = comp[spectype].nwaves;


  if (lambdamax > comp[spectype].xmod.w[nwav - 1] || lambdamin < comp[spectype].xmod.w[0])
  {
    Error ("emittance_continuum: freqmin %e freqmax %e for spectrum type  %d\n", freqmin, freqmax, spectype);
    Error ("emittance_continuum: emin %e emax %e\n", HEV * freqmin, HEV * freqmax);

    Error ("emittance_continuum: Requested wavelengths extend beyond models wavelengths for list %s\n", comp[spectype].name);
    Error ("lambda %f %f  model %f %f\n", lambdamin, lambdamax, comp[spectype].xmod.w[0], comp[spectype].xmod.w[nwav - 1]);
  }


  /* Now integrate over the desired wavelength/frequency interval */
  x = 0.0;
  for (n = 0; n < nwav - 1; n++)
  {
    w1 = comp[spectype].xmod.w[n];
    w2 = comp[spectype].xmod.w[n + 1];

    if (w2 < lambdamin || w1 > lambdamax)       //This band is outside the required band
    {
      continue;
    }
    else if (w1 < lambdamin && w2 > lambdamin)  //The bin brackets the lower band boundary
    {
      linterp (lambdamin, comp[integ_spectype].xmod.w, comp[integ_spectype].xmod.f, comp[integ_spectype].nwaves, &f_interp, 1); //Obtain the interpolated flux at lambdamin
      x += (w2 - lambdamin) * (f_interp + comp[spectype].xmod.f[n + 1]) / 2.0;  //Add the mean flux in the bin x the wavelength range
    }
    else if (w1 > lambdamin && w2 > lambdamin && w1 < lambdamax && w2 < lambdamax)      //The whole interval in in the requested band
    {
      x += (w2 - w1) * (comp[spectype].xmod.f[n] + comp[spectype].xmod.f[n + 1]) / 2.0; //Add the mean flux in the bin x the wavelength range
    }
    else if (w1 < lambdamax && w2 > lambdamax)  //The upper wavelength is outside the requested range
    {
      linterp (lambdamax, comp[integ_spectype].xmod.w, comp[integ_spectype].xmod.f, comp[integ_spectype].nwaves, &f_interp, 1); //Obtain the interpolated flux at lambdamin
      x += (lambdamax - w1) * (comp[spectype].xmod.f[n + 1] + f_interp) / 2.0;  //Add the mean flux in the bin x the wavelength range
    }
    else
    {
      continue;
    }
  }

  x *= 4. * PI;

  return (x);
}



/**********************************************************/
/** 
 * @brief      Compute f_lambda from a model for a given wavelength
 *
 * @param [in] double  ;ambda - the wavelength of interest
 * @param [in] void  params   An extra (unused) variable to make it paletable for the gsl integrator
 * @return     The flux for a model for a given wavelength
 *
 * @details
 * This is the integrand for the Romburg integration to obtain
 * the surface flux of a star from a continuum model. It is simply
 * a case of interpolating on the model, or returning zero if we
 * request a wavelength outside the model.
 *
 * ### Notes ###
 *
 * The number of the model we are using is set as an external variable.
 *
 *
 **********************************************************/


double
model_int (double lambda, void *params)
{
  double answer;                //The interpolated answer
  if (lambda < comp[integ_spectype].xmod.w[0])  //Our wavelength is below where we have a model
    answer = 0.0;               //retuen zero
  else if (lambda > comp[integ_spectype].xmod.w[comp[integ_spectype].nwaves - 1])       //Our wavelength is above where we have a model
    answer = 0.0;               //return zero
  else
  {
    linterp (lambda, comp[integ_spectype].xmod.w, comp[integ_spectype].xmod.f, comp[integ_spectype].nwaves, &answer, 1);        //Interpolate (in log space)
  }
  return (answer);
}
