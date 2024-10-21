
/***********************************************************/
/** @file  agn.c
 * @author nsh
 * @date   2011
 *
 * @brief  Subroutines relating to power law radiation sources
 *
 * The subroutines and functions in this file are to do with
 * power law SEDs, which are mainly used for AGN, hence the
 * name of the file.
 ***********************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "atomic.h"
#include "sirocco.h"
#include "log.h"


/**********************************************************/
/** 
 * @brief      Calculates the total luminosity of an AGN type source
 *
 * @param [in] double  r   radius of emitting object
 * @param [in] double  lum   the luminosity of the AGN (2-10keV) - used in the continuum model case 
 * @param [in] double  alpha   the spectral index of the PL source - also sometimes used for a temperature if a BB source is required
 * @param [in] double  freqmin   minimum frequency to integrate over
 * @param [in] double  freqmax   maximum frequency
 * @param [in] int  ioniz_or_extract   flag to say if we are in the ionization cycles or spectral cycle
 * @param [out] double   f    the returned luminosity  
 * @return     f - the luminosity - seems to be returned as well as set implicitly.
 *
 * @details
 * This routine is used when one is working out how many photons will be made
 * from each of several sources. This computes the total luminosity of a source
 * of radius r, with spectral index alpha between band boundaries given by freqmin 
 * and freqmax
 * 
 *
 * ### Notes ###
 * The 2-10keV luminosity of the PL source) is only used in the emittance_continuum
 * whwere it is used to scale the model luminosity to obtain the required value.
 * A slight inconsistency is that the bremstrahlung code uses values for T and alpha that
 * are stored in the geo structure, whilst the power law uses data supplied in the code.
 *
 **********************************************************/

double
agn_init (r, lum, alpha, freqmin, freqmax, ioniz_or_extract, f)
     double r, lum, alpha, freqmin, freqmax;
     int ioniz_or_extract;
     double *f;
{

  double t;
  double scaling;               //The scaling factor to get from a model to get the correct 2-10keV luminosity
  double emit, emit_2_10;
  int spectype;

  if (ioniz_or_extract == CYCLE_EXTRACT)
    spectype = geo.agn_spectype;        /* type for final spectrum */
  else
    spectype = geo.agn_ion_spectype;    /*type for ionization calculation */
  if (spectype >= 0)
  {
    /* This calls models - we *should* have only one model in this case, so we call the emittance
       continuum with dummy variables instead of T and g (the trailing 0.0s). */
    /* First we compute the emittance from 2-10keV to compare with the required 2-10keV luminosity */
    emit_2_10 = emittance_continuum (spectype, 4.84e17, 2.42e18, 0.0, 0.0);
    if (emit_2_10 <= 0.0)
    {
      Error ("agn_init: Supplied model SED has no power in 2-10keV band - cannot continue\n");
      exit (0);
    }
    scaling = lum / emit_2_10;  //This gives us a scaling factor which can be applied to whatever band we are dealing with


    emit = emittance_continuum (spectype, freqmin, freqmax, 0.0, 0.0);  //Compute emittance (luminosity) for the actual band
    *f = emit * scaling;        //And scale it
  }
  else if (spectype == SPECTYPE_POW)    //Power law - uses constant computed elsewhere and spectral index
  {
    /* Emittance_pow actucally returns the specific luminosity directly */
    emit = emittance_pow (freqmin, freqmax, alpha);
    *f = emit;
  }
  else if (spectype == SPECTYPE_BB)     //Blackbody - uses the temperature which is stored in alpha...
  {

    t = alpha;                  /* For this case the second variable is t */
    /* Note that emittance_bb is really the energy radiated per unit area  so
     * this must be multiplied by the area of the object */
    emit = emittance_bb (freqmin, freqmax, t);
    *f = emit;
    *f *= (4. * PI * r * r);
  }
  else if (spectype == SPECTYPE_CL_TAB) //A special broken power law mode made to match cloudy - mainly for testing purposes
  {
    /* Emittance_pow actually returns the specific luminosity directly */
    emit = emittance_bpow (freqmin, freqmax, alpha);
    *f = emit;
  }
  else if (spectype == SPECTYPE_BREM)   //Bremstrahlung - uses T and alpha which are stored in geo. 
  {
//    emit = qromb (integ_brem, freqmin, freqmax, 1e-4);
    emit = num_int (integ_brem, freqmin, freqmax, 1e-4);

    *f = emit;
  }

  return (*f);                  /* Return the luminosity    */
}


/* This routine returns the specific luminosity, that is the luminosity within the frequency interval */

/**********************************************************/
/** 
 * @brief      The luminosity of a power law over a frerquency range
 *
 * @param [in] double  freqmin   lower frequency bound
 * @param [in] double  freqmax   upper frequency bound
 * @param [in] double  alpha  spectral index
 * @return   double emit   the luminosoty
 *
 * @details
 * A Simple routine which just integrates a power law between two limits. There is
 * a little bit of code that deals with the case where an advanced mode is used
 * that cuts off the low frequency tail - and a bit of clever work to check for
 * the pathalogical case of alpha=-1 where one cannot symbolically integrate.
 *
 * ### Notes ###
 * 
 *
 **********************************************************/

double
emittance_pow (freqmin, freqmax, alpha)
     double freqmin, freqmax, alpha;
{
  double emit, this_fmin;

  /* if we have a PL cutoff, then we need to either adjust the minimum frequency,
     or return 0 luminosity if the cutoff is above our band */
  /* in advanced mode, this should always be zero */
  this_fmin = freqmin;          // default is no cutoff
  emit = 0.0;
  if (freqmax < geo.pl_low_cutoff)      //The maximum frequency requested is *lower* than the cut off
  {
    return (emit);              //Return zero
  }
  else if (freqmin < geo.pl_low_cutoff) //The minimum frequency is lower than the cxut off
  {
    this_fmin = geo.pl_low_cutoff;      //Reset the min frequency to the cutoff
  }

  /* conservative error check */
  if ((modes.iadvanced == 0) && (geo.pl_low_cutoff > 0))
    Error ("PL cutoff frequency is non-zero out of advanced mode!");

  /* we need to work out the luminosity between our limited frequency range */
  /* we may need some checking routines to make sure that the requested frequency range is within the defined range,
     or it could default to zero outside the defined range */

  if (alpha == -1.0)            //deal with the pathological case
  {
    emit = geo.const_agn * (log (freqmax) - log (this_fmin));
  }
  else
  {
    emit = geo.const_agn * ((pow (freqmax, alpha + 1.0) - pow (this_fmin, alpha + 1.0)) / (alpha + 1.0));
  }


  return (emit);
}


/**********************************************************/
/** 
 * @brief      Works out the luminosity of a broken power law
 * 
 * @param [in, out] double  freqmin   min frequency that we are generating photons over
 * @param [in, out] double  freqmax   max frequency that we are generating photons over
 * @param [in, out] double  lum   2-10kev luminosity
 * @param [in, out] double  alpha   is the 2-10kev slope.
 * @return    emit  The luminosity of the whole power law.
 *
 * @details
 *computes the emittance of a broken power law coded to try to replicate the cloudy table type power law. 
 * There are three sections - a central bit, which should cover the 2-10kev portion of the spectrum, 
 * a low frequency tail and a high frequency tail .
 *
 *
 * ### Notes ###
 * 
 *
 **********************************************************/

double
emittance_bpow (freqmin, freqmax, alpha)
     double freqmin, freqmax, alpha;
{
  double constant_low, constant_hi, emit;
  double e1, e2, e3;
  double atemp, ctemp, f1, f2;
  double pl_low, pl_hi;         //The low and high frequency breaks in the power law spectrum

#define   XFREQMIN  4.84e17
#define   XFREQMAX  2.42e18
  f1 = freqmin;                 /* NSH - 130506 added to reomve 03 compile errors */
  f2 = freqmax;                 /* NSH - 130506 added to reomve 03 compile errors */
  e1 = e2 = e3 = 0.0;           /* NSH - 130506 added to reomve 03 compile errors */
  /* first we need to calculate the constant for the 2-10 kev power law function */

  /* convert broken power law bands to freq */

  pl_low = geo.agn_cltab_low / HEV;
  pl_hi = geo.agn_cltab_hi / HEV;


  constant_low = xband.pl_const[0];     //we have already worked out the constants in bands.c
  constant_hi = xband.pl_const[xband.nbands - 1];


  /* now we need to work out the luminosity between our limited frequency range */
  emit = 0.0;


  /* first deal with energy in low frequency part of the spectrum */
  if (freqmin >= pl_low)
  {
    e1 = 0.0;                   //Our frequencies start above the cutoff - so there is no energy in the low frequency band
    Log ("Broken power law: Lowest frequency is above low frequency break\n");
  }
  else if (freqmax < pl_low)    //wierd situation, where all power is in the low frequency tail
  {
    atemp = geo.agn_cltab_low_alpha;
    f1 = freqmin;
    f2 = freqmax;
    ctemp = constant_low;
    e1 = ctemp * ((pow (f2, atemp + 1.0) - pow (f1, atemp + 1.0)) / (atemp + 1.0));
    Log ("Broken power law: Emittance below low frequency break is %e\n", e1);
  }
  else                          //some of the power is in the low frequency bit, so we need to integrate from fmin up to the low frequency boundary 
  {
    atemp = geo.agn_cltab_low_alpha;
    f1 = freqmin;
    f2 = pl_low;
    ctemp = constant_low;
    e1 = ctemp * ((pow (f2, atemp + 1.0) - pow (f1, atemp + 1.0)) / (atemp + 1.0));
    Log ("Broken power law: Emittance below low frequency break is %e\n", e1);
  }

  /* energy in main frequency band */
  if (freqmax > pl_low && freqmin < pl_hi)      //If there is any part of the required range between our central range
  {
    if (freqmin <= pl_low && freqmax >= pl_hi)
    {
      f1 = pl_low;
      f2 = pl_hi;
    }
    else if (freqmin <= pl_low && freqmax < pl_hi)
    {
      f1 = pl_low;
      f2 = freqmax;
    }
    else if (freqmin > pl_low && freqmax < pl_hi)
    {
      f1 = freqmin;
      f2 = freqmax;
    }
    else if (freqmin > pl_low && freqmax >= pl_hi)
    {
      f1 = freqmin;
      f2 = pl_hi;
    }

    e2 = emittance_pow (f1, f2, alpha);

    Log ("Broken power law: Emittance in centre is %e\n", e2);
  }
  else
  {
    Log ("Broken power law: No luminosity in central range\n");
  }

  /* energy in high frequency tail */
  if (freqmax <= pl_hi)
  {
    e3 = 0.0;                   //Our frequencies stop below the cutoff
    Log ("Broken power law: Highest frequency is below high frequency break\n");
  }
  else if (freqmin > pl_hi)     //Odd situation where all power is in high frequency tail
  {
    atemp = geo.agn_cltab_hi_alpha;
    f1 = freqmin;
    f2 = freqmax;
    ctemp = constant_hi;
    e3 = ctemp * ((pow (f2, atemp + 1.0) - pow (f1, atemp + 1.0)) / (atemp + 1.0));
    Log ("Broken power law: Emittance above high frequency break is %e\n", e3);
  }
  else                          //normal situation where we run from boundary to the max frequency
  {
    atemp = geo.agn_cltab_hi_alpha;
    f1 = pl_hi;
    f2 = freqmax;
    ctemp = constant_hi;
    e3 = ctemp * ((pow (f2, atemp + 1.0) - pow (f1, atemp + 1.0)) / (atemp + 1.0));
    Log ("Broken power law: Emittance above high frequency break is %e\n", e3);
  }

  emit = e1 + e2 + e3;

  return (emit);
}



/**********************************************************/
/** 
 * @brief      Generate a photon from an AGN object
 *
 * @param [in, out] PhotPtr  p   pointer to the photon array
 * @param [in, out] double  r   radius of the central source
 * @param [in, out] double  alpha   spectral index of the PL sourece - also used to communicate temperature if a BB or model is used.
 * @param [in, out] double  weight   the weight of each photon to be made
 * @param [in, out] double  f1   lower frequency bound
 * @param [in, out] double  f2   upper frequency nound
 * @param [in, out] int  spectype   the type of SED
 * @param [in, out] int  istart   the number of the first photon to be made - index into the p array
 * @param [in, out] int  nphot   the number of photons to be made
 * @return     0 (if successful or unscuccessful)
 *
 * @details
 * This subroutine produces photons from an AGN type source. There are a range of 
 * different SEDs that can be used, these are communicated by the spectype and different
 * choices cause different surbroutines to be called. The only *real* work done here
 * is to compute the strarting location and direction of the photons.
 *
 * ### Notes ###
 *
 **********************************************************/

int
photo_gen_agn (p, r, alpha, weight, f1, f2, spectype, istart, nphot)
     PhotPtr p;
     double r, alpha, weight;
     double f1, f2;             /* The freqency mininimum and maximum if a uniform distribution is selected */
     int spectype;              /*The spectrum type to generate: 0 is bb, 1 (or in fact anything but 0)
                                   is uniform in frequency space */
     int istart, nphot;         /* Respecitively the starting point in p and the number of photons to generate */
{
  double freqmin, freqmax, t;
  int i, iend;
  int n;
  double ftest;

  t = alpha;                    /* NSH 130605, this is a slightly odd statmenent, put in to avoid 03 compilation 
                                   errors, but stemming from the odd behaviour that one can give the AGN a 
                                   themal spectral type, and the alpha parameter is then used to transmit 
                                   the temperature. */


  if ((iend = istart + nphot) > NPHOT)  //Consistency check - if it fails then we are being asked to make more photons than expected
  {
    Error ("photo_gen_agn: iend %d > NPHOT %d\n", iend, NPHOT);
    Exit (0);
  }
  if (f2 < f1)                  //Another consistency check - it is not sensible to have the upper frequency lower than the lower frequency
  {
    Error ("photo_gen_agn: Cannot generate photons if freqmax %g < freqmin %g\n", f2, f1);
  }
  Log_silent ("photo_gen_agn creates nphot %5d photons from %5d to %5d \n", nphot, istart, iend);
  freqmin = f1;
  freqmax = f2;

  r = (1. + EPSILON) * r;       /* Generate photons just outside the photosphere unnecessary for the AGN perhaps? */

  /* Generate photons just outside the photosphere unnecessary for the AGN perhaps? */
  /* note this is only used in the spherical source */

  if (spectype == SPECTYPE_CL_TAB)      /* If we have a broken PL, we need to work out which alpha to use */
  {
    ftest = (freqmin + freqmax) / 2.0;
    for (n = 0; n < xband.nbands; n++)  /* work out what alpha and constant to use */
    {
      if (xband.f1[n] < ftest && xband.f2[n] > ftest)
      {
        alpha = xband.alpha[n];
      }
    }
  }




  for (i = istart; i < iend; i++)       //Loop over the number of photons we are asked to make
  {
    p[i].origin = PTYPE_AGN;    // For BL photons this is corrected in photon_gen 
    p[i].w = weight;            //Set the weight
    p[i].istat = p[i].nscat = p[i].nrscat = p[i].nmacro = 0;    //Initialise status, number of scatters and number of resonant scatters
    p[i].grid = 0;              //Set the grid number to zero 
    p[i].tau = 0.0;             //Set the opacity seen by the photon to zero
    p[i].nres = -1;             // It's a continuum photon - so it is not made in a resonance
    p[i].nnscat = 1;            // Set to one scatter

    if (spectype == SPECTYPE_BB)        //Blackbody spectrum, we use the supplied temperature
    {
      p[i].freq = planck (t, freqmin, freqmax);
    }
    else if (spectype == SPECTYPE_UNIFORM)
    {
      /*Produce a uniform distribution of frequencies */
      p[i].freq = random_number (freqmin, freqmax);     //Just a random frequency 
    }
    else if (spectype == SPECTYPE_POW)  //Power law spectrum
    {
      p[i].freq = get_rand_pow (freqmin, freqmax, alpha);
    }
    else if (spectype == SPECTYPE_CL_TAB)       //Broken power law
    {
      p[i].freq = get_rand_pow (freqmin, freqmax, alpha);
    }
    else if (spectype == SPECTYPE_BREM) //Bremstrahlung SED - the bremstrahlung parameters (T and alpha) are hidden in the geo structure
    {
      p[i].freq = get_rand_brem (freqmin, freqmax);
    }
    else if (spectype == SPECTYPE_MONO)
    {
      p[i].w = 1. / geo.pcycles;
      p[i].freq = geo.mono_freq;
    }
    else
    {
      p[i].freq = one_continuum (spectype, -1., -1., freqmin, freqmax); //A continuum (model) photon - we use t=g=-1 to flag that this is not a normal model
    }


    if (p[i].freq < freqmin || freqmax < p[i].freq)     //A check to see that we havent made a photon out of range.
    {
      Error_silent ("photo_gen_agn: phot no. %d freq %g out of range %g %g\n", i, p[i].freq, freqmin, freqmax);
    }


/* Now we work out where the photon starts.
   * The first case is the specail searchlight case
   * which is purely diagnostic
 */

    if (geo.ioniz_or_extract == CYCLE_EXTRACT && modes.searchlight)
    {
      stuff_v (geo.searchlight_x, p[i].x);
      stuff_v (geo.searchlight_lmn, p[i].lmn);
    }
    /*  original spherical X-ray source as used in e.g. Higginbottom+ 2013 */

    else if (geo.pl_geometry == PL_GEOMETRY_SPHERE)
    {
      randvec (p[i].x, r);      //Simple random coordinate on the surface of a sphere

      /* Added by SS August 2004 for finite disk. */
      if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
      {
        while (fabs (p[i].x[2]) < zdisk (p[i].x[0]))    //We just need to make sure that the photon isn't submerged in the extneded disk
        {
          randvec (p[i].x, r);
        }


        if (fabs (p[i].x[2]) < zdisk (r))
        {
          Error ("Photon_agn: agn photon %d in disk %g %g %g %g %g\n", i, p[i].x[0], p[i].x[1], p[i].x[2], zdisk (r), r);
          Exit (0);
        }
      }
      randvcos (p[i].lmn, p[i].x);      //Random direction centred on the previously randmised vector
    }

    /* if we have a lamp post geometry then we should generate isotropic photons at a height
       above the disk plane */
    else if (geo.pl_geometry == PL_GEOMETRY_LAMP_POST)
    {
      /* x and y coordinates are 0 */
      p[i].x[0] = p[i].x[1] = 0.0;

      /* need to set the z coordinate to the lamp post height, but allow it to be above or below */
      if (random_number (-1.0, 1.0) > 0.0)

      {                         /* Then the photon emerges in the upper hemisphere */
        p[i].x[2] = geo.lamp_post_height;
      }
      else
      {
        p[i].x[2] = -geo.lamp_post_height;
      }

      randvec (p[i].lmn, 1.0);  // lamp-post geometry is isotropic, so completely random vector
    }

    /* if we have a bubble geometry then we should generate isotropic photons at a random position
       within a sphere. */
    else if (geo.pl_geometry == PL_GEOMETRY_BUBBLE)
    {
      /* This is esssentially copied from spherical_get_random_location */
      double rrr;

      rrr =
        (geo.rstar * geo.rstar * geo.rstar) + (geo.bubble_size * geo.bubble_size * geo.bubble_size -
                                               geo.rstar * geo.rstar * geo.rstar) * random_number (0.0, 1.0);


      rrr = pow (rrr, (1. / 3.));

      /* The direction can be anywhere on a sphere, so completely random directins for origin */
      randvec (p[i].x, rrr);
      randvec (p[i].lmn, 1.0);  // lamp-post geometry is isotropic, so completely random vector
    }

    /* This is set up for looking at photons in spectral cycles at present */
    //if (modes.save_photons && geo.ioniz_or_extract == CYCLE_EXTRACT)
    //  save_photons (&p[i], "AGN");
  }



  return (0);
}
