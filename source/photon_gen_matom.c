/***********************************************************/
/** @file  photo_gen_matom.c
 * @author ksl, ss, jm
 * @date   January, 2018
 *
 * @brief functions for calculating emissivities and generating photons from macro-atoms and k-packets.
 *   during the spectral cycles. The actual functions which do the jumps inside an activated
 *  macro-atom are in matom.c. This is partly done to prevent overly long files (JM1504)
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

/**********************************************************/
/**
 * @brief      returns the specific luminosity in the band needed for the computation of the
 *        spectrum. It gets the total energy radiated by the process k-packet -> r-packet in the
 *        required wavelength range.
 *
 * @return double lum  The energy radiated by the process k-packet -> r-packet in the wind in the
 *            wavelength range required for the specrum calculation.
 *
 * ### Notes ###
 * Literally just loops through the kpkt_emiss entries in plasmamain.
 *
 **********************************************************/

double
get_kpkt_f ()
{
  int n;
  double lum;

  lum = 0.0;

  for (n = 0; n < NPLASMA; n++)
  {
    lum += plasmamain[n].kpkt_emiss;
  }


  return (lum);
}

/**********************************************************/
/**
 * @brief returns the specific luminosity in kpkts from nonthermal ("shock")
 *        heating. This is used to generate kpkts in the ionization cycles.
 *        This also populates the cell-by-cell kpkt luminosities in the
 *        variable plasmamain[n].kpkt_emiss.
 *
 * @return double lum
 *         The energy created by non-radiative heating throughout the
 *         computational domain.
 *
 **********************************************************/

double
get_kpkt_heating_f ()
{
  int n, nwind;
  double lum, shock_kpkt_luminosity;
  WindPtr one;

  lum = 0.0;

  for (n = 0; n < NPLASMA; n++)
  {
    nwind = plasmamain[n].nwind;
    one = &wmain[nwind];

    /* what we do depends on how the "net heating mode" is defined */
    if (KPKT_NET_HEAT_MODE)
      shock_kpkt_luminosity = (shock_heating (one) - plasmamain[n].cool_adiabatic);
    else
      shock_kpkt_luminosity = shock_heating (one);

    if (shock_kpkt_luminosity > 0)
    {
      if (geo.ioniz_or_extract == CYCLE_IONIZ)
        plasmamain[n].kpkt_emiss = shock_kpkt_luminosity;
      else
        plasmamain[n].kpkt_abs += shock_kpkt_luminosity;

      lum += shock_kpkt_luminosity;
    }
    else
      plasmamain[n].kpkt_emiss = 0.0;
  }

  return (lum);
}



/**********************************************************/
/**
 * @brief produces photon packets to account for creating of r-packets by k-packets.

 *
 * @param [in, out] PhotPtr  p   the ptr to the entire structure for the photons
 * @param [in] double  weight   the photon weight
 * @param [in] int  photstart   The first element of the photon stucure to be populated
 * @param [in] int  nphot   the number of photons to be generated
 * @return int nphot  The number of photon packets that were generated from k-packet elliminations.
 *
 * @details produces photon packets to account for creating of r-packets by k-packets in the spectrum calculation.
 * It should only be used once the total energy emitted in this way in the wavelength range in question is well known
 * (calculated in the ionization cycles). This routine is closely related to photo_gen_wind from which much of the code
 * has been copied.
 *
 * Photons are generated at a position in the Observer frame.
 * The weight of the photon should is the weight expected in the local
 * frame since photon is first created in thea local
 * rest frame, and then Doppler shifted to the Observer frame
 **********************************************************/

int
photo_gen_kpkt (p, weight, photstart, nphot)
     PhotPtr p;
     double weight;
     int photstart, nphot;
{
  int photstop;
  int icell;
  double xlum, xlumsum;
  struct photon pp;
  int nres, esc_ptr, which_out;
  int n;
  double test;
  int nnscat;
//OLD  int nplasma, ndom;
  int nplasma;
  int kpkt_mode;
  double freq_min, freq_max;

  photstop = photstart + nphot;
  Log ("photo_gen_kpkt  creates nphot %5d photons from %5d to %5d, weight %8.4e \n", nphot, photstart, photstop, weight);

  if (geo.ioniz_or_extract == CYCLE_IONIZ)
  {
    /* we are in the ionization cycles, so use all frequencies. kpkt_mode should allow all processes */
    freq_min = xband.f1[0];
    freq_max = xband.f2[xband.nbands - 1];
    kpkt_mode = KPKT_MODE_ALL;
  }
  else
  {
    /* we are in the spectral cycles, so use all the required frequency range */
    freq_min = geo.sfmin;
    freq_max = geo.sfmax;
    /* we only want k->r processes */
    kpkt_mode = KPKT_MODE_CONTINUUM;
  }

  for (n = photstart; n < photstop; n++)
  {
    /* locate the wind_cell in which the photon bundle originates. */

    xlum = random_number (0.0, 1.0) * geo.f_kpkt;

    xlumsum = 0;
    icell = 0;
    while (xlumsum < xlum)
    {
      if (wmain[icell].inwind >= 0)
      {
        nplasma = wmain[icell].nplasma;
        xlumsum += plasmamain[nplasma].kpkt_emiss;
      }
      icell++;
    }
    icell--;                    /* This is the cell in which the photon must be generated */

    /* Now generate a single photon in this cell */
    p[n].w = weight;

    pp.freq = 0.0;
    pp.grid = icell;

    /* This following block is a bad way of doing it - kpkt could be modified to
       do what we want in a more elegant way. */

    test = pp.freq;

    while (test > freq_max || test < freq_min)
    {
      pp.w = p[n].w;
      kpkt (&pp, &nres, &esc_ptr, kpkt_mode);

      if (esc_ptr == 0 && kpkt_mode == KPKT_MODE_CONTINUUM)
      {
        Error ("Didn't escape from kpkt despite being in continuum mode\n");
      }
      else
      {
        if (esc_ptr == 0)
        {
          macro_gov (&pp, &nres, 1, &which_out);
        }
        test = pp.freq;
      }
    }

    p[n].freq = pp.freq;
    p[n].nres = p[n].line_res = nres;
    p[n].w = pp.w;


    /* The photon frequency is now known. */

    /* Determine the position of the photon in the moving frame */

    get_random_location (icell, p[n].x);

    p[n].grid = icell;
    p[n].frame = F_LOCAL;

    nnscat = 1;
    // Determine the direction of the photon
    // Need to allow for anisotropic emission here
    if (p[n].nres < 0 || p[n].nres > NLINES || geo.scatter_mode == SCATTER_MODE_ISOTROPIC)
    {
      /*  It was either an electron scatter so the  distribution is isotropic, or it
         was a resonant scatter but we want isotropic scattering anyway.  */
      randvec (p[n].lmn, 1.0);  /* The photon is emitted isotropically */
    }
    else if (geo.scatter_mode == SCATTER_MODE_THERMAL)
    {                           //It was a line photon and we want the thermal trapping anisotropic model

      randwind_thermal_trapping (&p[n], &nnscat);
    }

    p[n].nnscat = nnscat;


//OLD    ndom = wmain[icell].ndom;

    /* Make an in-place transformation to the observer frame */
    if (local_to_observer_frame (&p[n], &p[n]))
    {
      Error ("photo_gen_kpkt:Frame transformation error\n");
    }

    p[n].istat = 0;
    p[n].tau = p[n].nscat = p[n].nrscat = p[n].nmacro = 0;
    p[n].origin = PTYPE_WIND;   // Call it a wind photon

    switch (geo.reverb)
    {                           //0715 SWM - Added path generation
    case REV_MATOM:
      line_paths_gen_phot (&wmain[icell], &p[n], nres);
      break;
    case REV_WIND:
      wind_paths_gen_phot (&wmain[icell], &p[n]);
      break;
    case REV_PHOTON:
      simple_paths_gen_phot (&p[n]);
      break;
    case REV_NONE:
    default:
      break;
    }

  }



  return (nphot);               /* Return the number of photons generated */


  /* All done. */
}

/**********************************************************/
/**
 * @brief      produces photon packets to account for creation of r-packets
 *      by deactivation of macro atoms in the spectrum calculation. It should only be used
 *      once the total energy emitted in this way in the wavelength range in question is well known
 *      (calculated in the ionization cycles).
 *
 * @param [in, out] PhotPtr  p   the ptr to the structire for the photons
 * @param [in] double  weight   the photon weight
 * @param [in] int  photstart   The position in the photon structure of the first photon to generate
 * @param [in] int  nphot   the number of the photons to be generated
 * @return int nphot When it finishes it should have generated nphot photons from macro atom deactivations.
 *
 * @details
 * This routine is closely related to photo_gen_kpkt from which much of the code has been copied.
 *
 * ### Notes ###
 * Consult Matthews thesis.
 *
 **********************************************************/

int
photo_gen_matom (p, weight, photstart, nphot)
     PhotPtr p;
     double weight;
     int photstart, nphot;
{
  int photstop;
  int icell;
  double xlum, xlumsum;
  struct photon pp;
  int nres;
  int n;
  double dot ();
//OLD  double test;
  int upper;
  int nnscat;
  int nplasma;
//OLD  int ndom;



  photstop = photstart + nphot;
  Log ("photo_gen_matom creates nphot %5d photons from %5d to %5d \n", nphot, photstart, photstop);

  for (n = photstart; n < photstop; n++)
  {
    /* locate the wind_cell in which the photon bundle originates. And also decide which of the macro
       atom levels will be sampled (identify that level as "upper"). */
    xlum = random_number (0.0, 1.0) * geo.f_matom;


    xlumsum = 0;
    icell = 0;
    upper = 0;
    while (xlumsum < xlum)
    {
      if (wmain[icell].inwind >= 0)
      {
        nplasma = wmain[icell].nplasma;
        xlumsum += macromain[nplasma].matom_emiss[upper];
        upper++;
        if (upper == nlevels_macro)
        {
          upper = 0;
          icell++;
        }
      }
      else
      {
        icell++;
        upper = 0;
      }
    }

    if (upper == 0)
    {
      /* If upper = 0 at this point it means the loop above finished with an increment of
         icell and setting upper = 0. In such a case the process we want was the LAST process
         in the previous cell. */
      icell--;
      upper = nlevels_macro;
    }

    upper--;                    /* This leaves the macro atom level that deactivaties. */

    /* Now generate a single photon in this cell */
    p[n].w = weight;

    pp.freq = 0.0;
    pp.grid = icell;

    /* This following block is a bad way of doing it but it'll do as
       a quick and dirty test for now. Really kpkt should be modified to
       do what we want in a more elegant way. */

//OLD    test = pp.freq;

//OLD    while (test > geo.sfmax || test < geo.sfmin)
//OLD    {

    /* Call routine that will select an emission process and a frequencyfor the
       deactivating macro atom. */
//OLD    If it deactivates outside the frequency
//OLD         range of interest then ignore it and try again. SS June 04. */

    emit_matom (wmain, &pp, &nres, upper, geo.sfmin, geo.sfmax);

//OLD      test = pp.freq;
//OLD    }


    p[n].freq = pp.freq;
    p[n].nres = p[n].line_res = nres;
    p[n].frame = F_LOCAL;



    /* The photon frequency is now known. */

    /* Determine the position of the photon in the moving frame */


    get_random_location (icell, p[n].x);

    p[n].grid = icell;


    /* Determine the direction of the photon
       Need to allow for anisotropic emission here
     */

    nnscat = 1;
    if (p[n].nres < 0 || p[n].nres > NLINES || geo.scatter_mode == SCATTER_MODE_ISOTROPIC)
    {
      /*  It was either an electron scatter so the  distribution is isotropic, or it
         was a resonant scatter but we want isotropic scattering anyway.  */
      randvec (p[n].lmn, 1.0);  /* The photon is emitted isotropically */
    }
    else if (geo.scatter_mode == SCATTER_MODE_THERMAL)
    {                           //It was a line photon and we want the thermal trapping anisotropic model

      randwind_thermal_trapping (&p[n], &nnscat);
    }
    p[n].nnscat = nnscat;


    /* The next two lines correct the frequency to first order, but do not result in
       forward scattering of the distribution */

//OLD    ndom = wmain[icell].ndom;


    /* Make an in-place transformation to the observer frame */
    if (local_to_observer_frame (&p[n], &p[n]))
    {
      Error ("photo_gen_matom:Frame transformation error\n");
    }

    p[n].istat = 0;
    p[n].tau = p[n].nscat = p[n].nrscat = p[n].nmacro = 0;
    p[n].origin = PTYPE_WIND;   // Call it a wind photon

    switch (geo.reverb)
    {                           //0715 SWM - Added path generation
    case REV_MATOM:
      line_paths_gen_phot (&wmain[icell], &p[n], nres);
      break;
    case REV_WIND:
      wind_paths_gen_phot (&wmain[icell], &p[n]);
      break;
    case REV_PHOTON:
      simple_paths_gen_phot (&p[n]);
      break;
    case REV_NONE:
    default:
      break;
    }
  }


  return (nphot);               /* Return the number of photons generated */


  /* All done. */
}
