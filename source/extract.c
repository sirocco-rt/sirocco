
/***********************************************************/
/** @file  extract.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  Routines for extracting photons during the
 * spectral generation phase.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief      The main routine for building
 * 	detailed spectra in the normal (extract) mode.
 *
 * @param [in] WindPtr  w   The entire wind
 * @param [in] PhotPtr  p   The photon to extract
 * @param [in] int  itype   An integer representing the type of photon
 * for the purpose of being extracted.
 * @return     Always returns 0
 *
 * @details
 * extract is called when a photon begins its flight and every time that photon
 * scatters, unless the user has exercised the "live or die" option, in
 * which case it is not called.
 *
 * extract carries out all of the preparatory steps to create the
 * photon to be extracted, and then calls a routine extract_one
 * to reduce the weight of the photon as it passes out of the system
 *
 * itype takes on the following values:
 * * PTYPE_STAR->the photon came for the star
 * * PTYPE_BL->the photon came from the boundary layer
 * * PTYPE_DISK->the photon being redirected arose in the disk,
 * * PTYPE_WIND->the photon being redirected arose in the wind,
 *
 * extract uses the types to prepare the photon for extraction, including
 * doppler shifting the photon if is of PTYPE_WIND or PTYPE_DISK.
 *
 * ### Notes ###
 *
 * The logic behind the weighting of the (disk/star) photons is described 
 * in Christian Knigge's thesis in
 * section 2.3.4.  According to equation 2.19
 * 	Pc/Pw=12 cos(theta)*(1+b cos(theta)/(3+2b) where b=1.5 corresponds to the
 * Eddington approximation.
 *
 * Usually, Python constructs a spectrum of all photons, but there are
 * advanced options which allone to restrict the spectrum created to
 * those produced with a certain number of scatters or from photons
 * that arise from the above or below the disk.  extract enforces
 * those choices before calling extract_one.
 * The parameters for this option all come in through
 * sirocco.h, and are contained in the spectrum structure.
 *
 * The basic choices, some of which can be used in tandom are as follows:
 *
 * *	If s[n].nscat>999; then one obtains the entire spectrum
 * *	if s[n].nscat is a positive number, then the spectrum is composed just
 * 		of photons with that number of scatters
 * *	if s[n].nscat is a negative number, then the spectrum contains all photons
 * 		with >=  the absolute value of s[n].nscat
 * * 	s[n].top_bot=0 -> accept all photons
 * * 	s[n].top_bot    >0  -> accept photons last above the disk
 * * 	s[n].top_bot    <0  -> accept photons below the disk
 * *    One can also select photons that are to be extracted from a particular
 * spherical region.
 *
 *
 **********************************************************/



int
extract (w, p, itype)
     WindPtr w;
     PhotPtr p;
     int itype;
{
  int n, mscat, mtopbot;
  struct photon pp, p_in, p_dummy;
  double xdiff[3];
  double p_norm, tau_norm;
  double dvds_max;
  int ierr;
  double x[3];
  double tau;
  double zz;
  double dvds;
  int ishell;
  double vel[3];
  double weight_scale;
  double w_orig;

  tau = 0.0;



  ierr = check_frame (p, F_OBSERVER, "extract_start");
  if (ierr)
  {
    Error ("extract: check_frame failure at very start for itype %d\n", itype);
  }

  stuff_phot (p, &p_in);

  if (itype == PTYPE_WIND)
  {
    if (geo.scatter_mode == SCATTER_MODE_THERMAL && p_in.nres <= NLINES && p_in.nres > -1)
    {
      /* Bound - Bound Scattering 

       * We increase weight to account for the number of scatters. This is
       * done because in extract we multiply by the escape probability
       * along a given direction, but we also need to divide the weight by
       * the mean escape probability, which is equal to 1/nnscat.  See
       * issue #710 for a more extended explanation of how the weight is
       * renormalized stochastically.
       *
       * We normalised our rejection method by the escape
       * probability along the vector of maximum velocity
       * gradient. First find the sobolev optical depth
       * along that vector. The -1 enforces calculation of
       * the ion density
       */

      dvds_max = get_dvds_max (&p_in);
      tau_norm = sobolev (&wmain[p_in.grid], p_in.x, -1.0, lin_ptr[p_in.nres], dvds_max);
      p_norm = p_escape_from_tau (tau_norm);
      p_in.w *= p_in.nnscat / p_norm;
      if ((ierr = observer_to_local_frame (&p_in, &p_in)))
        Error ("extract: wind photon not in observer frame %d\n", ierr);

    }

    else if (p_in.nres == NRES_ES)
    {
      /* Compton Scattering */

      if ((ierr = observer_to_local_frame (&p_in, &p_in)))
        Error ("extract: wind photon not in observer frame %d\n", ierr);

      lorentz_transform (&p_in, &p_in, velocity_electron);
      rescale (velocity_electron, -1, vel);     // Only need to do this once
    }

    else
    {

      if (p_in.nnscat != 1)
        Error
          ("trans_phot: nnscat is %i for photon %i in scatter mode %i! nres %i NLINES %i\n",
           p_in.nnscat, p_in.np, geo.scatter_mode, p_in.nres, NLINES);
      if ((ierr = observer_to_local_frame (&p_in, &p_in)))
        Error ("extract: wind photon not in observer frame %d\n", ierr);
    }



  }
  else if (itype == PTYPE_DISK)
  {
    if ((ierr = observer_to_local_frame_disk (&p_in, &p_in)))
      Error ("extract: disk photon not in observer frame %d\n", ierr);
  }

  /*
   * At this point were are in a local frame for WIND and DISK photons,
   * but the global frame for the central source and boundary layer
   */

  w_orig = p_in.w;
  for (n = MSPEC; n < nspectra; n++)
  {
    /*
     * If statement allows one to choose whether to construct the
     * spectrum from all photons or just from photons that have
     * scattered a specific number of times or in specific
     * regions of the wind. A region is specified by a position
     * and a radius.
     */


    if ((mscat = xxspec[n].nscat) >= MAXSCAT || p_in.nscat == mscat || (mscat < 0 && p_in.nscat >= (-mscat)))
    {
    }
    else
      continue;

    if ((mtopbot = xxspec[n].top_bot) == 0)
    {
    }
    //Then there are no positional parameters and we are done
    else if (mtopbot == -1 && p_in.x[2] < 0)
    {
    }
    else if (mtopbot == 1 && p_in.x[2] > 0)
    {
    }
    else if (mtopbot == 2)
      //Then to count, the photom must originate within sn.r of sn.x
    {
      vsub (p_in.x, xxspec[n].x, xdiff);
      if (length (xdiff) > xxspec[n].r)
        continue;
    }
    else
      continue;

    /*
     * Create a photon pp to use here and in extract_one,
     * and send it in the correct direction.  This
     * assures we have not modified p_in as part of
     * extract.  Also, allow for aberration of photons to
     * assure that we are extracting at the correct angle
     * in the observer frame
     * 
     * Note that a wind photon is in the local frame, and
     * the only thing we need to do is to figure out the
     * extraction direction in the local frame that will
     * produce the directon we want
     * 
     * A disk photon is in the observer frame
     */


    /* At this point we need to calculate the desired direction of the photon in
       the local frame of the disk or wind.  Additionally, we have to 
       mske a correction to the weights of these photons, due to
       an an effect of special relativty which accounts for solid angle
       corrections.

       PTYPE_STAR, PTYPE_BL do not need to be modified.  The rather
       peculiar rel_mode != REL_MODE_LINEAR arises because there
       are two modes that REL_MODE_SR_FREQ and  REL_MODE_FULL 
       that need frame transformations.

     */

    if (itype == PTYPE_WIND && rel_mode != REL_MODE_LINEAR)
    {
      stuff_phot (&p_in, &p_dummy);
      p_dummy.frame = F_OBSERVER;
      stuff_v (xxspec[n].lmn, p_dummy.lmn);

      observer_to_local_frame (&p_dummy, &p_dummy);

      weight_scale = p_dummy.w / w_orig;

      stuff_phot (&p_in, &pp);

      pp.w = w_orig / weight_scale / weight_scale;
      stuff_v (p_dummy.lmn, pp.lmn);
    }
    else if (itype == PTYPE_DISK && rel_mode != REL_MODE_LINEAR)
    {
      stuff_phot (&p_in, &p_dummy);
      p_dummy.frame = F_OBSERVER;
      stuff_v (xxspec[n].lmn, p_dummy.lmn);

      observer_to_local_frame_disk (&p_dummy, &p_dummy);

      weight_scale = p_dummy.w / w_orig;

      stuff_phot (&p_in, &pp);

      pp.w = w_orig / weight_scale / weight_scale;
      stuff_v (p_dummy.lmn, pp.lmn);
    }
    else
    {
      stuff_phot (&p_in, &pp);
      stuff_v (xxspec[n].lmn, pp.lmn);
    }



    /* At this stage, we are in the local frame for
       photons which are from the wind or the star.

       pp is the photon we are going to extract


       * Re-weight the photons. Note that photons
       * have already been frequency shifted prior
       * to entering extract. For disk and central
       * object photons, see Eqn 2.19 Knigge's
       * thesis
     */

    if (itype == PTYPE_STAR || itype == PTYPE_BL || (itype == PTYPE_AGN && geo.pl_geometry == PL_GEOMETRY_SPHERE))
    {
      stuff_v (pp.x, x);
      renorm (x, 1.);
      zz = fabs (dot (x, pp.lmn));
      pp.w *= zz * (2.0 + 3.0 * zz);
    }
    else if (itype == PTYPE_DISK)
    {

      zz = fabs (pp.lmn[2]);
      pp.w *= zz * (2.0 + 3.0 * zz);

    }
    else if (pp.nres == NRES_ES)
    {
      /* Reweight for electron scattering */

      double x1;

      x1 = PLANCK * p_in.freq / MELEC / VLIGHT / VLIGHT;
      if (x1 < 0.0001)
      {
        zz = dot (pp.lmn, p_in.lmn);
        pp.w *= 0.75 * (1 + zz * zz);
      }

      else
      {

        compton_reweight (&p_in, &pp);

      }






    }
    else if (pp.nres > -1 && pp.nres < NLINES)
    {

      /*
       * It was a wind photon.  In this
       * case, what we do depends on
       * whether it is a photon which arose
       * via line radiation or some other
       * process.
       * 
       * If
       * geo.scatter_mode==SCATTER_MODE_ISOTR
       * OPIC then there is no need to
       * reweight.  This is the isotropic
       * assumption.  Otherwise, one needs
       * to reweight
       * 
       */

      if (geo.scatter_mode == SCATTER_MODE_THERMAL)
      {
        dvds = dvwind_ds_cmf (&pp);
        ishell = pp.grid;
        tau = sobolev (&w[ishell], pp.x, -1.0, lin_ptr[pp.nres], dvds);
        if (tau > 0.0)
          pp.w *= p_escape_from_tau (tau);
        tau = 0.0;
      }
    }



    /*
     * Now , we need to transform the
     * photon back into the observer frame
     */

    if (itype == PTYPE_DISK)
    {
      if ((ierr = local_to_observer_frame_disk (&pp, &pp)))
        Error ("extract_one: disk photon not in local frame");
    }
    if (itype == PTYPE_WIND)
    {
      if (pp.nres == NRES_ES)
      {
        lorentz_transform (&pp, &pp, vel);
      }

      if ((ierr = local_to_observer_frame (&pp, &pp)))
        Error ("extract_one: wind photon not in local frame\n");
    }

    /* Make some final chacks before extracting the photon */
    if (tau > TAU_MAX)
    {
      Error ("extract: tau (%e) for photon %d should not be large \n", tau, pp.np);
      continue;
    }
    else if (geo.binary == TRUE && hit_secondary (&pp))
    {
      continue;
    }

    /* If one has reached this point, we extract the photon and increment the spectrum */


    extract_one (w, &pp, n);

  }


  return (0);
}


/**********************************************************/
/**
 * @brief      Reduce the weight of a single photon along a single line of sight.
 *
 * @param [in] WindPtr  w   The entire wind
 * @param [in] PhotPtr  pp  The photon to be extracted (in the observer frame)
 * @param [in] int  nspec   the spectrum which will be incremented
 * @return     The photon status after translation
 *
 * @details
 * This routine just extracts the photon, which is assumed to be in the
 * observer frame
 *
 * extract_one is analogous to the detailed portion of transphot except here the
 * basic point is to calculate the optical depth through the plasma in a certain
 * direction, and to increment the appropriate spectrum.
 *
 * ### Notes ###
 *
 * In Python, and in extract and transphot in particular, tau generally refers to the tau associated
 * with scattering processes, and the weight contains the effect of dimunition of the energy of
 * the photon bundle due to pure absorption processes.  So, in extract, we add pp->w * exp(-tau)
 * to the spectrum.
 *
 * Note that both linearly and logarithmically spaced spectra are produced.
 *
 **********************************************************/


int
extract_one (w, pp, nspec)
     WindPtr w;
     PhotPtr pp;
     int nspec;

{

  int istat, nres;
  struct photon pstart;
  struct photon pdummy, pdummy_orig;
  double weight_min;
  int icell;
  int k, k1;
  double tau;
  double lfreqmin, lfreqmax, ldfreq;
  double normal[3];

  /*
   * Preserve the starting position of the photon so one can use this
   * to determine whether the photon encountered the disk or star as it
   * tried to exist the wind.
   */

  weight_min = EPSILON * pp->w;
  tau = 0;
  pp->ds = 0;
  icell = 0;

  stuff_phot (pp, &pstart);
  stuff_phot (pp, &pdummy_orig);
  stuff_phot (pp, &pdummy);

  check_frame (pp, F_OBSERVER, "extract_one: photon not in observer frame at start");

  istat = P_INWIND;

  while (istat == P_INWIND)
  {
    istat = translate (w, pp, 20., &tau, &nres);
    icell++;

    stuff_phot (pp, &pdummy);
    istat = walls (pp, &pstart, normal);

    if (istat == -1)
    {

      Error ("Extract_one: Abnormal return from translate (icell %d) of phot no %5d\n", icell, pp->np);
      Error ("Extract_one: start %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
             pstart.x[0], pstart.x[1], pstart.x[2], pstart.lmn[0], pstart.lmn[1], pstart.lmn[2]);
      Error ("Extract_one:  orig %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
             pdummy_orig.x[0], pdummy_orig.x[1], pdummy_orig.x[2], pdummy_orig.lmn[0], pdummy_orig.lmn[1], pdummy_orig.lmn[2]);
      Error ("Extract_one:   was %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
             pdummy.x[0], pdummy.x[1], pdummy.x[2], pdummy.lmn[0], pdummy.lmn[1], pdummy.lmn[2]);
      Error ("Extract_one:    is %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
             pp->x[0], pp->x[1], pp->x[2], pp->lmn[0], pp->lmn[1], pp->lmn[2]);
      break;
    }
    if (pp->w < weight_min)
    {
      istat = P_ABSORB;         /* This photon was absorbed
                                 * within the wind */
      break;
    }
    if (istat == P_HIT_STAR)
    {                           /* It was absorbed in the
                                 * photosphere */
      break;
    }
    if (istat == P_HIT_DISK)
    {                           /* It was absorbed in the
                                 * disk */
      break;
    }
    if (istat == P_SCAT)
    {                           /* Cause the photon to scatter and
                                 * reinitialize */
      break;
    }
  }

//  if (modes.save_extract_photons)
//    save_photons (pp, "EXT");


  if (istat == P_ESCAPE)
  {

    if (!(0 <= tau && tau < 1.e4))
      Error_silent ("Warning: extract_one: ignoring very high tau  %8.2e at %g\n", tau, pp->freq);
    else
    {
      k = (int) ((pp->freq - xxspec[nspec].freqmin) / xxspec[nspec].dfreq);

      /*
       * Force the frequency to be in range of that
       * recorded in the spectrum
       */

      if (k < 0)
        k = 0;
      else if (k > NWAVE_EXTRACT - 1)
        k = NWAVE_EXTRACT - 1;


      lfreqmin = log10 (xxspec[nspec].freqmin);
      lfreqmax = log10 (xxspec[nspec].freqmax);
      ldfreq = (lfreqmax - lfreqmin) / NWAVE_EXTRACT;

      k1 = (int) ((log10 (pp->freq) - log10 (xxspec[nspec].freqmin)) / ldfreq);
      if (k1 < 0)
      {
        k1 = 0;
      }
      if (k1 > NWAVE_EXTRACT - 1)
      {
        k1 = NWAVE_EXTRACT - 1;
      }
      /*
       * Increment the spectrum.  Note that the photon
       * weight has not been diminished by its passage
       * through th wind, even though it may have
       * encounterd a number of resonance, and so the
       * weight must be reduced by tau
       */

      xxspec[nspec].f[k] += pp->w * exp (-(tau));
      xxspec[nspec].lf[k1] += pp->w * exp (-(tau));



      /*
       * If this photon was a wind photon, then also
       * increment the "reflected" spectrum
       */
      if (pp->origin == PTYPE_WIND || pp->origin == PTYPE_WIND_MATOM || pp->nscat > 0)
      {

        xxspec[nspec].f_wind[k] += pp->w * exp (-(tau));
        xxspec[nspec].lf_wind[k1] += pp->w * exp (-(tau));

      }
      /*
       * Records the total distance travelled by extracted
       * photon if in reverberation mode
       */
      if (geo.reverb != REV_NONE)
      {
        if (geo.reverb_filter_lines == -2 || pstart.nscat > 0 || pstart.origin > 9 || (pstart.nres > -1 && pstart.nres < nlines))
        {
          /*If this photon has scattered, been reprocessed, 
             or originated in the wind it 's important
           */
          pstart.w = pp->w * exp (-(tau));
          stuff_v (xxspec[nspec].lmn, pstart.lmn);
          delay_dump_single (&pstart, nspec);
        }
      }
    }

  }
  if (istat > -1 && istat < 9)
    xxspec[nspec].nphot[istat]++;
  else
    Error
      ("Extract: Abnormal photon %5d  %d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n",
       pp->np, istat, pp->x[0], pp->x[1], pp->x[2], pp->lmn[0], pp->lmn[1], pp->lmn[2]);

  return (istat);
}
