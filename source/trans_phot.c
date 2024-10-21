
/***********************************************************/
/** @file  trans_phot.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  This file contains high level routines which carry
 * out the propagation of photons through the wind/wind domains
 * both for ionization and spectrum (extract) cycles
 *
 * ### Notes ###
 *
 * The routines contained here are central to Python, and anyone
 * who wants to understand Python in general should spend time
 * understanding how they work
 *
 * There are two basic options associated with the routines.
 *
 * * live or die is the case in which one is only concerned
 * with how a photon goes through the wind and random numbers
 * are used to generate the scattering directions.  The photon is
 * followed until it leaves the system.  This is the only option
 * used in ionization cycles (currently).
 *
 * * extract is the mode where whenever a photon scatters we
 * also calculate what the photon weight would be if it scattered
 * in a set of certain directions.  When this happens the weight
 * of the photon is changed to reflect the fact that the direction
 * is not random.
 *
 * The extract option is used
 * normally during the spectral extraction cycles.
 * However, one can use the live or die
 * to construct the detailed spectrum.  One would not normally
 * want to do this, as many photons are "wasted" since they
 * don't scatter at the desired angle.  With sufficient numbers
 * of photons however the results of the two methods should
 * (by construction) be identical (or at least very very
 * similar).  As a result the live or die option is useful
 * if one has questions about whether the more complex 
 * extract option is functioning properly.  
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "atomic.h"
#include "sirocco.h"

long n_lost_to_dfudge = 0;


/**********************************************************/
/**
 * @brief      int (w,p,iextract) oversees the propagation of a "flight" of photons
 *
 * @param [in] WindPtr  w   The entire wind domain
 * @param [in, out] PhotPtr  p   A pointer to a "fligh" of photons
 * @param [in] int  iextract   An integer controlling whether we are to process the
 * flight in the live or die option (0) or whether we also need to extract photons in
 * specific directions (which is usually the case in constructing spectra
 * @return   Normally returns 0
 *
 * @details
 * This routine oversees the propagation of  individual photons.  The main loop
 * covers an entire "flight" of photons.   The routine generates the random
 * optical depth a photon can travel before scattering and monitors the
 * progress of the photon through the grid.
 *
 * The real physics is done elsewhere, in lower level routines.
 *
 * ### Notes ###
 *
 * At the end of the routine the position for each of the photons in p is the
 * last point where the photon was in the wind, * not the outer boundary of
 * the radiative transfer
 *
 **********************************************************/

int
trans_phot (WindPtr w, PhotPtr p, int iextract)
{
  int nphot;
  struct photon pp, pextract;
  int nreport;
  struct timeval timer_t0;

  xsignal (files.root, "%-20s Photon transport started\n", "NOK");

  nreport = NPHOT / 10;
  Log ("\n");

  timer_t0 = init_timer_t0 ();

  for (nphot = 0; nphot < NPHOT; nphot++)
  {
    p[nphot].np = nphot;
    check_frame (&p[nphot], F_OBSERVER, "trans_phot: photon not in observer frame as expeced\n");

    if (nphot % nreport == 0)
    {
      if (geo.ioniz_or_extract == CYCLE_IONIZ)
      {
        Log (" Ion. Cycle %d/%d of %s : Photon %10d of %10d or %6.1f per cent \n", geo.wcycle + 1, geo.wcycles, files.root, nphot, NPHOT,
             nphot * 100. / NPHOT);
      }
      else
      {
        Log ("Spec. Cycle %d/%d of %s : Photon %10d of %10d or %6.1f per cent \n", geo.pcycle + 1, geo.pcycles, files.root, nphot, NPHOT,
             nphot * 100. / NPHOT);
      }
    }

    Log_flush ();
    stuff_phot (&p[nphot], &pp);

    /* The next if statement is executed if we are calculating the detailed spectrum and
     * makes sure we always run extract on the original photon no matter where it
     * was generated */

    if (iextract)
    {
      stuff_phot (&p[nphot], &pextract);
      extract (w, &pextract, pextract.origin);
    }

    trans_phot_single (w, &p[nphot], iextract);
  }

  Log ("\n");

  print_timer_duration ("!!sirocco: photon transport completed in", timer_t0);
  //XXXX Delete when understand what is going on with state machines
  xsignal (files.root, "%-20s Photon transport completed\n", "NOK");

  /* Sometimes a photon will scatter near the edge of the wind and get pushed
   * out by DFUDGE. We record these. */

  if (n_lost_to_dfudge > 0)
  {
    Error
      ("trans_phot: %ld photons were lost due to DFUDGE (%8.4e) pushing them outside of the wind after scatter\n",
       n_lost_to_dfudge, DFUDGE);
  }

  n_lost_to_dfudge = 0;         // reset the counter

  return (0);
}

/**********************************************************/
/**
 * @brief      Transport a single photon photon through the wind.
 *
 * @param [in] WindPtr  w   The entire wind
 * @param [in, out] PhotPtr  p   A single photon (in the observer frame)
 * @param [in] int  iextract   If 0, then process this photon in the live or die option, without
 * calling extract
 *
 * @return     Always returns 0
 *
 * @details
 * This routine oversees the propagation of an individual photon through
 * the wind.  As the photon moves through the wind, its position, direction
 * and weight are updated.  In reverberation mode, the flight time is also
 * tracked.
 *
 * Basically what the routine does is generate a random number which is used to
 * determine the optical depth to a scatter, and then it calles translate
 * multiple times.   translate involves moving the photon only a single cell
 * (or alternatively a single transfer in the windless region), and returns
 * a status.  Depending on what this status is, trans_phot_single calls
 * trans_phot again doing nothing, but if the scattering depth has been
 * reached, then trans_phot_single causes the photon to scatter,
 * which changes its direction.  This process continues until the photon
 * exits the system or hits a barrier.  
 
 * If the photon hits a radiating
 * surface, the disk or star, then the photon may either be absorbed,
 * or scattered depending on the reflection/absorption mode. If the
 * reflection/absorption mode is set to reflection, then the 
 * program is redirected in the main loop, but if the surfaces are
 * set to absorb, hitting a surface will exit the routine.
 *
 *
 *
 * ### Notes ###
 *
 * This routine is called by trans_phot once for each photon in a flight of photons
 * Internally, there are two main PhotPtrs p, and pp.  pp is a place that
 * the photon has reached, and p is the location where it is going.  At
 * the end of the main loop before a new cycle, pp is updated.
 *
 * The photon starts in the observer frame (and remains that 
 * in this routine).  Comversions to the local frame occur
 * in the routines that are called, however.
 *
 **********************************************************/

int
trans_phot_single (WindPtr w, PhotPtr p, int iextract)
{
  double tau_scat, tau;
  int i, n_grid, ierr;
  enum istat_enum istat;
  int nnscat;
  int current_nres;
  double weight_min;
  struct photon pp, pextract;
  double normal[3];
  double rho, dz;


  /* Initialize parameters that are needed for the flight of the photon through the wind */

  stuff_phot (p, &pp);
  tau_scat = -log (1. - random_number (0.0, 1.0));
  weight_min = EPSILON * pp.w;
  istat = P_INWIND;
  tau = 0;

  /* This is the beginning of the loop for a single photon and executes until the photon leaves the wind */

  while (istat == P_INWIND)
  {

    /* The call to translate below involves only a single cell (or alternatively a single transfer 
       in the windless region). The returned value, istat, should either 1) be P_INWIND in which case the photon
       hit the other side of the cell without scattering, 2) P_SCAT in which case there was a scattering event
       in the cell, 3) P_ESCAPE in which case the photon reached the outside edge of the grid and escaped, 4)
       P_STAR in which case it reach the inner central object, etc. If the photon escapes then we leave the
       photon at the position of it's last scatter.  In most other cases though we store the final 
       position of the photon. */

    istat = translate (w, &pp, tau_scat, &tau, &current_nres);

    if (istat == P_ERROR)
    {
      Error ("trans_phot: abnormal return from translate on photon %d\n", p->np);
      break;
    }

    if (pp.w < weight_min)
    {
      pp.istat = P_ABSORB;
      pp.tau = VERY_BIG;
      stuff_phot (&pp, p);
      break;
    }

    /* Check boundary with walls - note that pp is the proposed new photon location
     * and p is the "original" location of the photon */

    istat = walls (&pp, p, normal);


    if (istat == P_HIT_STAR)
    {

      /*
       * The photon has hit the star. Reflect or absorb.
       */

      geo.lum_star_back += pp.w;
      spec_add_one (&pp, SPEC_HITSURF);

      /* The a new photon direction needs to be defined that will cause the photon to continue in the wind.
       * Since this is effectively a scattering event we also have to extract a photon to construct the
       * detailed spectrum.
       */

      if (geo.absorb_reflect == BACK_RAD_SCATTER)
      {
        randvcos (pp.lmn, normal);
        if (move_phot (&pp, DFUDGE))
        {
          Error ("trans_phot_single: photon not in correct frame when reflecting off of star\n");
        }

        p->ds = 0;
        tau_scat = -log (1. - random_number (0.0, 1.0));
        istat = pp.istat = P_INWIND;    /* Set the status back to P_INWIND so the photon will continue */
        tau = 0;
        stuff_phot (&pp, p);

        if (iextract)
        {
          stuff_phot (&pp, &pextract);
          extract (w, &pextract, PTYPE_STAR);   // Treat as stellar photon for purpose of extraction
        }
      }
      else                      /*Photons that hit the star are simply absorbed  */
      {
        stuff_phot (&pp, p);
        break;
      }
    }

    if (istat == P_HIT_DISK)
    {
      /*
       * The photon has hit the disk. Reflect or absorb.
       */

      /* Store the energy of the photon bundle into a disk structure so that one
         can determine later how much and where the disk was heated by photons.
         Note that the disk is defined from 0 to NRINGS-2. NRINGS-1 contains the
         position of the outer radius of the disk. */

      rho = sqrt (pp.x[0] * pp.x[0] + pp.x[1] * pp.x[1]);

      i = 0;
      while (rho > qdisk.r[i] && i < NRINGS - 1)
        i++;
      i--;                      /* So that the heating refers to the heating between i and i+1 */

      qdisk.nhit[i]++;
      geo.lum_disk_back = qdisk.heat[i] += pp.w;
      qdisk.ave_freq[i] += pp.w * pp.freq;

      if (geo.absorb_reflect == BACK_RAD_SCATTER)
      {
        /*
         * If the disk is vertically extended, then we need to move the photon
         * outside of the disk and push it by a little amount. It's unclear
         * to me why we haven't used dfudge here.
         */

        if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
        {
          dz = (zdisk (rho) - fabs (pp.x[2]));
          if (dz > 0)
          {
            if (pp.x[2] > 0)
            {
              pp.x[2] += (dz + 1000.);
            }
            else
            {
              pp.x[2] -= (dz + 1000.);
            }
          }
        }

        spec_add_one (&pp, SPEC_HITSURF);

        /* If we got here, a new photon direction needs to be defined that will cause the photon
         * to continue in the wind.  Since this is effectively a scattering event we also have to
         * extract a photon to construct the detailed spectrum.
         */

        randvcos (pp.lmn, normal);
        p->ds = 0;
        tau_scat = -log (1. - random_number (0.0, 1.0));
        istat = pp.istat = P_INWIND;
        tau = 0;
        stuff_phot (&pp, p);

        if (iextract)
        {
          stuff_phot (&pp, &pextract);
          extract (w, &pextract, PTYPE_DISK);
        }
      }
      else                      /* Photons that hit the disk are to be absorbed */
      {
        stuff_phot (&pp, p);
        break;
      }
    }

    if (istat == P_SCAT)
    {

      /*
       * The photon has scattered, as either a resonance or continuum scatter.
       */

      pp.grid = n_grid = where_in_grid (wmain[pp.grid].ndom, pp.x);

      if (n_grid < 0)
      {
        Error ("trans_phot: trying to scatter a photon which is not in the wind grid and the photon has been lost\n");
        Error ("trans_phot: %d grid %3d x %8.2e %8.2e %8.2e (%8.2e)\n", pp.np, pp.grid, pp.x[0], pp.x[1], pp.x[2],
               sqrt (pp.x[0] * pp.x[0] + pp.x[1] * pp.x[1] + pp.x[2] * pp.x[2]));
        pp.istat = P_ERROR;
        stuff_phot (&pp, p);
        break;
      }

      if (wmain[n_grid].nplasma == NPLASMA)     /* If the next error reoccurs, see Issue #154 (on GitHub) for discussion */
      {
        Error ("trans_phot: Trying to scatter a photon which is not in a cell in the plasma structure\n");
        Error ("trans_phot: %d grid %3d x %8.2e %8.2e %8.2e\n", pp.np, pp.grid, pp.x[0], pp.x[1], pp.x[2]);
        Error ("trans_phot: This photon is effectively lost!\n");
        pp.istat = P_ERROR;
        stuff_phot (&pp, p);
        break;
      }

      if (wmain[n_grid].inwind < 0)
      {
        Error ("trans_phot: Trying to scatter a photon in a cell with no wind volume and the photon has been lost\n");
        Error ("trans_phot: istat %d %d grid %3d x %8.2e %8.2e %8.2e\n", istat, pp.np, pp.grid, pp.x[0], pp.x[1], pp.x[2]);
        pp.istat = P_ERROR;
        stuff_phot (&pp, p);
        break;
      }

      /* Add path lengths for reverberation mapping */

      if ((geo.reverb == REV_WIND || geo.reverb == REV_MATOM) && geo.ioniz_or_extract == CYCLE_IONIZ && geo.wcycle == geo.wcycles - 1)
      {
        wind_paths_add_phot (&wmain[n_grid], &pp);
      }

      nnscat = 1;
      pp.nscat++;

      if (current_nres == NRES_ES)
      {
        stuff_phot (&pp, &pextract);
      }

      if ((ierr = scatter (&pp, &current_nres, &nnscat)))       // pp is modified
      {
        Error ("trans_phot_single: photon %d returned error code %d whilst scattering \n", pp.np, ierr);
      }

      if (geo.matom_radiation == 1 && geo.rt_mode == RT_MODE_MACRO && pp.w < weight_min)
      {
        pp.istat = P_ABSORB;
        pp.tau = VERY_BIG;
        stuff_phot (&pp, p);
        break;
      }

      /* If this is a BB interaction, calculate the line heating
       * and break the transport loop if it was absorbed */

      if (current_nres > -1 && current_nres < nlines)
      {
        pp.nrscat++;

        if (modes.track_resonant_scatters)
          track_scatters (&pp, wmain[n_grid].nplasma, "Resonant");

        plasmamain[wmain[n_grid].nplasma].scatters[line[current_nres].nion] += 1;

        if (geo.rt_mode == RT_MODE_2LEVEL)
        {
          line_heat (&plasmamain[wmain[n_grid].nplasma], &pp, current_nres);
        }

        if (pp.w < weight_min)
        {
          pp.istat = P_ABSORB;
          pp.tau = VERY_BIG;
          stuff_phot (&pp, p);
          break;
        }
      }

      if (pp.w < weight_min)
      {
        pp.istat = P_ABSORB;
        pp.tau = VERY_BIG;
        stuff_phot (&pp, p);
        break;
      }

      if (pp.istat == P_ERROR_MATOM || pp.istat == P_LOFREQ_FF || pp.istat == P_ADIABATIC)
      {
        p->istat = pp.istat;
        pp.tau = VERY_BIG;
        stuff_phot (&pp, p);
        break;
      }

      /* Now extract photons if we are in detailed the detailed spectrum portion of the program
       * N.B. To use the anisotropic scattering option, extract needs to follow scatter.
       * This is because the re-weighting which occurs in extract needs the pdf for scattering
       * to have been initialized
       *
       * For BB photons the photon we pass to extract is the one that has been scattered, but
       * for ES we pass the photon prior to scattering.
       */

      if (iextract)
      {
        if (current_nres != NRES_ES)
        {
          stuff_phot (&pp, &pextract);
        }
        pextract.nnscat = nnscat;
        extract (w, &pextract, PTYPE_WIND);
      }

      /* Reinitialize parameters for the scattered photon so it can can continue through the wind
       */

      tau = 0;
      tau_scat = -log (1. - random_number (0.0, 1.0));
      pp.istat = P_INWIND;
      pp.ds = 0;

      stuff_phot (&pp, p);
      istat = p->istat;
    }

    /*
     * Now we check if a photon has gotten stuck scattering in the wind, this
     * is mostly done for speed concerns as it is pointless to track photons which
     * are probably low weight and not contributing. However, in some cases MAXSCAT
     * should be increased to stop the code from throwing away too many photons
     */

    if (pp.nscat == MAXSCAT)
    {
      pp.istat = P_TOO_MANY_SCATTERS;
      stuff_phot (&pp, p);
      break;
    }

    if (pp.istat == P_ERROR_MATOM || pp.istat == P_LOFREQ_FF || pp.istat == P_ADIABATIC)
    {
      p->istat = pp.istat;
      stuff_phot (&pp, p);
      break;
    }

    /* This is an insurance policy but it is not obvious that, for example nscat
     * and nrscat, need to be updated */
    p->istat = istat;
    p->nscat = pp.nscat;
    p->nrscat = pp.nrscat;
    p->w = pp.w;
  }

  /* This is set up for looking at photons in spectral cycles at present */
  // if (modes.save_photons && geo.ioniz_or_extract == CYCLE_EXTRACT)
  //   save_photons (&pp, "End");

  return (0);
}
