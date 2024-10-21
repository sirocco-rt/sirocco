/***********************************************************/
/** @file  radiation.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  Routines to update radiation field parameters and to
 * calculate various opacities (free-free, bound-free, etc.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

#define COLMIN  0.01

int iicount = 0;


/**********************************************************/
/**
 * @brief      updates the  field parameters in the wind and reduces
 * the weight of the photon as a result of the effects of  free free, bound-free and
 * Compton scattering. The routine
 * also keeps track of the number of photoionizations for H and He in the
 * cell.
 *
 * @param [in,out] PhotPtr  p   the photon
 * @param [in] double  ds   the distance the photon has travelled in the cell
 * @return     Usually regurns kappa_tot 
 *
 * @details
 *
 * The pieces of the wind structure which are updated in each call to radiation
 * are:
 *   j,ave_freq,ntot, heat_photo, heat_ff, heat_h, heat_he1, heat_he2, heat_z,
 *   nioniz, and ioniz[].
 *
 * @details
 *
 * ### Notes ###
 * The # of ionizations of a specific ion = (w(0)-w(s))*n_i sigma_i/ (h nu * kappa_tot).  (The # of ionizations
 * is just given by the integral of n_i sigma_i w(s) / h nu ds, but if the density is assumed to
 * be constant and sigma is also constant [therefy ignoring velocity shifts in a grid cell], then
 * n_i sigma and h nu can be brought outside the integral.  Hence the integral is just over w(s),
 * but that is just given by (w(0)-w(smax))/kappa_tot.)  The routine calculates the number of ionizations per
 * unit volume.
 *
 * Inputs to radiation are assumed to be in the observer frame.  kappas are calculated
 * in the CMF frame, as elsewhere.  Where tau is calculted from kappa ds, one needs to
 * account for the difference in length in the two frames
 *
 **********************************************************/

double
radiation (PhotPtr p, double ds)
{
  TopPhotPtr x_top_ptr;

  WindPtr one;
  PlasmaPtr xplasma;


  double freq;
  double kappa_tot, kappa_tot_obs, frac_tot, frac_ff;
  /* Variables named frac are very pooorly named; they are actually opacities due to individual
     processes, not fractions of anything
   */
  double frac_z, frac_comp;     /* frac_comp - the heating in the cell due to Compton heating */
  double frac_ind_comp;         /* frac_ind_comp - the heating due to induced Compton heating */
  double frac_auger;
  double frac_tot_abs, frac_auger_abs, z_abs;
  double kappa_ion[NIONS];
  double frac_ion[NIONS];
  double kappa_inner_ion[n_inner_tot];
  double frac_inner_ion[n_inner_tot];
  double density, ft, tau, tau2;
  double energy_abs_obs, energy_abs_cmf;
  int n, nion;
  double q, x, z, z_obs;
  double w_ave_obs, w_in, w_out;
  int nconf;
  double p_in[3];               //The initial and final momentum.
  double freq_inner, freq_outer;
  double freq_min, freq_max;
  double frac_path, freq_xs;
  struct photon phot, phot_cmf, phot_mid, phot_mid_cmf, p_cmf;
  int ndom;
  double ds_cmf, w_ave_cmf;



  z = frac_path = freq_xs = 0;


  one = &wmain[p->grid];

  ndom = one->ndom;
  xplasma = &plasmamain[one->nplasma];

  /* See #870 if this next error is triggered many times ksl 210614 */
  if (check_plasma (xplasma, "radiation"))
  {
    Error ("Radiation: Photon %d is in wind cell %d which has no volume in the wind\n", p->np, p->grid);
    return (0);
  }

  /* JM 140321 -- #73 Bugfix
     In previous version we were comparing these frequencies in the global rest frame
     The photon frequency should be shifted to the rest frame of the cell in question
     We currently take the average of this frequency along ds. In principle
     this could be improved, so we throw an error if the difference between v1 and v2 is large */


  /* compute the initial momentum of the photon */

  stuff_v (p->lmn, p_in);       //Get the direction
  renorm (p_in, p->w / VLIGHT); //Renormalise to momentum

  /* Create phot, a photon at the position we are moving to
   *  note that the actual movement of the photon gets done after
   *  the call to radiation
   */

  stuff_phot (p, &phot);
  move_phot (&phot, ds);
  stuff_phot (p, &phot_mid);
  move_phot (&phot_mid, ds / 2.);

  /* calculate photon frequencies in rest frame of cell */

  if (observer_to_local_frame (&phot, &phot_cmf))
  {
    Error ("radiation: observer to local frame error\n");
  }

  if (observer_to_local_frame (p, &p_cmf))
  {
    Error ("radiation: observer to local frame error\n");
  }

  if (observer_to_local_frame (&phot_mid, &phot_mid_cmf))
  {
    Error ("radiation: observer to local frame error\n");
  }


  freq_inner = p_cmf.freq;
  freq_outer = phot_cmf.freq;

  phot_mid_cmf.freq = freq = 0.5 * (freq_inner + freq_outer);
  phot_mid.freq = 0.5 * (p->freq + phot.freq);

  /* calculate free-free, Compton and induced-Compton opacities
     note that we also call these with the average frequency along ds */

  kappa_tot = frac_ff = kappa_ff (xplasma, freq);
  kappa_tot += frac_comp = kappa_comp (xplasma, freq);
  kappa_tot += frac_ind_comp = kappa_ind_comp (xplasma, freq);

  frac_tot = frac_z = 0;
  frac_auger = 0;
  frac_tot_abs = frac_auger_abs = 0.0;

  /* Check which of the frequencies is larger.  */

  if (freq_outer > freq_inner)
  {
    freq_max = freq_outer;
    freq_min = freq_inner;
  }
  else
  {
    freq_max = freq_inner;
    freq_min = freq_outer;
  }

  if (freq > phot_freq_min)

  {
    if (geo.ioniz_or_extract == CYCLE_IONIZ)
    {
      for (nion = 0; nion < nions; nion++)
      {
        kappa_ion[nion] = 0;
        frac_ion[nion] = 0;
      }
      for (n = 0; n < n_inner_tot; n++)
      {
        kappa_inner_ion[n] = 0;
        frac_inner_ion[n] = 0;
      }
    }

    /* Next section is for photoionization with Topbase.  There may be more
       than one x-section associated with an ion, and so one has to keep track
       of the energy that goes into heating electrons carefully.  */

    /* JM 1405 -- I've added a check here that checks if a photoionization edge has been crossed.
       If it has, then we multiply sigma*density by a factor frac_path, which is equal to the how far along
       ds the edge occurs in frequency space  [(ft - freq_min) / (freq_max - freq_min)] */


    /* Next steps are a way to avoid the loop over photoionization x sections when it should not matter */
    if (DENSITY_PHOT_MIN > 0)
    {

      /* Loop over all photoionization xsections */
      for (n = 0; n < nphot_total; n++)
      {
        x_top_ptr = phot_top_ptr[n];
        ft = x_top_ptr->freq[0];
        if (ft > freq_min && ft < freq_max)
        {
          /* then the shifting of the photon causes it to cross an edge.
             Find out where between fmin and fmax the edge would be in freq space.
             frac_path is the fraction of the total path length above the absorption edge
             freq_xs is freq halfway between the edge and the max freq if an edge gets crossed */
          frac_path = (freq_max - ft) / (freq_max - freq_min);
          freq_xs = 0.5 * (ft + freq_max);
        }

        else if (ft > freq_max)
          break;                // The remaining transitions will have higher thresholds

        else if (ft < freq_min)
        {
          frac_path = 1.0;      // then the frequency of the photon is above the threshold all along the path
          freq_xs = freq;       // use the average frequency
        }

        if (freq_xs < x_top_ptr->freq[x_top_ptr->np - 1])
        {
          /* Need the appropriate density at this point.
             how we get this depends if we have a topbase (level by level)
             or vfky cross-section (ion by ion) */

          nion = x_top_ptr->nion;
          if (ion[nion].phot_info > 0)  // topbase or hybrid
          {
            nconf = x_top_ptr->nlev;
            density = den_config (xplasma, nconf);
          }

          else if (ion[nion].phot_info == 0)    // verner
            density = xplasma->density[nion];

          else
          {
            Error ("radiation.c: No type (%i) for xsection!\n");
            density = 0.0;
          }

          if (density > DENSITY_PHOT_MIN)
          {

            /* Note that this includes a filling factor  */
            kappa_tot += x = sigma_phot (x_top_ptr, freq_xs) * density * frac_path * zdom[ndom].fill;

            if (geo.ioniz_or_extract == CYCLE_IONIZ)
            {

              //This is the heating effect - i.e. the absorbed photon energy less the binding energy of the lost electron
              frac_tot += z = x * (freq_xs - ft) / freq_xs;
              //This is the absorbed energy fraction
              frac_tot_abs += z_abs = x;

              if (nion > 3)
              {
                frac_z += z;
              }

              frac_ion[nion] += z;
              kappa_ion[nion] += x;
            }

          }

        }
      }

      /* Loop over all inner shell cross sections as well! But only for VFKY ions - topbase has those edges in */

      if (freq > inner_freq_min)
      {
        for (n = 0; n < n_inner_tot; n++)
        {
          if (ion[inner_cross_ptr[n]->nion].phot_info != 1)
          {
            x_top_ptr = inner_cross_ptr[n];
            if (x_top_ptr->n_elec_yield != -1)
            {
              ft = x_top_ptr->freq[0];

              if (ft > freq_min && ft < freq_max)
              {
                frac_path = (freq_max - ft) / (freq_max - freq_min);
                freq_xs = 0.5 * (ft + freq_max);
              }
              else if (ft > freq_max)
                break;          // The remaining transitions will have higher thresholds
              else if (ft < freq_min)
              {
                frac_path = 1.0;        // then all frequency along ds are above edge
                freq_xs = freq; // use the average frequency
              }
              if (freq_xs < x_top_ptr->freq[x_top_ptr->np - 1])
              {
                nion = x_top_ptr->nion;
                if (ion[nion].phot_info == 0)   // verner only ion
                {
                  density = xplasma->density[nion];     //All these rates are from the ground state, so we just need the density of the ion.
                }
                else
                {
                  nconf = phot_top[ion[nion].ntop_ground].nlev; //The lower level of the ground state Pi cross section (should be GS!)
                  density = den_config (xplasma, nconf);
                }
                if (density > DENSITY_PHOT_MIN)
                {
                  kappa_tot += x = sigma_phot (x_top_ptr, freq_xs) * density * frac_path * zdom[ndom].fill;
//xxxx                  kappa_tot += x = exp (log_sigma_phot (x_top_ptr, log (freq_xs))) * density * frac_path * zdom[ndom].fill;

                  if (geo.ioniz_or_extract && x_top_ptr->n_elec_yield != -1)    // Calculate during ionization cycles only
                  {
                    frac_auger += z = x * (inner_elec_yield[x_top_ptr->n_elec_yield].Ea / EV2ERGS) / (freq_xs * HEV);
                    frac_auger_abs += z_abs = x;        //This is the absorbed energy fraction

                    if (nion > 3)
                    {
                      frac_z += z;
                    }
                    frac_inner_ion[n] += z;     //NSH We need to log the auger rate seperately - we do this by cross section
                    kappa_inner_ion[n] += x;    //NSH and we also og the opacity by ion
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  /* finished looping over cross-sections to calculate bf opacity
     we can now reduce weights and record certain estimators */

  kappa_tot_obs = kappa_tot / observer_to_local_frame_ds (p, 1);
  tau = kappa_tot_obs * ds;

  w_in = p->w;

  if (sane_check (tau))
  {
    Error ("Radiation:sane_check CHECKING ff=%e, comp=%e, ind_comp=%e\n", frac_ff, frac_comp, frac_ind_comp);
  }
/* Calculate the heating effect*/

  if (tau > 0.0001)
  {
    x = exp (-tau);
    energy_abs_obs = w_in * (1. - x);
  }
  else
  {
    tau2 = tau * tau;
    energy_abs_obs = w_in * (tau - 0.5 * tau2);
  }

  energy_abs_cmf = energy_abs_obs * phot_mid_cmf.freq / phot_mid.freq;

  /* Calculate the reduction in weight - Compton scattering is not included, it is now included at scattering
     however induced Compton heating is not implemented at scattering, so it should remain here for the time being
     to maimtain consistency. */

  tau = kappa_tot_obs * (1. - frac_comp / kappa_tot) * ds;

  if (sane_check (tau))
  {
    Error ("Radiation:sane_check CHECKING ff=%e, comp=%e, ind_comp=%e\n", frac_ff, frac_comp, frac_ind_comp);
  }
  /* Calculate the reduction in the w of the photon bundle along with the average
     weight in the cell */

  if (tau > 0.0001)
  {
    x = exp (-tau);
    p->w = w_out = w_in * x;
    w_ave_obs = (w_in - w_out) / tau;
  }
  else
  {
    tau2 = tau * tau;
    p->w = w_out = w_in * (1. - tau + 0.5 * tau2);
    w_ave_obs = w_in * (1. - 0.5 * tau + 0.1666667 * tau2);
  }

  phot_mid.w = w_ave_obs;
  phot_mid_cmf.w = w_ave_cmf = w_ave_obs * phot_mid_cmf.freq / phot_mid.freq;


  if (sane_check (p->w))
  {
    Error ("Radiation:sane_check photon weight is %e for photon %d tau %e\n", p->w, p->np, tau);
    p->w = 0.0;
    Error ("Radiation: after same check, set weight of photont %d to zerro\n", p->np);
  }

  if (geo.ioniz_or_extract == CYCLE_EXTRACT)
    return kappa_tot;

/* Everything after this point is only needed for ionization calculations */
/* Update the radiation parameters used ultimately in calculating t_r */

  if (freq > xplasma->max_freq) // check if photon frequency exceeds maximum frequency - use doppler shifted frequency
    xplasma->max_freq = freq;   // set maximum frequency sen in the cell to the mean doppler shifted freq - see bug #391

  if (modes.save_cell_stats && ncstat > 0)
  {
    save_photon_stats (one, p, ds, w_ave_obs);  // save photon statistics (extra diagnostics)
  }

  /* JM 1402 -- the banded versions of j, ave_freq etc. are now updated in update_banded_estimators,
   * which also updates the ionization parameters and scattered and direct contributions 
   *
   * Following bug #391, we now wish to use the mean, doppler shifted freqiency in the cell.
   * Update_banded_estimators requires photon freq and ds and w_ave in cmf frame */

  ds_cmf = observer_to_local_frame_ds (&phot_mid, ds);
  update_banded_estimators (xplasma, &phot_mid_cmf, ds_cmf, w_ave_cmf, ndom);
  update_flux_estimators (xplasma, &phot_mid, ds, w_ave_obs, ndom);


  if (sane_check (xplasma->j) || sane_check (xplasma->ave_freq))
  {
    Error ("radiation:sane_check Problem with j %g or ave_freq %g\n", xplasma->j, xplasma->ave_freq);
  }

  if (kappa_tot > 0)
  {

    // Use the cmf value of the energy aborbed 

    z = (energy_abs_cmf) / kappa_tot;
    xplasma->heat_ff += z * frac_ff;
    xplasma->heat_tot += z * frac_ff;
    xplasma->abs_tot += z * frac_ff;    /* The energy absorbed from the photon field in this cell */

    xplasma->heat_comp += z * frac_comp;        /* Calculate the heating in the cell due to Compton heating */
    xplasma->heat_tot += z * frac_comp; /* Add the Compton heating to the total heating for the cell */
    xplasma->abs_tot += z * frac_comp;  /* The energy absorbed from the photon field in this cell */
    xplasma->abs_tot += z * frac_ind_comp;      /* The energy absorbed from the photon field in this cell */

    xplasma->heat_tot += z * frac_ind_comp;     /* Calculate the heating in the cell due to induced Compton heating */
    xplasma->heat_ind_comp += z * frac_ind_comp;        /* Increment the induced Compton heating counter for the cell */
    if (freq > phot_freq_min)
    {
      xplasma->abs_photo += z * frac_tot_abs;   //Here we store the energy absorbed from the photon flux - different from the heating by the binding energy
      xplasma->abs_auger += z * frac_auger_abs; //same for auger
      xplasma->abs_tot += z * frac_tot_abs;     /* The energy absorbed from the photon field in this cell */
      xplasma->abs_tot += z * frac_auger_abs;   /* The energy absorbed from the photon field in this cell */

      xplasma->heat_photo += z * frac_tot;
      xplasma->heat_z += z * frac_z;
      xplasma->heat_tot += z * frac_tot;        //All of the photoinization opacities
      xplasma->heat_auger += z * frac_auger;
      xplasma->heat_tot += z * frac_auger;      //All the inner shell opacities

      q = (z) / (PLANCK * freq * xplasma->vol);

      /* So xplasma->ioniz for each species is just
         (energy_abs)*kappa_h/kappa_tot / PLANCK*freq / volume
         or the number of photons absorbed in this bundle per unit volume by this ion
       */

      for (nion = 0; nion < nions; nion++)
      {
        xplasma->ioniz[nion] += kappa_ion[nion] * q;
        xplasma->heat_ion[nion] += frac_ion[nion] * z;
      }
      for (n = 0; n < n_inner_tot; n++)
      {
        xplasma->heat_inner_ion[inner_cross_ptr[n]->nion] += frac_inner_ion[n] * z;     //This quantity is per ion - the ion number comes from the freq ordered cross section
        xplasma->inner_ioniz[n] += kappa_inner_ion[n] * q;      //This is the number of ionizations from this innershell cross section - at this point, inner_ioniz is ordered by frequency
      }
    }
  }

  z_obs = z * p_cmf.freq / p->freq;
  update_force_estimators (xplasma, p, &phot_mid, ds, w_ave_obs, ndom, z_obs, frac_ff, frac_auger, frac_tot);

  return kappa_tot;
}



/**********************************************************/
/**
 * @brief      calculates the free-free opacity allowing for stimulated emission
 *
 * @param [in] PlasmaPtr  xplasma   The plasma cell where the free-free optacity is to be calculated
 * @param [in] double  freq   The frequency at which kappa_ff is to be calculated
 * @return     kappa
 *
 * @details
 * Uses the formula from Allen, Astrophysical Quantities 3rd Edition, p 102.
 *
 * ### Notes ###
 *
 * Early versions of this routine included only singly ionized H and doubly ionized He
 * in the calculation and did not include a gaunt factorim and the routine will 
 * fall back to this if Sutherland's temperature averaged gaunt factors are not
 * included as part of the atomic data.
 *
 * Normally, Surtherland's data is part of the atomic data, and all all ions included 
 * in the calculation of the ff opancity.  The relevant prefactor fo kappa_ff is calculated in
 * the routine pop_kappa_ff_array and stored in kappa_ff_factor. The gaunt factor as currently
 * implemented is a frequency averaged one due to Sutherland (see the routine gaunt_ff), and so is 
 * approximate (but better than
 * just using 1). A future upgrade would use a more complex implementation where we
 * use the frequency dependant gaunt factor.
 *
 * See issue #750 for a bit more discussion.
 **********************************************************/

double
kappa_ff (xplasma, freq)
     PlasmaPtr xplasma;
     double freq;
{
  double x;
  double exp ();
  double x1, x2, x3;
  int ndom;

  if (xplasma->nwind < 0)
  {
    Error ("kappa_ff: xplasma->nplasma %d xplasma->nwind %d is invalid\n", xplasma->nplasma, xplasma->nwind);
    return (0.0);
  }



  ndom = wmain[xplasma->nwind].ndom;

  if (gaunt_n_gsqrd == 0)       //Maintain old behaviour if Sutherland gaunt factors are unavailable
  {
    if (nelements > 1)
    {
      x = x1 = 3.692e8 * xplasma->ne * (xplasma->density[1] + 4. * xplasma->density[4]);
    }
    else
    {
      x = x1 = 3.692e8 * xplasma->ne * (xplasma->density[1]);
    }
  }
  else
  {
    x = x1 = xplasma->kappa_ff_factor;
  }
  x *= x2 = (1. - exp (-H_OVER_K * freq / xplasma->t_e));
  x /= x3 = (sqrt (xplasma->t_e) * freq * freq * freq);

  x *= zdom[ndom].fill;

  return (x);
}



/**********************************************************/
/**
 * @brief      calculates the
 *  photionization crossection due to a Topbase level associated with
 *  x_ptr at frequency freq
 *
 * @param [in,out] struct topbase_phot *  x_ptr   The structure that contains
 * TopBase information about the photoionization x-section
 * @param [in] double  freq   The frequency where the x-section is to be calculated
 *
 * @return     The x-section
 *
 * @details
 * sigma_phot uses the Topbase x-sections to calculate the bound free
 * (or photoionization) xsection.   The data must have been into the
 * photoionization structures xphot with get_atomic_data and the
 * densities of individual ions must have been calculated previously.
 *
 * ### Notes ###
 * The x-section and frequency are stored in the topbase_phot struture
 * so that if one requests the same xsection with the same frequency  
 * again, then the calculation of the x-section is avoided.
 *
 **********************************************************/

double
sigma_phot (x_ptr, freq)
     struct topbase_phot *x_ptr;
     double freq;
{
  int nmax;
  double xsection;
  double frac, fbot, ftop;
  int linterp ();
  int nlast;

  if (freq < x_ptr->freq[0])
    return (0.0);               // Since this was below threshold
  if (freq == x_ptr->f)
    return (x_ptr->sigma);      // Avoid recalculating xsection

  if (x_ptr->nlast > -1)
  {
    nlast = x_ptr->nlast;
    if ((fbot = x_ptr->freq[nlast]) < freq && freq < (ftop = x_ptr->freq[nlast + 1]))
    {
//  frac = (log (freq) - log (fbot)) / (log (ftop) - log (fbot));
// xsection = exp ((1. - frac) * log (x_ptr->x[nlast]) + frac * log (x_ptr->x[nlast + 1]));
      frac = (log (freq) - x_ptr->log_freq[nlast]) / (x_ptr->log_freq[nlast + 1] - x_ptr->log_freq[nlast]);
      xsection = exp ((1. - frac) * x_ptr->log_x[nlast] + frac * x_ptr->log_x[nlast + 1]);


      x_ptr->sigma = xsection;
      x_ptr->f = freq;

      return (xsection);
    }
  }

/* Calculate the x-section */
  nmax = x_ptr->np;
  x_ptr->nlast = linterp (freq, &x_ptr->freq[0], &x_ptr->x[0], nmax, &xsection, 1);     //call linterp in log space

  x_ptr->sigma = xsection;
  x_ptr->f = freq;

  return (xsection);

}


/**********************************************************/
/**
 * @brief      returns the precalculated density
 *  of a particular "nlte" level.
 *
 * @param [in] PlasmaPtr  xplasma   The Plasma structure containing information about
 * the cell of interest.
 * @param [in] int  nconf   The running index that describes which level we are interested in
 * @return     The density for a particular level of an ion in the cell
 *
 * @details
 * The first case is when the density of the particular level
 * is in the wind array The second case is when the density of the levels
 * for an ion are not tracked, and hence only the total density for the
 * ion is known.  In that case, we assume the ion is in it's ground state.
 *
 * If levels are not defined for an ion it
 * assumes that all ions of are in the ground state.
 *
 * ### Notes ###
 *
 **********************************************************/

double
den_config (xplasma, nconf)
     PlasmaPtr xplasma;
     int nconf;
{
  double density;
  int nnlev, nion;

  nnlev = xconfig[nconf].nden;
  nion = xconfig[nconf].nion;

  if (nnlev >= 0)
  {                             // Then a "non-lte" level with a density
    density = xplasma->levden[nnlev] * xplasma->density[nion];
  }
  else if (nconf == ion[nion].firstlevel)
  {
/* Then we are using a Topbase photoionization x-section for this ion,
but we are not storing any densities, and so we assume it is completely in the
in the ground state */
    density = xplasma->density[nion];
  }
  else
  {
    density = 0;
  }

  return (density);

}



/**********************************************************/
/**
 * @brief      populates the multiplicative
 * constant, including a gaunt factor, to be  used in the 
 * calculating free free  emission (and absorption). 
 *
 * @return  Always returns 0
 *
 * @details
 * The routine populates plasmamain[].kappa_ff_factor
 *
 * The free-free multiplicative constant depends only
 * on the densities of ions in the cell, and the electron
 * temperature (which feeds into the gaunt factor) so it
 * saves time to calculate all this just the once. 
 *
 * ### Notes ###
 * This routine
 * is called just before the photons are 
 * sent through the wind.
 *
 * The gaunt factor is calculated in gaunt_ff (and uses
 * Sutherlands interpolation formula)  
 *
 * See #750 for more explanation
 *
 **********************************************************/

double
pop_kappa_ff_array ()
{

  double gsqrd, gaunt, sum;
  int i, j;


  for (i = 0; i < NPLASMA; i++)
  {
    sum = 0.0;
    for (j = 0; j < nions; j++)
    {
      if (ion[j].istate != 1)   //The neutral ion does not contribute
      {
        gsqrd = ((ion[j].istate - 1) * (ion[j].istate - 1) * RYD2ERGS) / (BOLTZMANN * plasmamain[i].t_e);
        gaunt = gaunt_ff (gsqrd);
        sum += plasmamain[i].density[j] * (ion[j].istate - 1) * (ion[j].istate - 1) * gaunt;
        if (sane_check (sum))
        {
          Error ("pop_kappa_ff_array:sane_check sum is %e this is a problem, possible in gaunt %e\n", sum, gaunt);
        }
      }
      else
      {
        sum += 0.0;
      }

    }
    plasmamain[i].kappa_ff_factor = plasmamain[i].ne * sum * 3.692e8;
  }

  return (0);
}

/**********************************************************/
/** 
 * @brief      returns a value for the mean intensity for a specific cell
 *             at a specific frequency
 *
 * @param[in] PlasmaPtr xplasma  The plasma cell
 * @param[in] double freq        The frequency at which we want j_bar
 * @param[in] int mode           If mode == MEAN_INTENSITY_BB_MODEL (1) then we
 *                               are happy to use a BB model as a fallback. Otherwise,
 *                               j_bar = 0 can be returned when there is not
 *                               a decent model.
 *
 * @return     The mean intensity j_bar at a specific frequency for the cell
 *
 * @details
 * This subroutine returns a value for the mean intensity J at a 
 * given frequency, using either a dilute blackbody model
 * or a spectral model depending on the value of geo.spec_mod, which counts the
 * number of spectral models available.
 *
 * For ionization models that make use of the crude spectra accumulated
 * in crude spectral bands, the routine uses these bands to
 * get the mean intensity.  If however, one is using one of the other
 * (older) ionization modes, then a dilute blackbody estimate is used instead.
 *
 * ### Notes ###
 * This subroutine was produced when we started to use a spectral model to
 * populate the upper state of a two level atom, as well as to calculate induced
 * Compton heating.
 *
 **********************************************************/

double
mean_intensity (xplasma, freq, mode)
     PlasmaPtr xplasma;
     double freq;
     int mode;
{
  double j_bar;

  if (geo.ioniz_mode == IONMODE_MATRIX_SPECTRALMODEL || geo.ioniz_mode == IONMODE_MATRIX_ESTIMATORS)
  {
    j_bar = mean_intensity_from_models (xplasma, freq, mode);
  }
  else
  {
    j_bar = mean_intensity_bb_estimate (freq, xplasma->t_r, xplasma->w);
  }

  return j_bar;
}

/**********************************************************/
/**
 * @brief  Get the mean intensity using the power law descriptions of the
 *         radiation field
 *
 * @param[in] PlasmaPtr xplasma  The specific plasma cell in question
 * @param[in] double freq        The frequency to calculate j_bar at
 * @param[in] int mode           Indicates whether to use BB estimator as a
 *                               fallback or not, use mode == 1 ==
 *                               MEAN_INTENSITY_BB_MODEL for this.
 *
 * @return  double j_bar  The mean intensity at a specific frequency
 *
 * @details
 *
 **********************************************************/

double
mean_intensity_from_models (PlasmaPtr xplasma, double freq, int mode)
{
  int i;
  double j_bar = 0.0;

  if (geo.spec_mod)
  {
    /*
     * Thus we have some spectral models, since we have run some ionization
     * cycles.
     */

    for (i = 0; i < geo.nxfreq; i++)
    {
      // Check that the band has the correct frequency range
      if (geo.xfreq[i] < freq && freq <= geo.xfreq[i + 1])
      {
        // Check we have a model for this band
        if (xplasma->spec_mod_type[i] > 0)
        {
          // Check the spectral model is defined for the frequency in question
          // todo: this seems redundant, but possibly can be important if the model
          //       in question isn't complete due to photon statistics
          if (freq > xplasma->fmin_mod[i] && freq < xplasma->fmax_mod[i])
          {
            if (xplasma->spec_mod_type[i] == SPEC_MOD_PL)
            {
              j_bar = pow (10, (xplasma->pl_log_w[i] + log10 (freq) * xplasma->pl_alpha[i]));
            }
            else if (xplasma->spec_mod_type[i] == SPEC_MOD_EXP)
            {
              j_bar = xplasma->exp_w[i] * exp ((-1 * PLANCK * freq) / (BOLTZMANN * xplasma->exp_temp[i]));
            }
            else
            {
              Error ("mean_intensity: unknown spectral model (%i) in band %i\n", xplasma->spec_mod_type[i], i);
              j_bar = 0.0;
            }
          }
          else
          {
            /* We have a spectral model, but it doesnt apply to the frequency
               in question. clearly this is a slightly odd situation, where last
               time we didnt get a photon of this frequency, but this time we did.
               Still this should only happen in very sparse cells, so induced Compton
               is unlikely to be important in such cells. We generate a warning, just
               so we can see if this is happening a lot */
            j_bar = 0.0;
            nerr_Jmodel_wrong_freq++;
          }
        }
        else                    /* There is no model in this band - this should not happen very often  */
        {
          j_bar = 0.0;
          nerr_no_Jmodel++;
        }
      }
    }
  }
  else
  {
    /*
     * We have not run any ionization cycles, so we either return j_bar = 0
     * or we can use a BB estimate if that mode is enabled
     */

    if (mode == MEAN_INTENSITY_BB_MODEL)
    {
      j_bar = mean_intensity_bb_estimate (freq, xplasma->t_r, xplasma->w);
    }
    else
    {
      j_bar = 0.0;
    }
  }

  return j_bar;
}

/**********************************************************/
/**
 * @brief Calculate the mean intensity of the field at a specific frequency
 *        for a specific cell.
 *
 * @param[in] double freq  The frequency to calculate j_bar at
 * @param[in] double t_r   The (radiation) temperature of the black body
 * @param[in] double w     The blackbody dilution factor for the cell
 *
 * @return  double j_bar  The mean intensity at a specific frequency
 *
 * @details
 *
 * The temperature of the blackbody is the radiation temperature of the cell,
 * as such we also pass the dilution factor of the same cell.
 *
 **********************************************************/

double
mean_intensity_bb_estimate (double freq, double t_r, double w)
{
  double j_bar, exponent;

  exponent = PLANCK * freq / BOLTZMANN / t_r;
  j_bar = 2 * PLANCK * pow (freq, 3) / VLIGHT / VLIGHT;
  j_bar *= 1 / (exp (exponent) - 1);
  j_bar *= w;

  return j_bar;
}
