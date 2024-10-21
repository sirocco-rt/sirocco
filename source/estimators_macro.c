
/***********************************************************/
/** @file  estimators_macro.c
 * @author Stuart Sim, James Matthews
 * @date   January, 2018
 *
 * @brief  Routines for dealing with macro-atom and simple-atom 
 *         Monte Carlo estimators. Consult sections 3.3 and 3.4 of 
 *         Matthews Phd Thesis: 
 *         \htmlonly
 *         <a href="https://doi.org/10.5281/zenodo.1256805">
 *         <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1256805.svg" 
 *         alt="DOI"></a>
 *         \endhtmlonly
 *
 * #Notes
 * This file includes both the routines to accumulate the data 
 * for creating the estimator during ionization cycles, but
 * also the routine mormalize_macro_estimators, which does
 * the normalization
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/* A couple of external things for use in the routines for computing gamma's below. */
struct topbase_phot *cont_ext_ptr2;     //continuum pointer passed externally
double temp_ext2;               //temperature passed externally
double temp_ext_rad;            //radiation temperature passed externally 

#define ALPHA_SP_CONSTANT 5.79618e-36   //


/**********************************************************/
/**
 * @brief increment the matom bound-free estimators
 *
 * @param [in] WindPtr  one pointer to cell
 * @param [in] PhotPtr  p the packet
 * @param [in] double  ds the path length
 * @return 0
 *
 * @details
 * increment the bound-free estimators as needed for the macro atom calculations. Whenever a packet
 * travels a path length ds through a shell the estimator for all bound-free
 * processes where the packet frequency is above the threshold frequency
 * is incremented. The estimator is not normalised here. That is done later.
 *
 * ### Notes ###
 * The estimators are not normalised by this routine.
 *
 * This routine also computes the total heating rate due
 * to simple ions.  Therefore the routine update_banded_esitmators is called
 * even for the simple non-macro atom case.
 * 
 * The bound-free estimators are described in section 3.3.3.1 Matthews Phd Thesis.
 *
 * bf_esimators_increment expects p
 * to be a co-moving frame photon and ds to be a co-moving frame path length
 * 
 **********************************************************/

int
bf_estimators_increment (one, p, ds)
     WindPtr one;
     PhotPtr p;
     double ds;

{
  double freq_av;
  double weight_of_packet;
  double x, ft;
  double y, yy;
  double exponential, heat_contribution;
  int n, m, llvl, nn;
  double density;
  double abs_cont;
  int nplasma, ndom;
  PlasmaPtr xplasma;
  MacroPtr mplasma;

  llvl = 0;
  density = 0.0;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  mplasma = &macromain[xplasma->nplasma];
  ndom = one->ndom;

  freq_av = p->freq;

  if (p->freq > xplasma->max_freq)
    xplasma->max_freq = p->freq;

  if (modes.save_cell_stats && ncstat > 0)
  {
    save_photon_stats (one, p, ds, p->w);
  }


  update_banded_estimators (xplasma, p, ds, p->w, ndom);
  update_flux_estimators (xplasma, p, ds, p->w, ndom);


  /* check that j and ave freq give sensible numbers */
  if (sane_check (xplasma->j) || sane_check (xplasma->ave_freq))
  {
    Error ("bf_estimators_increment:sane_check Problem with j %g or ave_freq %g\n", xplasma->j, xplasma->ave_freq);
  }



  for (nn = 0; nn < xplasma->kbf_nuse; nn++)
  {
    n = xplasma->kbf_use[nn];
    ft = phot_top[n].freq[0];   //This is the edge frequency (SS)

    if (ion[phot_top[n].nion].phot_info > 0)    //topbase or hybrid
    {
      llvl = phot_top[n].nlev;  //Returning lower level = correct (SS)
      density = den_config (xplasma, llvl);
    }
    else if (ion[phot_top[n].nion].phot_info == 0)      //vfky
    {
      density = xplasma->density[phot_top[n].nion];
      llvl = 0;                 // shouldn't ever be used 
    }

    if (kap_bf[nn] > 0.0 && (freq_av > ft))     // does the photon cause bf heating?
    {

      if (phot_top[n].macro_info == TRUE && geo.macro_simple == FALSE)  // it is a macro atom
      {

        x = kap_bf[nn] / (density * zdom[ndom].fill);   //this is the cross section

        /* Now identify which of the BF processes from this level this is. */

        m = 0;
        while (m < xconfig[llvl].n_bfu_jump && xconfig[llvl].bfu_jump[m] != n)
          m++;

        // m should now be the label to identify which of the bf processes from llvl
        // this is. Check that it is reasonable

        if (m > xconfig[llvl].n_bfu_jump - 1)
        {
          Error ("bf_estimators_increment: could not identify bf transition. Abort. \n");
          Exit (0);
        }

        // Now calculate the contributions and add them on.
        weight_of_packet = p->w;
        y = weight_of_packet * x * ds;

        exponential = y * exp (-(freq_av - ft) / BOLTZMANN / xplasma->t_e);

        /* Increment the photoionization rate estimator */

        mplasma->gamma[xconfig[llvl].bfu_indx_first + m] += y / freq_av;

        mplasma->alpha_st[xconfig[llvl].bfu_indx_first + m] += exponential / freq_av;

        mplasma->gamma_e[xconfig[llvl].bfu_indx_first + m] += y / ft;

        mplasma->alpha_st_e[xconfig[llvl].bfu_indx_first + m] += exponential / ft;

        /* Now record the contribution to the energy absorbed by macro atoms allowing
           for the filling factor. */

        yy = y * den_config (xplasma, llvl) * zdom[ndom].fill;

        mplasma->matom_abs[phot_top[n].uplev] += abs_cont = yy * ft / freq_av;

        xplasma->kpkt_abs += yy - abs_cont;

        /* Check for packets that appear to travelling  
           suspiciously large optical depth in the continuum */

        if ((yy / weight_of_packet) > 50)
        {
          Log ("bf_estimator_increment: A packet (%d) in cell %d (%d) survived an optical depth of %g\n", p->np, p->grid, wmain[p->grid],
               yy / weight_of_packet);
          Log ("bf_estimator_increment: freq_av %g, ft %g   yy %g  weight %g\n", freq_av, ft, yy, weight_of_packet);
        }
      }

      else                      // it is a simple ion
      {
        /* Now we are dealing with the heating due to the bf continua of simple ions. No stimulated
           recombination is included here. (SS, Apr 04) */
        if (density > DENSITY_PHOT_MIN)
        {
          x = sigma_phot (&phot_top[n], freq_av);

          weight_of_packet = p->w;
          y = weight_of_packet * x * ds;


          /* JM1411 -- added filling factor - density enhancement cancels with zdom[ndom].fill */
          xplasma->heat_photo += heat_contribution = y * density * (1.0 - (ft / freq_av)) * zdom[ndom].fill;

          xplasma->heat_tot += heat_contribution;
          /* This heat contribution is also the contibution to making k-packets in this volume. So we record it. */

          xplasma->kpkt_abs += heat_contribution;
        }
      }
    }
  }

  /* JM1411 -- the below processes have the factor zdom[ndom].fill incorporated directly
     into the kappa subroutines */

  if (xplasma->nwind < 0)
  {
    Error ("bf_estimators: Bad plasma cell %d  %d\n", xplasma->nplasma, xplasma->nwind);
    return (0);
  }

  /* Now for contribution to heating due to ff processes. (SS, Apr 04) */

  weight_of_packet = p->w;

  y = weight_of_packet * kappa_ff (xplasma, freq_av) * ds;

  xplasma->heat_ff += heat_contribution = y;    // record ff hea   

  /* This heat contribution is also the contibution to making k-packets in this volume. So we record it. */
  /* JM 2402 note that previously we incorrectly included Compton processes in kpkt_abs, which could lead to large 
     amounts of radiation coming out incorrectly in other k->r channels in spectral cycles */
  xplasma->kpkt_abs += heat_contribution;
     


  /* Now for contribution to heating due to compton processes. (JM, Sep 013) */

  y = weight_of_packet * kappa_comp (xplasma, freq_av) * ds;

  xplasma->heat_comp += y;      // record the compton heating
  heat_contribution += y;       // add compton to the heat contribution


  /* Now for contribution to heating due to induced compton processes. (JM, Sep 013) */

  y = weight_of_packet * kappa_ind_comp (xplasma, freq_av) * ds;

  xplasma->heat_ind_comp += y;  // record the induced compton heating
  heat_contribution += y;       // add induced compton to the heat contribution



  xplasma->heat_tot += heat_contribution;       // heat contribution is the contribution from compton, ind comp and ff processes

  return (0);
}


/**********************************************************/
/**
 * @brief increment the matom bound-bound estimators 
 *
 * @param [in] WindPtr  one   pointer to cell
 * @param [in] PhotPtr  p   the packet
 * @param [in] double  tau_sobolev   optical depth of line
 * @param [in] double  dvds   velocity gradient
 * @param [in] int  nn   the label for the line in question
 * @return   0 on success
 *           1 if the line was a "simple" line.  The routine should not have
 * 	           been called if that was the case.
 *
 * @details
*   increment the matom bound-bound estimators as needed for the macro atom calculations. Whenever a packet
 *  comes into resonance with a line this routine is called WHETHER OR NOT
 *  the packet interacts with the line.
 *  The estimator is not normalised here. That is done later.
 *
 * ### Notes ###
 * The estimators are not normalised by this routine
 *
 **********************************************************/

int
bb_estimators_increment (one, p, tau_sobolev, dvds, nn)
     WindPtr one;
     PhotPtr p;
     double tau_sobolev;
     double dvds;
     int nn;

{
  int llvl;
  int n;
  struct lines *line_ptr;
  double weight_of_packet;
  double y;
  int nmax;

  PlasmaPtr xplasma;
  MacroPtr mplasma;

  xplasma = &plasmamain[one->nplasma];
  mplasma = &macromain[xplasma->nplasma];


  line_ptr = lin_ptr[nn];

  if (line_ptr->macro_info == FALSE || geo.macro_simple == TRUE)
    return (-1);

  /* Start by identifying which estimator we want. Get lower level of
     transition . */

  llvl = line_ptr->nconfigl;

  /* Now find which of the bb transition from level llvl this is. */


  nmax = xconfig[llvl].n_bbu_jump;

  n = 0;
  while (n < nmax && line[xconfig[llvl].bbu_jump[n]].where_in_list != nn)
    n++;

  if (n == nmax)
  {
    Error ("bb_estimators_increment: could not identify bb transition. Abort. \n");
    Exit (0);
  }


  weight_of_packet = p->w;
  dvds = fabs (dvds);

  if (tau_sobolev > 0.00001)
  {
    y = weight_of_packet * (1. - exp (-tau_sobolev)) / tau_sobolev / dvds;
  }
  else
  {
    y = weight_of_packet / dvds;
  }

  if (y >= 0)
  {
    mplasma->jbar[xconfig[llvl].bbu_indx_first + n] += y;
  }
  else
  {
    Error ("bb_estimators_increment: trying to add negative contribution to jbar. %e. See #436\n", y);
    return (0);
  }

  /* Record contribution to energy absorbed by macro atoms. */

  mplasma->matom_abs[line_ptr->nconfigu] += weight_of_packet * (1. - exp (-tau_sobolev));

  return (0);
}

/**********************************************************/
/**
 * @brief the routine to normalise the matom bb and bf estimators. 
 *
 * @param [in, out] PlasmaPtr xplasma  The plasma (and associated macro) cell to normalise
 *
 * @return 0
 *
 * @details
*  This routine normalises the bb and bf mc estimators needed
 * for the macro atom jumping probabilities. During the main mc
 * simulation these were stored unnormalised. They are now
 * normalised and moved into the "old" slots for use in the next
 * iteration. The current slots are reset to zero for use in the
 * next iteration.
 *
 * It is performed as part of the wind update stage
 * of Python. The estimators should only be changed during the ionisation
 * cycles - after that they should be fixed.
 * It also now computes the bf heating rate using the normalised
 * estimators. heat_tot and heat_photo are incremented but nothing else
 * is for now.
 *
 * ### Notes ###
 * 
 * ksl -- This routine loops over nlte_levels, which in principle
 * could include non-macro ions, but that should not matter 
 * since nbfu_jumps will be zero for these.
 *
 **********************************************************/

int
normalise_macro_estimators (PlasmaPtr xplasma)
{
  double invariant_volume_time;
  int i, j, nlev_upper;
  double stimfac, line_freq, stat_weight_ratio;
  double heat_contribution, lower_density, upper_density;
  WindPtr one;
  MacroPtr mplasma;

  mplasma = &macromain[xplasma->nplasma];
  one = &wmain[xplasma->nwind];

  /* one->vol contains the CMF volume. This converts it to
     observer frame volume, or perhaps more correctly, it calculates the
     Lorentz invariant quantity (Delta V Delta t). Because Delta t_obs == 1
     in the code, (Delta V Delta t) == (Delta V_obs * 1.0).
     Note, this is also equivalent to multiplying Delta V_cmf by Delta t_cmf,
     which is Delta t_obs / gamma, i.e. 1.0/gamma. All other quantities that
     have been incremented should have been done so with CMF values.
   */

  invariant_volume_time = one->vol / one->xgamma_cen;

  /* bf estimators. During the mc calculation the quantity stored
     was weight * cross-section * path-length / frequency.
     To normalise this we need to put in:
     1 / h  (Planck constant)
     1 / Volume
     1 / Time
     I think that the weights are chosen such that Time = 1.
     So the normalisation of gamma and gamma_e is easy.
   */

  /* For the alpha_st estimators (stimulated recombination) the
     estimator also needs to be multiplied by a temperature
     dependent factor which is almost, but not quite, given by
     the LTE population ratio. The multiplicative factor
     is given by: */

  stimfac = 0.5 * pow (PLANCK * PLANCK / 2. / PI / MELEC / BOLTZMANN / xplasma->t_e, 3. / 2.);

  for (i = 0; i < nlte_levels; i++)
  {
    for (j = 0; j < xconfig[i].n_bfu_jump; j++)
    {
      mplasma->gamma_old[xconfig[i].bfu_indx_first + j] = mplasma->gamma[xconfig[i].bfu_indx_first + j] / PLANCK / invariant_volume_time;       // normalize by invariant_volume_time

      mplasma->gamma_e_old[xconfig[i].bfu_indx_first + j] =
        mplasma->gamma_e[xconfig[i].bfu_indx_first + j] / PLANCK / invariant_volume_time;

      mplasma->gamma[xconfig[i].bfu_indx_first + j] = 0.0;      //re-initialise for next iteration
      mplasma->gamma_e[xconfig[i].bfu_indx_first + j] = 0.0;

      /* For the stimulated recombination parts we need the the
         ratio of statistical weights too. For free electron statistical
         weight = 2 is included in stimfac above. */

      stat_weight_ratio = xconfig[phot_top[xconfig[i].bfu_jump[j]].uplev].g / xconfig[i].g;

      mplasma->alpha_st_old[xconfig[i].bfu_indx_first + j] =
        mplasma->alpha_st[xconfig[i].bfu_indx_first + j] * stimfac * stat_weight_ratio / PLANCK / invariant_volume_time;
      mplasma->alpha_st[xconfig[i].bfu_indx_first + j] = 0.0;

      mplasma->alpha_st_e_old[xconfig[i].bfu_indx_first + j] =
        mplasma->alpha_st_e[xconfig[i].bfu_indx_first + j] * stimfac * stat_weight_ratio / PLANCK / invariant_volume_time;
      mplasma->alpha_st_e[xconfig[i].bfu_indx_first + j] = 0.0;

      /* For continuua whose edges lie beyond freqmin assume that gamma
         is given by a black body. */

      /* For now place the limit at 7.5e12 which is 400000AA */
      /* Try also doing it for very high energy ones - greater than 50eV: 1.2e16 since up there the statistics of the estimators are very poor at the moment.
         Ideally we don't want to have this so should probably switch this back sometime (SS August 05) !!!BUG */

      if (phot_top[xconfig[i].bfu_jump[j]].freq[0] < 7.5e12 || phot_top[xconfig[i].bfu_jump[j]].freq[0] > 5e18)
      {
        mplasma->gamma_old[xconfig[i].bfu_indx_first + j] = get_gamma (&phot_top[xconfig[i].bfu_jump[j]], xplasma);
        mplasma->gamma_e_old[xconfig[i].bfu_indx_first + j] = get_gamma_e (&phot_top[xconfig[i].bfu_jump[j]], xplasma);
        mplasma->alpha_st_e_old[xconfig[i].bfu_indx_first + j] = get_alpha_st_e (&phot_top[xconfig[i].bfu_jump[j]], xplasma);
        mplasma->alpha_st_old[xconfig[i].bfu_indx_first + j] = get_alpha_st (&phot_top[xconfig[i].bfu_jump[j]], xplasma);
      }
    }

    /* That deals with the bf jumps. Now need to sort out the bb jumps. */

    /* For bb jumps the normalisation requires:
       1/ 4 PI
       1/ Volume
       1/ Time
       c / nu  - to convert dvds^-1 to dnuds^-1

       Also, I'm putting the correction factor for stimulated emission in here - so that it doesn't have to be
       computed in the macro atom jumping probabilities.
     */

    for (j = 0; j < xconfig[i].n_bbu_jump; j++)
    {
      nlev_upper = line[xconfig[i].bbu_jump[j]].nconfigu;
      lower_density = den_config (xplasma, i);
      upper_density = den_config (xplasma, nlev_upper);

      /* The correction for stimulated emission is (1 - n_lower * g_upper / n_upper / g_lower) */

      stimfac = upper_density / lower_density;
      stimfac = stimfac * xconfig[i].g / xconfig[line[xconfig[i].bbu_jump[j]].nconfigu].g;

      if (stimfac < 1.0 && stimfac >= 0.0)
      {
        stimfac = 1. - stimfac;
      }
      else if (upper_density > DENSITY_PHOT_MIN && lower_density > DENSITY_PHOT_MIN
               && xplasma->levden[xconfig[nlev_upper].nden] > DENSITY_MIN)
      {
        /* check for population inversions. We don't worry about this if the densities are extremely low or if the
           upper level has hit the density floor - the lower level is still allowed to hit this floor because it
           should never cause an inversion */
        Error ("normalise_macro_estimators: bb stimulated correction factor is out of bounds, 0 <= stimfac < 1 but got %g\n", stimfac);
        Error ("normalise_macro_estimators: upper_density %g lower_density %g xplasma->levden[config[nlev_upper].nden] %g\n",
               upper_density, lower_density, xplasma->levden[xconfig[nlev_upper].nden]);
        stimfac = 0.0;
      }
      else
      {
        stimfac = 0.0;
      }

      /* normalise jbar. Note that this uses the cell volume rather than the filled volume */

      line_freq = line[xconfig[i].bbu_jump[j]].freq;
      mplasma->jbar_old[xconfig[i].bbu_indx_first + j] =
        mplasma->jbar[xconfig[i].bbu_indx_first + j] * VLIGHT * stimfac / 4. / PI / invariant_volume_time / line_freq;
      mplasma->jbar[xconfig[i].bbu_indx_first + j] = 0.0;
    }
  }

  /* bb and bf now normalised. Done. */
  /* Get the heating contribution from macro atom bb transitions (the
     line heating). */

  xplasma->heat_lines += heat_contribution = macro_bb_heating (xplasma, xplasma->t_e);
  xplasma->heat_lines_macro = heat_contribution;
  xplasma->heat_tot += heat_contribution;

  /* Get the bf heating contributions here too. (SS June 04) */

  xplasma->heat_photo += heat_contribution = macro_bf_heating (xplasma, xplasma->t_e);
  xplasma->heat_photo_macro = heat_contribution;
  xplasma->heat_tot += heat_contribution;

  /* finally, check if we have any places where stimulated recombination wins over
     photoionization */

  check_stimulated_recomb (xplasma);

  /* Now that we have estimators, set the flag to use them for the level populations */

  geo.macro_ioniz_mode = MACRO_IONIZ_MODE_ESTIMATORS;

  /* force recalculation of k-packet rates and matrices, if applicable */
  mplasma->kpkt_rates_known = FALSE;
  mplasma->matrix_rates_known = FALSE;

  return (0);
}


/**********************************************************/
/**
 * @brief      computes the cooling rate due to free-bound recombinations
 *
 * @param [in] PlasmaPtr  xplasma 
 * @param [in] double  t_e  electron temperature
 * @param [in] double  f1  lower frequency
 * @param [in] double  f2  upper frequency
 * @return total
 *
 * @details
 * computes the cooling rate due to free-bound recombinations 
 * using MC estimators. It is modelled on total_fb
 * but makes use of the mc estimators for stimulated recombination.
 *
 * ### Notes ###
 * 
 * This routine loops over nlte_levels, which in principle
 * could include non-macro ions, but that should not matter since
 * since nbfu_jumps will be zero for these.
 *
 * free bound emission in this function is calculated CMF. 
 * 
 **********************************************************/

double
total_fb_matoms (xplasma, t_e, f1, f2)
     PlasmaPtr xplasma;
     double t_e;
     double f1, f2;
{
  double cool_contribution;
  double t_e_store;
  struct topbase_phot *cont_ptr;
  double total, density;
  int i, j;
  double q_ioniz ();
  MacroPtr mplasma;

  mplasma = &macromain[xplasma->nplasma];

  t_e_store = xplasma->t_e;     //store the temperature - will put it back at the end
  xplasma->t_e = t_e;           //for use in calls to alpha_sp below

  total = 0;

  if (geo.macro_simple == FALSE)        //allow for "only-simple" calculations (SS May04)
  {
    for (i = 0; i < nlte_levels; i++)
    {
      for (j = 0; j < xconfig[i].n_bfu_jump; j++)
      {
        /* Need the density for the upper level in the recombination
           process. */
        cont_ptr = &phot_top[xconfig[i].bfu_jump[j]];
        density = den_config (xplasma, cont_ptr->uplev);

        /* the cooling contribution for each transition is given by 
           density * ne * [(alpha_sp_e - alpha_sp) + (alpha_st_e - alpha_st)]
           we call alpha_sp() with modes 1 and 0 to get the different
           versions of the sp. recombination rate coefficient.
           This is essentially equation (33) of Lucy (2003) */

        cool_contribution =
          (mplasma->alpha_st_e_old[xconfig[i].bfu_indx_first + j] +
           alpha_sp (cont_ptr, xplasma, 1)
           - mplasma->alpha_st_old[xconfig[i].bfu_indx_first + j]
           - alpha_sp (cont_ptr, xplasma, 0)) * PLANCK * phot_top[xconfig[i].bfu_jump[j]].freq[0] * density * xplasma->ne * xplasma->vol;

        /* Now add the collisional ionization term. */
        density = den_config (xplasma, cont_ptr->nlev);
        cool_contribution +=
          q_ioniz (cont_ptr, t_e) * density * xplasma->ne * PLANCK * phot_top[xconfig[i].bfu_jump[j]].freq[0] * xplasma->vol;

        /* That's the bf cooling contribution. */
        total += cool_contribution;
      }
    }

    xplasma->t_e = t_e_store;   //restore the original value

  }

  return (total);
}

/**********************************************************/
/**
 * @brief      computes the total cooling in matom bb transitions
 *
 * @param [in] PlasmaPtr xplasma 
 * @param [in] doublet_e electron temperature
 * @return total
 *
 * @details
 * computes the total cooling in bb transitions
 * (due to collisions) for both macro atoms and simple ions. 
 *  It is used by the heating/cooling calculation to get the temperature.
 *
 * Total bound-bound emission in this function is calculated CMF. If called in
 * the cooling routines then this is left in the CMF. If called in the photon
 * generation routines, this is converted to the observer frame.
 *
 **********************************************************/

double
total_bb_cooling (xplasma, t_e)
     PlasmaPtr xplasma;
     double t_e;
{
  double cool_contribution;
  struct lines *line_ptr;
  double total, upper_density, lower_density;
  int i;
  double coll_rate, rad_rate;

  total = 0;                    // initialise
  for (i = 0; i < nlines; i++)
  {
    line_ptr = &line[i];
    if (line_ptr->macro_info == TRUE && geo.macro_simple == FALSE)
    {                           //This is a line from a macro atom for which we know
      //the upper and lower level populations
      lower_density = den_config (xplasma, line_ptr->nconfigl);
      cool_contribution = (lower_density * q12 (line_ptr, t_e)) * xplasma->ne * xplasma->vol * line_ptr->freq * PLANCK;
    }
    else
    {                           //It's a simple line - don't know the level populations
      // - just use a two-level-atom approach

      //The cooling rate is computed using the scattering probability formalism in KSL's notes on Python.

      two_level_atom (line_ptr, xplasma, &lower_density, &upper_density);
      coll_rate = q21 (line_ptr, t_e) * xplasma->ne * (1. - exp (-H_OVER_K * line_ptr->freq / t_e));

      cool_contribution =
        (lower_density * line_ptr->gu / line_ptr->gl -
         upper_density) * coll_rate / (exp (H_OVER_K * line_ptr->freq / t_e) - 1.) * xplasma->vol * line_ptr->freq * PLANCK;


      rad_rate = a21 (line_ptr) * p_escape (line_ptr, xplasma);
      cool_contribution *= rad_rate / (rad_rate + coll_rate);
    }

    /* That's the bb cooling contribution. */

    total += cool_contribution;

  }

  return (total);
}

/**********************************************************/
/**
 * @brief computes the total heating due to matom bb transitions
 *
 * @param [in] PlasmaPtr  xplasma  
 * @param [in] double  t_e electron temperature
 * @return total
 *
 * @details
 * computes the total heating due to bb transitions
 * (due to collisions) for macro atoms. The heating in simple ions
 * is taken care of elsewhere. 
 * It is used by the heating/cooling calculation to get the temperature.
 *
 **********************************************************/

double
macro_bb_heating (xplasma, t_e)
     PlasmaPtr xplasma;
     double t_e;
{
  double heat_contribution;
  struct lines *line_ptr;
  double total, upper_density;
  int i;


  total = 0;                    // initialise

  for (i = 0; i < nlines; i++)
  {
    line_ptr = &line[i];
    if (line_ptr->macro_info == TRUE && geo.macro_simple == FALSE)
    {                           //This is a line from a macro atom for which we know
      //the upper and lower level populations
      upper_density = den_config (xplasma, line_ptr->nconfigu);
      heat_contribution = upper_density * q21 (line_ptr, t_e) * xplasma->ne * xplasma->vol * line_ptr->freq * PLANCK;
      total += heat_contribution;
    }
  }

  return (total);
}


/**********************************************************/
/**
 * @brief computes the total heating due to bf transitions for macro atoms.
 *
 * @param [in] PlasmaPtr  xplasma  
 * @param [in] double t_e   electron temperature
 * @return total 
 *
 * @details
 * computes the total heating due to bf transitions for macro atoms. 
 * The heating in simple ions
 * is taken care of elsewhere. 
 * It is used by the heating/cooling calculation to get the temperature.
 *
 **********************************************************/

double
macro_bf_heating (xplasma, t_e)
     PlasmaPtr xplasma;
     double t_e;
{
  double heat_contribution;
  double total, upper_density, lower_density;
  int i, j;
  double q_recomb ();
  MacroPtr mplasma;

  mplasma = &macromain[xplasma->nplasma];

  total = 0;                    // initialise

  for (i = 0; i < nlte_levels; i++)
  {
    for (j = 0; j < xconfig[i].n_bfu_jump; j++)
    {
      heat_contribution = 0.0;
      /* Photoionization part. */
      lower_density = den_config (xplasma, phot_top[xconfig[i].bfu_jump[j]].nlev);
      heat_contribution =
        (mplasma->gamma_e_old[xconfig[i].bfu_indx_first + j] -
         mplasma->gamma_old[xconfig[i].bfu_indx_first +
                            j]) * PLANCK * phot_top[xconfig[i].bfu_jump[j]].freq[0] * lower_density * xplasma->vol;

      /* Three body recombination part. */
      upper_density = den_config (xplasma, phot_top[xconfig[i].bfu_jump[j]].uplev);
      heat_contribution +=
        q_recomb (&phot_top[xconfig[i].bfu_jump[j]],
                  t_e) * xplasma->ne * xplasma->ne * PLANCK * upper_density * xplasma->vol * phot_top[xconfig[i].bfu_jump[j]].freq[0];

      total += heat_contribution;

    }
  }

  return (total);
}


/**********************************************************/
/**
 * @brief  records the heating contribution from lines of simple elements. 
 *
 * @param [out] PlasmaPtrxplasma 
 * @param [in] PhotPtr  p   the packet
 * @param [in] double  tau_sobolev   optical depth of line
 * @param [in] int  nn   the label for the line in question
 * @return  0 on success
 *
 * @details
 * records the heating contribution from lines of simple elements. It is called whenever a packet comes into 
 * resonance with the line and records the heating contribution 
 * from that line and packet
 *
 * #Notes
 *
 * The heating contribution is modelled on the macro atom bb estimator 
 * calculations for the radiative excitation rate. This (energy) excitation
 * rate is multiplied by the destruction probability to get the heating. 
 * The destruction probability is obtained following the discussion in KSL's notes
 * on Python.
 *
 **********************************************************/

int
bb_simple_heat (xplasma, p, tau_sobolev, nn)
     PlasmaPtr xplasma;
     PhotPtr p;
     double tau_sobolev;
     int nn;

{
  double heat_contribution;
  double weight_of_packet;
  struct lines *line_ptr;
  double electron_temperature;
  double rad_rate, coll_rate, normalisation;


  weight_of_packet = p->w;
  line_ptr = lin_ptr[nn];
  electron_temperature = xplasma->t_e;

  rad_rate = a21 (line_ptr) * p_escape (line_ptr, xplasma);

  coll_rate = q21 (line_ptr, electron_temperature) * xplasma->ne * (1. - exp (-H_OVER_K * line_ptr->freq / electron_temperature));

  normalisation = rad_rate + coll_rate;


  /* Now add the heating contribution. */

  xplasma->heat_lines += heat_contribution = weight_of_packet * (coll_rate / normalisation) * (1. - exp (-1. * tau_sobolev));

  xplasma->heat_tot += heat_contribution;
  xplasma->kpkt_abs += heat_contribution;

  return (0);

}



/**********************************************************/
/**
 * @brief sanity check the stimulated recombination rate
 *
 * @param [in] PlasmaPtr  xplasma 
 * @return 0
 *
 * @details check that stimulated recombination rate doesn't exceed the 
 * photoionization rate (check if any levels have 
 * \f$\alpha_{st}n_e*\frac{n_u}{n_l} > \gamma\f$
 *
 **********************************************************/

int
check_stimulated_recomb (xplasma)
     PlasmaPtr xplasma;
{
  int i, j;
  struct topbase_phot *cont_ptr;
  int st_recomb_err;
  double st_recomb, gamma, coll_ioniz;
  MacroPtr mplasma;
  mplasma = &macromain[xplasma->nplasma];

  st_recomb_err = 0;

  for (i = 0; i < nlte_levels; i++)
  {
    for (j = 0; j < xconfig[i].n_bfu_jump; j++)
    {
      cont_ptr = &phot_top[xconfig[i].bfu_jump[j]];
      gamma = mplasma->gamma_old[xconfig[i].bfu_indx_first + j];
      st_recomb = mplasma->alpha_st_old[xconfig[i].bfu_indx_first + j];
      st_recomb *= xplasma->ne * den_config (xplasma, cont_ptr->uplev) / den_config (xplasma, cont_ptr->nlev);
      coll_ioniz = q_ioniz (cont_ptr, xplasma->t_e) * xplasma->ne;

      if (st_recomb > (gamma + coll_ioniz))
        st_recomb_err += 1;
    }
  }

  if (st_recomb_err > 0)
    Error ("check_stimulated_recomb: cell %i had %i bf jumps where ne*n_u/n_l*alpha_st > gamma\n", xplasma->nplasma, st_recomb_err);

  return (0);
}

/**********************************************************/
/**
 * @brief compute dilute matom estimators
 *
 * @param [in] PlasmaPtr  xplasma 
 * @return 
 *
 * @details
 * get_dilute_estimators computes the bound free and 
 * bound bound estimators for a cell from a dilute blackbody.
 * This is the default if macro_pops fails to find a solution
 * for an ion. Modifies mplasma.
 *
 **********************************************************/

int
get_dilute_estimators (xplasma)
     PlasmaPtr xplasma;
{

  int i, j;
  struct lines *line_ptr;
  MacroPtr mplasma;
  mplasma = &macromain[xplasma->nplasma];

  for (i = 0; i < nlte_levels; i++)
  {
    for (j = 0; j < xconfig[i].n_bfu_jump; j++)
    {
      mplasma->gamma_old[xconfig[i].bfu_indx_first + j] = get_gamma (&phot_top[xconfig[i].bfu_jump[j]], xplasma);
      mplasma->gamma_e_old[xconfig[i].bfu_indx_first + j] = get_gamma_e (&phot_top[xconfig[i].bfu_jump[j]], xplasma);
      mplasma->alpha_st_e_old[xconfig[i].bfu_indx_first + j] = get_alpha_st_e (&phot_top[xconfig[i].bfu_jump[j]], xplasma);
      mplasma->alpha_st_old[xconfig[i].bfu_indx_first + j] = get_alpha_st (&phot_top[xconfig[i].bfu_jump[j]], xplasma);
    }
    for (j = 0; j < xconfig[i].n_bbu_jump; j++)
    {
      line_ptr = &line[xconfig[i].bbu_jump[j]];
      mplasma->jbar_old[xconfig[i].bbu_indx_first + j] = mean_intensity (xplasma, line_ptr->freq, MEAN_INTENSITY_BB_MODEL);
    }
  }

  return (0);
}

/**********************************************************/
/**
 * @brief the dilute photoionization rate estimator
 *
 * @param [in] struct topbase_phot *  cont_ptr   
 * @param [in] PlasmaPtr  xplasma   
 * @return gamma_value
 *
 * @details
 * the dilute photoionization rate estimator, for a black body 
 * radiation field with known dilution factor and temperature 
 * see equation 3.64 of Matthews Phd thesis, see also Lucy (2003). 
 *
 **********************************************************/

double
get_gamma (cont_ptr, xplasma)
     struct topbase_phot *cont_ptr;
     PlasmaPtr xplasma;
{
  double gamma_value;
  double fthresh, flast;
  double qromb ();
  double gamma_integrand ();

  temp_ext2 = xplasma->t_r;     //external temperature
  cont_ext_ptr2 = cont_ptr;     //external cont pointer
  fthresh = cont_ptr->freq[0];  //first frequency in list
  flast = cont_ptr->freq[cont_ptr->np - 1];     //last frequency in list
  if ((H_OVER_K * (flast - fthresh) / temp_ext2) > ALPHA_MATOM_NUMAX_LIMIT)
  {
    //flast is currently very far into the exponential tail: so reduce flast to limit value of h nu / k T.
    flast = fthresh + temp_ext2 * ALPHA_MATOM_NUMAX_LIMIT / H_OVER_K;
  }


  // gamma_value = qromb (gamma_integrand, fthresh, flast, 1e-4);
  gamma_value = num_int (gamma_integrand, fthresh, flast, 1e-4);

  gamma_value *= 8 * PI / VLIGHT / VLIGHT * xplasma->w;

  return (gamma_value);

}

/**********************************************************/
/**
 * @brief integrand for the dilute photoionization rate estimator
 *
 * @param [in] double freq
 * @param [in] void  params   An extra (unused) variable to make it paletable for the gsl integrator
 * @return integrand
 *
 * @details
 * integrand for the dilute photoionization rate estimator, for a black body 
 * radiation field with known dilution factor and temperature 
 * see equation 3.64 of Matthews Phd thesis, see also Lucy (2003). 
 *
 **********************************************************/

double
gamma_integrand (double freq, void *params)
{
  double fthresh;
  double x;
  double integrand;
  double tt;

  fthresh = cont_ext_ptr2->freq[0];
  tt = temp_ext2;

  if (freq < fthresh)
    return (0.0);               // No photoionization at frequencies lower than the threshold freq occur

  x = sigma_phot (cont_ext_ptr2, freq); //this is the cross-section

  integrand = x * freq * freq / (exp (H_OVER_K * freq / tt) - 1);

  return (integrand);
}

/**********************************************************/
/**
 * @brief the dilute energy weighted photoionization rate estimator
 *
 * @param [in, out] struct topbase_phot *  cont_ptr 
 * @param [in, out] PlasmaPtr  xplasma 
 * @return gamma_e_value
 *
 * @details
 * the energy-weighted photoionization rate estimator, for a black body 
 * radiation field with known dilution factor and temperature 
 * equation 3.66 of Matthews Phd thesis, see also Lucy (2003). 
 *
 **********************************************************/

double
get_gamma_e (cont_ptr, xplasma)
     struct topbase_phot *cont_ptr;
     PlasmaPtr xplasma;
{
  double gamma_e_value;
  double fthresh, flast;
  double qromb ();
  double gamma_e_integrand ();

  temp_ext2 = xplasma->t_r;     //external temperature
  cont_ext_ptr2 = cont_ptr;     //external cont pointer
  fthresh = cont_ptr->freq[0];  //first frequency in list
  flast = cont_ptr->freq[cont_ptr->np - 1];     //last frequency in list
  if ((H_OVER_K * (flast - fthresh) / temp_ext2) > ALPHA_MATOM_NUMAX_LIMIT)
  {
    //flast is currently very far into the exponential tail: so reduce flast to limit value of h nu / k T.
    flast = fthresh + temp_ext2 * ALPHA_MATOM_NUMAX_LIMIT / H_OVER_K;
  }


//  gamma_e_value = qromb (gamma_e_integrand, fthresh, flast, 1e-4);
  gamma_e_value = num_int (gamma_e_integrand, fthresh, flast, 1e-4);

  gamma_e_value *= 8 * PI / VLIGHT / VLIGHT * xplasma->w;

  return (gamma_e_value);

}

/**********************************************************/
/**
 * @brief the integrand for the energy-weighted photoionization rate estimator, 
 *
 * @param [in] double  freq 
 * @param [in] void  params   An extra (unused) variable to make it paletable for the gsl integrator
 * @return integrand 
 *
 * @details
 * the integrand for the energy-weighted photoionization rate estimator, 
 * equation 3.66 of Matthews Phd thesis, see also Lucy (2003). 
 *
 **********************************************************/

double
gamma_e_integrand (double freq, void *params)
{
  double fthresh;
  double x;
  double integrand;
  double tt;

  fthresh = cont_ext_ptr2->freq[0];
  tt = temp_ext2;

  if (freq < fthresh)
    return (0.0);               // No photoionization at frequencies lower than the threshold freq occur

  x = sigma_phot (cont_ext_ptr2, freq); //this is the cross-section
  integrand = x * freq * freq * freq / (exp (H_OVER_K * freq / tt) - 1) / fthresh;

  return (integrand);
}



/**********************************************************/
/**
 * @brief the stimulated recombination rate estimator
 *
 * @param [in] struct topbase_phot *  cont_ptr 
 * @param [in] PlasmaPtr  xplasma 
 * @return alpha_st_value
 *
 * @details
 * the stimulated recombination rate estimator, 
 * equation 3.65 of Matthews Phd thesis, see also Lucy (2003). 
 *
 **********************************************************/

double
get_alpha_st (cont_ptr, xplasma)
     struct topbase_phot *cont_ptr;
     PlasmaPtr xplasma;
{
  double alpha_st_value;
  double fthresh, flast;
  double qromb ();
  double alpha_st_integrand ();

  temp_ext2 = xplasma->t_e;     //external for use in integrand
  temp_ext_rad = xplasma->t_r;
  cont_ext_ptr2 = cont_ptr;     //"
  fthresh = cont_ptr->freq[0];  //first frequency in list
  flast = cont_ptr->freq[cont_ptr->np - 1];     //last frequency in list

  if ((H_OVER_K * (flast - fthresh) / temp_ext2) > ALPHA_MATOM_NUMAX_LIMIT)
  {
    //flast is currently very far into the exponential tail: so reduce flast to limit value of h nu / k T.
    flast = fthresh + temp_ext2 * ALPHA_MATOM_NUMAX_LIMIT / H_OVER_K;
  }


  alpha_st_value = num_int (alpha_st_integrand, fthresh, flast, 1e-4);


  /* The lines above evaluate the integral in alpha_sp. Now we just want to multiply 
     through by the appropriate constant. */
  if (cont_ptr->macro_info == 1 && geo.macro_simple == FALSE)
  {
    alpha_st_value = alpha_st_value * xconfig[cont_ptr->nlev].g / xconfig[cont_ptr->uplev].g * pow (xplasma->t_e, -1.5);
  }
  else                          //case for simple element
  {
    alpha_st_value = alpha_st_value * xconfig[cont_ptr->nlev].g / ion[cont_ptr->nion + 1].g * pow (xplasma->t_e, -1.5); //g for next ion up used
  }

  alpha_st_value = alpha_st_value * ALPHA_SP_CONSTANT * xplasma->w;

  return (alpha_st_value);
}

/**********************************************************/
/**
 * @brief returns the integrand for alpha_st at a chosen frequency
 *
 * @param [in] double  freq 
 * @param [in] void  params   An extra (unused) variable to make it paletable for the gsl integrator
 * @return integrand 
 *
 * @details
 * integrand for the stimulated recombination rate estimator, 
 * in equation 3.65 of Matthews Phd thesis, see also Lucy (2003). 
 *
 **********************************************************/

double
alpha_st_integrand (double freq, void *params)
{
  double fthresh;
  double x;
  double integrand;
  double tt;
  double ttrr;

  fthresh = cont_ext_ptr2->freq[0];
  tt = temp_ext2;               //this is the electron temperature
  /* Also need the radiation temperature here */
  ttrr = temp_ext_rad;          //will do for now

  if (freq < fthresh)
    return (0.0);               // No recombination at frequencies lower than the threshold freq occur

  x = sigma_phot (cont_ext_ptr2, freq); //this is the cross-section

  integrand = x * freq * freq * exp (H_OVER_K * (fthresh - freq) / tt) / (exp (H_OVER_K * freq / ttrr) - 1);

  return (integrand);
}




/**********************************************************/
/**
 * @brief the energy-weighted stimulated recombination estimator 
 *
 * @param [in] struct topbase_phot *  cont_ptr
 * @param [in] PlasmaPtr  xplasma 
 * @return alpha_st_e_value  
 *
 * @details
 * the energy-weighted stimulated recombination estimator 
 * see equation 3.68 of Matthews Phd thesis, see also Lucy (2003). 
 *
 **********************************************************/

double
get_alpha_st_e (cont_ptr, xplasma)
     struct topbase_phot *cont_ptr;
     PlasmaPtr xplasma;
{
  double alpha_st_e_value;
  double fthresh, flast;
  double qromb ();
  double alpha_st_e_integrand ();

  temp_ext2 = xplasma->t_e;     //external for use in integrand
  temp_ext_rad = xplasma->t_r;  //"
  cont_ext_ptr2 = cont_ptr;     //"
  fthresh = cont_ptr->freq[0];  //first frequency in list
  flast = cont_ptr->freq[cont_ptr->np - 1];     //last frequency in list

  if ((H_OVER_K * (flast - fthresh) / temp_ext2) > ALPHA_MATOM_NUMAX_LIMIT)
  {
    //flast is currently very far into the exponential tail: so reduce flast to limit value of h nu / k T.
    flast = fthresh + temp_ext2 * ALPHA_MATOM_NUMAX_LIMIT / H_OVER_K;
  }


//  alpha_st_e_value = qromb (alpha_st_e_integrand, fthresh, flast, 1e-4);
  alpha_st_e_value = num_int (alpha_st_e_integrand, fthresh, flast, 1e-4);

  /* The lines above evaluate the integral in alpha_sp. Now we just want to multiply 
     through by the appropriate constant. */
  if (cont_ptr->macro_info == TRUE && geo.macro_simple == FALSE)
  {
    alpha_st_e_value = alpha_st_e_value * xconfig[cont_ptr->nlev].g / xconfig[cont_ptr->uplev].g * pow (xplasma->t_e, -1.5);
  }
  else                          //case for simple element
  {
    alpha_st_e_value = alpha_st_e_value * xconfig[cont_ptr->nlev].g / ion[cont_ptr->nion + 1].g * pow (xplasma->t_e, -1.5);     //g for next ion up used
  }

  alpha_st_e_value = alpha_st_e_value * ALPHA_SP_CONSTANT * xplasma->w;

  return (alpha_st_e_value);
}


/**********************************************************/
/**
 * @brief the integrand for alpha_st_e at a chosen frequency
 *
 * @param [in, out] double freq 
 * @param [in] void  params   An extra (unused) variable to make it paletable for the gsl integrator
 * @return integrand 
 *
 * @details
 * the integrand for the energy-weighted stimulated recombination estimator 
 * see equation 3.68 of Matthews Phd thesis, see also Lucy (2003). 
 *
 **********************************************************/

double
alpha_st_e_integrand (double freq, void *params)
{
  double fthresh;
  double x;
  double integrand;
  double tt;
  double ttrr;

  fthresh = cont_ext_ptr2->freq[0];
  tt = temp_ext2;               //this is the electron temperature
  /* Also need the radiation temperature here */
  ttrr = temp_ext_rad;          //will do for now

  if (freq < fthresh)
    return (0.0);               // No recombination at frequencies lower than the threshold freq occur

  x = sigma_phot (cont_ext_ptr2, freq); //this is the cross-section  

  integrand = x * freq * freq * exp (H_OVER_K * (fthresh - freq) / tt) / (exp (H_OVER_K * freq / ttrr) - 1) * freq / fthresh;

  return (integrand);
}
