/***********************************************************/
/** @file  estimators_simple.c
 * @author ksl
 * @date   July, 2020   
 *
 * @brief  Routines related to updating estimators for the
 * case of simple atoms.  The routine is primarily intended
 * for simple atoms, but also may be called when a macro-atom
 * calculation is carried out.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"
#include <gsl/gsl_integration.h>


/**********************************************************/
/** 
 * @brief      updates the estimators required for determining crude
 * spectra in each Plasma cell
 *
 * @param [in,out] PlasmaPtr  xplasma   PlasmaPtr for the cell of interest
 * @param [in] PhotPtr  p   Photon pointer in CMF frame
 * @param [in] double  ds   ds travelled in CMF frame
 * @param [in] double  w_ave   the weight of the photon in CMF frame 
 *
 * @return  Always returns 0
 *
 *
 *
 * @details
 * 
 * Increments the estimators that allow one to construct a crude
 * spectrum in each cell of the wind.  The frequency intervals
 * in which the spectra are constructed are defined for each Plasama element. This information
 * is used in different ways (or not at all) depending on the ionization mode.
 *
 * It also records the various parameters intended to describe the radiation field, 
 * including the IP.
 *
 * ### Notes ###
 *
 * The term direct refers to photons that have not been scattered by the wind.
 * 
 * In non macro atom mode, w_ave
 * this is an average weight (passed as w_ave), but 
 * in macro atom mode weights are never reduced (so p->w 
 * is used).
 *
 * This routine is called from bf_estimators in macro_atom modes
 * and from radiation.  Although the historical documentation
 * suggests it is only called for certain ionization modes, it appears
 * to be called in all cases, though clearly it is only provides diagnostic
 * information in some of them.
 *
 * update_banded_estimators takes CMF inputs
 * The photon frequency and number are used.
 * 
 **********************************************************/


/* External variables to prevent double counting of photons
   photons coming into a cell, in total and if they 
   are ionizing photons
*/

int previous_nioniz_nplasma = -1;
int previous_nioniz_np = -1;
int previous_nplasma = -1;
int previous_np = -1;

int
update_banded_estimators (PlasmaPtr xplasma, PhotPtr p, double ds, double w_ave, int ndom)
{
  int i;
  double log_freq;

  if (ds < 0)
  {
    Error ("updated_banded_estimators: ds %e < 0 in plasma cell %d when it has to be positive\n", ds, xplasma->nplasma);
    Exit (1);
  }

  if (w_ave < 0)
  {
    Error ("updated_banded_estimators: w_ave %e < 0 in plasma cell %d when it has to be positive\n", w_ave, xplasma->nplasma);
    Exit (1);
  }

  /*photon weight times distance in the shell is proportional to the mean intensity */

  xplasma->est.j += w_ave * ds;

  if (p->nscat == 0)
  {
    xplasma->est.j_direct += w_ave * ds;
  }
  else
  {
    xplasma->est.j_scatt += w_ave * ds;
  }



/* frequency weighted by the weights and distance in the shell .  See eqn 2 ML93 */
  xplasma->est.mean_ds += ds;
  xplasma->est.n_ds++;
  xplasma->est.ave_freq += p->freq * w_ave * ds;

  /* The next loop updates the banded versions of j and ave_freq, analogously to routine inradiation
     nxfreq refers to how many frequencies we have defining the bands. So, if we have 5 bands, we have 6 frequencies, 
     going from xfreq[0] to xfreq[5] Since we are breaking out of the loop when i>=nxfreq, this means the last loop 
     will be i=nxfreq-1 */

  /* note that here we can use the photon weight and don't need to calculate anm attenuated average weight
     as energy packets are indisivible in macro atom mode */


  for (i = 0; i < xplasma->state.nbands; i++)
  {
    if (xplasma->state.f1[i] < p->freq && p->freq <= xplasma->state.f2[i])
    {
      xplasma->est.xave_freq[i] += p->freq * w_ave * ds;        /* frequency weighted by weight and distance */
      xplasma->est.xsd_freq[i] += p->freq * p->freq * w_ave * ds;       /* input to allow standard deviation to be calculated */
      xplasma->est.xj[i] += w_ave * ds; /* photon weight times distance travelled */
      xplasma->est.nxtot[i]++;  /* increment the frequency banded photon counter */
      /* work out the range of frequencies within a band where photons have been seen */
      if (p->freq < xplasma->est.fmin[i])
      {
        xplasma->est.fmin[i] = p->freq;
      }
      if (p->freq > xplasma->est.fmax[i])
      {
        xplasma->est.fmax[i] = p->freq;
      }

    }
  }

  /* Increment the cell spectra  */

  log_freq = log10 (p->freq);
  i = (log_freq - geo.cell_log_freq_min) / geo.cell_delta_lfreq;
  if (i < 0)
    i = 0;
  if (i > geo.nbins_in_cell_spec - 1)
  {
    i = geo.nbins_in_cell_spec - 1;
  }
  xplasma->est.cell_spec_flux[i] += w_ave * ds;



  /* Increment counters which show the number of times a photon bundle of a certain type
     has passed through the cell.  Makde sure the photon is not counted multiple
     times in a single passage */

  if (xplasma->nplasma != previous_nplasma || p->np != previous_np)
  {
    xplasma->est.ntot++;

    if (p->origin == PTYPE_STAR)
      xplasma->est.ntot_star++;
    else if (p->origin == PTYPE_BL)
      xplasma->est.ntot_bl++;
    else if (p->origin == PTYPE_DISK)
      xplasma->est.ntot_disk++;
    else if (p->origin == PTYPE_WIND)
      xplasma->est.ntot_wind++;
    else if (p->origin == PTYPE_AGN)
      xplasma->est.ntot_agn++;
    previous_nplasma = xplasma->nplasma;
    previous_np = p->np;
  }


  /*
   * Calculate the number of H ionizing photons, see #255
   * EP 11-19: moving the number of ionizing photons counter into this
   * function so it will be incremented for both macro and non-macro modes
   */

  if (HEV * p->freq > 13.6)     // only record if above H ionization edge
  {

    if (xplasma->nplasma != previous_nioniz_nplasma || p->np != previous_nioniz_np)
    {
      xplasma->est.nioniz++;
      previous_nioniz_nplasma = xplasma->nplasma;
      previous_nioniz_np = p->np;
    }

    /* IP needs to be radiation density in the cell. We sum contributions from
       each photon, then it is normalised in wind_update. */
    xplasma->est.ip += ((w_ave * ds) / (PLANCK * p->freq));

    if (HEV * p->freq < 13600)  //Tartar et al integrate up to 1000Ryd to define the ionization parameter
    {
      xplasma->derived.xi += (w_ave * ds);
    }

    if (p->nscat == 0)
    {
      xplasma->est.ip_direct += ((w_ave * ds) / (PLANCK * p->freq));
    }
    else
    {
      xplasma->est.ip_scatt += ((w_ave * ds) / (PLANCK * p->freq));
    }
  }



  return (0);
}

/**********************************************************/
/** 
 * @brief      updates the estimators for calculating
 * the radiative force on each Plasma Cell
 * spectra in each Plasma cell
 *
 * @param [in,out] PlasmaPtr  xplasma   PlasmaPtr for the cell of interest
 * @param [in] PhotPtr  phot_mid   Photon pointer at midpt in Observer frame
 * @param [in] double  ds   ds travelled in Observer frame
 * @param [in] double  w_ave   the weight of the photon in Observer frame
 *
 * @return  Always returns 0
 *
 *
 *
 * @details
 *
 * 
 *
 * ### Notes ###
 *
 * In non macro atom mode, w_ave is an average weight (passed as w_ave), but
 * in macro atom mode weights are never reduced (so p->w is used).
 *
 * This routine is called from bf_estimators in macro_atom modes
 * and from radiation.  Although the historical documentation
 * suggests it is only called for certain ionization modes, it appears
 * to be called in all cases, though clearly it only provides diagnostic
 * information in some of them.
 *
 * The flux estimators are constructed in the Observer frame;
 * inputs are expected to be in the Observer frame.
 *
 * The z-component sign flip (for lower-hemisphere photons contributing to
 * the vertical flux) is applied only for CYLVAR, which still uses bilateral
 * symmetry.  CYLIND and RTHETA have per-hemisphere cells with correctly-signed
 * coordinates and need no flip.
 *
 **********************************************************/



int
update_flux_estimators (PlasmaPtr xplasma, PhotPtr phot_mid, double ds_obs, double w_ave, int ndom)
{
  double flux[3];
  double flux_orig[3];
  int iangle;
  double binw, angle;

/* Compute the weighted flux for this photon */

  stuff_v (phot_mid->lmn, flux_orig);   //Get the direction of the photon packet
  renorm (flux_orig, w_ave * ds_obs);   //Renormalise the direction into a flux vector 

  project_from_xyz_cyl (phot_mid->x, flux_orig, flux);  //Transform to the frame ine which fluxes are summed 

  double theta;
  double r = sqrt (pow (phot_mid->x[0], 2) + pow (phot_mid->x[1], 2));

  if (zdom[ndom].coord_type == CYLVAR && phot_mid->x[2] < 0)
  {                             //If the photon is in the lower hemisphere in CYLVAR - we need to reverse the sense of the z flux
    flux[2] *= (-1);
    theta = atan2 (r, -phot_mid->x[2]);
  }
  else
  {
    theta = atan2 (r, phot_mid->x[2]);
  }

  angle = 0.0;

  if (flux[2] > 0 && flux[0] > 0)
    angle = (atan (flux[0] / flux[2]) * RADIAN);
  else if (flux[2] < 0 && flux[0] > 0)
    angle = (atan (flux[0] / flux[2]) * RADIAN) + 180;
  else if (flux[2] < 0 && flux[0] < 0)
    angle = (atan (flux[0] / flux[2]) * RADIAN) + 180;
  else if (flux[2] > 0 && flux[0] < 0)
    angle = 360. + (atan (flux[0] / flux[2]) * RADIAN);

  binw = 360. / NFLUX_ANGLES;


  iangle = (angle) / binw;      //Turn the angle into an integer to pass into the flux array

  //xplasma->est.F_UV_ang_theta[iangle] += flux[0];
  //xplasma->est.F_UV_ang_phi[iangle] += flux[1];
  //xplasma->est.F_UV_ang_r[iangle] += flux[2];

  xplasma->est.F_UV_ang_r[iangle] += flux[0] * sin (theta) + flux[2] * cos (theta);
  xplasma->est.F_UV_ang_phi[iangle] += flux[1];
  xplasma->est.F_UV_ang_theta[iangle] += flux[0] * cos (theta) - flux[2] * sin (theta);




  /* This ends the calculation of the overall intensity as a funciton of solid angle 
   * The quantities calculated below are simply the weighted fluxes at the posiiton
   * of the photons
   */


  if (zdom[ndom].coord_type == SPHERICAL)
  {
    renorm (phot_mid->x, 1);    //Create a unit vector in the direction of the photon from the origin
    flux[0] = dot (flux_orig, phot_mid->x);     //In the spherical geometry, the first comonent is the radial flux
    flux[1] = flux[2] = 0.0;    //In the spherical geomerry, the theta and phi compnents are zero    
  }



  if (phot_mid->freq < UV_low)
  {
    vadd (xplasma->est.F_vis, flux, xplasma->est.F_vis);
    xplasma->est.F_vis[3] += length (flux);
  }
  else if (phot_mid->freq > UV_hi)
  {
    vadd (xplasma->est.F_Xray, flux, xplasma->est.F_Xray);
    xplasma->est.F_Xray[3] += length (flux);
  }
  else
  {
    vadd (xplasma->est.F_UV, flux, xplasma->est.F_UV);
    xplasma->est.F_UV[3] += length (flux);
  }


  return (0);
}




/**********************************************************/
/**
 * @brief Updates the estimators for calculating the radiation force.
 *
 * @param[in] PlasmaPtr xplasma  PlasmaPtr for the cell of interest
 * @param[in] PhotPtr p  Photon pointer at the start of its path in the observer frame
 * @param[in] PhotPtr phot_mid  Photon pointer at the midpoint in the observer frame
 * @param[in] double ds  The path length travelled in the observer frame
 * @param[in] int ndom  The domain of interest
 * @param[in] double w_ave  The average weight along the photon's path
 * @param[in] double z  The absorbed energy fraction
 * @param[in] double frac_ff  The fractional contribution for free free opacity
 * @param[in] double frac_auger The fractional contribution to opacity from Auger ionization
 * @param[in] double frac_tot The fractional contribution
 *
 * @return Always return 0
 *
 * @details
 *
 * In non-macro mode, w_ave is an average weight but in macro mode weights are
 * never reduced, so p->w is used instead.
 *
 * ### Notes ###
 *
 * We calculate the force estimators in the observer frame and thus require
 * input quantities in the observer frame.
 *
 * The z-component sign flip of the cylindrical momentum increment (for
 * lower-hemisphere photons) is applied only for CYLVAR, which still uses
 * bilateral symmetry.  CYLIND and RTHETA have per-hemisphere cells and need
 * no flip.
 *
 **********************************************************/


int
update_force_estimators (PlasmaPtr xplasma, PhotPtr p, PhotPtr phot_mid, double ds, double w_ave, int ndom, double z, double frac_ff,
                         double frac_auger, double frac_tot)
{
  int i;
  double p_out[3], dp_cyl[3];

  /*Deal with the special case of a spherical geometry */

  if (zdom[ndom].coord_type == SPHERICAL)
  {
    renorm (phot_mid->x, 1);    //Create a unit vector in the direction of the photon from the origin
  }

  stuff_v (p->lmn, p_out);
  renorm (p_out, z * frac_ff / VLIGHT);
  if (zdom[ndom].coord_type == SPHERICAL)
  {
    dp_cyl[0] = dot (p_out, phot_mid->x);       //In the spherical geometry, the first component is the radial component
    dp_cyl[1] = dp_cyl[2] = 0.0;
  }
  else
  {
    project_from_xyz_cyl (phot_mid->x, p_out, dp_cyl);
    /* CYLIND and RTHETA have per-hemisphere cells; project_from_xyz_cyl gives the
       correct signed dp_cyl[2] directly.  Only CYLVAR still uses bilateral symmetry
       and needs the z-component flipped for lower-hemisphere photons. */
    if (zdom[ndom].coord_type == CYLVAR && p->x[2] < 0)
      dp_cyl[2] *= (-1);
  }
  for (i = 0; i < 3; i++)
  {
    xplasma->est.rad_force_ff[i] += dp_cyl[i];
  }
  xplasma->est.rad_force_ff[3] += length (dp_cyl);

  stuff_v (p->lmn, p_out);
  renorm (p_out, (z * (frac_tot + frac_auger)) / VLIGHT);
  if (zdom[ndom].coord_type == SPHERICAL)
  {
    dp_cyl[0] = dot (p_out, phot_mid->x);       //In the spherical geometry, the first component is the radial component
    dp_cyl[1] = dp_cyl[2] = 0.0;
  }
  else
  {
    project_from_xyz_cyl (phot_mid->x, p_out, dp_cyl);
    if (zdom[ndom].coord_type == CYLVAR && p->x[2] < 0)
      dp_cyl[2] *= (-1);
  }
  for (i = 0; i < 3; i++)
  {
    xplasma->est.rad_force_bf[i] += dp_cyl[i];
  }

  xplasma->est.rad_force_bf[3] += length (dp_cyl);

  stuff_v (p->lmn, p_out);
  renorm (p_out, w_ave * ds * klein_nishina (p->freq));
  if (zdom[ndom].coord_type == SPHERICAL)
  {
    dp_cyl[0] = dot (p_out, phot_mid->x);       //In the spherical geometry, the first component is the radial component
    dp_cyl[1] = dp_cyl[2] = 0.0;
  }
  else
  {
    project_from_xyz_cyl (phot_mid->x, p_out, dp_cyl);
    if (zdom[ndom].coord_type == CYLVAR && p->x[2] < 0)
      dp_cyl[2] *= (-1);
  }
  for (i = 0; i < 3; i++)
  {
    xplasma->est.rad_force_es[i] += dp_cyl[i];
  }
  xplasma->est.rad_force_es[3] += length (dp_cyl);

  return 0;
}


typedef struct
{
  double T;                     // Temperature in Kelvin
} planck_params;

// Planck spectral radiance function
double
planck_spectral_radiance (double nu, void *params)
{
  planck_params *p = (planck_params *) params;
  double T = p->T;
  double exponent = PLANCK * nu / (BOLTZMANN * T);
  return (2 * PLANCK * pow (nu, 3) / pow (VLIGHT, 2)) / (exp (exponent) - 1);
}

// Function to calculate nu * J(nu)
double
nu_times_radiance (double nu, void *params)
{
  return nu * planck_spectral_radiance (nu, params);
}

// Function to calculate mean frequency over a specified range
double
mean_frequency_nu_range (double T, double nu_min, double nu_max)
{
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

  double numerator, denominator, error;
  planck_params params = { T };

  gsl_function F;
  F.function = &nu_times_radiance;
  F.params = &params;

  gsl_integration_qags (&F, nu_min, nu_max, 0, 1e-7, 1000, w, &numerator, &error);

  F.function = &planck_spectral_radiance;

  gsl_integration_qags (&F, nu_min, nu_max, 0, 1e-7, 1000, w, &denominator, &error);

  gsl_integration_workspace_free (w);

  return numerator / denominator;
}

// Function to estimate temperature given a mean frequency and a frequency range
double
estimate_temperature_from_mean_frequency (double mean_nu_target, double nu_min, double nu_max, double initial_guess)
{
  double T_guess = initial_guess;
  double tol = 1e-6;            // Convergence tolerance
  int max_iter = 100;           // Maximum number of iterations
  double mean_nu;


  for (int i = 0; i < max_iter; i++)
  {
    // Calculate the mean frequency for the current temperature guess
    mean_nu = mean_frequency_nu_range (T_guess, nu_min, nu_max);

    // Adjust the temperature guess
    double T_new = T_guess * (mean_nu_target / mean_nu);

    // Check for convergence
    if (fabs (T_new - T_guess) < tol)
    {
      return T_new;
    }

    T_guess = T_new;
  }

  return T_guess;               // Return the estimated temperature
}


/**********************************************************/
/**
 * @brief Normalise simple estimators and cell spectra
 *
 * @param[in] PlasmaPtr xplasma  PlasmaPtr for the cell of interest
 *
 * @return Always return 0
 *
 * @details 
 * This routine normalises "simple" estimators and in a few cases compleets
 * the determination of the quantity itself by performaing calculations
 * of various quantities that have been accumulated during photon transfer
 *
 * The quantities vary - some of them are banded estimators used for the ionization 
 * state spectral models, others are diagnostics such as xi, others are used
 * in dilute approximations (w and t_r). 
 *
 * ### Notes ###
 *
 * Most of the quantities that are renormalized here are quantities that
 * are measured in the co-moving frame, but the radiation force is
 * an observer frame quantity and so is notmralized differently
 *
 **********************************************************/

/* switch to choose how to calculate T_r. This could be made to depend on the ionization mode
   if desired, or it may be that TRUE is the desired behaviour long-term */
#define BAND_CORRECTED_TRAD FALSE

int
normalise_simple_estimators (PlasmaPtr xplasma)
{
  int i, j, nwind;
  double radiation_temperature, nh, wtest;
  double volume_obs, invariant_volume_time;
  double electron_density_obs;
  double freq_min, freq_max, dfreq;

//OLD  double wedge_volume, rmin, rmax, thetamin, thetamax;  //things needed for directional fluxes

  nwind = xplasma->nwind;

  invariant_volume_time = wmain[nwind].vol / wmain[nwind].xgamma_cen;

  if (xplasma->est.ntot > 0)
  {
    wtest = xplasma->est.ave_freq;
    xplasma->est.ave_freq /= xplasma->est.j;    /* Normalization to frequency moment */
    if (sane_check (xplasma->est.ave_freq))
    {
      Error ("normalise_simple_estimators:sane_check nwind %d nplasma %d ave_freq %e j %e ntot %d\n", nwind, xplasma->nplasma, wtest,
             xplasma->est.j, xplasma->est.ntot);
    }

    xplasma->est.j /= (4. * PI * invariant_volume_time);
    xplasma->est.j_direct /= (4. * PI * invariant_volume_time);
    xplasma->est.j_scatt /= (4. * PI * invariant_volume_time);

    xplasma->state.t_r_old = xplasma->state.t_r;        // Store the previous t_r in t_r_old immediately before recalculating

    /* the method of calculation of the band corrected radiation temperature depends on 
       the flag BAND_CORRECTED_TRAD -- see issue #1097 */
    if (BAND_CORRECTED_TRAD == FALSE)
    {
      radiation_temperature = xplasma->state.t_r = PLANCK * xplasma->est.ave_freq / (BOLTZMANN * 3.832);
    }
    else
    {
      radiation_temperature = xplasma->state.t_r =
        estimate_temperature_from_mean_frequency (xplasma->est.ave_freq, xband.f1[0], xband.f2[xband.nbands - 1], radiation_temperature);
    }

    xplasma->state.w =
      PI * xplasma->est.j / (STEFAN_BOLTZMANN * radiation_temperature * radiation_temperature * radiation_temperature *
                             radiation_temperature);

    if (xplasma->state.w > 1e10)
    {
      Error ("normalise_simple_estimators: Huge w %8.2e in cell %d trad %10.2e j %8.2e\n", xplasma->state.w, xplasma->nplasma,
             radiation_temperature, xplasma->est.j);
    }
    if (sane_check (radiation_temperature) || sane_check (xplasma->state.w))
    {
      Error ("normalise_simple_estimators:sane_check %d trad %8.2e w %8.2g\n", xplasma->nplasma, radiation_temperature, xplasma->state.w);
      Error ("normalise_simple_estimators: ave_freq %8.2e j %8.2e\n", xplasma->est.ave_freq, xplasma->est.j);
      Exit (0);
    }
  }
  else
  {
    xplasma->est.j = xplasma->est.j_direct = xplasma->est.j_scatt = 0;
    if (modes.fixed_temp != 1)
      xplasma->state.t_e *= 0.7;
    if (xplasma->state.t_e < MIN_TEMP)
      xplasma->state.t_e = MIN_TEMP;
    xplasma->state.w = 0;
  }

  /* Temporary location for constructing estimators from speectra */


  if (xdev == TRUE)
  {

    j = 0;
    for (i = 0; i < xplasma->state.nbands; i++)
    {
      xplasma->est.xave_freq[i] = 0;
      xplasma->est.xsd_freq[i] = 0;
      xplasma->est.xj[i] = 0;
      xplasma->est.nxtot[i] = 0;
      xplasma->est.fmin[i] = geo.cell_freq[geo.nbins_in_cell_spec];
      xplasma->est.fmax[i] = geo.cell_freq[0];

      while (geo.cell_freq[j] < xplasma->state.f2[i] && j < geo.nbins_in_cell_spec)
      {
        double ave_freq;
        if (xplasma->est.cell_spec_flux[j] > 0)
        {
          ave_freq = 0.5 * (geo.cell_freq[j + 1] + geo.cell_freq[j]);
          xplasma->est.xave_freq[i] += ave_freq * xplasma->est.cell_spec_flux[j];
          xplasma->est.xsd_freq[i] += ave_freq * ave_freq * xplasma->est.cell_spec_flux[j];
          xplasma->est.xj[i] += xplasma->est.cell_spec_flux[j];
          xplasma->est.nxtot[i]++;
          if (ave_freq < xplasma->est.fmin[i])
          {
            xplasma->est.fmin[i] = ave_freq;
          }
          if (ave_freq > xplasma->est.fmax[i])
          {
            xplasma->est.fmax[i] = ave_freq;
          }
        }
        j++;
      }

    }

  }




  /* End of temporary work on constructing estimators */

  /* Normalize and otherwise complete the information gathered about the 
     coarse spectra used for generation of spectral models */

  for (i = 0; i < xplasma->state.nbands; i++)
  {
    if (xplasma->est.nxtot[i] > 0)
    {
      xplasma->est.xave_freq[i] /= xplasma->est.xj[i];
      xplasma->est.xsd_freq[i] /= xplasma->est.xj[i];
      xplasma->est.xsd_freq[i] = sqrt (xplasma->est.xsd_freq[i] - (xplasma->est.xave_freq[i] * xplasma->est.xave_freq[i]));     /*Compute standard deviation */
      xplasma->est.xj[i] /= (4 * PI * invariant_volume_time);   /*Convert to radiation density */
    }
    else
    {
      xplasma->est.xj[i] = 0;
      xplasma->est.xave_freq[i] = 0;
      xplasma->est.xsd_freq[i] = 0;
    }
  }


  /* Normalize the cell spectra to J_nu - erg/s/cm^3/Sr/Hz */

  for (i = 0; i < geo.nbins_in_cell_spec; i++)
  {
    freq_min = geo.cell_log_freq_min + (i) * geo.cell_delta_lfreq;
    freq_max = freq_min + geo.cell_delta_lfreq;
    dfreq = pow (10., freq_max) - pow (10., freq_min);

    xplasma->est.cell_spec_flux[i] /= (4 * PI * invariant_volume_time * dfreq);

  }


  /* Normalise IP, which at this point should be
   * the number of photons in a cell by dividing by volume
   * and number density of hydrogen in the cell
   */

  nh = xplasma->state.rho * rho2nh;
  xplasma->est.ip /= (VLIGHT * invariant_volume_time * nh);
  xplasma->est.ip_direct /= (VLIGHT * invariant_volume_time * nh);
  xplasma->est.ip_scatt /= (VLIGHT * invariant_volume_time * nh);

  /* Normalise xi, which at this point should be the luminosity of
   * ionizing photons in a cell (just the sum of photon weights)
   */

  xplasma->derived.xi *= 4. * PI;
  xplasma->derived.xi /= (invariant_volume_time * nh);

  /*
   * The radiation force and flux estimators are all observer frame
   * quantities and hence are normalised with observer frame volumes
   * and densities
   */

  /* there was an error here before, see #1030. The observer frame density is gamma times the CMF one */
  electron_density_obs = xplasma->state.ne * wmain[nwind].xgamma_cen;   // Mihalas & Mihalas p146
  volume_obs = wmain[nwind].vol / wmain[nwind].xgamma_cen;

  for (i = 0; i < NFORCE_DIRECTIONS; i++)
  {
    xplasma->est.rad_force_es[i] *= (volume_obs * electron_density_obs) / (volume_obs * VLIGHT);
    xplasma->est.F_vis[i] /= volume_obs;
    xplasma->est.F_UV[i] /= volume_obs;
    xplasma->est.F_Xray[i] /= volume_obs;
  }

  for (i = 0; i < NFLUX_ANGLES; i++)
  {
    xplasma->est.F_UV_ang_theta[i] /= volume_obs;
    xplasma->est.F_UV_ang_phi[i] /= volume_obs;
    xplasma->est.F_UV_ang_r[i] /= volume_obs;
  }

  return (0);
}

/**********************************************************/
/**
 * @brief
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

void
update_persistent_directional_flux_estimators (int nplasma, double flux_persist_scale)
{
  int n;
  double flux_helper[3];

  if (geo.wcycle == 0)          //If this is the first time through, then the persistent flux is empty.
  {
    vadd (plasmamain[nplasma].derived.F_vis_persistent, plasmamain[nplasma].est.F_vis, plasmamain[nplasma].derived.F_vis_persistent);
    vadd (plasmamain[nplasma].derived.F_UV_persistent, plasmamain[nplasma].est.F_UV, plasmamain[nplasma].derived.F_UV_persistent);
    vadd (plasmamain[nplasma].derived.F_Xray_persistent, plasmamain[nplasma].est.F_Xray, plasmamain[nplasma].derived.F_Xray_persistent);
    vadd (plasmamain[nplasma].derived.rad_force_bf_persist, plasmamain[nplasma].est.rad_force_bf,
          plasmamain[nplasma].derived.rad_force_bf_persist);

    for (n = 0; n < NFLUX_ANGLES; n++)
    {
      plasmamain[nplasma].derived.F_UV_ang_theta_persist[n] =
        plasmamain[nplasma].derived.F_UV_ang_theta_persist[n] + plasmamain[nplasma].est.F_UV_ang_theta[n];
      plasmamain[nplasma].derived.F_UV_ang_phi_persist[n] =
        plasmamain[nplasma].derived.F_UV_ang_phi_persist[n] + plasmamain[nplasma].est.F_UV_ang_phi[n];
      plasmamain[nplasma].derived.F_UV_ang_r_persist[n] =
        plasmamain[nplasma].derived.F_UV_ang_r_persist[n] + plasmamain[nplasma].est.F_UV_ang_r[n];
    }
  }
  else
  {
    rescale (plasmamain[nplasma].derived.F_vis_persistent, (1 - flux_persist_scale), plasmamain[nplasma].derived.F_vis_persistent);
    rescale (plasmamain[nplasma].est.F_vis, flux_persist_scale, flux_helper);
    vadd (plasmamain[nplasma].derived.F_vis_persistent, flux_helper, plasmamain[nplasma].derived.F_vis_persistent);

    rescale (plasmamain[nplasma].derived.F_UV_persistent, (1 - flux_persist_scale), plasmamain[nplasma].derived.F_UV_persistent);
    rescale (plasmamain[nplasma].est.F_UV, flux_persist_scale, flux_helper);
    vadd (plasmamain[nplasma].derived.F_UV_persistent, flux_helper, plasmamain[nplasma].derived.F_UV_persistent);

    rescale (plasmamain[nplasma].derived.F_Xray_persistent, (1 - flux_persist_scale), plasmamain[nplasma].derived.F_Xray_persistent);
    rescale (plasmamain[nplasma].est.F_Xray, flux_persist_scale, flux_helper);
    vadd (plasmamain[nplasma].derived.F_Xray_persistent, flux_helper, plasmamain[nplasma].derived.F_Xray_persistent);

    rescale (plasmamain[nplasma].derived.rad_force_bf_persist, (1 - flux_persist_scale), plasmamain[nplasma].derived.rad_force_bf_persist);
    rescale (plasmamain[nplasma].est.rad_force_bf, flux_persist_scale, flux_helper);
    vadd (plasmamain[nplasma].derived.rad_force_bf_persist, flux_helper, plasmamain[nplasma].derived.rad_force_bf_persist);

    for (n = 0; n < NFLUX_ANGLES; n++)
    {
      plasmamain[nplasma].derived.F_UV_ang_theta_persist[n] =
        plasmamain[nplasma].derived.F_UV_ang_theta_persist[n] * (1 - flux_persist_scale);
      plasmamain[nplasma].derived.F_UV_ang_phi_persist[n] = plasmamain[nplasma].derived.F_UV_ang_phi_persist[n] * (1 - flux_persist_scale);
      plasmamain[nplasma].derived.F_UV_ang_r_persist[n] = plasmamain[nplasma].derived.F_UV_ang_r_persist[n] * (1 - flux_persist_scale);
      plasmamain[nplasma].derived.F_UV_ang_theta_persist[n] =
        plasmamain[nplasma].derived.F_UV_ang_theta_persist[n] + plasmamain[nplasma].est.F_UV_ang_theta[n] * flux_persist_scale;
      plasmamain[nplasma].derived.F_UV_ang_phi_persist[n] =
        plasmamain[nplasma].derived.F_UV_ang_phi_persist[n] + plasmamain[nplasma].est.F_UV_ang_phi[n] * flux_persist_scale;
      plasmamain[nplasma].derived.F_UV_ang_r_persist[n] =
        plasmamain[nplasma].derived.F_UV_ang_r_persist[n] + plasmamain[nplasma].est.F_UV_ang_r[n] * flux_persist_scale;
    }
  }
  plasmamain[nplasma].derived.F_vis_persistent[3] = length (plasmamain[nplasma].derived.F_vis_persistent);
  plasmamain[nplasma].derived.F_UV_persistent[3] = length (plasmamain[nplasma].derived.F_UV_persistent);
  plasmamain[nplasma].derived.F_Xray_persistent[3] = length (plasmamain[nplasma].derived.F_Xray_persistent);
  plasmamain[nplasma].derived.rad_force_bf_persist[3] = length (plasmamain[nplasma].derived.rad_force_bf_persist);

}
