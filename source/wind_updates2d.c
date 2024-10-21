/***********************************************************/
/** @file  wind_updates2d.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  This file contains the main routines for updating
 * and then reinitializing the wind after an ionization cycle
 *
 * The routines in this file are generic.  There is no dependence on a
 * particular wind model or any coordinate system dependencies.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

/**********************************************************/
/**
 * @brief      updates the parameters in the wind that are
 * 	affected by radiation, including ion densities.
 *
 * @param [in] WindPtr  The entire wind
 * @return     Always returns 0
 *
 * @details
 * This is the main routine used to update the wind at the end of
 * an ionization cycle (in preparation for a new cycle).  The routine
 * is parallelized to save time
 *
 * ### Notes ###
 * At the time wind_update is called the various quantities that are accumulated
 * during the photon transfer part of the cycle has been accumulated
 * and shared between the threads.  Here certain plasma cells are assigned to
 * each thread so that the ionization structure can be updated, and each thread
 * is responsible for calculaing the updates for a certain set of cells.  At
 * the end of the routine the updates are collected and reshared.
 *
 * The real need for prallelising the routine is the work done in ion_abundances
 *
 * Once this is done, various checks are made to determined what happened as a
 * function of the updates, various variables in geo are updated,  and for
 * a hydro model the results are written to a file
 *
 **********************************************************/

int
wind_update (WindPtr w)
{
  int n_plasma, i, j;
  double xsum, psum, fsum, lsum, csum, icsum, ausum, chexsum;
  double cool_sum, lum_sum, radiated_luminosity_sum;    //1706 - the total cooling and luminosity of the wind
  double apsum, aausum, abstot; //Absorbed photon energy from PI and auger
  double flux_persist_scale;
  double volume;
  double dt_r, dt_e;
  double t_r_ave_old, t_r_ave, t_e_ave_old, t_e_ave;
  int nmax_r, nmax_e;
  int nwind;
  int my_nmin, my_nmax;         //Note that these variables are still used even without MPI on
  int ndom;
  int n_cells_rank;

  dt_r = 0.0;
  dt_e = 0.0;
  nmax_r = -1;
  nmax_e = -1;
  t_r_ave_old = 0.0;
  t_r_ave = 0.0;
  t_e_ave_old = 0.0;
  t_e_ave = 0.0;

  xsignal (files.root, "%-20s Start wind update\n", "NOK");

#ifdef MPI_ON
  n_cells_rank = get_parallel_nrange (rank_global, NPLASMA, np_mpi_global, &my_nmin, &my_nmax);
#else
  my_nmin = 0;
  my_nmax = NPLASMA;
  n_cells_rank = NPLASMA;
#endif

  flux_persist_scale = 0.5;     //The amount of the latest flux that gets added into the persistent flux

  for (n_plasma = my_nmin; n_plasma < my_nmax; ++n_plasma)
  {
    nwind = plasmamain[n_plasma].nwind;
    volume = w[nwind].vol;

    /* Skip cells that are partially in the wind these are not to be included
       in the calculation */

    if (modes.partial_cells == PC_EXTEND && wmain[nwind].inwind == W_PART_INWIND)
    {
      continue;
    }

    if (plasmamain[n_plasma].ntot < 100)
    {
      Log
        ("!!wind_update: Cell %4d Dom %d  Vol. %8.2e r %8.2e theta %8.2e has only %4d photons\n",
         n_plasma, w[nwind].ndom, volume, w[nwind].rcen, w[nwind].thetacen, plasmamain[n_plasma].ntot);
    }

    /* Start with a call to the routine which normalises all the macro atom
       monte carlo radiation field estimators. It's best to do this first since
       some of the estimators include temperature terms (stimulated correction
       terms) which were included during the monte carlo simulation so we want
       to be sure that the SAME temperatures are used here. */

    if (geo.rt_mode == RT_MODE_MACRO && geo.macro_simple == FALSE)
    {
      normalise_macro_estimators (&plasmamain[n_plasma]);
    }

    /* this routine normalises the unbanded and banded estimators for simple atoms */
    normalise_simple_estimators (&plasmamain[n_plasma]);

    /* update the persistent fluxes */
    update_persistent_directional_flux_estimators (n_plasma, flux_persist_scale);

    /* If geo.adiabatic is true, then calculate the adiabatic cooling using the current, i.e
     * previous value of t_e.  Note that this may not be the best way to determine the cooling.
     * Changes made here should also be reflected in wind2d.c. At present, adiabatic cooling
     * is not included in updates to the temperature, even if the adiabatic cooling is calculated
     * here.
     */

    if (geo.adiabatic)
    {
      plasmamain[n_plasma].cool_adiabatic = adiabatic_cooling (&w[nwind], plasmamain[n_plasma].t_e);
    }
    else
    {
      plasmamain[n_plasma].cool_adiabatic = 0.0;
    }

    if (geo.nonthermal)
    {
      plasmamain[n_plasma].heat_shock = shock_heating (&w[nwind]);
    }
    else
    {
      plasmamain[n_plasma].heat_shock = 0.0;
    }

    /* Calculate the densities in various ways depending on the ioniz_mode */
    ion_abundances (&plasmamain[n_plasma], geo.ioniz_mode);
  }

  /*This is the end of the update loop that is parallised. We now need to exchange data between the tasks. */

  broadcast_updated_plasma_properties (my_nmin, my_nmax, n_cells_rank);
  if (geo.rt_mode == RT_MODE_MACRO && geo.macro_simple == FALSE)
  {
    broadcast_updated_macro_atom_properties (my_nmin, my_nmax, n_cells_rank);
  }

  /* Now we need to updated the densities immediately outside the wind so that the density interpolation in resonate will work.
     In this case all we have done is to copy the densities from the cell which is just in the wind (as one goes outward) to the
     cell that is just inside (or outside) the wind.

     SS asked whether we should also be extending the wind for other parameters, especially ne.  At present we do not interpolate
     on ne so this is not necessary.  If we did do that it would be required.

     In cylindrical coordinates, the fast dimension is z; grid positions increase up in z, and then out in r.
     In spherical polar coordinates, the fast dimension is theta; the grid increases in theta (measured)
     from the z axis), and then in r.
     In spherical coordinates, the grid increases as one might expect in r..
     *
   */

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    if (zdom[ndom].coord_type == CYLIND)
    {
      cylind_extend_density (ndom, w);
    }
    else if (zdom[ndom].coord_type == RTHETA)
    {
      rtheta_extend_density (ndom, w);
    }
    else if (zdom[ndom].coord_type == SPHERICAL)
    {
      spherical_extend_density (ndom, w);
    }
    else if (zdom[ndom].coord_type == CYLVAR)
    {
      cylvar_extend_density (ndom, w);
    }
    else
    {
      Error ("Wind_update2d: Unknown coordinate type %d for domain %d \n", zdom[ndom].coord_type, ndom);
      Exit (EXIT_FAILURE);
    }
  }
  /* Finished updating region outside of wind */

  /* Check the balance between the absorbed and the emitted flux */
  /* NSH 0717 - ensure the cooling and luminosities reflect the current temperature */

  /* We now need to re-calculate the wind cooling rates, as well as calculate the total luminosity of the wind.
   * Both of these calculations are parallelised, although are done in separation from another. There may be some
   * scope with a (extensive) re-factor to combine them into one parallelised section. */
  cool_sum = wind_cooling ();
  lum_sum = wind_luminosity (0.0, VERY_BIG, MODE_CMF_TIME);

  xsum = 0.0;
  psum = 0.0;
  ausum = 0.0;
  lsum = 0.0;
  fsum = 0.0;
  csum = 0.0;
  icsum = 0.0;
  apsum = 0.0;
  aausum = 0.0;
  abstot = 0.0;
  chexsum = 0.0;

  /* Each rank now has updated plasma cells (temperature, ion abundances, heat/cool rates, etc.), so we can now find
   * out what the max d_t is in the wind and also sum up properties to find the total/global values */
  for (n_plasma = 0; n_plasma < NPLASMA; ++n_plasma)
  {
    /* First we want to find the maximum change in temperature, which we will
     * use for reporting and to calculate the convergence */
    if ((fabs (plasmamain[n_plasma].t_r_old - plasmamain[n_plasma].t_r)) > fabs (dt_r))
    {
      dt_r = plasmamain[n_plasma].t_r - plasmamain[n_plasma].t_r_old;
      nmax_r = n_plasma;
    }
    if ((fabs (plasmamain[n_plasma].t_e_old - plasmamain[n_plasma].t_e)) > fabs (dt_e))
    {
      dt_e = plasmamain[n_plasma].t_e - plasmamain[n_plasma].t_e_old;
      nmax_e = n_plasma;
    }

    t_r_ave += plasmamain[n_plasma].t_r;
    t_e_ave += plasmamain[n_plasma].t_e;
    t_r_ave_old += plasmamain[n_plasma].t_r_old;
    t_e_ave_old += plasmamain[n_plasma].t_e_old;

    check_heating_rates_for_plasma_cell (n_plasma);

    plasmamain[n_plasma].cool_tot_ioniz = plasmamain[n_plasma].cool_tot;
    plasmamain[n_plasma].lum_ff_ioniz = plasmamain[n_plasma].lum_ff;
    plasmamain[n_plasma].cool_rr_ioniz = plasmamain[n_plasma].cool_rr;
    plasmamain[n_plasma].lum_rr_ioniz = plasmamain[n_plasma].lum_rr;
    plasmamain[n_plasma].cool_rr_metals_ioniz = plasmamain[n_plasma].cool_rr_metals;
    plasmamain[n_plasma].lum_lines_ioniz = plasmamain[n_plasma].lum_lines;
    plasmamain[n_plasma].cool_comp_ioniz = plasmamain[n_plasma].cool_comp;
    plasmamain[n_plasma].cool_dr_ioniz = plasmamain[n_plasma].cool_dr;
    plasmamain[n_plasma].cool_di_ioniz = plasmamain[n_plasma].cool_di;
    plasmamain[n_plasma].lum_tot_ioniz = plasmamain[n_plasma].lum_tot;
    plasmamain[n_plasma].cool_adiabatic_ioniz = plasmamain[n_plasma].cool_adiabatic;

    abstot += plasmamain[n_plasma].abs_tot;
    xsum += plasmamain[n_plasma].heat_tot;
    psum += plasmamain[n_plasma].heat_photo;
    ausum += plasmamain[n_plasma].heat_auger;
    fsum += plasmamain[n_plasma].heat_ff;
    lsum += plasmamain[n_plasma].heat_lines;
    csum += plasmamain[n_plasma].heat_comp;
    icsum += plasmamain[n_plasma].heat_ind_comp;
    apsum += plasmamain[n_plasma].abs_photo;
    aausum += plasmamain[n_plasma].abs_auger;
    chexsum += plasmamain[n_plasma].heat_ch_ex;
  }

  /* We can now calculate the average of the t */
  t_r_ave /= NPLASMA;
  t_e_ave /= NPLASMA;
  t_r_ave_old /= NPLASMA;
  t_e_ave_old /= NPLASMA;

  /* Update global ionization properties */
  geo.lum_ff_ioniz = geo.lum_ff;
  geo.cool_rr_ioniz = geo.cool_rr;
  geo.lum_rr_ioniz = geo.lum_rr;
  geo.lum_lines_ioniz = geo.lum_lines;
  geo.cool_comp_ioniz = geo.cool_comp;
  geo.cool_dr_ioniz = geo.cool_dr;
  geo.cool_di_ioniz = geo.cool_di;
  geo.cool_adiabatic_ioniz = geo.cool_adiabatic;
  geo.lum_disk_ioniz = geo.lum_disk;
  geo.lum_star_ioniz = geo.lum_star;
  geo.lum_bl_ioniz = geo.lum_bl;
  geo.lum_wind_ioniz = geo.lum_wind;
  geo.lum_tot_ioniz = geo.lum_tot;

  if (modes.zeus_connect == 1 && geo.hydro_domain_number > -1)  //If we are running in zeus connect mode, we output heating and cooling rates.
  {
    create_hydro_output_files ();
  }

  /* Added this system which counts number of times two situations occur (See #91)
     We only report these every 100,000 times (one can typically get ) */
  Log ("wind_update: note, errors from mean intensity can be high in a working model\n");
  Log
    ("wind_update: can be a problem with photon numbers if there are also errors from spectral_estimators and low photon number warnings\n");
  Log ("wind_update: mean_intensity: %8.4e occurrences, this cycle, this thread of 'no model exists in a band'\n", nerr_no_Jmodel);
  Log
    ("wind_update: mean intensity: %8.4e occurrences, this cycle, this thread of 'photon freq is outside frequency range of spectral model'\n",
     nerr_Jmodel_wrong_freq);

  /* zero the counters which record diagnostics from the function mean_intensity */
  nerr_Jmodel_wrong_freq = 0;
  nerr_no_Jmodel = 0;

  /* The lines differ only in that Wind_heating adds mechanical heating, that is adiabatic heating */
  Log
    ("!!wind_update: Absorbed flux    %8.2e  (photo %8.2e ff %8.2e compton %8.2e auger %8.2e induced_compton %8.2e lines %8.2e)\n",
     abstot, apsum, fsum, csum, aausum, icsum, lsum);
  Log
    ("!!wind_update: Wind heating     %8.2e  (photo %8.2e ff %8.2e compton %8.2e auger %8.2e induced_compton %8.2e lines %8.2e adiabatic %8.2e)\n",
     xsum + geo.heat_adiabatic, psum, fsum, csum, ausum, icsum, lsum, geo.heat_adiabatic);

  /* As was the case above, there are two almost identical lines.  Wind_cooling includes processes that do not produce photons,
   * not-only adiabatic cooling, but also goe.cool_comp, geo_cool_dr and geo.cool_di */
  Log
    ("!!wind_update: Wind luminosity  %8.2e (recomb %8.2e ff %8.2e lines %8.2e) after update\n",
     lum_sum, geo.lum_rr, geo.lum_ff, geo.lum_lines);
  radiated_luminosity_sum = wind_luminosity (xband.f1[0], xband.f2[xband.nbands - 1], MODE_CMF_TIME);
  Log
    ("!!wind_update: Rad luminosity  %8.2e (recomb %8.2e ff %8.2e lines %8.2e) after update\n",
     radiated_luminosity_sum, geo.lum_rr, geo.lum_ff, geo.lum_lines);
  Log
    ("!!wind_update: Wind cooling     %8.2e (recomb %8.2e ff %8.2e compton %8.2e DR %8.2e DI %8.2e lines %8.2e adiabatic %8.2e) after update\n",
     cool_sum, geo.cool_rr, geo.lum_ff, geo.cool_comp, geo.cool_dr, geo.cool_di, geo.lum_lines, geo.cool_adiabatic);

  if (modes.use_upweighting_of_simple_macro_atoms)
  {
    /* If we have "indivisible packet" mode on but are using the
       upweighting scheme for simple atoms then we report the flows into and out of the ion pool */
    if (geo.rt_mode == RT_MODE_MACRO)
    {
      report_bf_simple_ionpool ();
    }
  }

  /* report a warning if the induced Compton heating is greater than 10% of the heating, see #1016 */
  if (icsum >= (0.1 * xsum))
  {
    Error ("!!wind_update: Induced Compton is responsible for %3.1f percent of radiative heating. Could cause problems.\n",
           icsum / xsum * 100.0);
  }


  if (modes.zeus_connect == 1 || modes.fixed_temp == 1) //There is no point in computing temperature changes, because we have fixed them!
  {
    Log ("!!wind_update: We are running in fixed temperature mode - no temperature report\n");
  }

  if (modes.zeus_connect != TRUE || modes.fixed_temp != TRUE)
  {
    if (nmax_r != -1)
    {
      wind_n_to_ij (wmain[nmax_r].ndom, nmax_r, &i, &j);
      Log ("!!wind_update: Max change in t_r %6.0f at cell %4d (%d,%d)\n", dt_r, nmax_r, i, j);
      Log ("!!wind_update: Ave change in t_r %6.0f from %6.0f to %6.0f\n", (t_r_ave - t_r_ave_old), t_r_ave_old, t_r_ave);
    }
    else
    {
      Log ("!!wind_update: t_r did not change in any cells this cycle\n");
    }
    if (nmax_e != -1)
    {
      wind_n_to_ij (wmain[nmax_e].ndom, nmax_e, &i, &j);
      Log ("!!wind_update: Max change in t_e %6.0f at cell %4d (%d,%d)\n", dt_e, nmax_e, i, j);
      Log ("!!wind_update: Ave change in t_e %6.0f from %6.0f to %6.0f\n", (t_e_ave - t_e_ave_old), t_e_ave_old, t_e_ave);
    }
    else
    {
      Log ("!!wind_update: t_e did not change in any cells this cycle\n");
    }
    Log ("Summary  t_r  %6.0f   %6.0f  #t_r and dt_r on this update\n", t_r_ave, (t_r_ave - t_r_ave_old));
    Log ("Summary  t_e  %6.0f   %6.0f  #t_e and dt_e on this update\n", t_e_ave, (t_e_ave - t_e_ave_old));
  }

  check_convergence ();

  /* Summarize the radiative temperatures (ksl 04 mar) */
  xtemp_rad (w);

/* This next block is to allow the output of data relating to the abundances of ions when sirocco is being tested
 * with thin shell mode.
 */
  shell_output_wind_update_diagnostics (xsum, psum, fsum, csum, icsum, lsum, ausum, chexsum, cool_sum, lum_sum);

  xsignal (files.root, "%-20s Finished wind update\n", "NOK");

  return (0);
}

/**********************************************************/
/**
 * @brief This summarises the flows into and out of the ionization pool for
 *        simple ions in RT_MODE_MACRO
 *
 * @return    Always returns 0
 *
 **********************************************************/

int
report_bf_simple_ionpool (void)
{
  int n, m;
  int in_tot, out_tot;
  double total_in = 0.0;
  double total_out = 0.0;

  for (n = 0; n < NPLASMA; n++)
  {
    total_in += plasmamain[n].bf_simple_ionpool_in;
    total_out += plasmamain[n].bf_simple_ionpool_out;

    if (plasmamain[n].bf_simple_ionpool_out > plasmamain[n].bf_simple_ionpool_in)
    {
      Error ("The net flow out of simple ion pool (%8.4e) > than the net flow in (%8.4e) in cell %d\n",
             plasmamain[n].bf_simple_ionpool_out, plasmamain[n].bf_simple_ionpool_in, n);
    }
  }

  Log ("!! report_bf_simple_ionpool: Total flow into: %8.4e and out of: %8.4e bf_simple ion pool\n", total_in, total_out);

  total_in = total_out = 0;
  for (m = 0; m < nphot_total; m++)
  {
    in_tot = out_tot = 0;
    for (n = 0; n < NPLASMA; n++)
    {
      in_tot += plasmamain[n].n_bf_in[m];
      out_tot += plasmamain[n].n_bf_out[m];
    }

    Log ("!! report_bf:  %3d   %3d %3d %7d  %7d\n", m, phot_top[m].z, phot_top[m].istate, in_tot, out_tot);

    total_in += in_tot;
    total_out += out_tot;
  }

  Log ("!! report_bf tots:   %10.0f  %10.0f\n", total_in, total_out);

  return (0);
}

/**********************************************************/
/**
 * @brief      zeros those portions of the wind which contain the radiation properties
 * 	of the wind, i.e those portions which should be set to zeroed when the structure of the
 * 	wind has been changed or when you simply want to start off a calculation in a known state
 *
 * @details
 * The routine is called at the beginning of each ionization calculation
 * cycle.  It should zero all heating and radiation induced cooling in the Plasma structure.  Since
 * cooling is recalculated in wind_update, one needs to be sure that all of the appropriate
 * cooling terms are also rezeroed there as well.
 *
 * ### Notes ###
 *
 **********************************************************/

void
wind_rad_init ()
{
  init_plasma_rad_properties ();
  init_macro_rad_properties ();
}

/**********************************************************/
/**
 * @brief Check for inf or NaN for a plasma cell's heating rates
 *
 * @params [in]  int n_plasma  the index of the plasma cell
 *
 * @details
 *
 * This function uses `sane_check()` to check for NaN, inf and etc. for the
 * heating rates heat_tot, heat_photo, heat_auger, heat_photo_macro, heat_ff, heat_lines,
 * heat_lines_macro and heat_comp.
 *
 * When a heating rate iis not "sane", an error is printed.
 *
 * @TODO: this should really also include the cooling rates
 *
 **********************************************************/

void
check_heating_rates_for_plasma_cell (const int n_plasma)
{
  if (sane_check (plasmamain[n_plasma].heat_tot))
  {
    Error ("wind_update:sane_check w(%d).heat_tot is %e\n", n_plasma, plasmamain[n_plasma].heat_tot);
  }
  if (sane_check (plasmamain[n_plasma].heat_photo))
  {
    Error ("wind_update:sane_check w(%d).heat_photo is %e\n", n_plasma, plasmamain[n_plasma].heat_photo);
  }
  if (sane_check (plasmamain[n_plasma].heat_auger))
  {
    Error ("wind_update:sane_check w(%d).heat_auger is %e\n", n_plasma, plasmamain[n_plasma].heat_auger);
  }
  if (sane_check (plasmamain[n_plasma].heat_photo_macro))
  {
    Error ("wind_update:sane_check w(%d).heat_photo_macro is %e\n", n_plasma, plasmamain[n_plasma].heat_photo_macro);
  }
  if (sane_check (plasmamain[n_plasma].heat_ff))
  {
    Error ("wind_update:sane_check w(%d).heat_ff is %e\n", n_plasma, plasmamain[n_plasma].heat_ff);
  }
  if (sane_check (plasmamain[n_plasma].heat_lines))
  {
    Error ("wind_update:sane_check w(%d).heat_lines is %e\n", n_plasma, plasmamain[n_plasma].heat_lines);
  }
  if (sane_check (plasmamain[n_plasma].heat_lines_macro))
  {
    Error ("wind_update:sane_check w(%d).heat_lines_macro is %e\n", n_plasma, plasmamain[n_plasma].heat_lines_macro);
  }
  if (sane_check (plasmamain[n_plasma].heat_comp))
  {
    Error ("wind_update:sane_check w(%d).heat_comp is %e\n", n_plasma, plasmamain[n_plasma].heat_comp);
  }
}

/**********************************************************/
/**
 * @brief Initialises the radiative properties for all plasma cells
 *
 * @details
 *
 * This works over all plasma cells and does so in a serial loop. At the moment
 * there is no scope for making this parallelised, as NPLASMA is going to be
 * small enough for the far foreseeable future that setting the values of various
 * fields in each cell to zero is not a bottleneck.
 *
 **********************************************************/

void
init_plasma_rad_properties (void)
{
  int i;
  int j;

  for (i = 0; i < NPLASMA; ++i)
  {
    /* Start by initialising integer fields */
    plasmamain[i].j = 0;
    plasmamain[i].ave_freq = 0;
    plasmamain[i].ntot = 0;
    plasmamain[i].n_ds = 0;
    plasmamain[i].ntot_disk = 0;
    plasmamain[i].ntot_agn = 0;
    plasmamain[i].ntot_star = 0;
    plasmamain[i].ntot_bl = 0;
    plasmamain[i].nscat_es = 0;
    plasmamain[i].nscat_res = 0;
    plasmamain[i].ntot_wind = 0;
    plasmamain[i].nrad = 0;
    plasmamain[i].nioniz = 0;
    for (j = 0; j < nphot_total; j++)
    {
      plasmamain[i].n_bf_in[j] = 0;
      plasmamain[i].n_bf_out[j] = 0;
    }

    /* Next we'll initialise the rest of the fields, which are doubles */
    plasmamain[i].j_direct = 0.0;
    plasmamain[i].j_scatt = 0.0;
    plasmamain[i].ip = 0.0;
    plasmamain[i].xi = 0.0;
    plasmamain[i].ip_direct = 0.0;
    plasmamain[i].ip_scatt = 0.0;
    plasmamain[i].mean_ds = 0.0;
    plasmamain[i].heat_tot = 0.0;
    plasmamain[i].heat_ff = 0.0;
    plasmamain[i].heat_photo = 0.0;
    plasmamain[i].heat_lines = 0.0;
    plasmamain[i].abs_tot = 0.0;
    plasmamain[i].abs_auger = 0.0;
    plasmamain[i].abs_photo = 0.0;
    plasmamain[i].heat_z = 0.0;
    plasmamain[i].max_freq = 0.0;
    plasmamain[i].cool_tot = 0.0;
    plasmamain[i].lum_tot = 0.0;
    plasmamain[i].lum_lines = 0.0;
    plasmamain[i].lum_ff = 0.0;
    plasmamain[i].cool_rr = 0.0;
    plasmamain[i].cool_rr_metals = 0.0;
    plasmamain[i].lum_rr = 0.0;
    plasmamain[i].comp_nujnu = -1e99;
    plasmamain[i].cool_comp = 0.0;
    plasmamain[i].heat_comp = 0.0;
    plasmamain[i].heat_ind_comp = 0.0;
    plasmamain[i].heat_auger = 0.0;
    plasmamain[i].heat_ch_ex = 0.0;
    plasmamain[i].bf_simple_ionpool_out = 0.0;
    plasmamain[i].bf_simple_ionpool_in = 0.0;

    for (j = 0; j < N_DMO_DT_DIRECTIONS; j++)
    {
      plasmamain[i].dmo_dt[j] = 0.0;
    }
    for (j = 0; j < NFORCE_DIRECTIONS; j++)
    {
      plasmamain[i].rad_force_es[j] = 0.0;
      plasmamain[i].rad_force_ff[j] = 0.0;
      plasmamain[i].rad_force_bf[j] = 0.0;
      plasmamain[i].F_vis[j] = 0.0;
      plasmamain[i].F_UV[j] = 0.0;
      plasmamain[i].F_Xray[j] = 0.0;
      if (geo.wcycle == 0)      // Persistent values, so only initialise for first ionisation cycle
      {
        plasmamain[i].F_vis_persistent[j] = 0.0;
        plasmamain[i].F_UV_persistent[j] = 0.0;
        plasmamain[i].F_Xray_persistent[j] = 0.0;
        plasmamain[i].rad_force_bf_persist[j] = 0.0;
      }
    }
    for (j = 0; j < NFLUX_ANGLES; j++)
    {
      if (geo.wcycle == 0)      // Persistent values, so only initialise for first ionisation cycle
      {
        plasmamain[i].F_UV_ang_theta_persist[j] = 0.0;
        plasmamain[i].F_UV_ang_phi_persist[j] = 0.0;
        plasmamain[i].F_UV_ang_r_persist[j] = 0.0;
      }
      plasmamain[i].F_UV_ang_theta[j] = 0.0;
      plasmamain[i].F_UV_ang_phi[j] = 0.0;
      plasmamain[i].F_UV_ang_r[j] = 0.0;
    }

    /* Initialise  the frequency banded radiation estimators used for estimating the coarse spectra in each i */
    for (j = 0; j < NXBANDS; j++)
    {
      plasmamain[i].nxtot[j] = 0;
      plasmamain[i].xj[j] = 0.0;
      plasmamain[i].xave_freq[j] = 0.0;
      plasmamain[i].xsd_freq[j] = 0.0;
      plasmamain[i].fmin[j] = geo.xfreq[j + 1]; /* Set the minium frequency to the max frequency in the band */
      plasmamain[i].fmax[j] = geo.xfreq[j];     /* Set the maximum frequency to the min frequency in the band */
    }
    for (j = 0; j < NBINS_IN_CELL_SPEC; ++j)
    {
      plasmamain[i].cell_spec_flux[j] = 0.0;
    }

    for (j = 0; j < nions; j++)
    {
      plasmamain[i].ioniz[j] = 0.0;
      plasmamain[i].recomb[j] = 0.0;
      plasmamain[i].heat_ion[j] = 0.0;
      plasmamain[i].cool_rr_ion[j] = 0.0;
      plasmamain[i].lum_rr_ion[j] = 0.0;
      plasmamain[i].heat_inner_ion[j] = 0.0;

    }
    for (j = 0; j < n_inner_tot; j++)
    {
      plasmamain[i].inner_ioniz[j] = 0.0;
    }
  }
}

/**********************************************************/
/**
 * @brief Initialise the radiative properties for all the macro cells.
 *
 * @details
 *
 * Parts of this are parallelised, as calculating the stimulated recombination
 * rate alpha_sp is expensive. Especially for "high" resolution grids. The other
 * parts of this function are inexpensive due to NPLASMA being quite small for
 * the foreseeable future.
 *
 * The reason to not parallelise everything is to avoid communication overheads,
 * as the latency of communication (and packing and unpacking data) will
 * comparable or larger than the time it takes to do it in serial.
 *
 * @TODO: blocking broadcast should be replaced with non-blocking
 *
 **********************************************************/

void
init_macro_rad_properties (void)
{
  int n_plasma;
  int k;
  int macro_level;
  int n_start;
  int n_stop;
  int n_cells;

  /* Initialising these properties are inexpensive, so is not parallelised */
  for (n_plasma = 0; n_plasma < NPLASMA; ++n_plasma)
  {
    if (geo.rt_mode == RT_MODE_MACRO)
    {
      macromain[n_plasma].kpkt_rates_known = FALSE;
    }

    plasmamain[n_plasma].kpkt_emiss = 0.0;
    plasmamain[n_plasma].kpkt_abs = 0.0;

    for (macro_level = 0; macro_level < nlevels_macro; ++macro_level)
    {
      macromain[n_plasma].matom_abs[macro_level] = 0.0;
      macromain[n_plasma].matom_emiss[macro_level] = 0.0;

      for (k = 0; k < xconfig[macro_level].n_bbu_jump; ++k)
      {
        macromain[n_plasma].jbar[xconfig[macro_level].bbu_indx_first + k] = 0.0;
      }
      for (k = 0; k < xconfig[macro_level].n_bfu_jump; ++k)
      {
        macromain[n_plasma].gamma[xconfig[macro_level].bfu_indx_first + k] = 0.0;
        macromain[n_plasma].gamma_e[xconfig[macro_level].bfu_indx_first + k] = 0.0;
        macromain[n_plasma].alpha_st[xconfig[macro_level].bfd_indx_first + k] = 0.0;
        macromain[n_plasma].alpha_st_e[xconfig[macro_level].bfd_indx_first + k] = 0.0;
      }
    }
  }

  /* calculating recomb_sp and recomb_simple is expensive due to calls to
   * `alpha_sp()` , so we do this part of the initialisation in parallel */

#ifdef MPI_ON
  n_cells = get_parallel_nrange (rank_global, NPLASMA, np_mpi_global, &n_start, &n_stop);
#else
  n_start = 0;
  n_stop = NPLASMA;
  n_cells = NPLASMA;
#endif

  for (n_plasma = n_start; n_plasma < n_stop; ++n_plasma)
  {
    for (macro_level = 0; macro_level < nlevels_macro; ++macro_level)
    {
      for (k = 0; k < xconfig[macro_level].n_bfd_jump; ++k)
      {
        if (plasmamain[n_plasma].t_e > 1.0)
        {
          macromain[n_plasma].recomb_sp[xconfig[macro_level].bfd_indx_first + k] =
            alpha_sp (&phot_top[xconfig[macro_level].bfd_jump[k]], &plasmamain[n_plasma], 0);
          macromain[n_plasma].recomb_sp_e[xconfig[macro_level].bfd_indx_first + k] =
            alpha_sp (&phot_top[xconfig[macro_level].bfd_jump[k]], &plasmamain[n_plasma], 2);
        }
        else
        {
          macromain[n_plasma].recomb_sp[xconfig[macro_level].bfd_indx_first + k] = 0.0;
          macromain[n_plasma].recomb_sp_e[xconfig[macro_level].bfd_indx_first + k] = 0.0;
        }
      }
    }
    for (macro_level = 0; macro_level < ntop_phot; ++macro_level)
    {
      if ((geo.macro_simple == FALSE && phot_top[macro_level].macro_info == TRUE) || geo.rt_mode == RT_MODE_2LEVEL)
      {
        plasmamain[n_plasma].recomb_simple[macro_level] = 0.0;
        plasmamain[n_plasma].recomb_simple_upweight[macro_level] = 1.0;
      }
      else                      // we want a macro approach, but not for this ion so need recomb_simple instead
      {
        const double alpha_store = plasmamain[n_plasma].recomb_simple[macro_level] =
          alpha_sp (&phot_top[macro_level], &plasmamain[n_plasma], 2);
        plasmamain[n_plasma].recomb_simple_upweight[macro_level] =
          alpha_sp (&phot_top[macro_level], &plasmamain[n_plasma], 1) / alpha_store;
      }
    }
  }

  broadcast_macro_atom_recomb (n_start, n_stop, n_cells);
}

/**********************************************************/
/**
 * @brief  Print extra diagnostics for a shell wind
 *
 * @details
 *
 * This diagnostic information is for shell wind models, which should only
 * be used for testing purposes. This function is used to output data relating
 * to the abundances of ions. Note that this is very dependent on the peculiar
 * structure of the single shell model, which has only a single element in the
 * wind.
 *
 **********************************************************/

void
shell_output_wind_update_diagnostics (double xsum, double psum, double fsum, double csum, double icsum, double lsum, double ausum,
                                      double chexsum, double cool_sum, double lum_sum)
{
  int ndom;
  int i;
  int n;
  int nn;
  int m;
  int nshell;
  int first_ion_index;
  int last_ion_index;
  double total_density;
  double lum_h_line;
  double lum_he_line;
  double lum_c_line;
  double lum_n_line;
  double lum_o_line;
  double lum_fe_line;
  double agn_ip;

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    if (zdom[ndom].wind_type == SHELL)
    {
      /* nshell is the plasma cell that correspond to the second wind cell for the shell_wind model */
      nshell = wmain[zdom[ndom].nstart + 1].nplasma;
      n = plasmamain[nshell].nwind;
      WindPtr w = &wmain[n];
      for (i = 0; i < geo.nxfreq; i++)
      {                         /*loop over number of bands */
        Log
          ("Band %i f1 %e f2 %e model %i pl_alpha %f pl_log_w %e exp_t %e exp_w %e\n",
           i, geo.xfreq[i], geo.xfreq[i + 1],
           plasmamain[nshell].spec_mod_type[i],
           plasmamain[nshell].pl_alpha[i], plasmamain[nshell].pl_log_w[i], plasmamain[nshell].exp_temp[i], plasmamain[nshell].exp_w[i]);
      }
      /* Get some line diagnostics */

      lum_h_line = 0.0;
      lum_he_line = 0.0;
      lum_c_line = 0.0;
      lum_n_line = 0.0;
      lum_o_line = 0.0;
      lum_fe_line = 0.0;

      for (i = 0; i < nlines; i++)
      {
        if (lin_ptr[i]->z == 1)
          lum_h_line = lum_h_line + lin_ptr[i]->pow;
        else if (lin_ptr[i]->z == 2)
          lum_he_line = lum_he_line + lin_ptr[i]->pow;
        else if (lin_ptr[i]->z == 6)
          lum_c_line = lum_c_line + lin_ptr[i]->pow;
        else if (lin_ptr[i]->z == 7)
          lum_n_line = lum_n_line + lin_ptr[i]->pow;
        else if (lin_ptr[i]->z == 8)
          lum_o_line = lum_o_line + lin_ptr[i]->pow;
        else if (lin_ptr[i]->z == 26)
          lum_fe_line = lum_fe_line + lin_ptr[i]->pow;
      }
      agn_ip = geo.const_agn * (((pow (50000 / HEV, geo.alpha_agn + 1.0)) - pow (100 / HEV, geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));
      agn_ip /= (w[n].r * w[n].r);
      agn_ip /= plasmamain[nshell].rho * rho2nh;
      /* Report luminosities, IP and other diagnositic quantities */
      Log
        ("OUTPUT Lum_agn= %e T_e= %e N_h= %e N_e= %e alpha= %f IP(sim_2010)= %e Measured_IP(cloudy)= %e Measured_Xi= %e distance= %e volume= %e mean_ds=%e\n",
         geo.lum_agn, plasmamain[nshell].t_e,
         plasmamain[nshell].rho * rho2nh, plasmamain[nshell].ne,
         geo.alpha_agn, agn_ip, plasmamain[nshell].ip,
         plasmamain[nshell].xi, w[n].r, w[n].vol, plasmamain[nshell].mean_ds / plasmamain[nshell].n_ds);
      Log
        ("OUTPUT Absorbed_flux(ergs-1cm-3)    %8.2e  (photo %8.2e ff %8.2e compton %8.2e induced_compton %8.2e lines %8.2e auger %8.2e charge_ex %8.2e )\n",
         xsum / w[n].vol, psum / w[n].vol, fsum / w[n].vol, csum / w[n].vol, icsum / w[n].vol, lsum / w[n].vol, ausum / w[n].vol,
         chexsum / w[n].vol);
      /* Report all cooling mechanisms, including those that do not generate photons. */
      Log
        ("OUTPUT Wind_cooling(ergs-1cm-3)     %8.2e (recomb %8.2e ff %8.2e compton %8.2e DR %8.2e DI %8.2e adiabatic %8.2e lines %8.2e ) after update\n",
         cool_sum / w[n].vol, geo.cool_rr / w[n].vol,
         geo.lum_ff / w[n].vol, geo.cool_comp / w[n].vol,
         geo.cool_dr / w[n].vol, geo.cool_di / w[n].vol, geo.cool_adiabatic / w[n].vol, geo.lum_lines / w[n].vol);
      Log ("OUTPUT Wind_luminosity(ergs-1cm-3)     %8.2e (recomb %8.2e ff %8.2e lines %8.2e ) after update\n", lum_sum / w[n].vol,
           geo.lum_rr / w[n].vol, geo.lum_ff / w[n].vol, geo.lum_lines / w[n].vol);

      /* NSH 1701 calculate the recombination cooling for other elements */
      double c_rec = 0.0;
      double n_rec = 0.0;
      double o_rec = 0.0;
      double fe_rec = 0.0;
      double c_lum = 0.0;
      double n_lum = 0.0;
      double o_lum = 0.0;
      double fe_lum = 0.0;
      double c_dr = 0.0;
      double n_dr = 0.0;
      double o_dr = 0.0;
      double fe_dr = 0.0;
      double cool_dr_metals = 0.0;

      for (nn = 0; nn < nions; nn++)
      {
        if (ion[nn].z == 6)
        {
          c_dr = c_dr + plasmamain[nshell].cool_dr_ion[nn];
          c_rec = c_rec + plasmamain[nshell].cool_rr_ion[nn];
          c_lum = c_lum + plasmamain[nshell].lum_rr_ion[nn];
        }
        if (ion[nn].z == 7)
        {
          n_dr = n_dr + plasmamain[nshell].cool_dr_ion[nn];
          n_rec = n_rec + plasmamain[nshell].cool_rr_ion[nn];
          n_lum = n_lum + plasmamain[nshell].lum_rr_ion[nn];
        }
        if (ion[nn].z == 8)
        {
          o_dr = o_dr + plasmamain[nshell].cool_dr_ion[nn];
          o_rec = o_rec + plasmamain[nshell].cool_rr_ion[nn];
          o_lum = o_lum + plasmamain[nshell].lum_rr_ion[nn];
        }
        if (ion[nn].z == 26)
        {
          fe_dr = fe_dr + plasmamain[nshell].cool_dr_ion[nn];
          fe_rec = fe_rec + plasmamain[nshell].cool_rr_ion[nn];
          fe_lum = fe_lum + plasmamain[nshell].lum_rr_ion[nn];
        }
        if (ion[nn].z > 2)
        {
          cool_dr_metals = cool_dr_metals + plasmamain[nshell].cool_dr_ion[nn];
        }
      }

      Log ("Wind_line_cooling(ergs-1cm-3)  H %8.2e He %8.2e C %8.2e N %8.2e O %8.2e Fe %8.2e Metals %8.2e\n", lum_h_line / w[n].vol,
           lum_he_line / w[n].vol, lum_c_line / w[n].vol, lum_n_line / w[n].vol, lum_o_line / w[n].vol, lum_fe_line / w[n].vol);
      Log ("Wind_recomb_cooling(ergs-1cm-3)  H %8.2e He %8.2e C %8.2e N %8.2e O %8.2e Fe %8.2e Metals %8.2e\n",
           plasmamain[nshell].cool_rr_ion[0] / w[n].vol,
           (plasmamain[nshell].cool_rr_ion[2] + plasmamain[nshell].cool_rr_ion[3]) / w[n].vol, c_rec / w[n].vol, n_rec / w[n].vol,
           o_rec / w[n].vol, fe_rec / w[n].vol, plasmamain[nshell].cool_rr_metals / w[n].vol);
      Log ("Wind_recomb_lum(ergs-1cm-3)  H %8.2e He %8.2e C %8.2e N %8.2e O %8.2e Fe %8.2e Metals %8.2e\n",
           plasmamain[nshell].lum_rr_ion[0] / w[n].vol, (plasmamain[nshell].lum_rr_ion[2] + plasmamain[nshell].lum_rr_ion[3]) / w[n].vol,
           c_lum / w[n].vol, n_lum / w[n].vol, o_lum / w[n].vol, fe_lum / w[n].vol, plasmamain[nshell].lum_rr_metals / w[n].vol);
      Log ("Wind_dr_cooling(ergs-1cm-3)  H %8.2e He %8.2e C %8.2e N %8.2e O %8.2e Fe %8.2e Metals %8.2e\n",
           plasmamain[nshell].cool_dr_ion[0] / w[n].vol,
           (plasmamain[nshell].cool_dr_ion[2] + plasmamain[nshell].cool_dr_ion[3]) / w[n].vol, c_dr / w[n].vol, n_dr / w[n].vol,
           o_dr / w[n].vol, fe_dr / w[n].vol, cool_dr_metals / w[n].vol);
      /* 1110 NSH Added this line to report all cooling mechanisms, including those that do not generate photons. */
      Log ("Balance      Cooling=%8.2e Heating=%8.2e Lum=%8.2e T_e=%e after update\n", cool_sum, xsum, lum_sum, plasmamain[nshell].t_e);

      for (n = 0; n < nelements; n++)
      {
        first_ion_index = ele[n].firstion;
        last_ion_index = first_ion_index + ele[n].nions;
        Log ("OUTPUT %-5s ", ele[n].name);
        total_density = 0;
        for (m = first_ion_index; m < last_ion_index; m++)
        {
          total_density += plasmamain[nshell].density[m];
        }
        for (m = first_ion_index; m < last_ion_index; m++)
        {
          Log (" %8.2e", plasmamain[nshell].density[m] / total_density);
        }
        Log ("\n");
      }
      Log ("radial F_es %i %e \n", nshell, plasmamain[nshell].rad_force_es[0]);
      Log ("radial F_bf %i %e \n", nshell, plasmamain[nshell].rad_force_bf[0]);
      Log ("radial F_ff %i %e \n", nshell, plasmamain[nshell].rad_force_ff[0]);
      Log ("Radial Visible flux %e \n", plasmamain[nshell].F_vis[0]);
      Log ("Radial UV      flux %e \n", plasmamain[nshell].F_UV[0]);
      Log ("Radial Xray    flux %e \n", plasmamain[nshell].F_Xray[0]);
      Log ("Total Radial   flux %e \n", plasmamain[nshell].F_vis[0] + plasmamain[nshell].F_UV[0] + plasmamain[nshell].F_Xray[0]);
    }
  }
}
