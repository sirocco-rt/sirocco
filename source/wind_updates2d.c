/***********************************************************/
/** @file  wind_updates2d.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  This file contains the main routines for updating
 * and then reinitializing the wind after an ionization cycle
 *
 * The routines in this file are generic.  There is no dependence on a particlar wind model or
 * any coordinate system dependences.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

#define LINELEN 256

int num_updates = 0;


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
 *
 **********************************************************/
int
wind_update (w)
WindPtr (w);
{
  int n, i, j, ii;

  double xsum, psum, fsum, lsum, csum, icsum, ausum, chexsum;
  double cool_sum, lum_sum, rad_sum;    //1706 - the total cooling and luminosity of the wind
  double apsum, aausum, abstot; //Absorbed photon energy from PI and auger
  double c_rec, n_rec, o_rec, fe_rec;   //1701- NSH more outputs to show cooling from a few other elements
  double c_lum, n_lum, o_lum, fe_lum;   //1708- NSH and luminosities as well
  double cool_dr_metals;
  int nn;
  double flux_helper[3];
  double flux_persist_scale;

  double volume;
  double vol;
  char string[LINELEN];
  double t_r_old, t_e_old, dt_r, dt_e;
  double t_r_ave_old, t_r_ave, t_e_ave_old, t_e_ave;
  int iave, nmax_r, nmax_e;
  int nplasma, nshell;
  int nwind;
  int first, last, m;
  double tot, agn_ip;
  double lum_h_line, lum_he_line, lum_c_line, lum_n_line, lum_o_line, lum_fe_line;
  double c_dr, n_dr, o_dr, fe_dr;
  int my_nmin, my_nmax;         //Note that these variables are still used even without MPI on
  int ndom;
  FILE *fptr, *fptr2, *fptr3, *fptr4, *fptr5, *fopen ();        /*This is the file to communicate with zeus */
  double t_opt, t_UV, t_Xray, v_th, fhat[3];    /*This is the dimensionless optical depth parameter computed for communication to rad-hydro. */
  struct photon ptest;          //We need a test photon structure in order to compute t
  double kappa_es;              //The electron scattering opacity used for t

  fptr = fptr2 = fptr3 = fptr4 = fptr5 = NULL;

#ifdef MPI_ON
  int position, ndo, n_mpi, num_comm, n_mpi2;
  int size_of_commbuffer, size_of_specbuffer;
  char *commbuffer;

  /* JM 1409 -- Added for issue #110 to ensure correct reporting in parallel */
  int nmax_r_temp, nmax_e_temp;
  double dt_e_temp, dt_r_temp;


  /* The commbuffer needs to be larger enough to pack all variables in MPI_Pack and MPI_Unpack routines 
   * The cmombuffer is currently sized to be the minimum requred.  Therefore when variables are added, the
   * size must must be increased.
   *
   * The cell spectra are handled separately
   */

  size_of_commbuffer =
    8 * (n_inner_tot + 10 * nions + nlte_levels + 3 * nphot_total + 15 * NXBANDS + 133 * NFLUX_ANGLES) * (floor (NPLASMA / np_mpi_global) +
                                                                                                          1);

  size_of_specbuffer = 8 * NPLASMA * NBINS_IN_CELL_SPEC;

  size_of_commbuffer += size_of_specbuffer;

  commbuffer = (char *) malloc (size_of_commbuffer * sizeof (char));


  /* Initialise parallel only variables */
  nmax_r_temp = nmax_e_temp = -1;
  dt_e_temp = dt_r_temp = 0.0;

#endif
  dt_r = dt_e = 0.0;
  iave = 0;
  nmax_r = nmax_e = -1;
  t_r_ave_old = t_r_ave = t_e_ave_old = t_e_ave = 0.0;


  /* For MPI parallelisation, the following loop will be distributed over mutiple tasks.
     Note that the mynmim and mynmax variables are still used even without MPI on */
  my_nmin = 0;
  my_nmax = NPLASMA;
#ifdef MPI_ON
  ndo = get_parallel_nrange (rank_global, NPLASMA, np_mpi_global, &my_nmin, &my_nmax);
#endif

  /* Before we do anything let's record the average tr and te from the last cycle */
  /* JM 1409 -- Added for issue #110 to ensure correct reporting in parallel */
  for (n = 0; n < NPLASMA; n++)
  {
    t_r_ave_old += plasmamain[n].t_r;
    t_e_ave_old += plasmamain[n].t_e;
    /* macro-atom estimators need to be normalised for all cells. 
       Note they should have already been averaged over threads here */
    if (geo.rt_mode == RT_MODE_MACRO && geo.macro_simple == FALSE)
    {
      nwind = plasmamain[n].nwind;
      normalise_macro_estimators (nwind);

      /* force recalculation of kpacket rates and matrices, if applicable */
      macromain[n].kpkt_rates_known = FALSE;
      macromain[n].matrix_rates_known = FALSE;
    }
  }

  /* we now know how many cells this thread has to process - note this will be
     0-NPLASMA in serial mode */

  for (n = my_nmin; n < my_nmax; n++)
  {




    nwind = plasmamain[n].nwind;
    volume = w[nwind].vol;

    /* Skip cells that are partially in the wind these are not to be included
       in the calculation */

    if (modes.partial_cells == PC_EXTEND && wmain[nwind].inwind == W_PART_INWIND)
    {
      continue;
    }



    /* Start with a call to the routine which normalises all the macro atom
       monte carlo radiation field estimators. It's best to do this first since
       some of the estimators include temperature terms (stimulated correction
       terms) which were included during the monte carlo simulation so we want
       to be sure that the SAME temperatures are used here. (SS - Mar 2004). */

    /* Store some information so one can determine how much the temps are changing */
    t_r_old = plasmamain[n].t_r;
    t_e_old = plasmamain[n].t_e;
    iave++;

    if (plasmamain[n].ntot < 100)
    {
      Log
        ("!!wind_update: Cell %4d Dom %d  Vol. %8.2e r %8.2e theta %8.2e has only %4d photons\n",
         n, w[nwind].ndom, volume, w[nwind].rcen, w[nwind].thetacen, plasmamain[n].ntot);
    }

    /* this routine normalises the unbanded and banded estimators for simple atoms in this cell */
    normalise_simple_estimators (&plasmamain[n]);

    /* If geo.adiabatic is true, then calculate the adiabatic cooling using the current, i.e
     * previous value of t_e.  Note that this may not be  best way to determine the cooling.
     * Changes made here should also be reflected in wind2d.c.  At present, adiabatic cooling
     * is not included in updates to the temperature, even if the adiabatic cooling is calculated
     * here. 04nov -- ksl
     * 05apr -- ksl -- The index being used was incorrect.  This has been fixed now
     * 11sep -- nsh -- The index for the wind (&w) for adiabatic cooling was incorrect -
     * was being called with the plasma cell rather than the approriate wind cell - fixed
     * old: adiabatic_cooling (&w[n], plasmamain[n].t_e);
     */

    if (geo.adiabatic)
      plasmamain[n].cool_adiabatic = adiabatic_cooling (&w[nwind], plasmamain[n].t_e);
    else
      plasmamain[n].cool_adiabatic = 0.0;

    if (geo.nonthermal)
      plasmamain[n].heat_shock = shock_heating (&w[nwind]);
    else
      plasmamain[n].heat_shock = 0;

    /* Calculate the densities in various ways depending on the ioniz_mode */

    ion_abundances (&plasmamain[n], geo.ioniz_mode);



    /* Perform checks to see how much temperatures have changed in this iteration */

    if ((fabs (t_r_old - plasmamain[n].t_r)) > fabs (dt_r))
    {
      dt_r = plasmamain[n].t_r - t_r_old;
      nmax_r = n;
    }
    if ((fabs (t_e_old - plasmamain[n].t_e)) > fabs (dt_e))
    {
      dt_e = plasmamain[n].t_e - t_e_old;
      nmax_e = n;
    }
    t_r_ave += plasmamain[n].t_r;
    t_e_ave += plasmamain[n].t_e;
  }



  /*This is the end of the update loop that is parallised. We now need to exchange data between the tasks. */
#ifdef MPI_ON
  for (n_mpi = 0; n_mpi < np_mpi_global; n_mpi++)
  {
    position = 0;

    if (rank_global == n_mpi)
    {
      Log ("MPI task %d is working on cells %d to max %d (total size %d).\n", rank_global, my_nmin, my_nmax, NPLASMA);
      MPI_Pack (&ndo, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
      for (n = my_nmin; n < my_nmax; n++)
      {
        MPI_Pack (&n, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].nwind, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].nplasma, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ne, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].rho, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].vol, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].xgamma, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].density, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].partition, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].levden, nlte_levels, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].kappa_ff_factor, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].nscat_es, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].recomb_simple, nphot_total, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].recomb_simple_upweight, nphot_total, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].kpkt_abs, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].kbf_use, nphot_total, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].kbf_nuse, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].t_r, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].t_r_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].t_e, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].t_e_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].dt_e, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].dt_e_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_tot, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].abs_tot, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_tot_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_lines, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_ff, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_comp, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_ind_comp, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_lines_macro, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_photo_macro, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_photo, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_auger, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].abs_photo, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].abs_auger, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_z, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].w, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot_star, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot_bl, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot_disk, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot_wind, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot_agn, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].mean_ds, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].n_ds, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].nrad, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].nioniz, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].ioniz, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].inner_ioniz, n_inner_tot, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].recomb, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].scatters, nions, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].xscatters, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].heat_ion, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].heat_inner_ion, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].cool_rr_ion, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].lum_rr_ion, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].j, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].j_direct, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].j_scatt, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ave_freq, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_tot, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].xj, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].xave_freq, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].xsd_freq, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].nxtot, NXBANDS, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].F_vis, 4, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].F_UV, 4, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].F_Xray, 4, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].max_freq, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_lines, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_ff, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_adiabatic, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].comp_nujnu, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_comp, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_dr, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_ch_ex, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_di, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_rr, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_rr, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_rr_metals, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_rr_metals, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_tot, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_tot_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_tot_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_lines_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_ff_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_adiabatic_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_comp_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_dr_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_di_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_rr_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_rr_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_rr_metals_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_tot_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_shock, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].dmo_dt, 3, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].rad_force_es, 4, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].rad_force_ff, 4, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].rad_force_bf, 4, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].gain, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].converge_t_r, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].converge_t_e, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].converge_hc, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].trcheck, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].techeck, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].hccheck, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].converge_whole, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].converging, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].spec_mod_type, NXBANDS, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].pl_alpha, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].pl_log_w, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].exp_temp, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].exp_w, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].fmin_mod, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].fmax_mod, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ip, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ip_direct, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ip_scatt, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].xi, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].bf_simple_ionpool_in, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].bf_simple_ionpool_out, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].cell_spec_flux, NBINS_IN_CELL_SPEC, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].F_UV_ang_x, NFLUX_ANGLES, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].F_UV_ang_y, NFLUX_ANGLES, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].F_UV_ang_z, NFLUX_ANGLES, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&dt_e, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&dt_r, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&nmax_e, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&nmax_r, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
      }

      Log ("MPI task %d broadcasting plasma update information.\n", rank_global);
    }

    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Bcast (commbuffer, size_of_commbuffer, MPI_PACKED, n_mpi, MPI_COMM_WORLD);
    MPI_Barrier (MPI_COMM_WORLD);
    Log_silent ("MPI task %d survived broadcasting plasma update information.\n", rank_global);

    position = 0;

    if (rank_global != n_mpi)
    {
      MPI_Unpack (commbuffer, size_of_commbuffer, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
      for (n_mpi2 = 0; n_mpi2 < num_comm; n_mpi2++)
      {
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &n, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].nwind, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].nplasma, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ne, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].rho, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].vol, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].xgamma, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].density, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].partition, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].levden, nlte_levels, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].kappa_ff_factor, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].nscat_es, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].recomb_simple, nphot_total, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].recomb_simple_upweight, nphot_total, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].kpkt_abs, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].kbf_use, nphot_total, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].kbf_nuse, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_r_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].dt_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].dt_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].abs_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_tot_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_ind_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_lines_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_photo_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_photo, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_auger, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].abs_photo, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].abs_auger, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_z, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].w, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_star, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_bl, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_disk, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_wind, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_agn, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].mean_ds, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].n_ds, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].nrad, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].nioniz, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].ioniz, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].inner_ioniz, n_inner_tot, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].recomb, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].scatters, nions, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].xscatters, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].heat_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].heat_inner_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].cool_rr_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].lum_rr_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].j, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].j_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].j_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ave_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].xj, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].xave_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].xsd_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].nxtot, NXBANDS, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].F_vis, 4, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].F_UV, 4, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].F_Xray, 4, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].max_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_adiabatic, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].comp_nujnu, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_dr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_ch_ex, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_di, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_rr_metals, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_rr_metals, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_tot_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_tot_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_lines_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_ff_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_adiabatic_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_comp_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_dr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_di_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_rr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_rr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_rr_metals_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_tot_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_shock, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].dmo_dt, 3, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].rad_force_es, 4, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].rad_force_ff, 4, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].rad_force_bf, 4, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].gain, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_hc, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].trcheck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].techeck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].hccheck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_whole, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].converging, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].spec_mod_type, NXBANDS, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].pl_alpha, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].pl_log_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].exp_temp, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].exp_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].fmin_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].fmax_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ip, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ip_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ip_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].xi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].bf_simple_ionpool_in, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].bf_simple_ionpool_out, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].cell_spec_flux, NBINS_IN_CELL_SPEC, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].F_UV_ang_x, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].F_UV_ang_y, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].F_UV_ang_z, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &dt_e_temp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &dt_r_temp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &nmax_e_temp, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &nmax_r_temp, 1, MPI_INT, MPI_COMM_WORLD);

        /* JM 1409 -- Altered for issue #110 to ensure correct reporting in parallel */
        if (fabs (dt_e_temp) >= fabs (dt_e))
        {
          /* Check if any other threads found a higher maximum for te */
          dt_e = dt_e_temp;
          nmax_e = nmax_e_temp;
        }

        if (fabs (dt_r_temp) >= fabs (dt_r))
        {
          /* Check if any other threads found a higher maximum for tr */
          dt_r = dt_r_temp;
          nmax_r = nmax_r_temp;
        }

        t_r_ave += plasmamain[n].t_r;
        t_e_ave += plasmamain[n].t_e;
        iave++;

      }

    }

  }
  free (commbuffer);
#endif


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
      cylind_extend_density (ndom, w);
    else if (zdom[ndom].coord_type == RTHETA)
      rtheta_extend_density (ndom, w);
    else if (zdom[ndom].coord_type == SPHERICAL)
      spherical_extend_density (ndom, w);
    else if (zdom[ndom].coord_type == CYLVAR)
      cylvar_extend_density (ndom, w);
    else
    {
      Error ("Wind_update2d: Unknown coordinate type %d for domain %d \n", zdom[ndom].coord_type, ndom);
      Exit (0);
    }
  }
  flux_persist_scale = 0.5;     //The ammount of the latest flux that gets added into the persistent flux

  /* update the persistent fluxes */
  for (nplasma = 0; nplasma < NPLASMA; nplasma++)
  {
    if (geo.wcycle == 0)        //If this is the first time through, then the persistent flux is empty.
    {
      vadd (plasmamain[nplasma].F_vis_persistent, plasmamain[nplasma].F_vis, plasmamain[nplasma].F_vis_persistent);
      vadd (plasmamain[nplasma].F_UV_persistent, plasmamain[nplasma].F_UV, plasmamain[nplasma].F_UV_persistent);
      vadd (plasmamain[nplasma].F_Xray_persistent, plasmamain[nplasma].F_Xray, plasmamain[nplasma].F_Xray_persistent);
      vadd (plasmamain[nplasma].rad_force_bf_persist, plasmamain[nplasma].rad_force_bf, plasmamain[nplasma].rad_force_bf_persist);


      for (n = 0; n < NFLUX_ANGLES; n++)
      {
        plasmamain[nplasma].F_UV_ang_x_persist[n] = plasmamain[nplasma].F_UV_ang_x_persist[n] + plasmamain[nplasma].F_UV_ang_x[n];
        plasmamain[nplasma].F_UV_ang_y_persist[n] = plasmamain[nplasma].F_UV_ang_y_persist[n] + plasmamain[nplasma].F_UV_ang_y[n];
        plasmamain[nplasma].F_UV_ang_z_persist[n] = plasmamain[nplasma].F_UV_ang_z_persist[n] + plasmamain[nplasma].F_UV_ang_z[n];
      }
    }
    else
    {
      rescale (plasmamain[nplasma].F_vis_persistent, (1 - flux_persist_scale), plasmamain[nplasma].F_vis_persistent);
      rescale (plasmamain[nplasma].F_vis, flux_persist_scale, flux_helper);
      vadd (plasmamain[nplasma].F_vis_persistent, flux_helper, plasmamain[nplasma].F_vis_persistent);

      rescale (plasmamain[nplasma].F_UV_persistent, (1 - flux_persist_scale), plasmamain[nplasma].F_UV_persistent);
      rescale (plasmamain[nplasma].F_UV, flux_persist_scale, flux_helper);
      vadd (plasmamain[nplasma].F_UV_persistent, flux_helper, plasmamain[nplasma].F_UV_persistent);

      rescale (plasmamain[nplasma].F_Xray_persistent, (1 - flux_persist_scale), plasmamain[nplasma].F_Xray_persistent);
      rescale (plasmamain[nplasma].F_Xray, flux_persist_scale, flux_helper);
      vadd (plasmamain[nplasma].F_Xray_persistent, flux_helper, plasmamain[nplasma].F_Xray_persistent);

      rescale (plasmamain[nplasma].rad_force_bf_persist, (1 - flux_persist_scale), plasmamain[nplasma].rad_force_bf_persist);
      rescale (plasmamain[nplasma].rad_force_bf, flux_persist_scale, flux_helper);
      vadd (plasmamain[nplasma].rad_force_bf_persist, flux_helper, plasmamain[nplasma].rad_force_bf_persist);


      for (n = 0; n < NFLUX_ANGLES; n++)
      {

        plasmamain[nplasma].F_UV_ang_x_persist[n] = plasmamain[nplasma].F_UV_ang_x_persist[n] * (1 - flux_persist_scale);
        plasmamain[nplasma].F_UV_ang_y_persist[n] = plasmamain[nplasma].F_UV_ang_y_persist[n] * (1 - flux_persist_scale);
        plasmamain[nplasma].F_UV_ang_z_persist[n] = plasmamain[nplasma].F_UV_ang_z_persist[n] * (1 - flux_persist_scale);

        plasmamain[nplasma].F_UV_ang_x_persist[n] =
          plasmamain[nplasma].F_UV_ang_x_persist[n] + plasmamain[nplasma].F_UV_ang_x[n] * flux_persist_scale;
        plasmamain[nplasma].F_UV_ang_y_persist[n] =
          plasmamain[nplasma].F_UV_ang_y_persist[n] + plasmamain[nplasma].F_UV_ang_y[n] * flux_persist_scale;
        plasmamain[nplasma].F_UV_ang_z_persist[n] =
          plasmamain[nplasma].F_UV_ang_z_persist[n] + plasmamain[nplasma].F_UV_ang_z[n] * flux_persist_scale;


      }

    }
    plasmamain[nplasma].F_vis_persistent[3] = length (plasmamain[nplasma].F_vis_persistent);
    plasmamain[nplasma].F_UV_persistent[3] = length (plasmamain[nplasma].F_UV_persistent);
    plasmamain[nplasma].F_Xray_persistent[3] = length (plasmamain[nplasma].F_Xray_persistent);
    plasmamain[nplasma].rad_force_bf_persist[3] = length (plasmamain[nplasma].rad_force_bf_persist);

  }

  /*
     for (nplasma = 0; nplasma < NPLASMA; nplasma++)
     {
     if (geo.wcycle == 0)
     {
     vadd (plasmamain[nplasma].F_vis_persistent, plasmamain[nplasma].F_vis, plasmamain[nplasma].F_vis_persistent);
     vadd (plasmamain[nplasma].F_UV_persistent, plasmamain[nplasma].F_UV, plasmamain[nplasma].F_UV_persistent);
     vadd (plasmamain[nplasma].F_Xray_persistent, plasmamain[nplasma].F_Xray, plasmamain[nplasma].F_Xray_persistent);
     for (n = 0; n < NFLUX_ANGLES; n++)
     {
     plasmamain[nplasma].F_UV_ang_x_persist[n] = plasmamain[nplasma].F_UV_ang_x_persist[n] + plasmamain[nplasma].F_UV_ang_x[n];
     plasmamain[nplasma].F_UV_ang_y_persist[n] = plasmamain[nplasma].F_UV_ang_y_persist[n] + plasmamain[nplasma].F_UV_ang_y[n];
     plasmamain[nplasma].F_UV_ang_z_persist[n] = plasmamain[nplasma].F_UV_ang_z_persist[n] + plasmamain[nplasma].F_UV_ang_z[n];
     }
     }
     else
     {
     rescale (plasmamain[nplasma].F_vis_persistent, geo.wcycle, flux_helper);
     vadd (flux_helper, plasmamain[nplasma].F_vis, plasmamain[nplasma].F_vis_persistent);
     rescale (plasmamain[nplasma].F_vis_persistent, 1. / (geo.wcycle + 1), plasmamain[nplasma].F_vis_persistent);

     rescale (plasmamain[nplasma].F_UV_persistent, geo.wcycle, flux_helper);
     vadd (flux_helper, plasmamain[nplasma].F_UV, plasmamain[nplasma].F_UV_persistent);
     rescale (plasmamain[nplasma].F_UV_persistent, 1. / (geo.wcycle + 1), plasmamain[nplasma].F_UV_persistent);

     rescale (plasmamain[nplasma].F_Xray_persistent, geo.wcycle, flux_helper);
     vadd (flux_helper, plasmamain[nplasma].F_Xray, plasmamain[nplasma].F_Xray_persistent);
     rescale (plasmamain[nplasma].F_Xray_persistent, 1. / (geo.wcycle + 1), plasmamain[nplasma].F_Xray_persistent);


     for (n = 0; n < NFLUX_ANGLES; n++)
     {
     plasmamain[nplasma].F_UV_ang_x_persist[n] = plasmamain[nplasma].F_UV_ang_x_persist[n] * geo.wcycle;
     plasmamain[nplasma].F_UV_ang_y_persist[n] = plasmamain[nplasma].F_UV_ang_y_persist[n] * geo.wcycle;
     plasmamain[nplasma].F_UV_ang_z_persist[n] = plasmamain[nplasma].F_UV_ang_z_persist[n] * geo.wcycle;

     plasmamain[nplasma].F_UV_ang_x_persist[n] = plasmamain[nplasma].F_UV_ang_x_persist[n] + plasmamain[nplasma].F_UV_ang_x[n];
     plasmamain[nplasma].F_UV_ang_y_persist[n] = plasmamain[nplasma].F_UV_ang_y_persist[n] + plasmamain[nplasma].F_UV_ang_y[n];
     plasmamain[nplasma].F_UV_ang_z_persist[n] = plasmamain[nplasma].F_UV_ang_z_persist[n] + plasmamain[nplasma].F_UV_ang_z[n];

     plasmamain[nplasma].F_UV_ang_x_persist[n] = plasmamain[nplasma].F_UV_ang_x_persist[n] / (geo.wcycle + 1);
     plasmamain[nplasma].F_UV_ang_y_persist[n] = plasmamain[nplasma].F_UV_ang_y_persist[n] / (geo.wcycle + 1);
     plasmamain[nplasma].F_UV_ang_z_persist[n] = plasmamain[nplasma].F_UV_ang_z_persist[n] / (geo.wcycle + 1);
     }

     }
     plasmamain[nplasma].F_vis_persistent[3] = length (plasmamain[nplasma].F_vis_persistent);
     plasmamain[nplasma].F_UV_persistent[3] = length (plasmamain[nplasma].F_UV_persistent);
     plasmamain[nplasma].F_Xray_persistent[3] = length (plasmamain[nplasma].F_Xray_persistent);

     }
   */

  /* Finished updating region outside of wind */

  num_updates++;
  strcpy (string, "");
  sprintf (string, "# Wind update: Number %d", num_updates);

  if (modes.zeus_connect == 1 && geo.hydro_domain_number > -1)  //If we are running in zeus connect mode - we open a file for heatcool rates
  {
    Log ("Outputting heatcool file for connecting to zeus\n");
    fptr = fopen ("py_heatcool.dat", "w");
    fptr2 = fopen ("py_flux.dat", "w");
    fptr3 = fopen ("py_ion_data.dat", "w");
    fptr4 = fopen ("py_spec_data.dat", "w");
    fptr5 = fopen ("py_pcon_data.dat", "w");

    fprintf (fptr,
             "i j rcen thetacen vol temp xi ne heat_xray heat_comp heat_lines heat_ff cool_comp cool_lines cool_ff rho n_h rad_f_w rad_f_phi rad_f_z bf_f_w bf_f_phi bf_f_z\n");
    fprintf (fptr2, "i j F_vis_x F_vis_y F_vis_z F_UV_x F_UV_y F_UV_z F_Xray_x F_Xray_y F_Xray_z\n");   //directional flux by band

    fprintf (fptr3, "nions %i\n", nions);
    for (i = 0; i < nions; i++)
    {
      fprintf (fptr3, "ion %i %s %i %i\n", i, ele[ion[i].nelem].name, ion[i].z, ion[i].istate);
    }
    fprintf (fptr3, "nplasma %i\n", NPLASMA);

    fprintf (fptr4, "model %i\n", geo.ioniz_mode);
    fprintf (fptr4, "nbands %i\n", geo.nxfreq);
    fprintf (fptr4, "nplasma %i\n", NPLASMA);
    for (i = 0; i < geo.nxfreq + 1; i++)
      fprintf (fptr4, "%e ", geo.xfreq[i]);     //hard wired band edges
    fprintf (fptr4, "\n ");

    fprintf (fptr5, "nplasma %i\n", NPLASMA);

  }

  /* Check the balance between the absorbed and the emitted flux */

  /* NSH 0717 - ensure the cooling and luminosities reflect the current temperature */

  cool_sum = wind_cooling ();   /* We call wind_cooling here to obtain an up to date set of cooling rates */
  lum_sum = wind_luminosity (0.0, VERY_BIG, MODE_CMF_TIME);     /* and we also call wind_luminosity to get the luminosities */



  xsum = psum = ausum = lsum = fsum = csum = icsum = apsum = aausum = abstot = chexsum = 0;     //1108 NSH zero the new csum counter for compton heating

  for (nplasma = 0; nplasma < NPLASMA; nplasma++)
  {
    if (sane_check (plasmamain[nplasma].heat_tot))
      Error ("wind_update:sane_check w(%d).heat_tot is %e\n", nplasma, plasmamain[nplasma].heat_tot);
    if (sane_check (plasmamain[nplasma].heat_photo))
      Error ("wind_update:sane_check w(%d).heat_photo is %e\n", nplasma, plasmamain[nplasma].heat_photo);
    if (sane_check (plasmamain[nplasma].heat_auger))
      Error ("wind_update:sane_check w(%d).heat_auger is %e\n", nplasma, plasmamain[nplasma].heat_auger);
    if (sane_check (plasmamain[nplasma].heat_photo_macro))
      Error ("wind_update:sane_check w(%d).heat_photo_macro is %e\n", nplasma, plasmamain[nplasma].heat_photo_macro);
    if (sane_check (plasmamain[nplasma].heat_ff))
      Error ("wind_update:sane_check w(%d).heat_ff is %e\n", nplasma, plasmamain[nplasma].heat_ff);
    if (sane_check (plasmamain[nplasma].heat_lines))
      Error ("wind_update:sane_check w(%d).heat_lines is %e\n", nplasma, plasmamain[nplasma].heat_lines);
    if (sane_check (plasmamain[nplasma].heat_lines_macro))
      Error ("wind_update:sane_check w(%d).heat_lines_macro is %e\n", nplasma, plasmamain[nplasma].heat_lines_macro);
    /* 1108 NSH extra Sane check for compton heating */
    if (sane_check (plasmamain[nplasma].heat_comp))
      Error ("wind_update:sane_check w(%d).heat_comp is %e\n", nplasma, plasmamain[nplasma].heat_comp);

    abstot += plasmamain[nplasma].abs_tot;
    xsum += plasmamain[nplasma].heat_tot;
    psum += plasmamain[nplasma].heat_photo;
    ausum += plasmamain[nplasma].heat_auger;
    fsum += plasmamain[nplasma].heat_ff;
    lsum += plasmamain[nplasma].heat_lines;
    csum += plasmamain[nplasma].heat_comp;
    icsum += plasmamain[nplasma].heat_ind_comp;
    apsum += plasmamain[nplasma].abs_photo;
    aausum += plasmamain[nplasma].abs_auger;
    chexsum += plasmamain[nplasma].heat_ch_ex;


    plasmamain[nplasma].cool_tot_ioniz = plasmamain[nplasma].cool_tot;
    plasmamain[nplasma].lum_ff_ioniz = plasmamain[nplasma].lum_ff;
    plasmamain[nplasma].cool_rr_ioniz = plasmamain[nplasma].cool_rr;
    plasmamain[nplasma].lum_rr_ioniz = plasmamain[nplasma].lum_rr;
    plasmamain[nplasma].cool_rr_metals_ioniz = plasmamain[nplasma].cool_rr_metals;
    plasmamain[nplasma].lum_lines_ioniz = plasmamain[nplasma].lum_lines;
    plasmamain[nplasma].cool_comp_ioniz = plasmamain[nplasma].cool_comp;
    plasmamain[nplasma].cool_dr_ioniz = plasmamain[nplasma].cool_dr;
    plasmamain[nplasma].cool_di_ioniz = plasmamain[nplasma].cool_di;
    plasmamain[nplasma].lum_tot_ioniz = plasmamain[nplasma].lum_tot;
    plasmamain[nplasma].cool_adiabatic_ioniz = plasmamain[nplasma].cool_adiabatic;

  }



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

  /* Added this system which counts number of times two situations occur (See #91)
     We only report these every 100,000 times (one can typically get ) */
  Log ("wind_update: note, errors from mean intensity can be high in a working model\n");
  Log
    ("wind_update: can be a problem with photon numbers if there are also errors from spectral_estimators and low photon number warnings\n");
  Log ("wind_update: mean_intensity: %8.4e occurrences, this cycle, this thread of 'no model exists in a band'\n", (1.0 * nerr_no_Jmodel));
  Log
    ("wind_update: mean intensity: %8.4e occurrences, this cycle, this thread of 'photon freq is outside frequency range of spectral model'\n",
     (1.0 * nerr_Jmodel_wrong_freq));


  /* zero the counters which record diagnositics from mean_intensity */
  nerr_Jmodel_wrong_freq = 0;
  nerr_no_Jmodel = 0;






  if (modes.zeus_connect == 1 && geo.hydro_domain_number > -1)  //If we are running in zeus connect mode, we output heating and cooling rates.
  {
    for (nwind = zdom[geo.hydro_domain_number].nstart; nwind < zdom[geo.hydro_domain_number].nstop; nwind++)
    {
      if (wmain[nwind].inwind >= 0)
      {
        nplasma = wmain[nwind].nplasma;
        wind_n_to_ij (geo.hydro_domain_number, plasmamain[nplasma].nwind, &i, &j);
        i = i - 1;              //There is a radial 'ghost zone' in python, we need to make our i,j agree with zeus
        vol = w[plasmamain[nplasma].nwind].vol;
        fprintf (fptr, "%d %d %e %e %e ", i, j, w[plasmamain[nplasma].nwind].rcen, w[plasmamain[nplasma].nwind].thetacen / RADIAN, vol);        //output geometric things
        fprintf (fptr, "%e %e %e ", plasmamain[nplasma].t_e, plasmamain[nplasma].xi, plasmamain[nplasma].ne);   //output temp, xi and ne to ease plotting of heating rates
        fprintf (fptr, "%e ", (plasmamain[nplasma].heat_photo + plasmamain[nplasma].heat_auger) / vol); //Xray heating - or photoionization
        fprintf (fptr, "%e ", (plasmamain[nplasma].heat_comp) / vol);   //Compton heating
        fprintf (fptr, "%e ", (plasmamain[nplasma].heat_lines) / vol);  //Line heating 28/10/15 - not currently used in zeus
        fprintf (fptr, "%e ", (plasmamain[nplasma].heat_ff) / vol);     //FF heating 28/10/15 - not currently used in zeus
        fprintf (fptr, "%e ", (plasmamain[nplasma].cool_comp) / vol);   //Compton cooling
        fprintf (fptr, "%e ", (plasmamain[nplasma].lum_lines + plasmamain[nplasma].cool_rr + plasmamain[nplasma].cool_dr) / vol);       //Line cooling must include all recombination cooling
        fprintf (fptr, "%e ", (plasmamain[nplasma].lum_ff) / vol);      //ff cooling
        fprintf (fptr, "%e ", plasmamain[nplasma].rho); //density
        fprintf (fptr, "%e ", plasmamain[nplasma].rho * rho2nh);        //hydrogen number density
        fprintf (fptr, "%e ", plasmamain[nplasma].rad_force_es[0]);     //electron scattering radiation force in the w(x) direction
        fprintf (fptr, "%e ", plasmamain[nplasma].rad_force_es[1]);     //electron scattering radiation force in the phi(rotational) directionz direction
        fprintf (fptr, "%e ", plasmamain[nplasma].rad_force_es[2]);     //electron scattering radiation force in the z direction
        fprintf (fptr, "%e ", plasmamain[nplasma].rad_force_bf[0]);     //bound free scattering radiation force in the w(x) direction
        fprintf (fptr, "%e ", plasmamain[nplasma].rad_force_bf[1]);     //bound free scattering radiation force in the phi(rotational) direction
        fprintf (fptr, "%e \n", plasmamain[nplasma].rad_force_bf[2]);   //bound free scattering radiation force in the z direction
        fprintf (fptr2, "%d %d ", i, j);        //output geometric things               
        fprintf (fptr2, "%e %e %e ", plasmamain[nplasma].F_vis[0], plasmamain[nplasma].F_vis[1], plasmamain[nplasma].F_vis[2]); //directional flux by band
        fprintf (fptr2, "%e %e %e ", plasmamain[nplasma].F_UV[0], plasmamain[nplasma].F_UV[1], plasmamain[nplasma].F_UV[2]);    //directional flux by band
        fprintf (fptr2, "%e %e %e ", plasmamain[nplasma].F_Xray[0], plasmamain[nplasma].F_Xray[1], plasmamain[nplasma].F_Xray[2]);      //directional flux by band

        fprintf (fptr2, "\n");
        fprintf (fptr3, "%d %d ", i, j);        //output geometric things               
        for (ii = 0; ii < nions; ii++)
          fprintf (fptr3, "%e ", plasmamain[nplasma].density[ii]);
        fprintf (fptr3, "\n");

        fprintf (fptr4, "%d %d ", i, j);        //output geometric things       
        for (ii = 0; ii < geo.nxfreq; ii++)
          fprintf (fptr4, "%e %e %i %e %e %e %e ",
                   plasmamain[nplasma].fmin_mod[ii], plasmamain[nplasma].fmax_mod[ii], plasmamain[nplasma].spec_mod_type[ii],
                   plasmamain[nplasma].pl_log_w[ii], plasmamain[nplasma].pl_alpha[ii], plasmamain[nplasma].exp_w[ii],
                   plasmamain[nplasma].exp_temp[ii]);
        fprintf (fptr4, "\n ");


        //We need to compute the g factor for this cell and output it.


        v_th = pow ((2. * BOLTZMANN * plasmamain[nplasma].t_e / MPROT), 0.5);   //We need the thermal velocity for hydrogen
        stuff_v (w[plasmamain[nplasma].nwind].xcen, ptest.x);   //place our test photon at the centre of the cell
        ptest.grid = nwind;     //We need our test photon to know where it is 
        kappa_es = THOMPSON * plasmamain[nplasma].ne / plasmamain[nplasma].rho;

        //First for the optical band (up to 4000AA)     
        if (length (plasmamain[nplasma].F_vis) > 0.0)   //Only makes sense if flux in this band is non-zero
        {
          stuff_v (plasmamain[nplasma].F_vis, fhat);
          renorm (fhat, 1.);    //A unit vector in the direction of the flux - this can be treated as the lmn vector of a pretend photon
          stuff_v (fhat, ptest.lmn);    //place our test photon at the centre of the cell            
          t_opt = kappa_es * plasmamain[nplasma].rho * v_th / fabs (dvwind_ds_cmf (&ptest));
        }
        else
          t_opt = 0.0;          //Essentually a flag that there is no way of computing t (and hence M) in this cell.

        //Now for the UV band (up to 4000AA->100AA)                                             
        if (length (plasmamain[nplasma].F_UV) > 0.0)    //Only makes sense if flux in this band is non-zero
        {
          stuff_v (plasmamain[nplasma].F_UV, fhat);
          renorm (fhat, 1.);    //A unit vector in the direction of the flux - this can be treated as the lmn vector of a pretend photon
          stuff_v (fhat, ptest.lmn);    //place our test photon at the centre of the cell            
          t_UV = kappa_es * plasmamain[nplasma].rho * v_th / fabs (dvwind_ds_cmf (&ptest));
        }
        else
          t_UV = 0.0;           //Essentually a flag that there is no way of computing t (and hence M) in this cell.


        //And finally for the Xray band (up to 100AA and up)
        if (length (plasmamain[nplasma].F_Xray) > 0.0)  //Only makes sense if flux in this band is non-zero
        {
          stuff_v (plasmamain[nplasma].F_Xray, fhat);
          renorm (fhat, 1.);    //A unit vector in the direction of the flux - this can be treated as the lmn vector of a pretend photon
          stuff_v (fhat, ptest.lmn);    //place our test photon at the centre of the cell            
          t_Xray = kappa_es * plasmamain[nplasma].rho * v_th / fabs (dvwind_ds_cmf (&ptest));
        }
        else
          t_Xray = 0.0;         //Essentually a flag that there is no way of computing t (and hence M) in this cell.                

        fprintf (fptr5, "%i %i %e %e %e %e %e %e %e\n", i, j, plasmamain[nplasma].t_e, plasmamain[nplasma].rho,
                 plasmamain[nplasma].rho * rho2nh, plasmamain[nplasma].ne, t_opt, t_UV, t_Xray);
      }
    }
    fclose (fptr);
    fclose (fptr2);
    fclose (fptr3);
    fclose (fptr4);
    fclose (fptr5);
  }
  else if (modes.zeus_connect == 1 && geo.hydro_domain_number < 0)
  {
    Error ("wind_updates2d:  Attempting to access a hydro domain in a non hydro run - not writing out hydro file\n");
  }

  /* The lines differ only in that Wind_heating adds mechanical heating, that is adiabatic heating */

  Log
    ("!!wind_update: Absorbed flux    %8.2e  (photo %8.2e ff %8.2e compton %8.2e auger %8.2e induced_compton %8.2e lines %8.2e)\n",
     abstot, apsum, fsum, csum, aausum, icsum, lsum);

  Log
    ("!!wind_update: Wind heating     %8.2e  (photo %8.2e ff %8.2e compton %8.2e auger %8.2e induced_compton %8.2e lines %8.2e adiabatic %8.2e)\n",
     xsum + geo.heat_adiabatic, psum, fsum, csum, ausum, icsum, lsum, geo.heat_adiabatic);

  /* 1108 NSH added commands to report compton cooling 1110 removed,
   * As was the case above, there are two almost identical lines.  Wind_cooling includes processes that do not produce photons,
   * not-only adiabatic cooling, but also goe.cool_comp, geo_cool_dr and geo.cool_di */
  Log
    ("!!wind_update: Wind luminosity  %8.2e (recomb %8.2e ff %8.2e lines %8.2e) after update\n",
     lum_sum, geo.lum_rr, geo.lum_ff, geo.lum_lines);


  rad_sum = wind_luminosity (xband.f1[0], xband.f2[xband.nbands - 1], MODE_CMF_TIME);

  Log
    ("!!wind_update: Rad luminosity  %8.2e (recomb %8.2e ff %8.2e lines %8.2e) after update\n",
     rad_sum, geo.lum_rr, geo.lum_ff, geo.lum_lines);

  Log
    ("!!wind_update: Wind cooling     %8.2e (recomb %8.2e ff %8.2e compton %8.2e DR %8.2e DI %8.2e lines %8.2e adiabatic %8.2e) after update\n",
     cool_sum, geo.cool_rr, geo.lum_ff, geo.cool_comp, geo.cool_dr, geo.cool_di, geo.lum_lines, geo.cool_adiabatic);

  if (!modes.turn_off_upweighting_of_simple_macro_atoms)
  {
    /* If we have "indivisible packet" mode on but are using the 
       new BF_SIMPLE_EMISSIVITY_APPROACH then we report the flows into and out of the ion pool */
    if (geo.rt_mode == RT_MODE_MACRO)
      report_bf_simple_ionpool ();
  }



  if (modes.zeus_connect == 1 || modes.fixed_temp == 1) //There is no point in computing temperature changes, because we have fixed them!
  {
    Log ("!!wind_update: We are running in fixed temperature mode - no temperature report\n");
  }
  else
  {
    t_r_ave_old /= iave;
    t_e_ave_old /= iave;
    t_r_ave /= iave;
    t_e_ave /= iave;

    if (nmax_r != -1)
    {
      wind_n_to_ij (wmain[nmax_r].ndom, nmax_r, &i, &j);
      Log ("!!wind_update: Max change in t_r %6.0f at cell %4d (%d,%d)\n", dt_r, nmax_r, i, j);
      Log ("!!wind_update: Ave change in t_r %6.0f from %6.0f to %6.0f\n", (t_r_ave - t_r_ave_old), t_r_ave_old, t_r_ave);
    }
    else
      Log ("!!wind_update: t_r did not change in any cells this cycle\n");

    if (nmax_e != -1)
    {
      wind_n_to_ij (wmain[nmax_e].ndom, nmax_e, &i, &j);
      Log ("!!wind_update: Max change in t_e %6.0f at cell %4d (%d,%d)\n", dt_e, nmax_e, i, j);
      Log ("!!wind_update: Ave change in t_e %6.0f from %6.0f to %6.0f\n", (t_e_ave - t_e_ave_old), t_e_ave_old, t_e_ave);
    }
    else
      Log ("!!wind_update: t_e did not change in any cells this cycle\n");


    Log ("Summary  t_r  %6.0f   %6.0f  #t_r and dt_r on this update\n", t_r_ave, (t_r_ave - t_r_ave_old));
    Log ("Summary  t_e  %6.0f   %6.0f  #t_e and dt_e on this update\n", t_e_ave, (t_e_ave - t_e_ave_old));
  }

  check_convergence ();

  /* Summarize the radiative temperatures (ksl 04 mar) */

  xtemp_rad (w);

/* This next block is to allow the output of data relating to the abundances of ions when python is being tested
 * with thin shell mode.We will only want this to run if the wind mode is 9, for test or thin shell mode.
 *
 * Note that this section is very dependent on the peculiar structure of the single shell model, which has only
 * one element (namely element 2) in the wind.  nshell below correspond to that particular plasma shell. Recognizing
 * this was a key element to solving bug #412.
 */
  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {

    if (zdom[ndom].wind_type == SHELL)
    {

      /* nshell is the plasma cell that correspond to the second wind cell for the shell_wind model */
      nshell = wmain[zdom[ndom].nstart + 1].nplasma;
      n = plasmamain[nshell].nwind;
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

      c_rec = n_rec = o_rec = fe_rec = 0.0;
      c_lum = n_lum = o_lum = fe_lum = 0.0;
//OLD      h_dr = he_dr = c_dr = n_dr = o_dr = fe_dr = 0.0;
      c_dr = n_dr = o_dr = fe_dr = 0.0;
      cool_dr_metals = 0.0;

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
          cool_dr_metals = cool_dr_metals + plasmamain[nshell].cool_dr_ion[nn];
      }

      Log ("OUTPUT Wind_line_cooling(ergs-1cm-3)  H %8.2e He %8.2e C %8.2e N %8.2e O %8.2e Fe %8.2e Metals %8.2e\n", lum_h_line / w[n].vol,
           lum_he_line / w[n].vol, lum_c_line / w[n].vol, lum_n_line / w[n].vol, lum_o_line / w[n].vol, lum_fe_line / w[n].vol);
      Log ("OUTPUT Wind_recomb_cooling(ergs-1cm-3)  H %8.2e He %8.2e C %8.2e N %8.2e O %8.2e Fe %8.2e Metals %8.2e\n",
           plasmamain[nshell].cool_rr_ion[0] / w[n].vol, (plasmamain[nshell].cool_rr_ion[2] + plasmamain[nshell].cool_rr_ion[3]) / w[n].vol,
           c_rec / w[n].vol, n_rec / w[n].vol, o_rec / w[n].vol, fe_rec / w[n].vol, plasmamain[nshell].cool_rr_metals / w[n].vol);
      Log ("OUTPUT Wind_recomb_lum(ergs-1cm-3)  H %8.2e He %8.2e C %8.2e N %8.2e O %8.2e Fe %8.2e Metals %8.2e\n",
           plasmamain[nshell].lum_rr_ion[0] / w[n].vol, (plasmamain[nshell].lum_rr_ion[2] + plasmamain[nshell].lum_rr_ion[3]) / w[n].vol,
           c_lum / w[n].vol, n_lum / w[n].vol, o_lum / w[n].vol, fe_lum / w[n].vol, plasmamain[nshell].lum_rr_metals / w[n].vol);
      Log ("OUTPUT Wind_dr_cooling(ergs-1cm-3)  H %8.2e He %8.2e C %8.2e N %8.2e O %8.2e Fe %8.2e Metals %8.2e\n",
           plasmamain[nshell].cool_dr_ion[0] / w[n].vol, (plasmamain[nshell].cool_dr_ion[2] + plasmamain[nshell].cool_dr_ion[3]) / w[n].vol,
           c_dr / w[n].vol, n_dr / w[n].vol, o_dr / w[n].vol, fe_dr / w[n].vol, cool_dr_metals / w[n].vol);
      /* 1110 NSH Added this line to report all cooling mechanisms, including those that do not generate photons. */
      Log
        ("OUTPUT Balance      Cooling=%8.2e Heating=%8.2e Lum=%8.2e T_e=%e after update\n",
         cool_sum, xsum, lum_sum, plasmamain[nshell].t_e);

      for (n = 0; n < nelements; n++)
      {
        first = ele[n].firstion;
        last = first + ele[n].nions;
        Log ("OUTPUT %-5s ", ele[n].name);
        tot = 0;
        for (m = first; m < last; m++)
          tot += plasmamain[nshell].density[m];
        for (m = first; m < last; m++)
        {
          Log (" %8.2e", plasmamain[nshell].density[m] / tot);
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


  return (0);
}




/**********************************************************/
/**
 * @brief      zeros those portions of the wind which contain the radiation properties
 * 	of the wind, i.e those portions which should be set to zeroed when the structure of the
 * 	wind has been changed or when you simply want to start off a calculation in a known state
 *
 * @return    Always returns 0
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

int
wind_rad_init ()
{
  int n, i;
  int njump;
  double alpha_store;


  for (n = 0; n < NPLASMA; n++)
  {
    plasmamain[n].j = plasmamain[n].ave_freq = plasmamain[n].ntot = 0;
    plasmamain[n].j_direct = plasmamain[n].j_scatt = 0.0;
    plasmamain[n].ip = 0.0;
    plasmamain[n].xi = 0.0;

    plasmamain[n].ip_direct = plasmamain[n].ip_scatt = 0.0;
    plasmamain[n].mean_ds = 0.0;
    plasmamain[n].n_ds = 0;
    plasmamain[n].ntot_disk = plasmamain[n].ntot_agn = 0;
    plasmamain[n].ntot_star = plasmamain[n].ntot_bl = plasmamain[n].ntot_wind = 0;
    plasmamain[n].heat_tot = plasmamain[n].heat_ff = plasmamain[n].heat_photo = plasmamain[n].heat_lines = 0.0;
    plasmamain[n].abs_tot = plasmamain[n].abs_auger = plasmamain[n].abs_photo = 0.0;

    plasmamain[n].heat_z = 0.0;
    plasmamain[n].max_freq = 0.0;
    plasmamain[n].cool_tot = plasmamain[n].lum_tot = plasmamain[n].lum_lines = plasmamain[n].lum_ff = 0.0;
    plasmamain[n].cool_rr = plasmamain[n].cool_rr_metals = plasmamain[n].lum_rr = 0.0;
    plasmamain[n].nrad = plasmamain[n].nioniz = 0;
    plasmamain[n].comp_nujnu = -1e99;
    plasmamain[n].cool_comp = 0.0;
    plasmamain[n].heat_comp = 0.0;
    plasmamain[n].heat_ind_comp = 0.0;
    plasmamain[n].heat_auger = 0.0;
    plasmamain[n].heat_ch_ex = 0.0;

    /* zero the counters that record the flow into and out of the 
       ionization pool in indivisible packet mode */
    plasmamain[n].bf_simple_ionpool_out = 0.0;
    plasmamain[n].bf_simple_ionpool_in = 0.0;

    for (i = 0; i < 3; i++)
      plasmamain[n].dmo_dt[i] = 0.0;
    for (i = 0; i < 4; i++)
      plasmamain[n].rad_force_es[i] = 0.0;
    for (i = 0; i < 4; i++)
      plasmamain[n].rad_force_ff[i] = 0.0;
    for (i = 0; i < 4; i++)
      plasmamain[n].rad_force_bf[i] = 0.0;


    if (geo.rt_mode == RT_MODE_MACRO)
    {
      macromain[n].kpkt_rates_known = FALSE;
    }

/* Initialise  the frequency banded radiation estimators used for estimating the coarse spectra in each cell*/

    for (i = 0; i < geo.nxfreq; i++)
    {
      plasmamain[n].xj[i] = plasmamain[n].xave_freq[i] = plasmamain[n].nxtot[i] = 0;
      plasmamain[n].xsd_freq[i] = 0.0;
      plasmamain[n].fmin[i] = geo.xfreq[i + 1]; /* Set the minium frequency to the max frequency in the band */
      plasmamain[n].fmax[i] = geo.xfreq[i];     /* Set the maximum frequency to the min frequency in the band */
    }

    for (i = 0; i < NBINS_IN_CELL_SPEC; i++)
    {
      plasmamain[n].cell_spec_flux[i] = 0;
    }

    for (i = 0; i < 4; i++)
    {
      plasmamain[n].F_vis[i] = plasmamain[n].F_UV[i] = plasmamain[n].F_Xray[i] = 0.0;
      if (geo.wcycle == 0)
      {
        plasmamain[n].F_vis_persistent[i] = plasmamain[n].F_UV_persistent[i] = plasmamain[n].F_Xray_persistent[i] =
          plasmamain[n].rad_force_bf_persist[i] = 0.0;
      }
    }

    for (i = 0; i < NFLUX_ANGLES; i++)
    {
      if (geo.wcycle == 0)
      {
        plasmamain[n].F_UV_ang_x_persist[i] = 0.0;
        plasmamain[n].F_UV_ang_y_persist[i] = 0.0;
        plasmamain[n].F_UV_ang_z_persist[i] = 0.0;
      }
      if (i == 0 && n == 0)
        plasmamain[n].F_UV_ang_x[i] = 0.0;
      plasmamain[n].F_UV_ang_x[i] = 0.0;
      plasmamain[n].F_UV_ang_x[i] = 0.0;

    }





    for (i = 0; i < nions; i++)
    {
      plasmamain[n].ioniz[i] = plasmamain[n].recomb[i] = plasmamain[n].heat_ion[i] = plasmamain[n].cool_rr_ion[i] =
        plasmamain[n].lum_rr_ion[i] = plasmamain[n].heat_inner_ion[i] = 0.0;

    }
    for (i = 0; i < n_inner_tot; i++)
    {
      plasmamain[n].inner_ioniz[i] = 0.0;

    }
    /*Block added (Dec 08) to zero the auger rate estimators */
    /* commented out by NSH 2018 - removed code */
//    for (i = 0; i < nauger; i++)
//    {
//      plasmamain[n].gamma_inshl[i] = 0.0;
//    }

    /* Zero the Macro Atom estimators. */


    for (i = 0; i < nlevels_macro; i++)
    {
      for (njump = 0; njump < config[i].n_bbu_jump; njump++)
      {
        macromain[n].jbar[config[i].bbu_indx_first + njump] = 0.0;
      }
      for (njump = 0; njump < config[i].n_bfu_jump; njump++)
      {
        macromain[n].gamma[config[i].bfu_indx_first + njump] = 0.0;
        macromain[n].gamma_e[config[i].bfu_indx_first + njump] = 0.0;
        macromain[n].alpha_st[config[i].bfd_indx_first + njump] = 0.0;
        macromain[n].alpha_st_e[config[i].bfd_indx_first + njump] = 0.0;
      }


      for (njump = 0; njump < config[i].n_bfd_jump; njump++)
      {
        if (plasmamain[n].t_e > 1.0)
        {
          macromain[n].recomb_sp[config[i].bfd_indx_first + njump] = alpha_sp (&phot_top[config[i].bfd_jump[njump]], &plasmamain[n], 0);
          macromain[n].recomb_sp_e[config[i].bfd_indx_first + njump] = alpha_sp (&phot_top[config[i].bfd_jump[njump]], &plasmamain[n], 2);
        }
        else
        {
          macromain[n].recomb_sp[config[i].bfd_indx_first + njump] = 0.0;
          macromain[n].recomb_sp_e[config[i].bfd_indx_first + njump] = 0.0;
        }

      }
    }


    for (i = 0; i < ntop_phot; i++)
    {
      if ((geo.macro_simple == FALSE && phot_top[i].macro_info == TRUE) || geo.rt_mode == RT_MODE_2LEVEL)
      {
        plasmamain[n].recomb_simple[i] = 0.0;
        plasmamain[n].recomb_simple_upweight[i] = 1.0;
      }
      else
      {                         //we want a macro approach, but not for this ion so need recomb_simple
        plasmamain[n].recomb_simple[i] = alpha_store = alpha_sp (&phot_top[i], &plasmamain[n], 2);
        plasmamain[n].recomb_simple_upweight[i] = alpha_sp (&phot_top[i], &plasmamain[n], 1) / alpha_store;
      }
    }


    plasmamain[n].kpkt_emiss = 0.0;
    plasmamain[n].kpkt_abs = 0.0;
    for (i = 0; i < nlevels_macro; i++)
    {
      macromain[n].matom_abs[i] = 0.0;

      macromain[n].matom_emiss[i] = 0.0;

    }

  }


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
report_bf_simple_ionpool ()
{
  int n;
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

  return (0);
}
