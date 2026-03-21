/***********************************************************/
/** @file  communicate_plasma.c
 * @author EJP
 * @date   December 2023
 *
 * @brief Functions for communicating plasma properties
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

/* these next two numbers are the number of basic doubles and integers 
   in the plasma structure to be packed (i.e. things that are not arrays or pointers and
   have a size of MPI_INT or MPI_DOUBLE)*/
#define N_BASIC_DOUBLES 73
#define N_BASIC_INTS 22

/**********************************************************/
/**
 * @brief Broadcast the (initialised) plasma grid to all ranks
 *
 * @param [in] int n_start       The index of the first cell updated by this rank
 * @param [in] int n_stop        The index of the last cell updated by this rank
 * @param [in] int n_cells_rank  The number of cells this rank updated
 *
 * @details
 *
 * This should only be called once, after grid initialisation.
 *
 * The communication pattern and how the size of the communication buffer is
 * determined is documented in
 * $SIROCCO/docs/sphinx/source/developer/mpi_comms.rst.
 * To communicate a new variable, the communication buffer needs to be made
 * bigger and a new `MPI_Pack` and `MPI_Unpack` call need to be added. See the
 * developer documentation for more details.
 *
 **********************************************************/

void
broadcast_plasma_grid (const int n_start, const int n_stop, const int n_cells_rank)
{
#ifdef MPI_ON
  int i;
  int n_plasma;
  int current_rank;
  int num_comm;
  int position;

  PlasmaPtr cell;

  d_xsignal (files.root, "%-20s Begin communicating plasma grid\n", "NOK");
  const int n_cells_max = get_max_cells_per_rank (NPLASMA);

  const int num_ints = 1 + n_cells_max * (1 + N_BASIC_INTS + nphot_total + nions + NXBANDS + 2 * N_PHOT_PROC);
  const int num_doubles = n_cells_max * (N_BASIC_DOUBLES + 11 * nions + nlte_levels + 2 * nphot_total + n_inner_tot +
                                         11 * NXBANDS + NBINS_IN_CELL_SPEC + 6 * NFLUX_ANGLES +
                                         N_DMO_DT_DIRECTIONS + 12 * NFORCE_DIRECTIONS);


  const int comm_buffer_size = calculate_comm_buffer_size (num_ints, num_doubles);

  char *comm_buffer = malloc (comm_buffer_size);
  if (comm_buffer == NULL)
  {
    Error ("broadcast_plasma_grid: Error in allocating memory for comm_buffer\n");
    Exit (EXIT_FAILURE);
  }

  for (current_rank = 0; current_rank < np_mpi_global; current_rank++)
  {
    position = 0;
    if (rank_global == current_rank)
    {
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
      for (n_plasma = n_start; n_plasma < n_stop; ++n_plasma)
      {
        MPI_Pack (&n_plasma, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        cell = &plasmamain[n_plasma];
        MPI_Pack (&cell->nwind, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->nplasma, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->state.ne, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->state.rho, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->state.vol, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->state.xgamma, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->state.density, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->state.partition, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->state.levden, nlte_levels, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->state.kappa_ff_factor, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->state.recomb_simple, nphot_total, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->state.recomb_simple_upweight, nphot_total, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.kpkt_emiss, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.kpkt_abs, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->state.kbf_use, nphot_total, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->state.kbf_nuse, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->state.t_r, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->state.t_r_old, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->state.t_e, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->state.t_e_old, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.dt_e, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.dt_e_old, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.heat_tot, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.heat_tot_old, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.abs_tot, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.heat_lines, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.heat_ff, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.heat_comp, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.heat_ind_comp, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.heat_lines_macro, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.heat_photo_macro, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.cool_lines_macro, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.cool_bf_macro, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.heat_photo, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.heat_z, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.heat_auger, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.heat_ch_ex, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.abs_photo, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.abs_auger, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->state.w, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.ntot, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.ntot_star, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.ntot_bl, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.ntot_disk, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.ntot_wind, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.ntot_agn, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.nscat_es, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.nscat_res, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.nscat_bf, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.nscat_ff, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.mean_ds, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.n_ds, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.nrad, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.nioniz, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.ioniz, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.recomb, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.inner_ioniz, n_inner_tot, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.inner_recomb, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.scatters, nions, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.xscatters, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.heat_ion, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.heat_inner_ion, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.cool_rr_ion, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.lum_rr_ion, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.cool_dr_ion, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.j, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.ave_freq, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.xj, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.xave_freq, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.fmin, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.fmax, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->state.fmin_mod, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->state.fmax_mod, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.xsd_freq, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.nxtot, NXBANDS, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->state.spec_mod_type, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->state.pl_alpha, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->state.pl_log_w, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->state.exp_temp, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->state.exp_w, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.cell_spec_flux, NBINS_IN_CELL_SPEC, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.F_vis, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.F_UV, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.F_Xray, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.F_vis_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.F_UV_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.F_Xray_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.F_UV_ang_theta, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.F_UV_ang_phi, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.F_UV_ang_r, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.F_UV_ang_theta_persist, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.F_UV_ang_phi_persist, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.F_UV_ang_r_persist, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.j_direct, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.j_scatt, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.ip_direct, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.ip_scatt, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.max_freq, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.cool_tot, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.lum_lines, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.lum_ff, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.cool_adiabatic, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.lum_rr, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.lum_rr_metals, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.cool_comp, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.cool_di, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.cool_dr, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.cool_rr, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.cool_rr_metals, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.lum_tot, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.lum_tot_old, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.cool_tot_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.lum_lines_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.lum_ff_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.cool_adiabatic_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.lum_rr_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.cool_comp_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.cool_di_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.cool_dr_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.cool_rr_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.cool_rr_metals_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.lum_tot_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.heat_shock, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.bf_simple_ionpool_in, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.bf_simple_ionpool_out, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.n_bf_in, N_PHOT_PROC, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.n_bf_out, N_PHOT_PROC, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.comp_nujnu, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.dmo_dt, N_DMO_DT_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.rad_force_es, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.rad_force_ff, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->est.rad_force_bf, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->derived.rad_force_es_persist, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (cell->derived.rad_force_ff_persist, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (cell->derived.rad_force_bf_persist, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.gain, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.converge_t_r, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.converge_t_e, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.converge_hc, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.trcheck, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.techeck, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.hccheck, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.converge_whole, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.converging, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->est.ip, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->derived.xi, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
      }
    }

    MPI_Bcast (comm_buffer, comm_buffer_size, MPI_PACKED, current_rank, MPI_COMM_WORLD);
    position = 0;

    if (rank_global != current_rank)
    {
      MPI_Unpack (comm_buffer, comm_buffer_size, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
      for (i = 0; i < num_comm; i++)
      {
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &n_plasma, 1, MPI_INT, MPI_COMM_WORLD);
        cell = &plasmamain[n_plasma];
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->nwind, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->nplasma, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->state.ne, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->state.rho, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->state.vol, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->state.xgamma, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->state.density, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->state.partition, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->state.levden, nlte_levels, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->state.kappa_ff_factor, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->state.recomb_simple, nphot_total, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->state.recomb_simple_upweight, nphot_total, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.kpkt_emiss, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.kpkt_abs, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->state.kbf_use, nphot_total, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->state.kbf_nuse, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->state.t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->state.t_r_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->state.t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->state.t_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.dt_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.dt_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.heat_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.heat_tot_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.abs_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.heat_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.heat_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.heat_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.heat_ind_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.heat_lines_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.heat_photo_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.cool_lines_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.cool_bf_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.heat_photo, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.heat_z, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.heat_auger, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.heat_ch_ex, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.abs_photo, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.abs_auger, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->state.w, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.ntot, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.ntot_star, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.ntot_bl, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.ntot_disk, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.ntot_wind, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.ntot_agn, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.nscat_es, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.nscat_res, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.nscat_bf, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.nscat_ff, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.mean_ds, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.n_ds, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.nrad, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.nioniz, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.ioniz, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.recomb, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.inner_ioniz, n_inner_tot, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.inner_recomb, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.scatters, nions, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.xscatters, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.heat_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.heat_inner_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.cool_rr_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.lum_rr_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.cool_dr_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.j, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.ave_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.xj, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.xave_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.fmin, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.fmax, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->state.fmin_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->state.fmax_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.xsd_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.nxtot, NXBANDS, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->state.spec_mod_type, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->state.pl_alpha, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->state.pl_log_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->state.exp_temp, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->state.exp_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.cell_spec_flux, NBINS_IN_CELL_SPEC, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.F_vis, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.F_UV, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.F_Xray, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.F_vis_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.F_UV_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.F_Xray_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.F_UV_ang_theta, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.F_UV_ang_phi, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.F_UV_ang_r, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.F_UV_ang_theta_persist, NFLUX_ANGLES, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.F_UV_ang_phi_persist, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.F_UV_ang_r_persist, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.j_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.j_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.ip_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.ip_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.max_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.cool_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.lum_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.lum_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.cool_adiabatic, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.lum_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.lum_rr_metals, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.cool_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.cool_di, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.cool_dr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.cool_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.cool_rr_metals, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.lum_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.lum_tot_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.cool_tot_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.lum_lines_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.lum_ff_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.cool_adiabatic_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.lum_rr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.cool_comp_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.cool_di_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.cool_dr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.cool_rr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.cool_rr_metals_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.lum_tot_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.heat_shock, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.bf_simple_ionpool_in, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.bf_simple_ionpool_out, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.n_bf_in, N_PHOT_PROC, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.n_bf_out, N_PHOT_PROC, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.comp_nujnu, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.dmo_dt, N_DMO_DT_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.rad_force_es, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.rad_force_ff, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->est.rad_force_bf, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.rad_force_es_persist, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.rad_force_ff_persist, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->derived.rad_force_bf_persist, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.gain, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.converge_t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.converge_t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.converge_hc, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.trcheck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.techeck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.hccheck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.converge_whole, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.converging, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->est.ip, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->derived.xi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
      }
    }
  }

  free (comm_buffer);

  /* Barrier to ensure shared memory writes are visible to all node-local ranks */
  MPI_Barrier (node_comm);

  d_xsignal (files.root, "%-20s Finished communicating plasma grid\n", "OK");
#endif
}

/**********************************************************/
/**
 * @brief  Communicate the luminosity properties for plasma cells
 *
 * @param [in] int n_start       The index of the first cell updated by this rank
 * @param [in] int n_stop        The index of the last cell updated by this rank
 * @param [in] int n_cells_rank  The number of cells this rank updated
 *
 * @details
 *
 * The communication pattern and how the size of the communication buffer is
 * determined is documented in
 * $SIROCCO/docs/sphinx/source/developer/mpi_comms.rst.
 * To communicate a new variable, the communication buffer needs to be made
 * bigger and a new `MPI_Pack` and `MPI_Unpack` call need to be added. See the
 * developer documentation for more details.
 *
 * ### Notes ###
 *
 * When this is called in wind update, there is redundant information being
 * communicated in `broadcast_updated_plasma_properties` which communicates the exact (but
 * probably incorrect) data this function does. A refactor to clean this up could
 * be done in the future to avoid the extra communication latency from
 * communicating the data twice.
 *
 **********************************************************/

void
broadcast_wind_luminosity (const int n_start, const int n_stop, const int n_cells_rank)
{
#ifdef MPI_ON
  int n_plasma;
  int current_rank;
  int position;
  int num_comm;

  d_xsignal (files.root, "%-20s Begin communicating wind luminosity\n", "NOK");
  const int n_cells_max = get_max_cells_per_rank (NPLASMA);
  const int comm_buffer_size = calculate_comm_buffer_size (1 + n_cells_max, 4 * n_cells_max);
  char *const comm_buffer = malloc (comm_buffer_size);
  if (comm_buffer == NULL)
  {
    Error ("broadcast_wind_luminosity: Error in allocating memory for comm_buffer\n");
    Exit (EXIT_FAILURE);
  }

  for (current_rank = 0; current_rank < np_mpi_global; ++current_rank)
  {
    position = 0;

    if (rank_global == current_rank)
    {
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
      for (n_plasma = n_start; n_plasma < n_stop; ++n_plasma)
      {
        MPI_Pack (&n_plasma, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.lum_tot, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.lum_ff, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.lum_rr, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.lum_lines, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
      }
    }

    MPI_Bcast (comm_buffer, comm_buffer_size, MPI_PACKED, current_rank, MPI_COMM_WORLD);

    position = 0;

    if (rank_global != current_rank)
    {
      MPI_Unpack (comm_buffer, comm_buffer_size, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
      for (n_plasma = 0; n_plasma < num_comm; ++n_plasma)
      {
        int cell;
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[cell].derived.lum_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[cell].derived.lum_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[cell].derived.lum_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[cell].derived.lum_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);

      }
    }
  }

  /* Barrier to ensure shared memory writes are visible to all node-local ranks */
  MPI_Barrier (node_comm);

  d_xsignal (files.root, "%-20s Finished communicating wind luminosity\n", "OK");
  free (comm_buffer);
#endif
}

/**********************************************************/
/**
 * @brief Communicate wind luminosity information between ranks
 *
 * @param [in] int n_start       The index of the first cell updated by this rank
 * @param [in] int n_stop        The index of the last cell updated by this rank
 * @param [in] int n_cells_rank  The number of cells this rank updated
 *
 * @details
 *
 * The communication pattern and how the size of the communication buffer is
 * determined is documented in
 * $SIROCCO/docs/sphinx/source/developer/mpi_comms.rst.
 * To communicate a new variable, the communication buffer needs to be made
 * bigger and a new `MPI_Pack` and `MPI_Unpack` call need to be added. See the
 * developer documentation for more details.
 *
 * ### Notes ###
 *
 * When this is called in wind update, there is redundant information being
 * communicated in `broadcast_updated_plasma_properties` which communicates the exact (but
 * probably incorrect) data this function does. A refactor to clean this up could
 * be done in the future to avoid the extra communication latency from
 * communicating the data twice.
 *
 **********************************************************/

void
broadcast_wind_cooling (const int n_start, const int n_stop, const int n_cells_rank)
{
#ifdef MPI_ON
  int i;
  int current_rank;
  int position;
  int num_comm;

  d_xsignal (files.root, "%-20s Begin communicating wind cooling\n", "NOK");
  const int n_cells_max = get_max_cells_per_rank (NPLASMA);
  const int comm_buffer_size = calculate_comm_buffer_size (1 + n_cells_max, 9 * n_cells_max);
  char *const comm_buffer = malloc (comm_buffer_size);
  if (comm_buffer == NULL)
  {
    Error ("broadcast_wind_cooling: Error in allocating memory for comm_buffer\n");
    Exit (EXIT_FAILURE);
  }

  for (current_rank = 0; current_rank < np_mpi_global; ++current_rank)
  {
    position = 0;

    if (rank_global == current_rank)
    {
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
      for (i = n_start; i < n_stop; ++i)
      {
        MPI_Pack (&i, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].est.cool_tot, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].derived.lum_ff, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].derived.lum_lines, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].derived.cool_rr, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].derived.cool_comp, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].derived.cool_di, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].derived.cool_dr, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].derived.cool_adiabatic, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].derived.heat_shock, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
      }
    }

    MPI_Bcast (comm_buffer, comm_buffer_size, MPI_PACKED, current_rank, MPI_COMM_WORLD);

    position = 0;

    if (rank_global != current_rank)
    {
      MPI_Unpack (comm_buffer, comm_buffer_size, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
      for (i = 0; i < num_comm; ++i)
      {
        int n;
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &n, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].est.cool_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].derived.lum_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].derived.lum_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].derived.cool_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].derived.cool_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].derived.cool_di, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].derived.cool_dr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].derived.cool_adiabatic, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].derived.heat_shock, 1, MPI_DOUBLE, MPI_COMM_WORLD);
      }
    }
  }

  /* Barrier to ensure shared memory writes are visible to all node-local ranks */
  MPI_Barrier (node_comm);

  d_xsignal (files.root, "%-20s Finished communicating wind cooling\n", "OK");
  free (comm_buffer);
#endif
}

/**********************************************************/
/**
 * @brief Communicate changing properties in the plasma cells between ranks.
 *
 * @param [in]  int n_start      the first cell this rank will communicate
 * @param [in]  int n_stop       the last cell this rank will communicate
 * @param [in]  int n_cells_rank the number of cells this rank will communicate
 *
 * @details
 *
 * The communication pattern and how the size of the communication buffer is
 * determined is documented in
 * $SIROCCO/docs/sphinx/source/developer/mpi_comms.rst.
 * To communicate a new variable, the communication buffer needs to be made
 * bigger and a new `MPI_Pack` and `MPI_Unpack` call need to be added. See the
 * developer documentation for more details.
 *
 **********************************************************/

int
broadcast_updated_plasma_properties (const int n_start_rank, const int n_stop_rank, const int n_cells_rank)
{
#ifdef MPI_ON
  int i;
  int n_plasma;
  int position;
  int n_mpi;
  int num_cells_communicated;

  d_xsignal (files.root, "%-20s Begin communicating updated plasma properties\n", "NOK");
  const int n_cells_max = get_max_cells_per_rank (NPLASMA);
  //OLD  const int num_ints = 1 + n_cells_max * (20 + nphot_total + 2 * NXBANDS + 2 * N_PHOT_PROC + nions);
  const int num_ints = 1 + n_cells_max * (N_BASIC_INTS + nphot_total + 2 * NXBANDS + 2 * N_PHOT_PROC + nions);
  const int num_doubles =
    n_cells_max * (N_BASIC_DOUBLES + 1 * 3 + 9 * 4 + 6 * NFLUX_ANGLES + 3 * NFORCE_DIRECTIONS + 9 * nions + 1 * nlte_levels +
                   3 * nphot_total + 1 * n_inner_tot + 9 * NXBANDS + 1 * NBINS_IN_CELL_SPEC);

  const int size_of_comm_buffer = calculate_comm_buffer_size (num_ints, num_doubles);
  char *const comm_buffer = malloc (size_of_comm_buffer);
  if (comm_buffer == NULL)
  {
    Error ("broadcast_updated_plasma_properties: Error in allocating memory for comm_buffer\n");
    Exit (EXIT_FAILURE);
  }

  for (n_mpi = 0; n_mpi < np_mpi_global; n_mpi++)
  {
    position = 0;

    if (rank_global == n_mpi)
    {
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
      for (n_plasma = n_start_rank; n_plasma < n_stop_rank; n_plasma++)
      {
        // cell number
        MPI_Pack (&n_plasma, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].nwind, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].nplasma, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].state.ne, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].state.rho, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].state.vol, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].state.xgamma, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].state.density, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].state.partition, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].state.levden, nlte_levels, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].state.kappa_ff_factor, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].state.recomb_simple, nphot_total, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].state.recomb_simple_upweight, nphot_total, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.kpkt_emiss, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.kpkt_abs, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].state.kbf_use, nphot_total, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].state.kbf_nuse, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].state.t_r, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].state.t_r_old, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].state.t_e, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].state.t_e_old, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.dt_e, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.dt_e_old, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.heat_tot, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.heat_tot_old, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.abs_tot, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.heat_lines, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.heat_ff, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.heat_comp, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.heat_ind_comp, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.heat_lines_macro, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.heat_photo_macro, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.cool_lines_macro, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.cool_bf_macro, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.heat_photo, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.heat_z, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.heat_auger, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.heat_ch_ex, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.abs_photo, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.abs_auger, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].state.w, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.ntot, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.ntot_star, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.ntot_bl, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.ntot_disk, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.ntot_wind, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.ntot_agn, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.nscat_es, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.nscat_bf, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.nscat_ff, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.mean_ds, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.n_ds, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.nrad, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.nioniz, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.ioniz, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].derived.recomb, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.inner_ioniz, n_inner_tot, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].derived.scatters, nions, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].derived.xscatters, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.heat_ion, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.heat_inner_ion, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].derived.cool_rr_ion, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].derived.lum_rr_ion, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.j, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.ave_freq, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.xj, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.xave_freq, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].state.fmin_mod, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].state.fmax_mod, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.xsd_freq, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.nxtot, NXBANDS, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].state.spec_mod_type, NXBANDS, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].state.pl_alpha, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].state.pl_log_w, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].state.exp_temp, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].state.exp_w, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.cell_spec_flux, NBINS_IN_CELL_SPEC, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.F_vis, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.F_UV, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.F_Xray, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].derived.F_vis_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].derived.F_UV_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].derived.F_Xray_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer,
                  &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.F_UV_ang_theta, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.F_UV_ang_phi, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.F_UV_ang_r, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].derived.F_UV_ang_theta_persist, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, size_of_comm_buffer,
                  &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].derived.F_UV_ang_phi_persist, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].derived.F_UV_ang_r_persist, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.j_direct, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.j_scatt, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.ip_direct, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.ip_scatt, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.max_freq, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.cool_tot, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.lum_lines, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.lum_ff, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.cool_adiabatic, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.lum_rr, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.lum_rr_metals, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.cool_comp, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.cool_di, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.cool_dr, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.cool_rr, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.cool_rr_metals, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.lum_tot, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.lum_tot_old, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.cool_tot_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.lum_lines_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.lum_ff_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.cool_adiabatic_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.lum_rr_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.cool_comp_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.cool_di_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.cool_dr_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.cool_rr_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.cool_rr_metals_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.lum_tot_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.heat_shock, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.bf_simple_ionpool_in, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.bf_simple_ionpool_out, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.n_bf_in, N_PHOT_PROC, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.n_bf_out, N_PHOT_PROC, MPI_INT, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.comp_nujnu, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].derived.dmo_dt, N_DMO_DT_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.rad_force_es, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.rad_force_ff, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].est.rad_force_bf, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].derived.rad_force_es_persist, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer,
                  &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].derived.rad_force_ff_persist, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer,
                  &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].derived.rad_force_bf_persist, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer,
                  &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.gain, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.converge_t_r, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.converge_t_e, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.converge_hc, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.trcheck, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.techeck, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.hccheck, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.converge_whole, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.converging, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].est.ip, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].derived.xi, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);

      }
    }

    MPI_Bcast (comm_buffer, size_of_comm_buffer, MPI_PACKED, n_mpi, MPI_COMM_WORLD);

    position = 0;

    if (rank_global != n_mpi)
    {
      MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &num_cells_communicated, 1, MPI_INT, MPI_COMM_WORLD);
      for (i = 0; i < num_cells_communicated; ++i)
      {
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &n_plasma, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].nwind, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].nplasma, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].state.ne, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].state.rho, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].state.vol, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].state.xgamma, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].state.density, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].state.partition, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].state.levden, nlte_levels, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].state.kappa_ff_factor, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].state.recomb_simple, nphot_total, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].state.recomb_simple_upweight, nphot_total, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.kpkt_emiss, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.kpkt_abs, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].state.kbf_use, nphot_total, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].state.kbf_nuse, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].state.t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].state.t_r_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].state.t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].state.t_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.dt_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.dt_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.heat_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.heat_tot_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.abs_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.heat_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.heat_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.heat_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.heat_ind_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.heat_lines_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.heat_photo_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.cool_lines_macro, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.cool_bf_macro, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.heat_photo, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.heat_z, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.heat_auger, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.heat_ch_ex, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.abs_photo, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.abs_auger, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].state.w, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.ntot, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.ntot_star, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.ntot_bl, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.ntot_disk, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.ntot_wind, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.ntot_agn, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.nscat_es, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.nscat_bf, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.nscat_ff, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.mean_ds, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.n_ds, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.nrad, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.nioniz, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.ioniz, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].derived.recomb, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.inner_ioniz, n_inner_tot, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].derived.scatters, nions, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].derived.xscatters, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.heat_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.heat_inner_ion, nions, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].derived.cool_rr_ion, nions, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].derived.lum_rr_ion, nions, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.j, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.ave_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.xj, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.xave_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].state.fmin_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].state.fmax_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.xsd_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.nxtot, NXBANDS, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].state.spec_mod_type, NXBANDS, MPI_INT,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].state.pl_alpha, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].state.pl_log_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].state.exp_temp, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].state.exp_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.cell_spec_flux, NBINS_IN_CELL_SPEC, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.F_vis, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.F_UV, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.F_Xray, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].derived.F_vis_persistent, NFORCE_DIRECTIONS,
                    MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].derived.F_UV_persistent, NFORCE_DIRECTIONS,
                    MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].derived.F_Xray_persistent, NFORCE_DIRECTIONS,
                    MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.F_UV_ang_theta, NFLUX_ANGLES, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.F_UV_ang_phi, NFLUX_ANGLES, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.F_UV_ang_r, NFLUX_ANGLES, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].derived.F_UV_ang_theta_persist, NFLUX_ANGLES,
                    MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].derived.F_UV_ang_phi_persist, NFLUX_ANGLES,
                    MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].derived.F_UV_ang_r_persist, NFLUX_ANGLES, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.j_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.j_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.ip_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.ip_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.max_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.cool_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.lum_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.lum_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.cool_adiabatic, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.lum_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.lum_rr_metals, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.cool_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.cool_di, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.cool_dr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.cool_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.cool_rr_metals, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.lum_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.lum_tot_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.cool_tot_ioniz, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.lum_lines_ioniz, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.lum_ff_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.cool_adiabatic_ioniz, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.lum_rr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.cool_comp_ioniz, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.cool_di_ioniz, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.cool_dr_ioniz, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.cool_rr_ioniz, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.cool_rr_metals_ioniz, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.lum_tot_ioniz, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.heat_shock, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.bf_simple_ionpool_in, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.bf_simple_ionpool_out, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.n_bf_in, N_PHOT_PROC, MPI_INT,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.n_bf_out, N_PHOT_PROC, MPI_INT,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.comp_nujnu, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].derived.dmo_dt, N_DMO_DT_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.rad_force_es, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.rad_force_ff, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].est.rad_force_bf, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].derived.rad_force_es_persist, NFORCE_DIRECTIONS,
                    MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].derived.rad_force_ff_persist, NFORCE_DIRECTIONS,
                    MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].derived.rad_force_bf_persist, NFORCE_DIRECTIONS,
                    MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.gain, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.converge_t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.converge_t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.converge_hc, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.trcheck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.techeck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.hccheck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.converge_whole, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.converging, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].est.ip, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].derived.xi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
      }
    }
  }

  free (comm_buffer);

  /* Barrier to ensure shared memory writes are visible to all node-local ranks */
  MPI_Barrier (node_comm);

  d_xsignal (files.root, "%-20s Finished communicating updated plasma properties\n", "OK");
#endif
  return EXIT_SUCCESS;
}

/**********************************************************/
/**
 * @brief      communicates the MC estimators between tasks
 *
 * @details
 * communicates the MC estimators between tasks relating to
 * spectral models, heating and cooling and cell diagnostics like IP.
 * In the case of some variables, the quantities are maxima and minima so the
 * flag MPI_MAX or MPI_MIN is used in MPI_Reduce. For summed
 * quantities like heating we use MPI_SUM.
 *
 * This routine should only do anything if the MPI_ON flag was present
 * in compilation. It communicates all the information
 * required for the spectral model ionization scheme, and
 * also heating and cooling quantities in cells.
 **********************************************************/

int
reduce_simple_estimators (void)
{
#ifdef MPI_ON                   // these routines should only be called anyway in parallel but we need these to compile

  int mpi_i, mpi_j;
  double *maxfreqhelper, *maxfreqhelper2;
  /*NSH 131213 the next line introduces new helper arrays for the max and min frequencies in bands */
  double *maxbandfreqhelper, *maxbandfreqhelper2, *minbandfreqhelper, *minbandfreqhelper2;
  double *redhelper, *redhelper2, *qdisk_helper, *qdisk_helper2;
  double *ion_helper, *ion_helper2;
  double *inner_ion_helper, *inner_ion_helper2;
  double *flux_helper, *flux_helper2;

//OLD  int nspec, size_of_commbuffer;
  int size_of_commbuffer;

  int *iredhelper, *iredhelper2, *iqdisk_helper, *iqdisk_helper2;
  // int size_of_helpers;
  int plasma_double_helpers, plasma_int_helpers;

  d_xsignal (files.root, "%-20s Begin reduction of simple estimators\n", "NOK");

  /* The size of the helper array for doubles. We transmit 10 numbers
     for each cell, plus three arrays, each of length NXBANDS */

  plasma_double_helpers = (39 + 3 * NXBANDS) * NPLASMA;

  /* The size of the helper array for integers. We transmit 7 numbers
     for each cell, plus one array of length NXBANDS */
  plasma_int_helpers = (7 + NXBANDS) * NPLASMA;


  maxfreqhelper = calloc (sizeof (double), NPLASMA);
  maxfreqhelper2 = calloc (sizeof (double), NPLASMA);
  /* NSH 131213 - allocate memory for the band limited max and min frequencies */
  maxbandfreqhelper = calloc (sizeof (double), NPLASMA * NXBANDS);
  maxbandfreqhelper2 = calloc (sizeof (double), NPLASMA * NXBANDS);
  minbandfreqhelper = calloc (sizeof (double), NPLASMA * NXBANDS);
  minbandfreqhelper2 = calloc (sizeof (double), NPLASMA * NXBANDS);
  redhelper = calloc (sizeof (double), plasma_double_helpers);
  redhelper2 = calloc (sizeof (double), plasma_double_helpers);

  ion_helper = calloc (sizeof (double), NPLASMA * nions);
  ion_helper2 = calloc (sizeof (double), NPLASMA * nions);
  inner_ion_helper = calloc (sizeof (double), NPLASMA * n_inner_tot);
  inner_ion_helper2 = calloc (sizeof (double), NPLASMA * n_inner_tot);
  /* Routine to average the qdisk quantities. The 3 is because
     we have three doubles to worry about (emit, heat and ave_freq) and
     two integers (nhit and nphot) */
  qdisk_helper = calloc (sizeof (double), NRINGS * 3);
  qdisk_helper2 = calloc (sizeof (double), NRINGS * 3);

  flux_helper = calloc (sizeof (double), NPLASMA * NFLUX_ANGLES * 3);
  flux_helper2 = calloc (sizeof (double), NPLASMA * NFLUX_ANGLES * 3);

  // the following blocks gather all the estimators to the zeroth (Master) thread
  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    maxfreqhelper[mpi_i] = plasmamain[mpi_i].est.max_freq;
    redhelper[mpi_i] = plasmamain[mpi_i].est.j / np_mpi_global;
    redhelper[mpi_i + NPLASMA] = plasmamain[mpi_i].est.ave_freq / np_mpi_global;
    redhelper[mpi_i + 2 * NPLASMA] = plasmamain[mpi_i].est.cool_tot / np_mpi_global;
    redhelper[mpi_i + 3 * NPLASMA] = plasmamain[mpi_i].est.heat_tot / np_mpi_global;
    redhelper[mpi_i + 4 * NPLASMA] = plasmamain[mpi_i].est.heat_lines / np_mpi_global;
    redhelper[mpi_i + 5 * NPLASMA] = plasmamain[mpi_i].est.heat_ff / np_mpi_global;
    redhelper[mpi_i + 6 * NPLASMA] = plasmamain[mpi_i].est.heat_comp / np_mpi_global;
    redhelper[mpi_i + 7 * NPLASMA] = plasmamain[mpi_i].est.heat_ind_comp / np_mpi_global;
    redhelper[mpi_i + 8 * NPLASMA] = plasmamain[mpi_i].est.heat_photo / np_mpi_global;
    redhelper[mpi_i + 9 * NPLASMA] = plasmamain[mpi_i].est.ip / np_mpi_global;
    redhelper[mpi_i + 10 * NPLASMA] = plasmamain[mpi_i].est.j_direct / np_mpi_global;
    redhelper[mpi_i + 11 * NPLASMA] = plasmamain[mpi_i].est.j_scatt / np_mpi_global;
    redhelper[mpi_i + 12 * NPLASMA] = plasmamain[mpi_i].est.ip_direct / np_mpi_global;
    redhelper[mpi_i + 13 * NPLASMA] = plasmamain[mpi_i].est.ip_scatt / np_mpi_global;
    redhelper[mpi_i + 14 * NPLASMA] = plasmamain[mpi_i].est.heat_auger / np_mpi_global;
    redhelper[mpi_i + 15 * NPLASMA] = plasmamain[mpi_i].est.rad_force_es[0] / np_mpi_global;
    redhelper[mpi_i + 16 * NPLASMA] = plasmamain[mpi_i].est.rad_force_es[1] / np_mpi_global;
    redhelper[mpi_i + 17 * NPLASMA] = plasmamain[mpi_i].est.rad_force_es[2] / np_mpi_global;
    redhelper[mpi_i + 18 * NPLASMA] = plasmamain[mpi_i].est.rad_force_es[3] / np_mpi_global;
    redhelper[mpi_i + 19 * NPLASMA] = plasmamain[mpi_i].est.F_vis[0] / np_mpi_global;
    redhelper[mpi_i + 20 * NPLASMA] = plasmamain[mpi_i].est.F_vis[1] / np_mpi_global;
    redhelper[mpi_i + 21 * NPLASMA] = plasmamain[mpi_i].est.F_vis[2] / np_mpi_global;
    redhelper[mpi_i + 22 * NPLASMA] = plasmamain[mpi_i].est.F_vis[3] / np_mpi_global;
    redhelper[mpi_i + 23 * NPLASMA] = plasmamain[mpi_i].est.F_UV[0] / np_mpi_global;
    redhelper[mpi_i + 24 * NPLASMA] = plasmamain[mpi_i].est.F_UV[1] / np_mpi_global;
    redhelper[mpi_i + 25 * NPLASMA] = plasmamain[mpi_i].est.F_UV[2] / np_mpi_global;
    redhelper[mpi_i + 26 * NPLASMA] = plasmamain[mpi_i].est.F_UV[3] / np_mpi_global;
    redhelper[mpi_i + 27 * NPLASMA] = plasmamain[mpi_i].est.F_Xray[0] / np_mpi_global;
    redhelper[mpi_i + 28 * NPLASMA] = plasmamain[mpi_i].est.F_Xray[1] / np_mpi_global;
    redhelper[mpi_i + 29 * NPLASMA] = plasmamain[mpi_i].est.F_Xray[2] / np_mpi_global;
    redhelper[mpi_i + 30 * NPLASMA] = plasmamain[mpi_i].est.F_Xray[3] / np_mpi_global;
    redhelper[mpi_i + 31 * NPLASMA] = plasmamain[mpi_i].est.rad_force_bf[0] / np_mpi_global;
    redhelper[mpi_i + 32 * NPLASMA] = plasmamain[mpi_i].est.rad_force_bf[1] / np_mpi_global;
    redhelper[mpi_i + 33 * NPLASMA] = plasmamain[mpi_i].est.rad_force_bf[2] / np_mpi_global;
    redhelper[mpi_i + 34 * NPLASMA] = plasmamain[mpi_i].est.rad_force_bf[3] / np_mpi_global;
    redhelper[mpi_i + 35 * NPLASMA] = plasmamain[mpi_i].est.rad_force_ff[0] / np_mpi_global;
    redhelper[mpi_i + 36 * NPLASMA] = plasmamain[mpi_i].est.rad_force_ff[1] / np_mpi_global;
    redhelper[mpi_i + 37 * NPLASMA] = plasmamain[mpi_i].est.rad_force_ff[2] / np_mpi_global;
    redhelper[mpi_i + 38 * NPLASMA] = plasmamain[mpi_i].est.rad_force_ff[3] / np_mpi_global;

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      redhelper[mpi_i + (39 + mpi_j) * NPLASMA] = plasmamain[mpi_i].est.xj[mpi_j] / np_mpi_global;
      redhelper[mpi_i + (39 + NXBANDS + mpi_j) * NPLASMA] = plasmamain[mpi_i].est.xave_freq[mpi_j] / np_mpi_global;
      redhelper[mpi_i + (39 + 2 * NXBANDS + mpi_j) * NPLASMA] = plasmamain[mpi_i].est.xsd_freq[mpi_j] / np_mpi_global;
      /* 131213 NSH populate the band limited min and max frequency arrays */
      maxbandfreqhelper[mpi_i * NXBANDS + mpi_j] = plasmamain[mpi_i].est.fmax[mpi_j];
      minbandfreqhelper[mpi_i * NXBANDS + mpi_j] = plasmamain[mpi_i].est.fmin[mpi_j];
    }
    for (mpi_j = 0; mpi_j < nions; mpi_j++)
    {
      ion_helper[mpi_i * nions + mpi_j] = plasmamain[mpi_i].est.ioniz[mpi_j] / np_mpi_global;
    }
    for (mpi_j = 0; mpi_j < n_inner_tot; mpi_j++)
    {
      inner_ion_helper[mpi_i * n_inner_tot + mpi_j] = plasmamain[mpi_i].est.inner_ioniz[mpi_j] / np_mpi_global;
    }
    for (mpi_j = 0; mpi_j < NFLUX_ANGLES; mpi_j++)
    {
      flux_helper[mpi_i * (3 * NFLUX_ANGLES) + mpi_j] = plasmamain[mpi_i].est.F_UV_ang_theta[mpi_j] / np_mpi_global;
      flux_helper[mpi_i * (3 * NFLUX_ANGLES) + NFLUX_ANGLES + mpi_j] = plasmamain[mpi_i].est.F_UV_ang_phi[mpi_j] / np_mpi_global;
      flux_helper[mpi_i * (3 * NFLUX_ANGLES) + 2 * NFLUX_ANGLES + mpi_j] = plasmamain[mpi_i].est.F_UV_ang_r[mpi_j] / np_mpi_global;
    }
  }

  for (mpi_i = 0; mpi_i < NRINGS; mpi_i++)
  {
    qdisk_helper[mpi_i] = qdisk.heat[mpi_i] / np_mpi_global;
    qdisk_helper[mpi_i + NRINGS] = qdisk.ave_freq[mpi_i] / np_mpi_global;
    qdisk_helper[mpi_i + 2 * NRINGS] = qdisk.emit[mpi_i] / np_mpi_global;
  }


  /* Reduced and communicate all ranks. Some operations are MIN/MAX but most are
   * sums to compute the average across ranks */
  MPI_Allreduce (minbandfreqhelper, minbandfreqhelper2, NPLASMA * NXBANDS, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce (maxbandfreqhelper, maxbandfreqhelper2, NPLASMA * NXBANDS, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (maxfreqhelper, maxfreqhelper2, NPLASMA, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (redhelper, redhelper2, plasma_double_helpers, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (flux_helper, flux_helper2, NPLASMA * 3 * NFLUX_ANGLES, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (ion_helper, ion_helper2, NPLASMA * nions, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (inner_ion_helper, inner_ion_helper2, NPLASMA * n_inner_tot, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (qdisk_helper, qdisk_helper2, 3 * NRINGS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /* Unpacking stuff */
  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    plasmamain[mpi_i].est.max_freq = maxfreqhelper2[mpi_i];
    plasmamain[mpi_i].est.j = redhelper2[mpi_i];
    plasmamain[mpi_i].est.ave_freq = redhelper2[mpi_i + NPLASMA];
    plasmamain[mpi_i].est.cool_tot = redhelper2[mpi_i + 2 * NPLASMA];
    plasmamain[mpi_i].est.heat_tot = redhelper2[mpi_i + 3 * NPLASMA];
    plasmamain[mpi_i].est.heat_lines = redhelper2[mpi_i + 4 * NPLASMA];
    plasmamain[mpi_i].est.heat_ff = redhelper2[mpi_i + 5 * NPLASMA];
    plasmamain[mpi_i].est.heat_comp = redhelper2[mpi_i + 6 * NPLASMA];
    plasmamain[mpi_i].est.heat_ind_comp = redhelper2[mpi_i + 7 * NPLASMA];
    plasmamain[mpi_i].est.heat_photo = redhelper2[mpi_i + 8 * NPLASMA];
    plasmamain[mpi_i].est.ip = redhelper2[mpi_i + 9 * NPLASMA];
    plasmamain[mpi_i].est.j_direct = redhelper2[mpi_i + 10 * NPLASMA];
    plasmamain[mpi_i].est.j_scatt = redhelper2[mpi_i + 11 * NPLASMA];
    plasmamain[mpi_i].est.ip_direct = redhelper2[mpi_i + 12 * NPLASMA];
    plasmamain[mpi_i].est.ip_scatt = redhelper2[mpi_i + 13 * NPLASMA];
    plasmamain[mpi_i].est.heat_auger = redhelper2[mpi_i + 14 * NPLASMA];
    plasmamain[mpi_i].est.rad_force_es[0] = redhelper2[mpi_i + 15 * NPLASMA];
    plasmamain[mpi_i].est.rad_force_es[1] = redhelper2[mpi_i + 16 * NPLASMA];
    plasmamain[mpi_i].est.rad_force_es[2] = redhelper2[mpi_i + 17 * NPLASMA];
    plasmamain[mpi_i].est.rad_force_es[3] = redhelper2[mpi_i + 18 * NPLASMA];
    plasmamain[mpi_i].est.F_vis[0] = redhelper2[mpi_i + 19 * NPLASMA];
    plasmamain[mpi_i].est.F_vis[1] = redhelper2[mpi_i + 20 * NPLASMA];
    plasmamain[mpi_i].est.F_vis[2] = redhelper2[mpi_i + 21 * NPLASMA];
    plasmamain[mpi_i].est.F_vis[3] = redhelper2[mpi_i + 22 * NPLASMA];
    plasmamain[mpi_i].est.F_UV[0] = redhelper2[mpi_i + 23 * NPLASMA];
    plasmamain[mpi_i].est.F_UV[1] = redhelper2[mpi_i + 24 * NPLASMA];
    plasmamain[mpi_i].est.F_UV[2] = redhelper2[mpi_i + 25 * NPLASMA];
    plasmamain[mpi_i].est.F_UV[3] = redhelper2[mpi_i + 26 * NPLASMA];
    plasmamain[mpi_i].est.F_Xray[0] = redhelper2[mpi_i + 27 * NPLASMA];
    plasmamain[mpi_i].est.F_Xray[1] = redhelper2[mpi_i + 28 * NPLASMA];
    plasmamain[mpi_i].est.F_Xray[2] = redhelper2[mpi_i + 29 * NPLASMA];
    plasmamain[mpi_i].est.F_Xray[3] = redhelper2[mpi_i + 30 * NPLASMA];
    plasmamain[mpi_i].est.rad_force_bf[0] = redhelper2[mpi_i + 31 * NPLASMA];
    plasmamain[mpi_i].est.rad_force_bf[1] = redhelper2[mpi_i + 32 * NPLASMA];
    plasmamain[mpi_i].est.rad_force_bf[2] = redhelper2[mpi_i + 33 * NPLASMA];
    plasmamain[mpi_i].est.rad_force_bf[3] = redhelper2[mpi_i + 34 * NPLASMA];
    plasmamain[mpi_i].est.rad_force_ff[0] = redhelper2[mpi_i + 35 * NPLASMA];
    plasmamain[mpi_i].est.rad_force_ff[1] = redhelper2[mpi_i + 36 * NPLASMA];
    plasmamain[mpi_i].est.rad_force_ff[2] = redhelper2[mpi_i + 37 * NPLASMA];
    plasmamain[mpi_i].est.rad_force_ff[3] = redhelper2[mpi_i + 38 * NPLASMA];

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      plasmamain[mpi_i].est.xj[mpi_j] = redhelper2[mpi_i + (39 + mpi_j) * NPLASMA];
      plasmamain[mpi_i].est.xave_freq[mpi_j] = redhelper2[mpi_i + (39 + NXBANDS + mpi_j) * NPLASMA];
      plasmamain[mpi_i].est.xsd_freq[mpi_j] = redhelper2[mpi_i + (39 + NXBANDS * 2 + mpi_j) * NPLASMA];

      /* 131213 NSH And unpack the min and max banded frequencies to the plasma array */
      plasmamain[mpi_i].est.fmax[mpi_j] = maxbandfreqhelper2[mpi_i * NXBANDS + mpi_j];
      plasmamain[mpi_i].est.fmin[mpi_j] = minbandfreqhelper2[mpi_i * NXBANDS + mpi_j];
    }
    for (mpi_j = 0; mpi_j < nions; mpi_j++)
    {
      plasmamain[mpi_i].est.ioniz[mpi_j] = ion_helper2[mpi_i * nions + mpi_j];
    }
    for (mpi_j = 0; mpi_j < n_inner_tot; mpi_j++)
    {
      plasmamain[mpi_i].est.inner_ioniz[mpi_j] = inner_ion_helper2[mpi_i * n_inner_tot + mpi_j];
    }
    for (mpi_j = 0; mpi_j < NFLUX_ANGLES; mpi_j++)
    {
      plasmamain[mpi_i].est.F_UV_ang_theta[mpi_j] = flux_helper2[mpi_i * (3 * NFLUX_ANGLES) + mpi_j];
      plasmamain[mpi_i].est.F_UV_ang_phi[mpi_j] = flux_helper2[mpi_i * (3 * NFLUX_ANGLES) + NFLUX_ANGLES + mpi_j];
      plasmamain[mpi_i].est.F_UV_ang_r[mpi_j] = flux_helper2[mpi_i * (3 * NFLUX_ANGLES) + 2 * NFLUX_ANGLES + mpi_j];
    }
  }


  for (mpi_i = 0; mpi_i < NRINGS; mpi_i++)
  {
    qdisk.heat[mpi_i] = qdisk_helper2[mpi_i];
    qdisk.ave_freq[mpi_i] = qdisk_helper2[mpi_i + NRINGS];
    qdisk.emit[mpi_i] = qdisk_helper2[mpi_i + 2 * NRINGS];
  }

  /* now we've done all the doubles so we can free their helper arrays */
  free (qdisk_helper);
  free (qdisk_helper2);
  free (redhelper);
  free (redhelper2);
  free (maxfreqhelper);
  free (maxfreqhelper2);
  free (maxbandfreqhelper);
  free (maxbandfreqhelper2);
  free (minbandfreqhelper);
  free (minbandfreqhelper2);

  free (ion_helper);
  free (ion_helper2);

  free (inner_ion_helper);
  free (inner_ion_helper2);

  free (flux_helper);
  free (flux_helper2);

  /* allocate the integer helper arrays, set a barrier, then do all the integers. */
  iqdisk_helper = calloc (sizeof (int), NRINGS * 2);
  iqdisk_helper2 = calloc (sizeof (int), NRINGS * 2);
  iredhelper = calloc (sizeof (int), plasma_int_helpers);
  iredhelper2 = calloc (sizeof (int), plasma_int_helpers);

  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    iredhelper[mpi_i] = plasmamain[mpi_i].est.ntot;
    iredhelper[mpi_i + NPLASMA] = plasmamain[mpi_i].est.ntot_star;
    iredhelper[mpi_i + 2 * NPLASMA] = plasmamain[mpi_i].est.ntot_bl;
    iredhelper[mpi_i + 3 * NPLASMA] = plasmamain[mpi_i].est.ntot_disk;
    iredhelper[mpi_i + 4 * NPLASMA] = plasmamain[mpi_i].est.ntot_wind;
    iredhelper[mpi_i + 5 * NPLASMA] = plasmamain[mpi_i].est.ntot_agn;
    iredhelper[mpi_i + 6 * NPLASMA] = plasmamain[mpi_i].est.nioniz;

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      iredhelper[mpi_i + (7 + mpi_j) * NPLASMA] = plasmamain[mpi_i].est.nxtot[mpi_j];
    }
  }

  for (mpi_i = 0; mpi_i < NRINGS; mpi_i++)
  {
    iqdisk_helper[mpi_i] = qdisk.nphot[mpi_i];
    iqdisk_helper[mpi_i + NRINGS] = qdisk.nhit[mpi_i];
  }

  MPI_Allreduce (iredhelper, iredhelper2, plasma_int_helpers, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (iqdisk_helper, iqdisk_helper2, 2 * NRINGS, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    plasmamain[mpi_i].est.ntot = iredhelper2[mpi_i];
    plasmamain[mpi_i].est.ntot_star = iredhelper2[mpi_i + NPLASMA];
    plasmamain[mpi_i].est.ntot_bl = iredhelper2[mpi_i + 2 * NPLASMA];
    plasmamain[mpi_i].est.ntot_disk = iredhelper2[mpi_i + 3 * NPLASMA];
    plasmamain[mpi_i].est.ntot_wind = iredhelper2[mpi_i + 4 * NPLASMA];
    plasmamain[mpi_i].est.ntot_agn = iredhelper2[mpi_i + 5 * NPLASMA];
    plasmamain[mpi_i].est.nioniz = iredhelper2[mpi_i + 6 * NPLASMA];

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      plasmamain[mpi_i].est.nxtot[mpi_j] = iredhelper2[mpi_i + (7 + mpi_j) * NPLASMA];
    }
  }

  for (mpi_i = 0; mpi_i < NRINGS; mpi_i++)
  {
    qdisk.nphot[mpi_i] = iqdisk_helper2[mpi_i];
    qdisk.nhit[mpi_i] = iqdisk_helper2[mpi_i + NRINGS];
  }

  free (iredhelper);
  free (iredhelper2);
  free (iqdisk_helper);
  free (iqdisk_helper2);


/* Now during ionization cycles, process the cell spectra */

  /* The size of the commbuffers need to be the number of spectra x the length of each */

  if (geo.ioniz_or_extract == CYCLE_IONIZ)
  {
    size_of_commbuffer = NPLASMA * NBINS_IN_CELL_SPEC;

    redhelper = calloc (sizeof (double), size_of_commbuffer);
    redhelper2 = calloc (sizeof (double), size_of_commbuffer);

    for (mpi_i = 0; mpi_i < NBINS_IN_CELL_SPEC; mpi_i++)
    {
      for (mpi_j = 0; mpi_j < NPLASMA; mpi_j++)
      {
        redhelper[mpi_i * NPLASMA + mpi_j] = plasmamain[mpi_j].est.cell_spec_flux[mpi_i] / np_mpi_global;

      }
    }

    MPI_Allreduce (redhelper, redhelper2, size_of_commbuffer, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (mpi_i = 0; mpi_i < NBINS_IN_CELL_SPEC; mpi_i++)
    {
      for (mpi_j = 0; mpi_j < NPLASMA; mpi_j++)
      {
        plasmamain[mpi_j].est.cell_spec_flux[mpi_i] = redhelper2[mpi_i * NPLASMA + mpi_j];

      }
    }

    free (redhelper);
    free (redhelper2);
  }

  d_xsignal (files.root, "%-20s Finished reduction of simple estimators\n", "OK");

#endif
  return (0);
}
