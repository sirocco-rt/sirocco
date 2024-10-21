/***********************************************************/
/** @file  communicate_plasma.c
 * @author EJP
 * @date   December 2023
 *
 * @brief Functions for communicating plasma properties
 *
 * @TODO: as much as this as possible should use non-blocking communication
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
 * @brief
 *
 * @param [in] int n_start       The index of the first cell updated by this rank
 * @param [in] int n_stop        The index of the last cell updated by this rank
 * @param [in] int n_cells_rank  The number of cells this rank updated
 *
 * @details
 *
 * The communication pattern is as outlined in dated_plasma_properties.
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
  const int comm_buffer_size = calculate_comm_buffer_size (1 + n_cells_max * (1 + 20 + nphot_total + nions + NXBANDS + 2 * N_PHOT_PROC),
                                                           n_cells_max * (71 + 11 * nions + nlte_levels + 2 * nphot_total + n_inner_tot +
                                                                          11 * NXBANDS + NBINS_IN_CELL_SPEC + 6 * NFLUX_ANGLES +
                                                                          N_DMO_DT_DIRECTIONS + 12 * NFORCE_DIRECTIONS));
  char *comm_buffer = malloc (comm_buffer_size);
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
        MPI_Pack (&cell->ne, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->rho, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->vol, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->xgamma, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->density, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->partition, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->levden, nlte_levels, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->kappa_ff_factor, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->recomb_simple, nphot_total, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->recomb_simple_upweight, nphot_total, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->kpkt_emiss, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->kpkt_abs, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->kbf_use, nphot_total, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->kbf_nuse, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->t_r, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->t_r_old, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->t_e, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->t_e_old, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->dt_e, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->dt_e_old, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->heat_tot, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->heat_tot_old, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->abs_tot, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->heat_lines, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->heat_ff, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->heat_comp, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->heat_ind_comp, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->heat_lines_macro, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->heat_photo_macro, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->heat_photo, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->heat_z, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->heat_auger, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->heat_ch_ex, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->abs_photo, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->abs_auger, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->w, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->ntot, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->ntot_star, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->ntot_bl, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->ntot_disk, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->ntot_wind, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->ntot_agn, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->nscat_es, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->nscat_res, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->mean_ds, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->n_ds, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->nrad, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->nioniz, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->ioniz, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->recomb, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->inner_ioniz, n_inner_tot, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->inner_recomb, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->scatters, nions, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->xscatters, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->heat_ion, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->heat_inner_ion, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->cool_rr_ion, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->lum_rr_ion, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->cool_dr_ion, nions, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->j, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->ave_freq, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->xj, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->xave_freq, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->fmin, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->fmax, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->fmin_mod, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->fmax_mod, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->xsd_freq, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->nxtot, NXBANDS, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->spec_mod_type, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->pl_alpha, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->pl_log_w, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->exp_temp, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->exp_w, NXBANDS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->cell_spec_flux, NBINS_IN_CELL_SPEC, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->F_vis, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->F_UV, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->F_Xray, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->F_vis_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->F_UV_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->F_Xray_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->F_UV_ang_theta, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->F_UV_ang_phi, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->F_UV_ang_r, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->F_UV_ang_theta_persist, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->F_UV_ang_phi_persist, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->F_UV_ang_r_persist, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->j_direct, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->j_scatt, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->ip_direct, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->ip_scatt, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->max_freq, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->cool_tot, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->lum_lines, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->lum_ff, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->cool_adiabatic, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->lum_rr, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->lum_rr_metals, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->cool_comp, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->cool_di, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->cool_dr, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->cool_rr, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->cool_rr_metals, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->lum_tot, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->lum_tot_old, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->cool_tot_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->lum_lines_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->lum_ff_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->cool_adiabatic_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->lum_rr_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->cool_comp_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->cool_di_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->cool_dr_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->cool_rr_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->cool_rr_metals_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->lum_tot_ioniz, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->heat_shock, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->bf_simple_ionpool_in, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->bf_simple_ionpool_out, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->n_bf_in, N_PHOT_PROC, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->n_bf_out, N_PHOT_PROC, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->comp_nujnu, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->dmo_dt, N_DMO_DT_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->rad_force_es, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->rad_force_ff, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->rad_force_bf, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->rad_force_es_persist, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->rad_force_ff_persist, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (cell->rad_force_bf_persist, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->gain, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->converge_t_r, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->converge_t_e, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->converge_hc, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->trcheck, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->techeck, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->hccheck, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->converge_whole, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->converging, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->ip, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&cell->xi, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
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
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->ne, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->rho, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->vol, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->xgamma, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->density, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->partition, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->levden, nlte_levels, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->kappa_ff_factor, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->recomb_simple, nphot_total, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->recomb_simple_upweight, nphot_total, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->kpkt_emiss, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->kpkt_abs, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->kbf_use, nphot_total, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->kbf_nuse, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->t_r_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->t_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->dt_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->dt_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->heat_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->heat_tot_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->abs_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->heat_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->heat_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->heat_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->heat_ind_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->heat_lines_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->heat_photo_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->heat_photo, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->heat_z, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->heat_auger, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->heat_ch_ex, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->abs_photo, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->abs_auger, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->w, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->ntot, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->ntot_star, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->ntot_bl, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->ntot_disk, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->ntot_wind, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->ntot_agn, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->nscat_es, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->nscat_res, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->mean_ds, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->n_ds, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->nrad, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->nioniz, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->ioniz, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->recomb, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->inner_ioniz, n_inner_tot, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->inner_recomb, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->scatters, nions, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->xscatters, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->heat_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->heat_inner_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->cool_rr_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->lum_rr_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->cool_dr_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->j, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->ave_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->xj, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->xave_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->fmin, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->fmax, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->fmin_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->fmax_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->xsd_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->nxtot, NXBANDS, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->spec_mod_type, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->pl_alpha, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->pl_log_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->exp_temp, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->exp_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->cell_spec_flux, NBINS_IN_CELL_SPEC, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->F_vis, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->F_UV, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->F_Xray, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->F_vis_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->F_UV_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->F_Xray_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->F_UV_ang_theta, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->F_UV_ang_phi, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->F_UV_ang_r, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->F_UV_ang_theta_persist, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->F_UV_ang_phi_persist, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->F_UV_ang_r_persist, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->j_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->j_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->ip_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->ip_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->max_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->cool_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->lum_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->lum_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->cool_adiabatic, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->lum_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->lum_rr_metals, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->cool_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->cool_di, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->cool_dr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->cool_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->cool_rr_metals, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->lum_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->lum_tot_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->cool_tot_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->lum_lines_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->lum_ff_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->cool_adiabatic_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->lum_rr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->cool_comp_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->cool_di_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->cool_dr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->cool_rr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->cool_rr_metals_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->lum_tot_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->heat_shock, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->bf_simple_ionpool_in, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->bf_simple_ionpool_out, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->n_bf_in, N_PHOT_PROC, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->n_bf_out, N_PHOT_PROC, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->comp_nujnu, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->dmo_dt, N_DMO_DT_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->rad_force_es, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->rad_force_ff, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->rad_force_bf, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->rad_force_es_persist, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->rad_force_ff_persist, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, cell->rad_force_bf_persist, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->gain, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->converge_t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->converge_t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->converge_hc, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->trcheck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->techeck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->hccheck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->converge_whole, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->converging, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->ip, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &cell->xi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
      }
    }
  }

  free (comm_buffer);
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
 * The communication pattern is as outlined in dated_plasma_properties.
 *
 * ### Notes ###
 *
 * When this is called in wind update, there is redundant information being
 * communicated in `dated_plasma_properties` which communicates the exact (but
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

  for (current_rank = 0; current_rank < np_mpi_global; ++current_rank)
  {
    position = 0;

    if (rank_global == current_rank)
    {
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
      for (n_plasma = n_start; n_plasma < n_stop; ++n_plasma)
      {
        MPI_Pack (&n_plasma, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_tot, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_ff, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_rr, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_lines, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
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
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[cell].lum_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[cell].lum_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[cell].lum_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[cell].lum_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);

      }
    }
  }

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
 * The communication pattern is as outlined in dated_plasma_properties.
 *
 * ### Notes ###
 *
 * When this is called in wind update, there is redundant information being
 * communicated in `dated_plasma_properties` which communicates the exact (but
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

  for (current_rank = 0; current_rank < np_mpi_global; ++current_rank)
  {
    position = 0;

    if (rank_global == current_rank)
    {
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
      for (i = n_start; i < n_stop; ++i)
      {
        MPI_Pack (&i, 1, MPI_INT, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].cool_tot, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].lum_ff, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].lum_lines, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].cool_rr, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].cool_comp, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].cool_di, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].cool_dr, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].cool_adiabatic, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[i].heat_shock, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &position, MPI_COMM_WORLD);
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
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].cool_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].lum_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].lum_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].cool_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].cool_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].cool_di, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].cool_dr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].cool_adiabatic, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &position, &plasmamain[n].heat_shock, 1, MPI_DOUBLE, MPI_COMM_WORLD);
      }
    }
  }

  d_xsignal (files.root, "%-20s Finished communicating wind cooling\n", "OK");
  free (comm_buffer);
#endif
}

/**********************************************************/
/**
 * @brief Communicate changing properties in the plasma cells between ranks.
 *
 * @param [in]  int n_start  the first cell this rank will communicate
 * @param [in]  int n_stop   the last cell this rank will communicate
 *
 * @details
 *
 * This makes sure each rank has an updated plasma grid. The way the
 * communication is setup is as follows:
 *
 * - We create a loop over each MPI rank in the MPI_COMM_WORLD communicator
 * - If the loop variable is equal to the current rank, the subset of cells that
 *   rank worked on are packed into `comm_buffer` which is broadcast to all
 *   ranks.
 * - All other ranks unpack that data into the plasma cell.
 *
 * As well as the properties of the plasma cells, the number of cells
 * communicated and the cell numbers are also communicated. The size of the
 * comm buffer is currently the minimum size required. To communicate more data
 * you need to increase the size of the comm buffer.
 *
 * @TODO: we need to find out what data is not required
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
  const int num_ints = 1 + n_cells_max * (20 + nphot_total + 2 * NXBANDS + 2 * N_PHOT_PROC + nions);
  const int num_doubles =
    n_cells_max * (71 + 1 * 3 + 9 * 4 + 6 * NFLUX_ANGLES + 3 * NFORCE_DIRECTIONS + 9 * nions + 1 * nlte_levels + 3 * nphot_total +
                   1 * n_inner_tot + 9 * NXBANDS + 1 * NBINS_IN_CELL_SPEC);
  const int size_of_comm_buffer = calculate_comm_buffer_size (num_ints, num_doubles);
  char *const comm_buffer = malloc (size_of_comm_buffer);

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
        MPI_Pack (&plasmamain[n_plasma].ne, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].rho, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].vol, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].xgamma, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].density, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].partition, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].levden, nlte_levels, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].kappa_ff_factor, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].recomb_simple, nphot_total, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].recomb_simple_upweight, nphot_total, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].kpkt_emiss, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].kpkt_abs, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].kbf_use, nphot_total, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].kbf_nuse, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].t_r, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].t_r_old, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].t_e, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].t_e_old, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].dt_e, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].dt_e_old, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].heat_tot, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].heat_tot_old, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].abs_tot, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].heat_lines, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].heat_ff, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].heat_comp, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].heat_ind_comp, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].heat_lines_macro, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].heat_photo_macro, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].heat_photo, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].heat_z, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].heat_auger, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].heat_ch_ex, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].abs_photo, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].abs_auger, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].w, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].ntot, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].ntot_star, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].ntot_bl, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].ntot_disk, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].ntot_wind, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].ntot_agn, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].nscat_es, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].mean_ds, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].n_ds, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].nrad, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].nioniz, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].ioniz, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].recomb, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].inner_ioniz, n_inner_tot, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].scatters, nions, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].xscatters, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].heat_ion, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].heat_inner_ion, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].cool_rr_ion, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].lum_rr_ion, nions, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].j, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].ave_freq, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].xj, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].xave_freq, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].fmin_mod, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].fmax_mod, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].xsd_freq, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].nxtot, NXBANDS, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].spec_mod_type, NXBANDS, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].pl_alpha, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].pl_log_w, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].exp_temp, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].exp_w, NXBANDS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].cell_spec_flux, NBINS_IN_CELL_SPEC, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].F_vis, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].F_UV, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].F_Xray, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].F_vis_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].F_UV_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].F_Xray_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].F_UV_ang_theta, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].F_UV_ang_phi, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].F_UV_ang_r, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].F_UV_ang_theta_persist, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].F_UV_ang_phi_persist, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].F_UV_ang_r_persist, NFLUX_ANGLES, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].j_direct, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].j_scatt, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].ip_direct, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].ip_scatt, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].max_freq, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].cool_tot, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_lines, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_ff, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].cool_adiabatic, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_rr, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_rr_metals, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].cool_comp, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].cool_di, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].cool_dr, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].cool_rr, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].cool_rr_metals, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_tot, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_tot_old, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].cool_tot_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_lines_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_ff_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].cool_adiabatic_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_rr_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].cool_comp_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].cool_di_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].cool_dr_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].cool_rr_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].cool_rr_metals_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_tot_ioniz, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].heat_shock, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].bf_simple_ionpool_in, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].bf_simple_ionpool_out, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].n_bf_in, N_PHOT_PROC, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].n_bf_out, N_PHOT_PROC, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].comp_nujnu, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].dmo_dt, N_DMO_DT_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].rad_force_es, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].rad_force_ff, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].rad_force_bf, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].rad_force_es_persist, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].rad_force_ff_persist, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n_plasma].rad_force_bf_persist, NFORCE_DIRECTIONS, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position,
                  MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].gain, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].converge_t_r, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].converge_t_e, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].converge_hc, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].trcheck, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].techeck, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].hccheck, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].converge_whole, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].converging, 1, MPI_INT, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].ip, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].xi, 1, MPI_DOUBLE, comm_buffer, size_of_comm_buffer, &position, MPI_COMM_WORLD);

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
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].ne, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].rho, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].vol, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].xgamma, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].density, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].partition, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].levden, nlte_levels, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].kappa_ff_factor, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].recomb_simple, nphot_total, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].recomb_simple_upweight, nphot_total, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].kpkt_emiss, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].kpkt_abs, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].kbf_use, nphot_total, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].kbf_nuse, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].t_r_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].t_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].dt_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].dt_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].heat_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].heat_tot_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].abs_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].heat_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].heat_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].heat_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].heat_ind_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].heat_lines_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].heat_photo_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].heat_photo, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].heat_z, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].heat_auger, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].heat_ch_ex, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].abs_photo, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].abs_auger, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].w, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].ntot, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].ntot_star, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].ntot_bl, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].ntot_disk, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].ntot_wind, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].ntot_agn, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].nscat_es, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].mean_ds, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].n_ds, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].nrad, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].nioniz, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].ioniz, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].recomb, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].inner_ioniz, n_inner_tot, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].scatters, nions, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].xscatters, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].heat_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].heat_inner_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].cool_rr_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].lum_rr_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].j, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].ave_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].xj, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].xave_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].fmin_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].fmax_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].xsd_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].nxtot, NXBANDS, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].spec_mod_type, NXBANDS, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].pl_alpha, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].pl_log_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].exp_temp, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].exp_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].cell_spec_flux, NBINS_IN_CELL_SPEC, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].F_vis, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].F_UV, NFORCE_DIRECTIONS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].F_Xray, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].F_vis_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].F_UV_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].F_Xray_persistent, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].F_UV_ang_theta, NFLUX_ANGLES, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].F_UV_ang_phi, NFLUX_ANGLES, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].F_UV_ang_r, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].F_UV_ang_theta_persist, NFLUX_ANGLES, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].F_UV_ang_phi_persist, NFLUX_ANGLES, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].F_UV_ang_r_persist, NFLUX_ANGLES, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].j_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].j_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].ip_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].ip_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].max_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].cool_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].lum_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].lum_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].cool_adiabatic, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].lum_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].lum_rr_metals, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].cool_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].cool_di, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].cool_dr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].cool_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].cool_rr_metals, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].lum_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].lum_tot_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].cool_tot_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].lum_lines_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].lum_ff_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].cool_adiabatic_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].lum_rr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].cool_comp_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].cool_di_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].cool_dr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].cool_rr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].cool_rr_metals_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].lum_tot_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].heat_shock, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].bf_simple_ionpool_in, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].bf_simple_ionpool_out, 1, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].n_bf_in, N_PHOT_PROC, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].n_bf_out, N_PHOT_PROC, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].comp_nujnu, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].dmo_dt, N_DMO_DT_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].rad_force_es, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].rad_force_ff, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].rad_force_bf, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].rad_force_es_persist, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].rad_force_ff_persist, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, plasmamain[n_plasma].rad_force_bf_persist, NFORCE_DIRECTIONS, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].gain, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].converge_t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].converge_t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].converge_hc, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].trcheck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].techeck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].hccheck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].converge_whole, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].converging, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].ip, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, size_of_comm_buffer, &position, &plasmamain[n_plasma].xi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
      }
    }
  }

  free (comm_buffer);
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
    maxfreqhelper[mpi_i] = plasmamain[mpi_i].max_freq;
    redhelper[mpi_i] = plasmamain[mpi_i].j / np_mpi_global;
    redhelper[mpi_i + NPLASMA] = plasmamain[mpi_i].ave_freq / np_mpi_global;
    redhelper[mpi_i + 2 * NPLASMA] = plasmamain[mpi_i].cool_tot / np_mpi_global;
    redhelper[mpi_i + 3 * NPLASMA] = plasmamain[mpi_i].heat_tot / np_mpi_global;
    redhelper[mpi_i + 4 * NPLASMA] = plasmamain[mpi_i].heat_lines / np_mpi_global;
    redhelper[mpi_i + 5 * NPLASMA] = plasmamain[mpi_i].heat_ff / np_mpi_global;
    redhelper[mpi_i + 6 * NPLASMA] = plasmamain[mpi_i].heat_comp / np_mpi_global;
    redhelper[mpi_i + 7 * NPLASMA] = plasmamain[mpi_i].heat_ind_comp / np_mpi_global;
    redhelper[mpi_i + 8 * NPLASMA] = plasmamain[mpi_i].heat_photo / np_mpi_global;
    redhelper[mpi_i + 9 * NPLASMA] = plasmamain[mpi_i].ip / np_mpi_global;
    redhelper[mpi_i + 10 * NPLASMA] = plasmamain[mpi_i].j_direct / np_mpi_global;
    redhelper[mpi_i + 11 * NPLASMA] = plasmamain[mpi_i].j_scatt / np_mpi_global;
    redhelper[mpi_i + 12 * NPLASMA] = plasmamain[mpi_i].ip_direct / np_mpi_global;
    redhelper[mpi_i + 13 * NPLASMA] = plasmamain[mpi_i].ip_scatt / np_mpi_global;
    redhelper[mpi_i + 14 * NPLASMA] = plasmamain[mpi_i].heat_auger / np_mpi_global;
    redhelper[mpi_i + 15 * NPLASMA] = plasmamain[mpi_i].rad_force_es[0] / np_mpi_global;
    redhelper[mpi_i + 16 * NPLASMA] = plasmamain[mpi_i].rad_force_es[1] / np_mpi_global;
    redhelper[mpi_i + 17 * NPLASMA] = plasmamain[mpi_i].rad_force_es[2] / np_mpi_global;
    redhelper[mpi_i + 18 * NPLASMA] = plasmamain[mpi_i].rad_force_es[3] / np_mpi_global;
    redhelper[mpi_i + 19 * NPLASMA] = plasmamain[mpi_i].F_vis[0] / np_mpi_global;
    redhelper[mpi_i + 20 * NPLASMA] = plasmamain[mpi_i].F_vis[1] / np_mpi_global;
    redhelper[mpi_i + 21 * NPLASMA] = plasmamain[mpi_i].F_vis[2] / np_mpi_global;
    redhelper[mpi_i + 22 * NPLASMA] = plasmamain[mpi_i].F_vis[3] / np_mpi_global;
    redhelper[mpi_i + 23 * NPLASMA] = plasmamain[mpi_i].F_UV[0] / np_mpi_global;
    redhelper[mpi_i + 24 * NPLASMA] = plasmamain[mpi_i].F_UV[1] / np_mpi_global;
    redhelper[mpi_i + 25 * NPLASMA] = plasmamain[mpi_i].F_UV[2] / np_mpi_global;
    redhelper[mpi_i + 26 * NPLASMA] = plasmamain[mpi_i].F_UV[3] / np_mpi_global;
    redhelper[mpi_i + 27 * NPLASMA] = plasmamain[mpi_i].F_Xray[0] / np_mpi_global;
    redhelper[mpi_i + 28 * NPLASMA] = plasmamain[mpi_i].F_Xray[1] / np_mpi_global;
    redhelper[mpi_i + 29 * NPLASMA] = plasmamain[mpi_i].F_Xray[2] / np_mpi_global;
    redhelper[mpi_i + 30 * NPLASMA] = plasmamain[mpi_i].F_Xray[3] / np_mpi_global;
    redhelper[mpi_i + 31 * NPLASMA] = plasmamain[mpi_i].rad_force_bf[0] / np_mpi_global;
    redhelper[mpi_i + 32 * NPLASMA] = plasmamain[mpi_i].rad_force_bf[1] / np_mpi_global;
    redhelper[mpi_i + 33 * NPLASMA] = plasmamain[mpi_i].rad_force_bf[2] / np_mpi_global;
    redhelper[mpi_i + 34 * NPLASMA] = plasmamain[mpi_i].rad_force_bf[3] / np_mpi_global;
    redhelper[mpi_i + 35 * NPLASMA] = plasmamain[mpi_i].rad_force_ff[0] / np_mpi_global;
    redhelper[mpi_i + 36 * NPLASMA] = plasmamain[mpi_i].rad_force_ff[1] / np_mpi_global;
    redhelper[mpi_i + 37 * NPLASMA] = plasmamain[mpi_i].rad_force_ff[2] / np_mpi_global;
    redhelper[mpi_i + 38 * NPLASMA] = plasmamain[mpi_i].rad_force_ff[3] / np_mpi_global;

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      redhelper[mpi_i + (39 + mpi_j) * NPLASMA] = plasmamain[mpi_i].xj[mpi_j] / np_mpi_global;
      redhelper[mpi_i + (39 + NXBANDS + mpi_j) * NPLASMA] = plasmamain[mpi_i].xave_freq[mpi_j] / np_mpi_global;
      redhelper[mpi_i + (39 + 2 * NXBANDS + mpi_j) * NPLASMA] = plasmamain[mpi_i].xsd_freq[mpi_j] / np_mpi_global;
      /* 131213 NSH populate the band limited min and max frequency arrays */
      maxbandfreqhelper[mpi_i * NXBANDS + mpi_j] = plasmamain[mpi_i].fmax[mpi_j];
      minbandfreqhelper[mpi_i * NXBANDS + mpi_j] = plasmamain[mpi_i].fmin[mpi_j];
    }
    for (mpi_j = 0; mpi_j < nions; mpi_j++)
    {
      ion_helper[mpi_i * nions + mpi_j] = plasmamain[mpi_i].ioniz[mpi_j] / np_mpi_global;
    }
    for (mpi_j = 0; mpi_j < n_inner_tot; mpi_j++)
    {
      inner_ion_helper[mpi_i * n_inner_tot + mpi_j] = plasmamain[mpi_i].inner_ioniz[mpi_j] / np_mpi_global;
    }
    for (mpi_j = 0; mpi_j < NFLUX_ANGLES; mpi_j++)
    {
      flux_helper[mpi_i * (3 * NFLUX_ANGLES) + mpi_j] = plasmamain[mpi_i].F_UV_ang_theta[mpi_j] / np_mpi_global;
      flux_helper[mpi_i * (3 * NFLUX_ANGLES) + NFLUX_ANGLES + mpi_j] = plasmamain[mpi_i].F_UV_ang_phi[mpi_j] / np_mpi_global;
      flux_helper[mpi_i * (3 * NFLUX_ANGLES) + 2 * NFLUX_ANGLES + mpi_j] = plasmamain[mpi_i].F_UV_ang_r[mpi_j] / np_mpi_global;
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
    plasmamain[mpi_i].max_freq = maxfreqhelper2[mpi_i];
    plasmamain[mpi_i].j = redhelper2[mpi_i];
    plasmamain[mpi_i].ave_freq = redhelper2[mpi_i + NPLASMA];
    plasmamain[mpi_i].cool_tot = redhelper2[mpi_i + 2 * NPLASMA];
    plasmamain[mpi_i].heat_tot = redhelper2[mpi_i + 3 * NPLASMA];
    plasmamain[mpi_i].heat_lines = redhelper2[mpi_i + 4 * NPLASMA];
    plasmamain[mpi_i].heat_ff = redhelper2[mpi_i + 5 * NPLASMA];
    plasmamain[mpi_i].heat_comp = redhelper2[mpi_i + 6 * NPLASMA];
    plasmamain[mpi_i].heat_ind_comp = redhelper2[mpi_i + 7 * NPLASMA];
    plasmamain[mpi_i].heat_photo = redhelper2[mpi_i + 8 * NPLASMA];
    plasmamain[mpi_i].ip = redhelper2[mpi_i + 9 * NPLASMA];
    plasmamain[mpi_i].j_direct = redhelper2[mpi_i + 10 * NPLASMA];
    plasmamain[mpi_i].j_scatt = redhelper2[mpi_i + 11 * NPLASMA];
    plasmamain[mpi_i].ip_direct = redhelper2[mpi_i + 12 * NPLASMA];
    plasmamain[mpi_i].ip_scatt = redhelper2[mpi_i + 13 * NPLASMA];
    plasmamain[mpi_i].heat_auger = redhelper2[mpi_i + 14 * NPLASMA];
    plasmamain[mpi_i].rad_force_es[0] = redhelper2[mpi_i + 15 * NPLASMA];
    plasmamain[mpi_i].rad_force_es[1] = redhelper2[mpi_i + 16 * NPLASMA];
    plasmamain[mpi_i].rad_force_es[2] = redhelper2[mpi_i + 17 * NPLASMA];
    plasmamain[mpi_i].rad_force_es[3] = redhelper2[mpi_i + 18 * NPLASMA];
    plasmamain[mpi_i].F_vis[0] = redhelper2[mpi_i + 19 * NPLASMA];
    plasmamain[mpi_i].F_vis[1] = redhelper2[mpi_i + 20 * NPLASMA];
    plasmamain[mpi_i].F_vis[2] = redhelper2[mpi_i + 21 * NPLASMA];
    plasmamain[mpi_i].F_vis[3] = redhelper2[mpi_i + 22 * NPLASMA];
    plasmamain[mpi_i].F_UV[0] = redhelper2[mpi_i + 23 * NPLASMA];
    plasmamain[mpi_i].F_UV[1] = redhelper2[mpi_i + 24 * NPLASMA];
    plasmamain[mpi_i].F_UV[2] = redhelper2[mpi_i + 25 * NPLASMA];
    plasmamain[mpi_i].F_UV[3] = redhelper2[mpi_i + 26 * NPLASMA];
    plasmamain[mpi_i].F_Xray[0] = redhelper2[mpi_i + 27 * NPLASMA];
    plasmamain[mpi_i].F_Xray[1] = redhelper2[mpi_i + 28 * NPLASMA];
    plasmamain[mpi_i].F_Xray[2] = redhelper2[mpi_i + 29 * NPLASMA];
    plasmamain[mpi_i].F_Xray[3] = redhelper2[mpi_i + 30 * NPLASMA];
    plasmamain[mpi_i].rad_force_bf[0] = redhelper2[mpi_i + 31 * NPLASMA];
    plasmamain[mpi_i].rad_force_bf[1] = redhelper2[mpi_i + 32 * NPLASMA];
    plasmamain[mpi_i].rad_force_bf[2] = redhelper2[mpi_i + 33 * NPLASMA];
    plasmamain[mpi_i].rad_force_bf[3] = redhelper2[mpi_i + 34 * NPLASMA];
    plasmamain[mpi_i].rad_force_ff[0] = redhelper2[mpi_i + 35 * NPLASMA];
    plasmamain[mpi_i].rad_force_ff[1] = redhelper2[mpi_i + 36 * NPLASMA];
    plasmamain[mpi_i].rad_force_ff[2] = redhelper2[mpi_i + 37 * NPLASMA];
    plasmamain[mpi_i].rad_force_ff[3] = redhelper2[mpi_i + 38 * NPLASMA];

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      plasmamain[mpi_i].xj[mpi_j] = redhelper2[mpi_i + (39 + mpi_j) * NPLASMA];
      plasmamain[mpi_i].xave_freq[mpi_j] = redhelper2[mpi_i + (39 + NXBANDS + mpi_j) * NPLASMA];
      plasmamain[mpi_i].xsd_freq[mpi_j] = redhelper2[mpi_i + (39 + NXBANDS * 2 + mpi_j) * NPLASMA];

      /* 131213 NSH And unpack the min and max banded frequencies to the plasma array */
      plasmamain[mpi_i].fmax[mpi_j] = maxbandfreqhelper2[mpi_i * NXBANDS + mpi_j];
      plasmamain[mpi_i].fmin[mpi_j] = minbandfreqhelper2[mpi_i * NXBANDS + mpi_j];
    }
    for (mpi_j = 0; mpi_j < nions; mpi_j++)
    {
      plasmamain[mpi_i].ioniz[mpi_j] = ion_helper2[mpi_i * nions + mpi_j];
    }
    for (mpi_j = 0; mpi_j < n_inner_tot; mpi_j++)
    {
      plasmamain[mpi_i].inner_ioniz[mpi_j] = inner_ion_helper2[mpi_i * n_inner_tot + mpi_j];
    }
    for (mpi_j = 0; mpi_j < NFLUX_ANGLES; mpi_j++)
    {
      plasmamain[mpi_i].F_UV_ang_theta[mpi_j] = flux_helper2[mpi_i * (3 * NFLUX_ANGLES) + mpi_j];
      plasmamain[mpi_i].F_UV_ang_phi[mpi_j] = flux_helper2[mpi_i * (3 * NFLUX_ANGLES) + NFLUX_ANGLES + mpi_j];
      plasmamain[mpi_i].F_UV_ang_r[mpi_j] = flux_helper2[mpi_i * (3 * NFLUX_ANGLES) + 2 * NFLUX_ANGLES + mpi_j];
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

  /* allocate the integer helper arrays, set a barrier, then do all the integers. */
  iqdisk_helper = calloc (sizeof (int), NRINGS * 2);
  iqdisk_helper2 = calloc (sizeof (int), NRINGS * 2);
  iredhelper = calloc (sizeof (int), plasma_int_helpers);
  iredhelper2 = calloc (sizeof (int), plasma_int_helpers);

  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    iredhelper[mpi_i] = plasmamain[mpi_i].ntot;
    iredhelper[mpi_i + NPLASMA] = plasmamain[mpi_i].ntot_star;
    iredhelper[mpi_i + 2 * NPLASMA] = plasmamain[mpi_i].ntot_bl;
    iredhelper[mpi_i + 3 * NPLASMA] = plasmamain[mpi_i].ntot_disk;
    iredhelper[mpi_i + 4 * NPLASMA] = plasmamain[mpi_i].ntot_wind;
    iredhelper[mpi_i + 5 * NPLASMA] = plasmamain[mpi_i].ntot_agn;
    iredhelper[mpi_i + 6 * NPLASMA] = plasmamain[mpi_i].nioniz;

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      iredhelper[mpi_i + (7 + mpi_j) * NPLASMA] = plasmamain[mpi_i].nxtot[mpi_j];
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
    plasmamain[mpi_i].ntot = iredhelper2[mpi_i];
    plasmamain[mpi_i].ntot_star = iredhelper2[mpi_i + NPLASMA];
    plasmamain[mpi_i].ntot_bl = iredhelper2[mpi_i + 2 * NPLASMA];
    plasmamain[mpi_i].ntot_disk = iredhelper2[mpi_i + 3 * NPLASMA];
    plasmamain[mpi_i].ntot_wind = iredhelper2[mpi_i + 4 * NPLASMA];
    plasmamain[mpi_i].ntot_agn = iredhelper2[mpi_i + 5 * NPLASMA];
    plasmamain[mpi_i].nioniz = iredhelper2[mpi_i + 6 * NPLASMA];

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      plasmamain[mpi_i].nxtot[mpi_j] = iredhelper2[mpi_i + (7 + mpi_j) * NPLASMA];
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
        redhelper[mpi_i * NPLASMA + mpi_j] = plasmamain[mpi_j].cell_spec_flux[mpi_i] / np_mpi_global;

      }
    }

    MPI_Allreduce (redhelper, redhelper2, size_of_commbuffer, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (mpi_i = 0; mpi_i < NBINS_IN_CELL_SPEC; mpi_i++)
    {
      for (mpi_j = 0; mpi_j < NPLASMA; mpi_j++)
      {
        plasmamain[mpi_j].cell_spec_flux[mpi_i] = redhelper2[mpi_i * NPLASMA + mpi_j];

      }
    }

    free (redhelper);
    free (redhelper2);
  }

  d_xsignal (files.root, "%-20s Finished reduction of simple estimators\n", "OK");

#endif
  return (0);
}
