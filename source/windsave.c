
/***********************************************************/
/** @file  windsave.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  Routines to save and read in the structues the constitute
 * a model and and spectra which have been genrerated
 *
 * 
 * The first two routines in this file write and read the wind structure.  		
 * The second two routines do the same thing for the spectrum structure
 * 
 * ### Notes ###
 *
 * The files here are all written out as binary files.  They are
 * used for restars, and also by routines like swind and windsave2talbe
 * which inspect what is happening in the wind.
 *
 * There are separate ascii_writing 
 * routines for writing the spectra out for plotting.)
 * 
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/** 
 * @brief      Save all of the strutures associated with the 
 * wind to a file
 *
 * @param [in] char  filename[]   The name of the file to write to
 * @return     The number of successful writes
 *
 * @details
 *
 * ### Notes ###
 *
 * For the most part, adding a variable to the structures geo,
 * or plasma, does not require changes to this routine, unless
 * new variable length arrays are involved.
 *
 **********************************************************/

int
wind_save (char filename[])
{
  FILE *fptr;
  char header[LINELENGTH];
  int ndom;
  int m;
  int n;

  if ((fptr = fopen (filename, "w")) == NULL)
  {
    Error ("wind_save: Unable to open %s\n", filename);
    Exit (0);
  }

  memset (header, ' ', sizeof (header));
  header[sizeof (header) - 1] = '\0';   // ensure termination
  snprintf (header, sizeof (header), "Version %s\n", VERSION);
  n = fwrite (header, sizeof (header), 1, fptr);


//OLD  sprintf (header, "Version %s\n", VERSION);
//OLD  n = fwrite (header, sizeof (header), 1, fptr);

  n += fwrite (&geo, sizeof (geo), 1, fptr);

  n += fwrite (zdom, sizeof (domain_dummy), geo.ndomain, fptr);
  for (ndom = 0; ndom < geo.ndomain; ++ndom)
  {
    n += fwrite (zdom[ndom].wind_x, sizeof (double), zdom[ndom].ndim, fptr);
    n += fwrite (zdom[ndom].wind_z, sizeof (double), zdom[ndom].mdim, fptr);
    n += fwrite (zdom[ndom].wind_midx, sizeof (double), zdom[ndom].ndim, fptr);
    n += fwrite (zdom[ndom].wind_midz, sizeof (double), zdom[ndom].mdim, fptr);

    if (zdom[ndom].coord_type == CYLVAR)
    {
      n += fwrite (zdom[ndom].wind_z_var, sizeof (double), zdom[ndom].ndim * zdom[ndom].mdim, fptr);
      n += fwrite (zdom[ndom].wind_midz_var, sizeof (double), zdom[ndom].ndim * zdom[ndom].mdim, fptr);
    }
  }

  n += fwrite (wmain, sizeof (wind_dummy), NDIM2, fptr);
  n += fwrite (&disk, sizeof (disk), 1, fptr);
  n += fwrite (&qdisk, sizeof (disk), 1, fptr);
  n += fwrite (plasmamain, sizeof (plasma_dummy), NPLASMA, fptr);

/* Write out the variable length arrays
in the plasma structure */

  for (m = 0; m < NPLASMA; m++)
  {
    n += fwrite (plasmamain[m].state.density, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].state.partition, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].est.ioniz, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].derived.recomb, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].derived.inner_recomb, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].derived.scatters, sizeof (int), nions, fptr);
    n += fwrite (plasmamain[m].derived.xscatters, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].est.heat_ion, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].derived.cool_rr_ion, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].derived.cool_dr_ion, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].derived.lum_rr_ion, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].state.levden, sizeof (double), nlte_levels, fptr);
    n += fwrite (plasmamain[m].state.recomb_simple, sizeof (double), nphot_total, fptr);
    n += fwrite (plasmamain[m].state.recomb_simple_upweight, sizeof (double), nphot_total, fptr);
    n += fwrite (plasmamain[m].state.kbf_use, sizeof (double), nphot_total, fptr);

    /* Fixed-size arrays now in contiguous blocks */
    n += fwrite (plasmamain[m].state.f1, sizeof (double), NXBANDS + 1, fptr);
    n += fwrite (plasmamain[m].state.f2, sizeof (double), NXBANDS + 1, fptr);
    n += fwrite (plasmamain[m].state.spec_mod_type, sizeof (int), NXBANDS, fptr);
    n += fwrite (plasmamain[m].state.pl_alpha, sizeof (double), NXBANDS, fptr);
    n += fwrite (plasmamain[m].state.pl_log_w, sizeof (double), NXBANDS, fptr);
    n += fwrite (plasmamain[m].state.exp_temp, sizeof (double), NXBANDS, fptr);
    n += fwrite (plasmamain[m].state.exp_w, sizeof (double), NXBANDS, fptr);
    n += fwrite (plasmamain[m].state.fmin_mod, sizeof (double), NXBANDS, fptr);
    n += fwrite (plasmamain[m].state.fmax_mod, sizeof (double), NXBANDS, fptr);
    n += fwrite (plasmamain[m].derived.F_vis_persistent, sizeof (double), NFORCE_DIRECTIONS, fptr);
    n += fwrite (plasmamain[m].derived.F_UV_persistent, sizeof (double), NFORCE_DIRECTIONS, fptr);
    n += fwrite (plasmamain[m].derived.F_Xray_persistent, sizeof (double), NFORCE_DIRECTIONS, fptr);
    n += fwrite (plasmamain[m].derived.rad_force_es_persist, sizeof (double), NFORCE_DIRECTIONS, fptr);
    n += fwrite (plasmamain[m].derived.rad_force_ff_persist, sizeof (double), NFORCE_DIRECTIONS, fptr);
    n += fwrite (plasmamain[m].derived.rad_force_bf_persist, sizeof (double), NFORCE_DIRECTIONS, fptr);
    n += fwrite (plasmamain[m].derived.F_UV_ang_theta_persist, sizeof (double), NFLUX_ANGLES, fptr);
    n += fwrite (plasmamain[m].derived.F_UV_ang_phi_persist, sizeof (double), NFLUX_ANGLES, fptr);
    n += fwrite (plasmamain[m].derived.F_UV_ang_r_persist, sizeof (double), NFLUX_ANGLES, fptr);
    n += fwrite (plasmamain[m].derived.n_bf_in, sizeof (int), nphot_total, fptr);
    n += fwrite (plasmamain[m].derived.n_bf_out, sizeof (int), nphot_total, fptr);
  }

/* Now write out the macro atom info */

  if (geo.nmacro)
  {
    n += fwrite (macromain, sizeof (macro_dummy), NPLASMA, fptr);
    for (m = 0; m < NPLASMA; m++)
    {
      n += fwrite (macromain[m].est.jbar, sizeof (double), size_Jbar_est, fptr);
      n += fwrite (macromain[m].state.jbar_old, sizeof (double), size_Jbar_est, fptr);
      n += fwrite (macromain[m].est.gamma, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].state.gamma_old, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].est.gamma_e, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].state.gamma_e_old, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].est.alpha_st, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].state.alpha_st_old, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].est.alpha_st_e, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].state.alpha_st_e_old, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].est.recomb_sp, sizeof (double), size_alpha_est, fptr);
      n += fwrite (macromain[m].est.recomb_sp_e, sizeof (double), size_alpha_est, fptr);
      n += fwrite (macromain[m].derived.matom_emiss, sizeof (double), nlevels_macro, fptr);
      n += fwrite (macromain[m].est.matom_abs, sizeof (double), nlevels_macro, fptr);
    }
  }

  /* Per-wind-cell diagnostic counters: both hemispheres (2*NDIM2 entries).
   * Written last so that files produced before this field was added can still
   * be read — wind_read() silently leaves the array at zero on EOF. */
  n += fwrite (wind_counts_main, sizeof (wind_counts_dummy), 2 * NDIM2, fptr);

  fclose (fptr);

  Log_silent
    ("wind_write sizes: NPLASMA %d size_Jbar_est %d size_gamma_est %d size_alpha_est %d nlevels_macro %d\n",
     NPLASMA, size_Jbar_est, size_gamma_est, size_alpha_est, nlevels_macro);

  return (n);

}

/*

   wind_read (filename)

   History
	11dec	ksl	Updated so returns -1 if it cannot open the windsave file.  This
			was done to enable one to handle missing files differently in
			different cases
	14jul	nsh	Code added to read in variable length arrays in plasma structure
	15aug	ksl	Updated to read domain structure
	15oct	ksl	Updated to read disk and qdisk stuctures
*/


/**********************************************************/
/** 
 * @brief      Read back the windsavefile 
 *
 * @param [in] char  filename[]   The full name of the windsave file
 * @return     The number of successful reads, or -1 if the file cannot 
 * be opened
 *
 * @details
 * 
 * The routine reads in both the windsave file and the
 * associated atomic data files for a model. It also reads the
 * disk and qdisk structures.
 *
 *
 * ### Notes ###
 *
 * ### Programming Comment ### 
 * This routine calls wind_complete. This looks superfluous, since 
 * wind_complete and its subsidiary routines but it
 * also appears harmless.  ksl 
 *
 **********************************************************/

int
wind_read (char filename[])
{
  FILE *fptr;
  int ndom;
  int n, m;
  char header[LINELENGTH];
  char version[LINELENGTH];
  struct stat file_stat;        // Used to check the atomic data exists

  if ((fptr = fopen (filename, "r")) == NULL)
  {
    return (-1);
  }

  n = fread (header, sizeof (header), 1, fptr);
  sscanf (header, "%*s %s", version);
  Log ("Reading Windfile %s created with sirocco version %s with sirocco version %s\n", filename, version, VERSION);

  /* Now read in the geo structure */

  n += fread (&geo, sizeof (geo), 1, fptr);

  /* Null out pointer fields that were serialized as raw bytes — they will be
   * re-allocated when bands_init() runs.  Without this, the stale pointer
   * from the previous process could cause a double-free or corruption. */
  geo.cell_freq = NULL;

  /* Read the atomic data file.  This is necessary to do here in order to establish the 
   * values for the dimensionality of some of the variable length structures, associated 
   * with macro atoms, especially but likely to be a good idea ovrall
   */

  if (stat (geo.atomic_filename, &file_stat))
  {
    if (system ("Setup_Sirocco_Dir"))
    {
      Error ("Unable to open %s or create link for atomic data\n", geo.atomic_filename);
      Exit (1);
    }
  }

  get_atomic_data (geo.atomic_filename);

/* Now allocate space for the wind array */

  NDIM2 = geo.ndim2;
  NPLASMA = geo.nplasma;

  n += fread (zdom, sizeof (domain_dummy), geo.ndomain, fptr);
  for (ndom = 0; ndom < geo.ndomain; ++ndom)
  {
    allocate_domain_wind_coords (ndom);
    n += fread (zdom[ndom].wind_x, sizeof (double), zdom[ndom].ndim, fptr);
    n += fread (zdom[ndom].wind_z, sizeof (double), zdom[ndom].mdim, fptr);
    n += fread (zdom[ndom].wind_midx, sizeof (double), zdom[ndom].ndim, fptr);
    n += fread (zdom[ndom].wind_midz, sizeof (double), zdom[ndom].mdim, fptr);
    if (zdom[ndom].coord_type == CYLVAR)
    {
      cylvar_allocate_domain (ndom);
      n += fread (zdom[ndom].wind_z_var, sizeof (double), zdom[ndom].ndim * zdom[ndom].mdim, fptr);
      n += fread (zdom[ndom].wind_midz_var, sizeof (double), zdom[ndom].ndim * zdom[ndom].mdim, fptr);
    }
  }

  calloc_wind (2 * NDIM2);
#ifdef MPI_ON
  if (np_mpi_global > 1)
  {
    if (node_rank == 0)
    {
      n += fread (wmain, sizeof (wind_dummy), NDIM2, fptr);
    }
    else
    {
      /* Skip past the wind data in the file without reading into shared memory */
      fseek (fptr, (long) NDIM2 * sizeof (wind_dummy), SEEK_CUR);
    }
    MPI_Barrier (node_comm);
  }
  else
#endif
  {
    n += fread (wmain, sizeof (wind_dummy), NDIM2, fptr);
  }

  /* Reconstruct lower-hemisphere transport cells from the upper-hemisphere
   * cells that were just read from disk. */
  make_transport_grid ();

  /* Allocate per-wind-cell diagnostic counters; populated from the save file
   * below (backward-compatible: silently zero if the file predates this). */
  calloc_wind_counts (2 * NDIM2);

  /* Read the disk and qdisk structures */

  n += fread (&disk, sizeof (disk), 1, fptr);
  n += fread (&qdisk, sizeof (disk), 1, fptr);

  calloc_plasma (NPLASMA);

  n += fread (plasmamain, sizeof (plasma_dummy), NPLASMA, fptr);

  /*Allocate space for the dynamically allocated plasma arrays */

  calloc_dyn_plasma (NPLASMA);

  /* Read in the dynamically allocated plasma arrays */

  for (m = 0; m < NPLASMA; m++)
  {

    n += fread (plasmamain[m].state.density, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].state.partition, sizeof (double), nions, fptr);

    n += fread (plasmamain[m].est.ioniz, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].derived.recomb, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].derived.inner_recomb, sizeof (double), nions, fptr);

    n += fread (plasmamain[m].derived.scatters, sizeof (int), nions, fptr);
    n += fread (plasmamain[m].derived.xscatters, sizeof (double), nions, fptr);

    n += fread (plasmamain[m].est.heat_ion, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].derived.cool_rr_ion, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].derived.cool_dr_ion, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].derived.lum_rr_ion, sizeof (double), nions, fptr);

    n += fread (plasmamain[m].state.levden, sizeof (double), nlte_levels, fptr);
    n += fread (plasmamain[m].state.recomb_simple, sizeof (double), nphot_total, fptr);
    n += fread (plasmamain[m].state.recomb_simple_upweight, sizeof (double), nphot_total, fptr);
    n += fread (plasmamain[m].state.kbf_use, sizeof (double), nphot_total, fptr);

    /* Fixed-size arrays now in contiguous blocks */
    n += fread (plasmamain[m].state.f1, sizeof (double), NXBANDS + 1, fptr);
    n += fread (plasmamain[m].state.f2, sizeof (double), NXBANDS + 1, fptr);
    n += fread (plasmamain[m].state.spec_mod_type, sizeof (int), NXBANDS, fptr);
    n += fread (plasmamain[m].state.pl_alpha, sizeof (double), NXBANDS, fptr);
    n += fread (plasmamain[m].state.pl_log_w, sizeof (double), NXBANDS, fptr);
    n += fread (plasmamain[m].state.exp_temp, sizeof (double), NXBANDS, fptr);
    n += fread (plasmamain[m].state.exp_w, sizeof (double), NXBANDS, fptr);
    n += fread (plasmamain[m].state.fmin_mod, sizeof (double), NXBANDS, fptr);
    n += fread (plasmamain[m].state.fmax_mod, sizeof (double), NXBANDS, fptr);
    n += fread (plasmamain[m].derived.F_vis_persistent, sizeof (double), NFORCE_DIRECTIONS, fptr);
    n += fread (plasmamain[m].derived.F_UV_persistent, sizeof (double), NFORCE_DIRECTIONS, fptr);
    n += fread (plasmamain[m].derived.F_Xray_persistent, sizeof (double), NFORCE_DIRECTIONS, fptr);
    n += fread (plasmamain[m].derived.rad_force_es_persist, sizeof (double), NFORCE_DIRECTIONS, fptr);
    n += fread (plasmamain[m].derived.rad_force_ff_persist, sizeof (double), NFORCE_DIRECTIONS, fptr);
    n += fread (plasmamain[m].derived.rad_force_bf_persist, sizeof (double), NFORCE_DIRECTIONS, fptr);
    n += fread (plasmamain[m].derived.F_UV_ang_theta_persist, sizeof (double), NFLUX_ANGLES, fptr);
    n += fread (plasmamain[m].derived.F_UV_ang_phi_persist, sizeof (double), NFLUX_ANGLES, fptr);
    n += fread (plasmamain[m].derived.F_UV_ang_r_persist, sizeof (double), NFLUX_ANGLES, fptr);
    n += fread (plasmamain[m].derived.n_bf_in, sizeof (int), nphot_total, fptr);
    n += fread (plasmamain[m].derived.n_bf_out, sizeof (int), nphot_total, fptr);
  }


  /*Allocate space for macro-atoms and read in the data */

  if (geo.nmacro > 0)
  {
    calloc_macro (NPLASMA);
    n += fread (macromain, sizeof (macro_dummy), NPLASMA, fptr);
    calloc_estimators (NPLASMA);
    calloc_matom_matrix (NPLASMA);

    for (m = 0; m < NPLASMA; m++)
    {
      n += fread (macromain[m].est.jbar, sizeof (double), size_Jbar_est, fptr);
      n += fread (macromain[m].state.jbar_old, sizeof (double), size_Jbar_est, fptr);
      n += fread (macromain[m].est.gamma, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].state.gamma_old, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].est.gamma_e, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].state.gamma_e_old, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].est.alpha_st, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].state.alpha_st_old, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].est.alpha_st_e, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].state.alpha_st_e_old, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].est.recomb_sp, sizeof (double), size_alpha_est, fptr);
      n += fread (macromain[m].est.recomb_sp_e, sizeof (double), size_alpha_est, fptr);
      n += fread (macromain[m].derived.matom_emiss, sizeof (double), nlevels_macro, fptr);
      n += fread (macromain[m].est.matom_abs, sizeof (double), nlevels_macro, fptr);

      /* Force recalculation of kpkt_rates and matrix rates */

      macromain[m].derived.kpkt_rates_known = FALSE;
      macromain[m].derived.matrix_rates_known = FALSE;
    }

  }

  /* Per-wind-cell diagnostic counters — appended last for backward compat.
   * If the file was written before this field existed, fread returns 0 and
   * wind_counts_main stays zero (already calloc'd above). */
  {
    int nwc = fread (wind_counts_main, sizeof (wind_counts_dummy), 2 * NDIM2, fptr);
    if (nwc < 2 * NDIM2)
      Log_silent ("wind_read: wind_counts_main not in save file; counters zeroed\n");
  }

  fclose (fptr);

  wind_complete ();

  Log ("Read geometry and wind structures from windsavefile %s\n", filename);

  return (n);

}




/**********************************************************/
/** 
 * @brief      A driver routine that calls coordinate-system specific routines
 * that complete the description of the wind
 *
 * @return     Always returns 0
 *
 * @details
 *
 * For the most point, the various routines that are called
 * just recalculate some of the various arrays used for 
 * finding the boundaries of an individual cell.
 *
 * These basically are just 1-d versions of the coordinate
 * grids at the edge and mid-points of each grid cell
 *
 * ### Notes ###
 *
 * The need to call this routine whenever a windsave file is read
 * in could be obviated if all of the variable length arrays
 * in domains were actually written out.  But since these
 * variables are easy to calculate again, it is done here.
 * 
 *
 **********************************************************/

void
wind_complete ()
{
  int ndom;

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    if (zdom[ndom].coord_type == SPHERICAL)
    {
      spherical_wind_complete (ndom, wmain);
    }
    else if (zdom[ndom].coord_type == CYLIND)
    {
      cylind_wind_complete (ndom, wmain);
    }
    else if (zdom[ndom].coord_type == RTHETA)
    {
      rtheta_wind_complete (ndom, wmain);
    }
    else if (zdom[ndom].coord_type == CYLVAR)
    {
      cylvar_wind_complete (ndom, wmain);
    }
    else
    {
      Error ("wind_complete: Don't know how to complete coord_type %d\n", zdom[ndom].coord_type);
      Exit (0);
    }
  }
}


/**********************************************************/
/** 
 * @brief      Save all the spectra to a binary file
 *
 * @param [in] char  filename[]   The file to write to 
 * @return     The number of successful writes
 *
 * @details
 *
 * ### Notes ###
 * 
 * The program exits if one cannot write the file
 *
 **********************************************************/

int
spec_save (char filename[])
{

  FILE *fptr;
  char header[LINELENGTH];
  int count;
  int i;

  if ((fptr = fopen (filename, "w")) == NULL)
  {
    Error ("spec_save: Unable to open %s\n", filename);
    Exit (EXIT_FAILURE);
  }

  sprintf (header, "Version %s  nspectra %d NWAVE_IONIZ %d NWAVE_EXTRACT %d NWAVE_MAX %d\n", VERSION, nspectra, NWAVE_IONIZ,
           NWAVE_EXTRACT, NWAVE_MAX);

  count = (int) fwrite (header, sizeof (header), 1, fptr);
  count += (int) fwrite (xxspec, sizeof (spectrum_dummy), nspectra, fptr);

  for (i = 0; i < nspectra; ++i)
  {
    count += (int) fwrite (xxspec[i].f, sizeof (*xxspec[i].f), NWAVE_MAX, fptr);
    count += (int) fwrite (xxspec[i].lf, sizeof (*xxspec[i].lf), NWAVE_MAX, fptr);
    count += (int) fwrite (xxspec[i].f_wind, sizeof (*xxspec[i].f_wind), NWAVE_MAX, fptr);
    count += (int) fwrite (xxspec[i].lf_wind, sizeof (*xxspec[i].lf_wind), NWAVE_MAX, fptr);
  }

  fclose (fptr);

  return (count);
}



/**********************************************************/
/** 
 * @brief      Read the binary file containing all the spectra
 *
 * @param [in] char  filename[]   The name of the file
 * @return     The number of spectra that were read in
 *
 * @details
 *
 * The routine allocates space for the spectra and reads
 * then in.
 *
 * ### Notes ###
 *
 * The first line of the file contains the Sirocco version
 * and the number of spectra to be read in
 * 
 * The program exits if the file does not exist
 *
 **********************************************************/

int
spec_read (char filename[])
{
  FILE *fptr;
  int nhead, nwave_ioniz_check;
  int count;
  int i;
  char header[LINELENGTH];
  char version[LINELENGTH];

  if ((fptr = fopen (filename, "r")) == NULL)
  {
    Error ("spec_read: Unable to open %s\n", filename);
    Exit (1);
  }

  count = (int) fread (header, sizeof (header), 1, fptr);

  nhead = sscanf (header, "%*s %s %*s %d %*s %d %*s %d %*s %d", version, &nspectra, &nwave_ioniz_check, &NWAVE_EXTRACT, &NWAVE_MAX);
  if (nhead != 5)
  {
    Error ("Incorrect header format in %s\n", files.specsave);
    Exit (EXIT_FAILURE);
  }
  if (nwave_ioniz_check != (int) NWAVE_IONIZ)
  {
    Error ("The current NWAVE_IONIZ (%d) value is incompatible with the spec_save file which has NWAVE_IONIZ = %d\n", NWAVE_IONIZ,
           nwave_ioniz_check);
    Exit (EXIT_FAILURE);
  }

  Log
    ("Reading specfile %s with %d spectra and %d wavelength bins, created with sirocco version %s and currently using sirocco version %s\n",
     filename, nspectra, NWAVE_EXTRACT, version, VERSION);

  /* First allocate space */

  xxspec = calloc (nspectra, sizeof (spectrum_dummy));
  if (xxspec == NULL)
  {
    Error ("spectrum_init: Could not allocate memory for %d spectra\n", nspectra);
    Exit (EXIT_FAILURE);
  }
  count += (int) fread (xxspec, sizeof (spectrum_dummy), nspectra, fptr);

  /* Now read the rest of the file */

  spectrum_allocate (nspectra);
  for (i = 0; i < nspectra; ++i)
  {
    count += (int) fread (xxspec[i].f, sizeof (*xxspec[i].f), NWAVE_MAX, fptr);
    count += (int) fread (xxspec[i].lf, sizeof (*xxspec[i].lf), NWAVE_MAX, fptr);
    count += (int) fread (xxspec[i].f_wind, sizeof (*xxspec[i].f_wind), NWAVE_MAX, fptr);
    count += (int) fread (xxspec[i].lf_wind, sizeof (*xxspec[i].lf_wind), NWAVE_MAX, fptr);
  }

  fclose (fptr);

  Log ("Read spec structures from specfile %s\n", filename);

  return (count);

}
