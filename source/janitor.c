/***********************************************************/
/** @file  janitor.c
 * @author EJP
 * @date   December, 2023
 *
 * @brief  Contains functions for cleaning things up
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atomic.h"
#include "sirocco.h"

/**********************************************************/
/**
 * @brief  Free memory associated with the domains
 *
 * @details
 *
 **********************************************************/

void
free_domains (void)
{
  int i;

  for (i = 0; i < geo.ndomain; ++i)
  {
    free (zdom[i].wind_x);
    free (zdom[i].wind_midx);
    free (zdom[i].wind_z);
    free (zdom[i].wind_midz);

    if (zdom[i].coord_type == RTHETA)
    {
      free (zdom[i].cones_rtheta);
    }
    else if (zdom[i].coord_type == CYLVAR)
    {
      free (zdom[i].wind_z_var[0]);
      free (zdom[i].wind_z_var);
      free (zdom[i].wind_midz_var[0]);
      free (zdom[i].wind_midz_var);
    }
  }

  free (zdom);
}

/**********************************************************/
/**
 * @brief  Free memory associated with the wind grid
 *
 * @details
 *
 **********************************************************/

void
free_wind_grid (void)
{
  int n_wind;

  for (n_wind = 0; n_wind < NDIM2; ++n_wind)
  {
    free (wmain[n_wind].paths);
    free (wmain[n_wind].line_paths);
  }

  free (wmain);
}

/**********************************************************/
/**
 * @brief Free memory associated with the plasma grid
 *
 * @details
 *
 **********************************************************/

void
free_plasma_grid (void)
{
  int n_plasma;

  for (n_plasma = 0; n_plasma < NPLASMA + 1; ++n_plasma)
  {
    free (plasmamain[n_plasma].density);
    free (plasmamain[n_plasma].partition);
    free (plasmamain[n_plasma].ioniz);
    free (plasmamain[n_plasma].recomb);
    free (plasmamain[n_plasma].scatters);
    free (plasmamain[n_plasma].xscatters);
    free (plasmamain[n_plasma].heat_ion);
    free (plasmamain[n_plasma].heat_inner_ion);
    free (plasmamain[n_plasma].cool_rr_ion);
    free (plasmamain[n_plasma].lum_rr_ion);
    free (plasmamain[n_plasma].inner_recomb);
    free (plasmamain[n_plasma].inner_ioniz);
    free (plasmamain[n_plasma].cool_dr_ion);
    free (plasmamain[n_plasma].levden);
    free (plasmamain[n_plasma].recomb_simple);
    free (plasmamain[n_plasma].recomb_simple_upweight);
    free (plasmamain[n_plasma].kbf_use);
  }

  free (plasmamain);
}

/**********************************************************/
/**
 * @brief Free memory associated with the macro grid
 *
 * @details
 *
 **********************************************************/

void
free_macro_grid (void)
{
  int n_plasma;

  for (n_plasma = 0; n_plasma < NPLASMA + 1; n_plasma++)
  {
    free (macromain[n_plasma].jbar);
    free (macromain[n_plasma].jbar_old);
    free (macromain[n_plasma].gamma);
    free (macromain[n_plasma].gamma_old);
    free (macromain[n_plasma].gamma_e);
    free (macromain[n_plasma].gamma_e_old);
    free (macromain[n_plasma].alpha_st);
    free (macromain[n_plasma].alpha_st_old);
    free (macromain[n_plasma].alpha_st_e);
    free (macromain[n_plasma].alpha_st_e_old);
    free (macromain[n_plasma].recomb_sp);
    free (macromain[n_plasma].recomb_sp_e);
    free (macromain[n_plasma].matom_emiss);
    free (macromain[n_plasma].matom_abs);
    free (macromain[n_plasma].cooling_bf);
    free (macromain[n_plasma].cooling_bf_col);
    free (macromain[n_plasma].cooling_bb);
    if (macromain[n_plasma].store_matom_matrix == TRUE)
    {
      free (macromain[n_plasma].matom_matrix[0]);
      free (macromain[n_plasma].matom_matrix);
    }
  }

  free (macromain);
}

/**********************************************************/
/**
 * @brief Free memory associated with photon storage
 *
 * @details
 *
 **********************************************************/

void
free_photons (void)
{
  free (photstoremain);
  free (matomphotstoremain);
  free (photmain);
}

/**********************************************************/
/**
 * @brief Free memory associated with atomic data
 *
 * @details
 *
 **********************************************************/

void
free_atomic_data (void)
{
  free (ele);
  free (ion);
  free (xconfig);
  free (line);
  free (auger_macro);
}

/**********************************************************/
/**
 * @brief Free memory associated with the spectrum structure
 *
 * @details
 *
 **********************************************************/

void
free_spectra (void)
{
  int n_spec;

  for (n_spec = 0; n_spec < nspectra; ++n_spec)
  {
    free (xxspec[n_spec].f);
    free (xxspec[n_spec].lf);
    free (xxspec[n_spec].f_wind);
    free (xxspec[n_spec].lf_wind);
  }

  free (xxspec);
}

/**********************************************************/
/**
 * @brief Clean up memory usage at the end of the program
 *
 * @details
 *
 * This should be used to clean up allocated global structures
 * at the end of the program. Even though the OS will handle this
 * for us, cleaning up at the end will prevent false positives for
 * leak detection in tools like Valgrind.
 *
 **********************************************************/

void
clean_on_exit (void)
{
  free_domains ();
  free_wind_grid ();
  free_plasma_grid ();
  if (nlevels_macro > 0)
  {
    free_macro_grid ();
  }
  free_photons ();
  free_spectra ();
  free_atomic_data ();
}
