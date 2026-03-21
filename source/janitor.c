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
 * @brief  Free a contiguous block, using MPI_Win_free for shared blocks.
 *
 * @param [in,out]  ptr         Address of the block pointer (set to NULL after freeing)
 * @param [in]      is_shared   TRUE if the block was allocated with MPI shared memory
 *
 * @details
 * This is the counterpart to alloc_block_double/alloc_block_int in gridwind.c.
 * For shared blocks, we cannot call free() since the memory was allocated by
 * MPI_Win_allocate_shared; instead we would need MPI_Win_free on the corresponding
 * window.  However, at program exit the MPI runtime handles cleanup, so for
 * shared blocks we simply NULL the pointer without calling free.
 * For private (non-shared) blocks, we call free() as normal.
 **********************************************************/

static void
free_plasma_block (void **ptr, int is_shared)
{
  if (*ptr == NULL)
    return;

#ifdef MPI_ON
  if (is_shared && np_mpi_global > 1)
  {
    /* Shared memory is freed via MPI_Win_free, which is handled
     * by the cleanup in calloc_dyn_plasma/calloc_estimators when
     * re-allocating, or by MPI_Finalize at exit. */
    *ptr = NULL;
    return;
  }
#endif

  (void) is_shared;
  free (*ptr);
  *ptr = NULL;
}

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
  /* Free contiguous blocks instead of per-cell pointers, since all cell
   * pointers now index into contiguous blocks managed by plasma_block_ptrs.
   * For MPI shared memory blocks, use MPI_Win_free instead of free. */

  int is_shared = FALSE;
#ifdef MPI_ON
  is_shared = plasma_block_ptrs.shared_memory_active;
#endif

  /* state blocks (shared in MPI mode) */
  if (plasma_block_ptrs.density_block != NULL)
  {
    free_plasma_block ((void **) &plasma_block_ptrs.density_block, is_shared);
    free_plasma_block ((void **) &plasma_block_ptrs.partition_block, is_shared);
    free_plasma_block ((void **) &plasma_block_ptrs.levden_block, is_shared);
    free_plasma_block ((void **) &plasma_block_ptrs.recomb_simple_block, is_shared);
    free_plasma_block ((void **) &plasma_block_ptrs.recomb_simple_upweight_block, is_shared);
    free_plasma_block ((void **) &plasma_block_ptrs.kbf_use_block, is_shared);

    /* est blocks (always private) */
    free_plasma_block ((void **) &plasma_block_ptrs.ioniz_block, FALSE);
    free_plasma_block ((void **) &plasma_block_ptrs.heat_ion_block, FALSE);
    free_plasma_block ((void **) &plasma_block_ptrs.heat_inner_ion_block, FALSE);
    free_plasma_block ((void **) &plasma_block_ptrs.inner_ioniz_block, FALSE);

    /* derived blocks (shared in MPI mode, except scatters/xscatters which are private) */
    free_plasma_block ((void **) &plasma_block_ptrs.recomb_block, is_shared);
    free_plasma_block ((void **) &plasma_block_ptrs.scatters_block, FALSE);
    free_plasma_block ((void **) &plasma_block_ptrs.xscatters_block, FALSE);
    free_plasma_block ((void **) &plasma_block_ptrs.cool_rr_ion_block, is_shared);
    free_plasma_block ((void **) &plasma_block_ptrs.lum_rr_ion_block, is_shared);
    free_plasma_block ((void **) &plasma_block_ptrs.cool_dr_ion_block, is_shared);
    free_plasma_block ((void **) &plasma_block_ptrs.inner_recomb_block, is_shared);

    /* Fixed-size array blocks (shared in MPI mode) */
    free_plasma_block ((void **) &plasma_block_ptrs.state_xbands_dblock, is_shared);
    free_plasma_block ((void **) &plasma_block_ptrs.state_spec_mod_type_block, is_shared);
    free_plasma_block ((void **) &plasma_block_ptrs.derived_persist_force_block, is_shared);
    free_plasma_block ((void **) &plasma_block_ptrs.derived_persist_angle_block, is_shared);
    free_plasma_block ((void **) &plasma_block_ptrs.n_bf_in_block, FALSE);
    free_plasma_block ((void **) &plasma_block_ptrs.n_bf_out_block, FALSE);
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

  /* Free contiguous blocks instead of per-cell pointers */
  int is_shared = FALSE;
#ifdef MPI_ON
  is_shared = macro_block_ptrs.shared_memory_active;
#endif

  if (macro_block_ptrs.jbar_block != NULL)
  {
    /* state blocks (shared in MPI mode) */
    free_plasma_block ((void **) &macro_block_ptrs.jbar_old_block, is_shared);
    free_plasma_block ((void **) &macro_block_ptrs.gamma_old_block, is_shared);
    free_plasma_block ((void **) &macro_block_ptrs.gamma_e_old_block, is_shared);
    free_plasma_block ((void **) &macro_block_ptrs.alpha_st_old_block, is_shared);
    free_plasma_block ((void **) &macro_block_ptrs.alpha_st_e_old_block, is_shared);

    /* est blocks (always private) */
    free_plasma_block ((void **) &macro_block_ptrs.jbar_block, FALSE);
    free_plasma_block ((void **) &macro_block_ptrs.gamma_block, FALSE);
    free_plasma_block ((void **) &macro_block_ptrs.gamma_e_block, FALSE);
    free_plasma_block ((void **) &macro_block_ptrs.alpha_st_block, FALSE);
    free_plasma_block ((void **) &macro_block_ptrs.alpha_st_e_block, FALSE);
    free_plasma_block ((void **) &macro_block_ptrs.recomb_sp_block, FALSE);
    free_plasma_block ((void **) &macro_block_ptrs.recomb_sp_e_block, FALSE);
    free_plasma_block ((void **) &macro_block_ptrs.matom_abs_block, FALSE);
    free_plasma_block ((void **) &macro_block_ptrs.cooling_bf_block, FALSE);
    free_plasma_block ((void **) &macro_block_ptrs.cooling_bf_col_block, FALSE);
    free_plasma_block ((void **) &macro_block_ptrs.cooling_bb_block, FALSE);

    /* derived blocks (shared in MPI mode) */
    free_plasma_block ((void **) &macro_block_ptrs.matom_emiss_block, is_shared);
  }

  /* matom_matrix is still allocated per-cell (not part of contiguous blocks) */
  for (n_plasma = 0; n_plasma < NPLASMA + 1; n_plasma++)
  {
    if (macromain[n_plasma].state.store_matom_matrix == TRUE)
    {
      free (macromain[n_plasma].derived.matom_matrix[0]);
      free (macromain[n_plasma].derived.matom_matrix);
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
  Log_flush ();
//  if (photmain != NULL && photmain_allocated == TRUE)
//  {
//    Log ("CCC - Releasing photmain as part of janitor process\n");
//    free (photmain);
//    photmain_allocated = FALSE;
//  }
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
