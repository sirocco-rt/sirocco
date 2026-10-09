
/***********************************************************/
/** @file  wind_update_model.c
 * @author ksl
 * @date   September, 2026
 *
 * @brief  Update an existing wind with the density, velocity and
 * temperature of a new imported model
 *
 * ###Notes###
 *
 * This is used when Sirocco is coupled to a hydrodynamics code.
 * The hydro code advances the density, velocity and temperature,
 * and these are painted onto a windsave file from an earlier
 * Sirocco run, keeping the ionization state as a starting guess
 * for the next ionization cycles.
 *
 * The quantities are set following the same conventions as
 * define_wind, so that the updated wind is what Sirocco would have
 * produced from the new model, apart from the ionization state.
 * See issue #1190.
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
 * @brief      Update the wind in a domain from a new imported model
 *
 * @param [in] int  ndom   The domain to update
 * @param [in] char *filename   The model file, in the format used for
 * imported models
 * @return     0 on success.  The routine exits if the model does not
 * match the existing wind.
 *
 * @details
 *
 * The model must describe the same grid as the existing wind, and
 * the same cells must be in the wind, since the number of plasma
 * cells in a windsave file cannot be changed.
 *
 * For each cell the routine
 *
 * - sets the velocity from the model, and recalculates the velocity
 *   gradients, dv/ds and gamma factors
 * - sets rho to the clump density, model density / (fill * gamma),
 *   as in define_wind, and scales the ion densities by the change
 *   in rho, so that the ion fractions are unchanged
 * - recalculates the electron density from the ion densities
 * - if the model contains temperatures, sets t_e and t_r and
 *   recalculates the partition functions
 *
 **********************************************************/

int
update_wind_from_model (int ndom, char *filename)
{
  int n, nplasma, nion;
  int nstart, nstop;
  int nbad_x, nbad_inwind;
  double dx[3], old_xgamma_cen, old_rho, new_rho, fill;
  WindPtr scratch, cell;
  PlasmaPtr xplasma;

  if (zdom[ndom].wind_type != IMPORT)
  {
    Error ("update_wind_from_model: domain %d is not an imported model (wind_type %d)\n", ndom, zdom[ndom].wind_type);
    Exit (EXIT_FAILURE);
  }

  fill = zdom[ndom].fill;
  if (fill <= 0.0)
  {
    Error ("update_wind_from_model: invalid filling factor %e in domain %d\n", fill, ndom);
    Exit (EXIT_FAILURE);
  }

  import_wind2 (ndom, filename);

  if (imported_model[ndom].ndim != zdom[ndom].ndim || imported_model[ndom].mdim != zdom[ndom].mdim)
  {
    Error ("update_wind_from_model: model %s is %d x %d but the wind in domain %d is %d x %d\n",
           filename, imported_model[ndom].ndim, imported_model[ndom].mdim, ndom, zdom[ndom].ndim, zdom[ndom].mdim);
    Exit (EXIT_FAILURE);
  }

  /* Build the grid of the new model in a scratch copy of the wind, and check
   * that it matches the existing wind */

  if ((scratch = (WindPtr) calloc (NDIM2 + 1, sizeof (wind_dummy))) == NULL)
  {
    Error ("update_wind_from_model: could not allocate scratch wind\n");
    Exit (EXIT_FAILURE);
  }
  import_make_grid (ndom, scratch);

  /* import_make_grid also resets the coordinate arrays in zdom from the
   * model; recalculate them from wmain, as is done after the grid is
   * made in define_wind */
  wind_complete ();

  nstart = zdom[ndom].nstart;
  nstop = zdom[ndom].nstop;
  nbad_x = nbad_inwind = 0;

  for (n = nstart; n < nstop; n++)
  {
    vsub (scratch[n].x, wmain[n].x, dx);
    if (length (dx) > 1e-6 * length (wmain[n].x))
    {
      if (nbad_x < 10)
        Error ("update_wind_from_model: cell %d is at %e %e %e in the model but %e %e %e in the wind\n",
               n, scratch[n].x[0], scratch[n].x[1], scratch[n].x[2], wmain[n].x[0], wmain[n].x[1], wmain[n].x[2]);
      nbad_x++;
    }
    if (scratch[n].inwind != wmain[n].inwind)
    {
      if (nbad_inwind < 10)
        Error ("update_wind_from_model: cell %d has inwind %d in the model but %d in the wind\n", n, scratch[n].inwind, wmain[n].inwind);
      nbad_inwind++;
    }
  }

  if (nbad_x > 0 || nbad_inwind > 0)
  {
    Error ("update_wind_from_model: model %s does not match the wind: %d cells with different positions, %d with different inwind\n",
           filename, nbad_x, nbad_inwind);
    Exit (EXIT_FAILURE);
  }

  /* Update all the velocities first, since the velocities of imported
   * models are interpolated from those in wmain */

  for (n = nstart; n < nstop; n++)
  {
    stuff_v (scratch[n].v, wmain[n].v);
  }

  /* Now update everything that depends on the velocities.  The volume
   * includes the gamma factor at the cell centre, so it is rescaled
   * when that changes */

  for (n = nstart; n < nstop; n++)
  {
    set_cell_velocity_gradient (&wmain[n]);
  }

  /* dv/ds is found from the velocity gradients of neighbouring cells, so
   * it is calculated once all the gradients have been updated */

  for (n = nstart; n < nstop; n++)
  {
    cell = &wmain[n];
    old_xgamma_cen = cell->xgamma_cen;

    set_cell_dvds_and_gamma (cell);

    if (rel_mode == REL_MODE_FULL && old_xgamma_cen > 0.0)
    {
      cell->vol *= cell->xgamma_cen / old_xgamma_cen;
    }
  }

  free (scratch);

  /* Update the plasma cells, following create_plasma_grid */

  for (n = nstart; n < nstop; n++)
  {
    cell = &wmain[n];
    if (cell->vol <= 0.0)
    {
      continue;
    }

    nplasma = cell->nplasma;
    xplasma = &plasmamain[nplasma];

    xplasma->state.xgamma = cell->xgamma_cen;
    xplasma->state.vol = cell->vol * fill;

    old_rho = xplasma->state.rho;
    new_rho = model_rho (ndom, cell->xcen) / (fill * xplasma->state.xgamma);

    if (old_rho <= 0.0)
    {
      Error ("update_wind_from_model: plasma cell %d has density %e in the wind\n", nplasma, old_rho);
      Exit (EXIT_FAILURE);
    }

    /* Keep the ion fractions, but scale the ion densities to the new density */
    for (nion = 0; nion < nions; nion++)
    {
      xplasma->state.density[nion] *= new_rho / old_rho;
    }
    xplasma->state.rho = new_rho;
    xplasma->state.ne = get_ne (xplasma->state.density);

    if (imported_model[ndom].init_temperature == FALSE)
    {
      xplasma->state.t_r = import_temperature (ndom, cell->xcen, FALSE);
      xplasma->state.t_e = import_temperature (ndom, cell->xcen, TRUE);

      /* As in hydro_restart, set the level populations to those of the
       * ground state for the new temperature */
      partition_functions (xplasma, NEBULARMODE_LTE_GROUND);
    }
  }

  Log ("update_wind_from_model: updated domain %d from %s\n", ndom, filename);

  return (0);
}
