/***********************************************************/
/** @file  transport_grid.c
 * @brief  Build and maintain the transport grid.
 *
 * @details
 * The transport grid extends wmain beyond the plasma (upper-hemisphere)
 * cells to cover additional geometric regions that share a single plasma
 * cell.  The many-to-one mapping is expressed through wmain[n].nplasma,
 * which points to the plasma cell that supplies physical properties for
 * photons travelling through transport cell n.
 *
 * The first extension implemented here is the lower-hemisphere mirror:
 * for each upper-hemisphere plasma cell k (in nstart..nstop-1), a
 * corresponding lower-hemisphere transport cell is created at k + NDIM2.
 * Its geometry is identical except that x[2] and xcen[2] are negated.
 * Its nplasma points to the same plasma cell as the upper cell, so
 * photons in either hemisphere deposit estimators into the same plasma.
 *
 * Per-domain indices:
 *   Upper hemisphere: nstart .. nstop-1
 *   Lower hemisphere: nstart_2 .. nstop_2-1  (= nstart+NDIM2 .. nstop+NDIM2-1)
 *
 * ### Notes ###
 * Call make_transport_grid() after create_plasma_grid() (so that
 * wmain[k].nplasma is already populated for k = 0..NDIM2-1).
 * wmain must have been allocated with calloc_wind(2 * NDIM2).
 *
 **********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief  Populate lower-hemisphere transport cells in wmain.
 *
 * @return  0 on success
 *
 * @details
 * For every upper-hemisphere wind cell k, creates a mirror cell at
 * k + NDIM2 with negated z coordinates and the same nplasma index.
 *
 * All wind-model-derived quantities (velocity, velocity gradient,
 * divergence, dvds, and relativistic gamma factors) are computed
 * directly from the wind model at the reflected position.  No
 * assumption of z-symmetry is made beyond what is inherent in the
 * wind model itself.  Only purely geometric quantities (cell
 * boundaries, cone intercepts) are set by reflection.
 *
 * Also sets nstart_2 and nstop_2 in each domain struct.
 *
 **********************************************************/

int
make_transport_grid (void)
{
  int ndom, k, k_lower;
  double v_cen[3];
  WindPtr cell;

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    zdom[ndom].nstart_2 = NDIM2 + zdom[ndom].nstart;
    zdom[ndom].nstop_2 = NDIM2 + zdom[ndom].nstop;

    for (k = zdom[ndom].nstart; k < zdom[ndom].nstop; k++)
    {
      k_lower = k + NDIM2;
      cell = &wmain[k_lower];

      wmain[k_lower] = wmain[k];        /* copy all fields from upper cell */

      cell->x[2] *= -1.0;       /* reflect z position */
      cell->xcen[2] *= -1.0;    /* reflect z cell centre */
      cell->nwind = k_lower;    /* correct self-index */
      cell->nwind_dom = k_lower;        /* correct domain-local index */

      /* nplasma is already correct from the copy: both hemispheres share
       * the same plasma cell, so no change needed. */

      /* Reflect the cone boundary: the lower-hemisphere cone is the mirror
       * image of the upper one, so its z-intercept changes sign. dzdr
       * is unchanged (same slope). */
      cell->wcone.z *= -1.0;

      /* Reflect the explicit z cell boundaries for CYLIND cells.
       * The lower-hemisphere cell occupies [-cell_z_max, -cell_z_min]. */
      cell->cell_z_min = -wmain[k].cell_z_max;
      cell->cell_z_max = -wmain[k].cell_z_min;

      /* Compute all wind-model-derived quantities directly at the
       * lower-hemisphere position.  This avoids any assumption of
       * z-symmetry and remains correct for genuinely asymmetric winds. */

      if (zdom[ndom].wind_type != IMPORT)
        model_velocity (ndom, cell->x, cell->v);

      model_vgrad (ndom, cell->x, cell->v_grad);
      wind_div_v (ndom, cell);
      calculate_cell_dvds_ave (ndom, cell);
      calculate_cell_dvds_max (ndom, cell);

      /* Relativistic gamma factors.  cell->vol was copied from the upper
       * cell, where it already has the upper xgamma_cen baked in.  Recover
       * the geometric volume by dividing it out, then apply the lower
       * hemisphere's own xgamma_cen.  In non-relativistic mode both
       * xgamma_cen values are 1.0, so the formula is correct in all cases. */
      if (rel_mode == REL_MODE_FULL)
      {
        cell->xgamma = calculate_gamma_factor (cell->v);
        model_velocity (ndom, cell->xcen, v_cen);
        cell->xgamma_cen = calculate_gamma_factor (v_cen);
        cell->vol = (wmain[k].vol / wmain[k].xgamma_cen) * cell->xgamma_cen;
      }
      else
      {
        cell->xgamma = 1.0;
        cell->xgamma_cen = 1.0;
        /* vol is unchanged: geometric volume is the same for both hemispheres */
      }
    }

    Log ("make_transport_grid: domain %d lower-hemisphere cells %d to %d\n", ndom, zdom[ndom].nstart_2, zdom[ndom].nstop_2 - 1);
  }

  return 0;
}
