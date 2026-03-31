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
 * Also sets nstart_2 and nstop_2 in each domain struct.
 *
 **********************************************************/

int
make_transport_grid (void)
{
  int ndom, k, k_lower;

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    zdom[ndom].nstart_2 = NDIM2 + zdom[ndom].nstart;
    zdom[ndom].nstop_2 = NDIM2 + zdom[ndom].nstop;

    for (k = zdom[ndom].nstart; k < zdom[ndom].nstop; k++)
    {
      k_lower = k + NDIM2;

      wmain[k_lower] = wmain[k];        /* copy all fields from upper cell */

      wmain[k_lower].x[2] *= -1.0;     /* reflect z position */
      wmain[k_lower].xcen[2] *= -1.0;  /* reflect z cell centre */
      wmain[k_lower].nwind = k_lower;  /* correct self-index */
      wmain[k_lower].nwind_dom = k_lower;   /* correct domain-local index */

      /* nplasma is already correct from the copy: both hemispheres share
       * the same plasma cell, so no change needed. */
    }

    Log ("make_transport_grid: domain %d lower-hemisphere cells %d to %d\n",
         ndom, zdom[ndom].nstart_2, zdom[ndom].nstop_2 - 1);
  }

  return 0;
}
