
/***********************************************************/
/** @file  photon2d.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  This file contains the routines which translate
 * photons, both inside and outside of the domains.
 *
 * Except for translating the photon the routines do as
 * litle as possible w.r.t modifying any of the other
 * information in the photon structure.
 *
 * There is a steering routine translate which calle
 * either translate_in_wind or translate_in_space, depending
 * on whether a photon bundle is in a wind domain or in space.
 * The routines translate the photon until hits something,
 * either the edge of the wind, or the edge of a cell, or
 * a resonance and then return with a status noting what
 * has happened.
 *
 * ### Notes ###
 *
 * Some of these routines pass the WindPtr for the entire
 * wind which should not be necessary and had been removed
 * in favor of using wmain.  It would probably be a good idea
 * to that here.
 *
 * These routines are called from trans_phot.  One could argue
 * that the filename is misnamed for that reason.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief      a steering routine that either calls translate_in_space or translate_in_wind  depending upon the
 * 	current location of the photon.
 *
 * @param [in] WindPtr  w   A pointer to the wind
 * @param [in, out] PhotPtr  pp   A photon
 * @param [in] double  tau_scat   the depth at which the photon should scatter
 * @param [ out] double *  tau   The optical depth associated with a resonanant scatter when
 * this is the reason the photon has stopped
 * @param [in, out] int *  nres   The number of the resonance if the photon stopped.
 * due to reaching the scattering optical depth
 * @return     A status that states what caused the photon to stp as it did
 *
 * @details
 *
 * The routine takes the current position of the photon and determines
 * where it is.  Depending on this it calls either translate_in_space
 * or translate_in_wind
 *
 * On exit, the position of the photon will have been updated.
 *
 * ### Notes ###
 * Translate controls the flow of the photon through one grid cell.  In Python,
 * there are additional possibilities since there are wind free regions.
 *
 * The routine first calls where_in_wind which returns the correct domain for a given position, or domain 0 if the
 * position is not in the wind of a domain.  The photon structure does not have the domain
 * number directly encoded, but it can be obtained from the grid number, which where in_grid updates
 *
 *
 * This construction is
 * used so that the calling routines for translate (trans_phot and extract) can be the same
 * in the 1 and 2d versions.
 *
 **********************************************************/

int
translate (w, pp, tau_scat, tau, nres)
     WindPtr w;                 //w here refers to entire wind, not a single element
     PhotPtr pp;
     double tau_scat;
     double *tau;
     int *nres;
{
  int istat;
  int ndomain;
  int ichoose;


  ichoose = where_in_wind (pp->x, &ndomain);


  if (ichoose < 0)
  {
    istat = translate_in_space (pp);

  }
  else if ((pp->grid = where_in_grid (ndomain, pp->x)) >= 0)
  {

    istat = translate_in_wind (w, pp, tau_scat, tau, nres);
  }
  else
  {
    istat = pp->istat = P_ERROR;        /* It's not in the wind and it's not in the grid.  Bummer! */
    Error ("translate: Found photon that was not in wind or grid, istat %i\n", where_in_wind (pp->x, &ndomain));
  }

  return (istat);
}




/**********************************************************/
/**
 * @brief      translates the photon from its current position to the
 * 	edge of the wind.
 *
 * @param [in, out] PhotPtr  pp   A photon
 * @return     A status flag indication why the photon stopped
 *
 * @details
 * The routine translate the photon (which is assumed to be
 * in the wind, but not in an active wind cell) to the closest
 * boundary.
 *
 * The photon position is updated in the process.
 *
 *
 * ### Notes ###
 * @bug  There are questions in the comments about why an additonal
 * check is needed as to whther the photon as hit the star. It seems
 * superfluous so someone should check whether this addtional check
 * can be removed.
 *
 **********************************************************/

int
translate_in_space (pp)
     PhotPtr pp;
{
  double ds, delta, s, smax, prhosq;
  int ndom, ndom_next;
  struct photon ptest;

  ds = ds_to_wind (pp, &ndom);

  /* For IMPORT, although we have reached the edge of the wind, we may be in a cell that is
   * not really in the wind, so we have to address this situtation here.  The first problem
   * we have though is that we do not know what domain we are in.*/

  if (ndom >= 0 && (zdom[ndom].wind_type == IMPORT))
  {
    stuff_phot (pp, &ptest);
    move_phot (&ptest, ds + DFUDGE);    /* So now ptest is at the edge of the wind as defined by the boundary
                                           From here on we should be in the grid  */
    ds += DFUDGE;               //Fix for Bug #592 - we need to keep track of the little DFUDGE we moved the test photon


    /* Note there is a possibility that we reach the other side 
     * of the grid without actually encountering a
     * wind cell
     */

    prhosq = (pp->x[0] * pp->x[0]) + (pp->x[1] * pp->x[1]);

    if ((prhosq > (zdom[ndom].wind_rhomin_at_disk * zdom[ndom].wind_rhomin_at_disk)) &&
        (prhosq < (zdom[ndom].wind_rhomax_at_disk * zdom[ndom].wind_rhomax_at_disk)))
    {
      stuff_phot (pp, &ptest);
      ds = 0.0;
    }


    if (where_in_wind (ptest.x, &ndom_next) < 0)
    {
      smax = ds_to_wind (&ptest, &ndom_next);   // This is the maximum distance can go in this domain
      s = 0;
      while (s < smax && where_in_wind (ptest.x, &ndom_next) < 0)
      {
        if ((delta = ds_in_cell (ndom, &ptest)) > 0)
        {
          move_phot (&ptest, delta + DFUDGE);
          s += delta + DFUDGE;  // The distance the photon has moved
        }
        else
        {
          break;
        }
      }
      /* So at this point we have either gotten out of the domain or we have found a cell that
       * is actually in the wind or we encoutered the error above
       */
      if (s > 0)
      {
        ds += s - DFUDGE;       /* We are going to add DFUDGE back later */
      }
    }

  }

  move_phot (pp, ds + DFUDGE);

  return (pp->istat);
}




/**********************************************************/
/**
 * @brief      calculates the photon pathlength to the edge of the wind.
 *
 * @param [in, out] PhotPtr  pp   A photon
 * @param [in, out] int *  ndom_current   The current domain
 * @return     The distance to the nearest boundary of the wind and the domain for which
 * 	the boudary applies.
 *
 * @details
 * Python defines the boundaries of the wind in term of the intersection
 * of a biconical flow (defined by an inner and outer windcone) and
 * an inner and outer radius centered on the star.	 If the user is interested
 * in a biconical flow, he/she will most likely set the radii so that they
 * do not truncate the wind.  Similarly, by choosing parameters for the windcones
 * properly, one can assure that one has a completely spherical flow.  However,
 * the user may also perversely set the parameters to that both the conical
 * flow boundaries and the min and maximum radius come into play.
 *
 * Usually, one would
 * efine the two types of The two types of boundaries are
 * sually defined so that one includes the other and so the cones apply when one is
 * ealing with a CV type wind, while the inner and outer radius applies for a spherical
 * odel.  However this routine does not require this to be the case, since it just
 * alculates where the edges are.
 *
 * n any event, if you are inside the wind already ds_to_wind calculates the distance to the edge of the wind.
 * If you are outside, It will also be to the nearest edge.
 *
 * The routine distinguishes between  two basic cases.  Either the photon is already in the wind
 * in which case one looks for the nearest boundary, or the photon is not currently
 * in the wind and you have to translate it until it reaches the the boundary (or
 * VERY_BIG)
 *
 * ### Notes ###
 * There is no guarantee that you will still be in the region defined by the grid.
 *
 * @bug 1802 -ksl - At present this routine for imported models this routine only deals with
 * cylindrical models.  Additionally for imported models we skip all of the
 * of the wind_cones.  This is inefficient, and needs to be corrected for
 * rtheta and spherical models which can easily be handled using wind cones.
 *
 * ksl 2201 - It is unclear what the comment I made in 1802 actually means anymore.
 *
 **********************************************************/

double
ds_to_wind (pp, ndom_current)
     PhotPtr pp;
     int *ndom_current;
{
  struct photon ptest, qtest;
  double ds, x, rho, z;
  int ndom;

  stuff_phot (pp, &ptest);

  /* First calculated the distance to the edge of the of
     all of the "computatational domain */

  ds = ds_to_sphere (geo.rmax, &ptest);
  *ndom_current = (-1);
  xxxbound = BOUND_NONE;

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    if ((zdom[ndom].wind_type != IMPORT || (zdom[ndom].wind_type == IMPORT && zdom[ndom].coord_type != CYLIND)) &&
        zdom[ndom].wind_type != CORONA)
    {
      /* Check if the photon hits the inner or outer radius of the wind */
      if ((x = ds_to_sphere (zdom[ndom].rmax, &ptest)) < ds)
      {
        ds = x;
        *ndom_current = ndom;
        xxxbound = BOUND_RMIN;
      }

      if ((x = ds_to_sphere (zdom[ndom].rmin, &ptest)) < ds)
      {
        ds = x;
        *ndom_current = ndom;
        xxxbound = BOUND_RMAX;
      }

      /* Check if the photon hits the inner or outer windcone */

      if ((x = ds_to_cone (&zdom[ndom].windcone[0], &ptest)) < ds)
      {
        ds = x;
        *ndom_current = ndom;
        xxxbound = BOUND_INNER_CONE;
      }
      if ((x = ds_to_cone (&zdom[ndom].windcone[1], &ptest)) < ds)
      {
        ds = x;
        *ndom_current = ndom;
        xxxbound = BOUND_OUTER_CONE;
      }
    }

    /* For this rectangular region we check whether we are in side the grid,
     * which should effectively.  For * an imported region file we may not be
     * inside the wind, since some cells may be empty
     */

    else if (zdom[ndom].wind_type == CORONA || (zdom[ndom].wind_type == IMPORT && zdom[ndom].coord_type == CYLIND))
    {
      x = ds_to_plane (&zdom[ndom].windplane[0], &ptest, TRUE);
      if (x > 0 && x < ds)
      {
        stuff_phot (pp, &qtest);
        move_phot (&qtest, x);
        rho = sqrt (qtest.x[0] * qtest.x[0] + qtest.x[1] * qtest.x[1]);
        if (zdom[ndom].wind_rhomin_at_disk <= rho && rho <= zdom[ndom].wind_rhomax_at_disk)
        {
          ds = x;
          *ndom_current = ndom;
          xxxbound = BOUND_ZMIN;
        }
      }
      x = ds_to_plane (&zdom[ndom].windplane[1], &ptest, TRUE);
      if (x > 0 && x < ds)
      {
        stuff_phot (pp, &qtest);
        move_phot (&qtest, x);
        rho = sqrt (qtest.x[0] * qtest.x[0] + qtest.x[1] * qtest.x[1]);
        if (zdom[ndom].wind_rhomin_at_disk <= rho && rho <= zdom[ndom].wind_rhomax_at_disk)
        {
          ds = x;
          *ndom_current = ndom;
          xxxbound = BOUND_ZMAX;
        }
      }

      x = ds_to_cylinder (zdom[ndom].wind_rhomin_at_disk, &ptest);
      if (x > 0 && x < ds)
      {
        stuff_phot (pp, &qtest);
        move_phot (&qtest, x);
        z = fabs (qtest.x[2]);
        if (zdom[ndom].zmin <= z && z <= zdom[ndom].zmax)
        {
          ds = x;
          *ndom_current = ndom;
          xxxbound = BOUND_INNER_RHO;
        }
      }

      x = ds_to_cylinder (zdom[ndom].wind_rhomax_at_disk, &ptest);
      if (x > 0 && x < ds)
      {
        stuff_phot (pp, &qtest);
        move_phot (&qtest, x);
        z = fabs (qtest.x[2]);
        if (zdom[ndom].zmin <= z && z <= zdom[ndom].zmax)
        {
          ds = x;
          *ndom_current = ndom;
          xxxbound = BOUND_OUTER_RHO;
        }
      }
    }
    else
    {
      Error ("ds_to_wind:Do not know how to deal with this combination of coordinate type %d and wind_type %d\n", zdom[ndom].coord_type,
             zdom[ndom].wind_type);
      Exit (0);
    }
  }

  return (ds);
}

/**********************************************************/
/**
 * @brief      translates the photon within a single cell in the wind.
 *
 * @param [in] WindPtr  w   The entire wind
 * @param [in, out] PhotPtr  p   A photon
 * @param [in] double  tau_scat   The depth at which the photon will scatter
 * @param [out] double *  tau   The tau of a resonance
 * @param [out] int *  nres   The resonance which caused the photon to stop
 * @return     A status indicated whether the photon has stopped for a scattering
 * even of for some other reason
 *
 * @details
 * It calculates and updates the final position of the photon p, the optical depth tau after
 * having translated the photon, and, if the photon is scattered, the number of the resonance
 * responsible for scattering (or absorbing) the photon bundle.
 *
 * In the simple atom case, the weight of the photon is also reduced as a result
 * of continuum absorption, so that what is returned final weight of the photon
 * in the observer frame..
 *
 * ### Notes ###
 *
 * In addition to moving the photon, the routine also updates some values
 * having to do with the radiation field via calls to radiation (simple atoms)
 * or update_bf_estimators (macro-atoms)
 *
 **********************************************************/
int
translate_in_wind (w, p, tau_scat, tau, nres)
     WindPtr w;                 //w here refers to entire wind, not a single element
     PhotPtr p;
     double tau_scat, *tau;
     int *nres;
{
  int n;
  double smax, ds_current, ds_cmf;
  int istat;
  int nplasma;

  WindPtr one;
  PlasmaPtr xplasma;
  struct photon phot_mid, phot_mid_cmf; // Photon at the midpt of its path in the cell

  /* First verify that the photon is in the grid, and if not
     return and record an error */

  if ((p->grid = n = where_in_grid (wmain[p->grid].ndom, p->x)) < 0)
  {
    return (n);                 /* Photon was not in grid */
  }
  /* Assign the pointers for the cell containing the photon */

  one = &wmain[n];              /* one is the grid cell where the photon is */
  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  /* Calculate the maximum distance the photon can travel in the cell */

  smax = smax_in_cell (p);

  /* We now determine whether scattering prevents the photon from reaching the far edge of
     the cell.  calculate_ds calculates whether there are scatterings and makes use of the
     current position of the photon and the position of the photon at the far edge of the
     shell.  It needs a "trial photon at the maximum distance however...
     Note that ds_current does not alter p in any way */


  if ((modes.partial_cells == PC_EXTEND && one->inwind == W_PART_INWIND) || one->inwind == W_IGNORE)
  {
    ds_current = smax;
    istat = 0;
    if (smax < 0)
    {
      Error ("Houston, there is a problem %e\n", smax);
    }
  }
  else
  {
    ds_current = calculate_ds (w, p, tau_scat, tau, nres, smax, &istat);

    if (p->nres == NRES_ES)
      xplasma->nscat_es++;
    else if (p->nres > 0)
      xplasma->nscat_res++;

    /* We now increment the radiation field in the cell, translate the photon and wrap
     * things up.  For simple atoms, the routine radiation also reduces
     * the weight of the photon due to continuum absorption, e.g. free free.
     */

    if (geo.rt_mode == RT_MODE_MACRO)
    {
      /* In the macro-method, b-f and other continuum processes do not reduce the photon
         weight, but are treated as as scattering processes.  Therefore most of what was in
         subroutine radiation for the simple atom case can be avoided.  
       */
      if (geo.ioniz_or_extract == CYCLE_IONIZ)
      {
        /* Provide inputs to bf_estimators in the local frame  */
        stuff_phot (p, &phot_mid);
        move_phot (&phot_mid, 0.5 * ds_current);
        observer_to_local_frame (&phot_mid, &phot_mid_cmf);
        ds_cmf = observer_to_local_frame_ds (&phot_mid, ds_current);
        if (p->grid >= 0 && p->grid < geo.ndim2)
        {
          bf_estimators_increment (&w[p->grid], &phot_mid_cmf, ds_cmf);
        }
        else
        {
          Error ("translate_in_wind: Cannot call bf_estimators for photon not in wind: %d\n", p->np);
        }

      }
    }
    else
    {
      radiation (p, ds_current);
    }

    if (*nres > -1 && *nres <= NLINES && *nres == p->nres && istat == P_SCAT)
    {
      if (ds_current < wmain[p->grid].dfudge)
      {
        Error
          ("translate_in_wind: found repeated resonance scattering nres %5d after motion of %10.3e for photon %d in plasma cell %d)\n",
           *nres, ds_current, p->np, wmain[p->grid].nplasma);
      }
    }

    p->nres = *nres;
    if (p->nres > -1 && p->nres < nlines)
      p->line_res = p->nres;
  }

  p->istat = istat;



  move_phot (p, ds_current);
  return (p->istat);
}




/* ************************************************************************* */
/**
 * @brief           Calculate the maximum distance a photon can travel in its
 *                  current cell.
 *
 * @param[in]       PhotPtr   p       The currently transported photon packet
 *
 * @return          double    smax    The maximum distance the photon packet
 *                                    can move
 *
 * @details
 *
 * smax can be found in various ways depending on if the photon is in the wind
 * or not.
 *
 * This routine operates in the global/observer frame.
 *
 * This section of code was put into its own function in an attempt to avoid
 * code duplication in the optical depth diagnostic ray tracing routine.
 *
 * ************************************************************************** */

double
smax_in_cell (PhotPtr p)
{
  int n, ndom, ndom_current;
  double s, smax;
  int hit_disk;

  WindPtr one;

  n = p->grid;
  one = &wmain[n];              /* one is the grid cell where the photon is */
  ndom = one->ndom;

  /* Calculate the maximum distance the photon can travel in the cell */

  if ((smax = ds_in_cell (ndom, p)) < 0)
  {
    return ((int) smax);
  }
  if (one->inwind == W_PART_INWIND)
  {                             /* The cell is partially in the wind */
    s = ds_to_wind (p, &ndom_current);  /* smax is set to be the distance to edge of the wind */
    if (s < smax)
      smax = s;
    s = ds_to_disk (p, 0, &hit_disk);   /* the 0 implies ds_to_disk can not return a negative distance */
    if (s > 0 && s < smax)
      smax = s;
  }
  else if (one->inwind == W_IGNORE)
  {
    smax += one->dfudge;
    //move_phot (p, smax);
    return (smax);
  }
  else if (one->inwind == W_NOT_INWIND)
  {                             /* The cell is not in the wind at all */

    Error ("translate_in_wind: Grid cell %d of photon is not in wind, moving photon %.2e\n", n, smax);
    Error ("translate_in_wind: photon %d position: x %g y %g z %g\n", p->np, p->x[0], p->x[1], p->x[2]);
    move_phot (p, smax);
    return (smax);

  }

  /* At this point we now know how far the photon can travel in it's current grid cell */

  smax += one->dfudge;          /* dfudge is to force the photon through the cell boundaries. */

  /* Set limits the distance a photon can travel.  There are
     a good many photons which travel more than this distance without this
     limitation, at least in the standard 30 x 30 instantiation.  It does
     make small differences in the structure of lines in some cases.
     The choice of SMAX_FRAC can affect execution time. */

  if (smax > SMAX_FRAC * length (p->x))
  {
    smax = SMAX_FRAC * length (p->x);
  }

  return smax;
}




/**********************************************************/
/**
 * @brief      calculates the distance photon can travel within the cell
 * 	that it is currently in.
 *
 * @param [in, out] int  ndom   The current domain
 * @param [in, out] PhotPtr  p   A photon
 * @return     A distance indicating how far the photon can travel
 *
 * @details
 * The routine is basically just a steering routine
 * and calls ds_in_whatever for various coordinate systems
 *
 * ### Notes ###
 *
 * The other routines are contained in routines like cylindrical.c,
 * rtheta.c etc.
 *
 **********************************************************/

double
ds_in_cell (ndom, p)
     int ndom;
     PhotPtr p;

{
  int n;
  double smax;

  /* First verify that the photon is in the grid, and if not
     return and record an error */

  if ((p->grid = n = where_in_grid (ndom, p->x)) < 0)
  {
    return (n);
  }

  if (zdom[ndom].coord_type == CYLIND)
  {
    smax = cylind_ds_in_cell (ndom, p); // maximum distance the photon can travel in a cell
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {
    smax = rtheta_ds_in_cell (ndom, p);
  }
  else if (zdom[ndom].coord_type == SPHERICAL)
  {
    smax = spherical_ds_in_cell (ndom, p);
  }
  else if (zdom[ndom].coord_type == CYLVAR)
  {
    smax = cylvar_ds_in_cell (ndom, p);
  }
  else
  {
    smax = 0;
    Error ("ds_in_cell: Don't know how to find ds_in_cell in this coord system %d\n", zdom[ndom].coord_type);
    Exit (0);
  }

  return (smax);
}
