
/***********************************************************/
/** @file  gradv.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  Routines related to calculating velocity gradients
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief      find the gradient velocity vector v, dv_ds for a photon
 * at a certain position travelling in a certain direction in the cmf
 * frame.
 *
 * @param [in] PhotPtr  p   A photon
 * @return     dvds    on successful completion
 * the value of where in grid if the photon outside the wind grid.
 *
 * @details
 *
 * ### Notes ###
 * For spherical coordinates the routine calculates the gradient
 * on the fly, that is moves the photon a bit and calculates dvds
 *
 * For 2d systems, the velocity gradient is calculated using
 * the velocity gradient tensors, which contain the velocity
 * gradient in the xz plane.
 *
 * It's not clear how much faster the method used for 2d systems
 * actually is, and one might want to consider the on-the-
 * fly appoach for all systems.
 *
 * Note that as calculated here dvds can be negative, although
 * in many cases what one wants is the absolute value of dvds.
 * But this has to be calculated at the spot where dvds is used,
 * since oftend dvds is calculated at two positions and interperlated.
 * So it would be a mistake to calculate the absolute value here.
 *
 **********************************************************/

double
dvwind_ds_cmf (p)
     PhotPtr p;
{
  double v_grad[3][3];
  double lmn[3], dvel_ds[3], dvds;
  int j, k, nn;
  struct photon pp;
  int nnn[4], nelem;
  double frac[4];
  double x;
  int ndom;

  ndom = wmain[p->grid].ndom;

  /* We want the change in velocity along the line of sight, but we
     need to be careful because of the fact that we have elected to
     combine the upper and lower hemispheres in the wind array.  Since
     we are only concerned with the the scalar dv_ds, the safest thing
     to do is to create a new photon that is only in the upper hemisphere
     02jan ksl */

  stuff_phot (p, &pp);
  if (pp.x[2] < 0.0)
  {                             /*move the photon to the northern hemisphere */
    pp.x[2] = -pp.x[2];
    pp.lmn[2] = -pp.lmn[2];
  }

  /* JM 1411 -- ideally, we want to do an interpolation on v_grad here. However,
     the interpolation was incorrect in spherical coordinates (see issue #118).
     For the moment, I've adopted an on the fly method for spherical coordinates.
     This should be improved by figuring out out how the velocity gradient
     tensor ought to be rotated in order to give the right answer for spherical
     coordinates. */

  if (zdom[ndom].coord_type == SPHERICAL || USE_GRADIENTS == FALSE)
  {
    struct photon pnew;
    double v1[3], v2[3], dv[3], diff[3];
    double ds;
    /* choose a small distance which is dependent on the cell size */
    vsub (pp.x, wmain[pp.grid].x, diff);
    vsub (wmain[pp.grid].xcen, wmain[pp.grid].x, diff);
    ds = 0.000001 * length (diff);
    /* calculate the velocity at the position of the photon */
    /* note we use model velocity, which could potentially be slow,
       but avoids interpolating (see #118) */
    model_velocity (ndom, pp.x, v1);

    /* copy the photon and move it by ds, and evaluate the velocity
       at the new point */
    stuff_phot (&pp, &pnew);
    /*
     * Put the photon into the observer frame, this way move_phot won't throw
     * an error. Should be ok since we are only using move_photon to move a
     * photon some vector ds.
     * */
    pnew.frame = F_OBSERVER;
    move_phot (&pnew, ds);
    model_velocity (ndom, pnew.x, v2);

    /* calculate the relevant gradient */
    if (rel_mode == REL_MODE_FULL)
    {
      observer_to_local_frame_velocity (v2, v1, dv);
      dvds = length (dv) / ds;
    }
    else
    {
      dvds = fabs (dot (v1, pp.lmn) - dot (v2, pp.lmn)) / ds;
    }
  }

  else                          // for non spherical coords we interpolate on v_grad
  {

    coord_fraction (ndom, 0, pp.x, nnn, frac, &nelem);


    for (j = 0; j < 3; j++)
    {
      for (k = 0; k < 3; k++)
      {
        x = 0;
        for (nn = 0; nn < nelem; nn++)
          x += wmain[nnn[nn]].v_grad[j][k] * frac[nn];

        v_grad[j][k] = x;

      }
    }

    /* v_grad is in cylindrical coordinates, or more precisely intended
       to be azimuthally symmetric.  One could either
       (a) rotate  v_grad to be correct at the position of the photon or
       (b) rotate the direction of photon travel so that is is correct
       (assuming azimuthal symmetry) in the xz plane.

       Possibility b is more straightforward and that is what is done
     */

    project_from_xyz_cyl (pp.x, pp.lmn, lmn);

    dvds = dot_tensor_vec (v_grad, lmn, dvel_ds);

    /* Note that the vector dvel_ds is also in an azimuthally symmetric system in the
     * xx plane, and could be rotated back if it were needed
     */
  }

  if (sane_check (dvds))
  {
    Error ("dvwind_ds_cmf: sane_check %f\n", dvds);
  }

  return (dvds);

}



#define N_DVDS_AVE	10000

/**********************************************************/
/**
 * @brief      Calculate the direction averaged (and maximum)
 * dv_ds in each grid cell of the wind
 *
 * @param [in]  int ndom  the domain of the cell
 * @param [in,out] WindPtr cell the cell to find dvds_ave for
 *
 * @return     Always returns 0
 *
 * @details
 * The routine calculates the average value of |dv_ds| at the center of the wind
 * cell by randomly generating directions and then calculating dv_ds in these
 * directions. If the RNG is uniform, then this should be accurate.
 *
 * ### Notes ###
 * The routine is called during the initialisation process and fills
 * the following elements of wmain
 *
 *  * calculate_cell_dvds_ave - the average the absolute value of dvds
 *
 * There is an advanced mode which prints this information to
 * file.
 *
 **********************************************************/


int
calculate_cell_dvds_ave (int ndom, WindPtr cell)
{
  int n;
  int icell;
  static int file_init = FALSE;
  struct photon p, pp;
  double v_zero[3], delta[3], vdelta[3], diff[3];
  double dvds_sum, dvds, ds;
  double dvds_max, lmn[3];
  double dvds_min, lmn_min[3];
  char filename[LINELENGTH];

  icell = cell->nwind;
  dvds_max = 0.0;
  dvds_min = 1.e30;

  /* Find the center of the cell */
  stuff_v (cell->xcen, p.x);

  /* Define a small length */
  vsub (p.x, cell->x, diff);
  ds = 0.001 * length (diff);

  /* Find the velocity at the center of the cell */
  model_velocity (ndom, p.x, v_zero);

  dvds_sum = 0.0;
  for (n = 0; n < N_DVDS_AVE; n++)
  {
    randvec (delta, ds);
    if (p.x[2] + delta[2] < 0)
    {                           // Then the new position would punch through the disk
      delta[0] = (-delta[0]);   // so we reverse the direction of the vector
      delta[1] = (-delta[1]);
      delta[2] = (-delta[2]);
    }
    vadd (p.x, delta, pp.x);
    model_velocity (ndom, pp.x, vdelta);
    vsub (vdelta, v_zero, diff);
    dvds = length (diff);

    /* Find the maximum and minimum values of dvds and the direction
     * for this. This is only required for additional diagnostics */
    if (modes.print_dvds_info)
    {
      if (dvds > dvds_max)
      {
        dvds_max = dvds;
        renorm (delta, 1.0);
        stuff_v (delta, lmn);
      }
      if (dvds < dvds_min)
      {
        dvds_min = dvds;
        renorm (delta, 1.0);
        stuff_v (delta, lmn_min);
      }
    }

    dvds_sum += dvds;
  }

  /* Store the result for the cell */
  cell->dvds_ave = dvds_sum / (N_DVDS_AVE * ds);

  if (modes.print_dvds_info)
  {
    sprintf (filename, "%.100s.dvds.diag", files.root);
    optr = fopen (filename, file_init == FALSE ? "w" : "a");
    file_init = TRUE;
    fprintf (optr,
             "%d %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e \n",
             icell, p.x[0], p.x[1], p.x[2], dvds_max / ds,
             dvds_min / ds, lmn[0], lmn[1], lmn[2], lmn_min[0], lmn_min[1], lmn_min[2], dot (lmn, lmn_min));
    fclose (optr);
  }

  return (0);
}

/**********************************************************/
/**
 * @brief      Calculate the maximum
 * dv_ds in each grid cell of the wind
 *
 * @param [in]  int ndom  the domain of the cell
 * @param [in,out] WindPtr cell the cell to find dvds_max for
 *
 * @return     Always returns 0
 *
 * @details
 * The routine calculates the maximum value of the absolute value of dv_ds at
 * the corner of the wind cell by randomly generating directions and then
 * calculating dv_ds in these directions. If the RNG is uniform, then this
 * should be accurate.
 * *
 * ### Notes ###
 * The routine is called during the initialisation process and fills
 * the following elements of wmain
 *
 *  * dvds_max - the maximum value of |dvds|
 *
 * Unlike calculate_cell_dvds_ave, these values are at the corners of cells,
 * and are intended to be interpolated.  
 *
 * There is an advanced mode which prints this information to
 * file.
 * 
 * See #888 for a discussion of why the absolute value is used
 *
 **********************************************************/

int
calculate_cell_dvds_max (int ndom, WindPtr cell)
{
  int n;
  int icell;
  static int file_init = FALSE;
  struct photon p;
  double dvds;
  double dvds_max, lmn[3];
  double dvds_min, lmn_min[3];
  char filename[LINELENGTH];

  (void) ndom;                  /* ndom is an argument to maintain consistency with other functions */
  icell = cell->nwind;
  dvds_max = 0.0;
  dvds_min = 1.e30;

  /*  x is the corner of the cell */
  stuff_v (cell->x, p.x);

  /* Cannot calculate the velocity gradient along the z axis so fudge this */
  if (p.x[0] == 0)
  {
    p.x[0] = 0.1 * cell->xcen[0];
  }

  p.grid = icell;

  for (n = 0; n < N_DVDS_AVE; n++)
  {
    randvec (p.lmn, 1);
    dvds = fabs (dvwind_ds_cmf (&p));

    /* Find the maximum and minimum values of dvds and the direction
     * for this
     */
    if (dvds > dvds_max)
    {
      dvds_max = dvds;
      if (modes.print_dvds_info)
      {
        stuff_v (p.lmn, lmn);
      }
    }
    if (dvds < dvds_min)
    {
      dvds_min = dvds;
      if (modes.print_dvds_info)
      {
        stuff_v (p.lmn, lmn_min);
      }
    }
  }

  /* Store the result for the cell */
  cell->dvds_max = dvds_max;

  if (modes.print_dvds_info)
  {
    sprintf (filename, "%.100s.dvds.diag", files.root);
    optr = fopen (filename, file_init == FALSE ? "w" : "a");
    file_init = TRUE;
    fprintf (optr,
             "%d %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e \n",
             icell, p.x[0], p.x[1], p.x[2], dvds_max,
             dvds_min, lmn[0], lmn[1], lmn[2], lmn_min[0], lmn_min[1], lmn_min[2], dot (lmn, lmn_min));
    fclose (optr);
  }

  return (0);
}




/**********************************************************/
/**
 * @brief      Calculate the maximum
 * dv_ds at a particular position in a grid
 *
 * @param [in] PhotPtr  p   A photon

 * @return     Returns dvds_max at the position of the photon
 *
 * @details
 * The routine interpolates dvds_max given the position of
 * a photon in a cell
 *
 * dvds_max at the vertex points of cells must have been
 * initialized using the routine dvds_max
 *
 * ### Notes ###
 *
 * The routine uses both the position of the photon and 
 * the grid cell in which the photon exists, so this
 * must be acurrate.
 **********************************************************/


double
get_dvds_max (p)
     PhotPtr p;
{
  int ndom, nn, nnn[4], nelem;
  double frac[4];
  double dvds;

  ndom = wmain[p->grid].ndom;

  coord_fraction (ndom, 0, p->x, nnn, frac, &nelem);

  dvds = 0;

  for (nn = 0; nn < nelem; nn++)
  {
    dvds += wmain[nnn[nn]].dvds_max;
  }

  return dvds;

}
