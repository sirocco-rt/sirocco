
/***********************************************************/
/** @file  wind_util.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  These are very simple utilities within Sirocco that do not
 *   fall into a larger category like those in vector.
 *
 * ### Notes ###
 *
 * These routines all have to do with wind cells. 
 *
 * It may be useful to rethink routins at some point, as they were all
 * written prior to the existence of domains, and domains are
 * just grafted on.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "log.h"
#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief      calculate the fractional
 * 	contributions in interpolations of values of quantities,
 * 	such as electron density, by various cells in the
 * 	grid in a coordinate system independent way.
 *
 * @param [in] int  ndom   The domain in where the interpolation will take place
 * @param [in] int  ichoice  interpolate on vertices (0), interpolate on
 * centers (1)
 * @param [in] double  x[]   the 3-vector position for which you want
 * the fractinal position
 * @param [out] int  ii[]    an array that contains the 1-d element
 * numbers that must be summed.  These are absolute indices into wmain,
 * that is they include zdom[ndom].nstart, for every coordinate system
 * @param [out] double  frac[]   the array that contains the fractional
 * contribution of the corresponding
 * element
 * @param [out] int *  nelem   the number of elements that must be
 * summed, nominally 2 for a spherical grid
 * and 4 for a two dimensional grid.  (For
 * a 3 d grid it would be 6, but we have not
 * implemented this.
 *
 * @return
 * 	 1 if  inside the grid
 * 	-2 if outside the grid
 * 	-1 if inside the grid
 *
 *
 * @details
 *
 * ### Notes ###
 * There are numerous times when one wants the value
 * 	of an interpolated  variable in the wind.  There
 * 	is no easy way to interpolate the variable easily.
 * 	What this routine does is calculate the fractional
 * 	contributions of elements in the array to that
 * 	position.  Then one must sum up the actual values
 * 	elsewhere.
 *
 * 	If positions are outside the grid, coord_fraction
 * 	attempts to give you the value at the edge of the
 * 	grid.
 *
 * 	It's possible that coord_fraction could be used
 * 	to interpolate beyond the edge of the grid where
 * 	a variable is defined, although this is not done
 * 	at present.
 *
 * 	Coordinate conversion for each grid type:
 * 	- CYLIND:  r = rho (cylindrical radius), z = x[2] (signed; negative for
 * 	  lower hemisphere).  wind_z[] stores signed z values covering both
 * 	  hemispheres.
 * 	- RTHETA:  r = |x| (spherical radius), z = acos(x[2]/r)*RADIAN giving
 * 	  theta in [0, 180] deg over the full doubled grid.  No fabs is applied
 * 	  to x[2], so lower-hemisphere positions (x[2] < 0) give theta > 90 deg
 * 	  and are correctly mapped to lower-hemisphere cells.
 * 	- SPHERICAL: r = |x|, z unused (1d grid).
 *
 **********************************************************/

int
coord_fraction (int ndom, int ichoice, double x[], int ii[], double frac[], int *nelem)
{
  double r, z;
  double *xx, *zz;
  int ix, iz;
  double dr, dz;


  /* Assign pointers to the xx and zz depending on whether
   * one wants to interpolate on vertex points (0) or
   * midpoints (1)
   */
  if (ichoice == 0)
  {
    xx = zdom[ndom].wind_x;
    zz = zdom[ndom].wind_z;
  }
  else
  {
    xx = zdom[ndom].wind_midx;
    zz = zdom[ndom].wind_midz;
  }

  /* Now convert x to the appropriate coordinate system */
  r = z = 0.0;
  if (zdom[ndom].coord_type == CYLIND)
  {
    r = sqrt (x[0] * x[0] + x[1] * x[1]);
    /* Data-based hemisphere detection: if wind_z[0] >= 0 the grid only
       covers the upper hemisphere (imported single-hemisphere file), so fold
       lower-hemisphere positions via fabs.  If wind_z[0] < 0 the grid covers
       both hemispheres (parameterized doubled grid or two-hemisphere import),
       so use signed z directly. */
    if (zdom[ndom].wind_z[0] >= 0.0)
      z = fabs (x[2]);
    else
      z = x[2];                 /* signed: lower-hemisphere cells have wind_z[j] < 0 */
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {
    r = length (x);
    /* Data-based hemisphere detection: if wind_z[mdim-1] <= 90 the grid only
       covers the upper hemisphere (imported single-hemisphere file), so fold
       lower-hemisphere positions.  If wind_z[mdim-1] > 90 the grid covers
       both hemispheres (parameterized doubled grid or two-hemisphere import),
       so use signed theta over the full 0-180 deg range. */
    if (zdom[ndom].wind_z[zdom[ndom].mdim - 1] <= 90.0)
      z = acos (fabs (x[2]) / r) * RADIAN;
    else
      z = acos (x[2] / r) * RADIAN;     /* signed: gives 0-180 deg for full doubled grid */
  }
  else if (zdom[ndom].coord_type == SPHERICAL)
  {
    r = length (x);
    z = 0;                      // To avoid -O3 warning
  }
  else
  {
    Error ("coord_fraction: Unknown coordinate type %d for doman\n", zdom[ndom].coord_type, ndom);
    Exit (0);
  }

  if (zdom[ndom].coord_type == SPHERICAL)
  {                             /* We are dealing with a 1d system */
    fraction (r, xx, zdom[ndom].ndim, &ix, &dr, 0);     //linear space
    ii[0] = zdom[ndom].nstart + ix;
    frac[0] = (1. - dr);
    ii[1] = zdom[ndom].nstart + ix + 1;
    frac[1] = dr;
    *nelem = 2;
//    if (sane_check (dr))
//    {
//      Error ("coord_frac:sane_check dr=%f for spherical coords. \n", dr);
//    }
  }
  else
  {                             /* We are dealing with a 2d system */
    fraction (r, xx, zdom[ndom].ndim, &ix, &dr, 0);
    fraction (z, zz, zdom[ndom].mdim, &iz, &dz, 0);

    ii[0] = zdom[ndom].nstart + ix * zdom[ndom].mdim + iz;
    frac[0] = (1. - dz) * (1. - dr);

    ii[1] = zdom[ndom].nstart + (ix + 1) * zdom[ndom].mdim + iz;
    frac[1] = (1. - dz) * dr;

    ii[2] = zdom[ndom].nstart + ix * zdom[ndom].mdim + iz + 1;
    frac[2] = (dz) * (1. - dr);

    ii[3] = zdom[ndom].nstart + (ix + 1) * zdom[ndom].mdim + iz + 1;
    frac[3] = (dz) * (dr);
    *nelem = 4;

//    if (sane_check (dr) || sane_check (dz))
//    {
//      Error ("coord_frac:sane_check dr=%f dz=%f for 2d coords\n", dr, dz);
//    }
  }
  /* At this point i,j are just outside the x position */
  /* Check to see if x is outside the region of the calculation */
  /* Note that this is a very incoplethe check in the sneste that
   * the posision could be out of the grid in other directions */

  if (r > xx[zdom[ndom].ndim - 1])
  {
    return (-2);                /* x is outside grid */
  }
  else if (r < xx[0])
  {
    return (-1);                /*x is inside grid */
  }

  return (1);

}




/**********************************************************/
/**
 * @brief      Translate a cell number in the wind to the two-d grid
 * position in a specific domain
 *
 * @param [in] int  ndom   The domain of interest
 * @param [in] int  n   The 1-d position in the wind domain
 * @param [out] int *  i   The first index for a 2d position in a given domain
 * @param [out] int *  j   The second index for a 2d postion in a given domain
 * @return     i,j position in a given domain
 *
 * @details
 * wmain contains elements for the entire winds.  This routine translates
 * from a wmain index to the two-d index of a specific domain so that one
 * can interpolate.
 *
 * ### Notes ###
 * For 2d, the underlying 1-d grid is organized so that as n increases,
 * one goes up in z, and then steps out in rho., or in theta and then
 * z as the case may be.
 *
 * With domains, this routine has to be used carefully, and a better
 * statement would be that these routines map to and from the element
 * in wmain to the normally 2d element in an individual domain.
 *
 * So this means that if you want to get the ith and jth element of
 * domain 1, then one needs to give nstart+whatever to wind_n_to_ij
 *
 **********************************************************/

int
wind_n_to_ij (int ndom, int n, int *i, int *j)
{
  int n_use;
  if (zdom[ndom].coord_type == SPHERICAL)
  {
    *i = n - zdom[ndom].nstart;
    *j = 0;
  }
  n_use = n - zdom[ndom].nstart;
  *i = n_use / zdom[ndom].mdim;
  *j = n_use - (*i) * zdom[ndom].mdim;
  return (0);
}


/**********************************************************/
/**
 * @brief      Translate from the 2d element for an individual domain to the wind
 * element number
 *
 * @param [in] int  ndom   The domain of interest
 * @param [in] int  i   The element number for the first (x) dimension
 * @param [in] int  j   The element number for the second (z or theta) dimension
 * @param [in, out] int *  n   The element number in the wind
 * @return     
 * * 0 normally
 * * 1 if i or j are out of bounds, or if
 *      wind_ij_to_n is called for a SPHERICAL coordinate 
 *
 * @details
 * This is just a simple translation to find the wind domain element
 * number
 *
 * ### Notes ###
 * If i or j are out of bounds (high), then an element that is
 * in the grid is returned, but this is a serious error
 *
 **********************************************************/

int
wind_ij_to_n (int ndom, int i, int j, int *n)
{
  int ierror = 0;
  int ii, jj;

  if (zdom[ndom].coord_type == SPHERICAL)
  {
    Error ("wind_ij_to_n: wind_ij_to_n being called for spherical coordinates %d %d\n", i, j);
    *n = i;
    return (1);
  }

  ii = i;
  jj = j;

  if (ii >= zdom[ndom].ndim)
  {
    ii = zdom[ndom].ndim - 1;
    ierror = 1;
  }
  if (jj >= zdom[ndom].mdim)
  {
    jj = zdom[ndom].mdim - 1;
    ierror = 1;
  }

  if (ierror)
  {
    Error ("wind_ij_to_n: i %d >= %d j %d >= %d being requested\n", i, zdom[ndom].ndim, j, zdom[ndom].mdim);
  }

  *n = zdom[ndom].nstart + ii * zdom[ndom].mdim + jj;   // MDIM because the array is in z order
  return (ierror);
}


/**********************************************************/
/**
 * @brief      Determine the wind domain element number from a postion
 *
 * @param [in] double  x[]   A three vector describing a position
 * @param [out] int *  n   The element number in the wind domain
 * @return     Also retruns the element number in the wind domain
 *
 * @details
 *
 * The routine simply cycles through the various domains until
 * if finds one in where x is in that domain.
 *
 * ### Notes ###
 *
 * @bug  There are no checks for whether x is within the wind.
 * This seems very dangerous.
 *
 **********************************************************/

int
wind_x_to_n (double x[], int *n)
{
  int ndom, i, j;
  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    for (i = 0; i < zdom[ndom].ndim - 1; i++)
    {
      if (x[0] > zdom[ndom].wind_x[i] && x[0] < zdom[ndom].wind_x[i + 1])
      {
        for (j = 0; j < zdom[ndom].mdim - 1; j++)
        {
          if (x[2] > zdom[ndom].wind_z[j] && x[2] < zdom[ndom].wind_z[j + 1])
          {
            wind_ij_to_n (ndom, i, j, n);
            return (*n);
          }
        }
      }
    }
  }
  return (*n);
}
