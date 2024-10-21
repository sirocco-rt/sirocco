
/***********************************************************/
/** @file  spherical.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  These are generic routines for setting up 1d spherical
 * coordinate system, as well as a special case of this the so-called
 * shell-wind model
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

/* Notes on spherical coordinates

   In spherical coordinates, we set NDIM = to whatever and MDIM = 1.

   Some of the variables in the Wind_ptr array, like the divergence
   have to be calculated at a certain place.  Logically, one would
   do this along one of the principle axes.  But at present we have
   calculated this along a 45 degree angle.  The reason we did this
   was because I want the code to run for normal bipolar winds.  But
   this can be an issue for anisotropic scattering!!  It should not
   affect items like the divergence, but it certainly affects the
   gradients.

*/


/**********************************************************/
/**
 * @brief      calculates the distance to the far
 *         boundary of the cell in which the photon bundle resides.
 *
 * @param [in] int  ndom   The domain in which the photon resides
 * @param [in, out] PhotPtr  p   Photon pointer
 * @return     Distance to the far boundary of the cell in which the photon
 * 	currently resides.  Negative numbers (and zero) should be
 * 	regarded as errors.
 *
 * @details
 *
 * The routine simply determines which grid cell the photon is
 * in and then solves two quadradic equations for the distance
 * to the inner and outer boundary of the cell.  The routine
 * returns the smallest positive distance.
 *
 * ### Notes ###
 *
 **********************************************************/

double
spherical_ds_in_cell (ndom, p)
     int ndom;
     PhotPtr p;

{

  int n, ix;
  double s, smax;


  if ((p->grid = n = where_in_grid (ndom, p->x)) < 0)
  {
    Error ("spherical_ds_in_cell: Photon %d not in grid when routine entered\n", p->np);
    return (n);
  }

  ix = n;

  /* Set up the quadratic equations in the radial  direction */

  smax = ds_to_sphere (zdom[ndom].wind_x[ix], p);
  s = ds_to_sphere (zdom[ndom].wind_x[ix + 1], p);

  if (smax == VERY_BIG && s == VERY_BIG)
  {
    Error ("spherical: ds_in_cell: s and smax returning VERY_BIG in cell %i nudging photon %d by DFUDGE\n", p->grid, p->np);
    return (DFUDGE);
  }


  if (s < smax)
    smax = s;

  if (smax <= 0)
  {
    Error ("spherical: ds_in_cell: smax %f negative for photon %d\n", smax, p->np);
  }
  return (smax);
}





/**********************************************************/
/**
 * @brief      defines the cells in a spherical grid
 *
 * @param [in] WindPtr  w   The structure which defines the wind in Python
 * @param [in] int  ndom   The domain number
 * @return     Always returns 0
 *
 * @details
 *
 * The routine defines the boundaries of cells for a spherical
 * grid in the appropriate parts of the wind domain
 * In spherical coordinates w runs from 0 to n.
 *
 * ### Notes ###
 *
 * The centers of the grid cells are defined in the xz plane at a 45
 * degree angle.  This was done so that one would be in a plausible
 * region of a biconical wind.
 *
 *
 **********************************************************/

int
spherical_make_grid (int ndom, WindPtr w)
{
  double dr, dlogr;
  int j, n;
  int ndim;
  double xfudge;


  ndim = zdom[ndom].ndim;


  for (j = 0; j < ndim; j++)
  {
    n = j + zdom[ndom].nstart;  // This is the element in wmain

    /*Define the grid points */
    if (zdom[ndom].log_linear == 1)
    {                           // linear intervals

      dr = (zdom[ndom].rmax - zdom[ndom].rmin) / (ndim - 3);
      w[n].r = zdom[ndom].rmin + j * dr;
      w[n].rcen = w[n].r + 0.5 * dr;
    }
    else
    {                           //logarithmic intervals

      dlogr = (log10 (zdom[ndom].rmax / zdom[ndom].rmin)) / (ndim - 3);
      w[n].r = zdom[ndom].rmin * pow (10., dlogr * (j - 1));
      w[n].rcen = 0.5 * zdom[ndom].rmin * (pow (10., dlogr * (j)) + pow (10., dlogr * (j - 1)));
    }

    /* Now calculate the positions of these points in the xz plane.
       There is a choice about how one does this.  Here we  have elected
       to calculate this at a 45 degree anglein the hopes this will
       be a reasonable place in the wind in a biconical flow.
     */

    w[n].x[1] = w[n].xcen[1] = 0.0;
    w[n].x[0] = w[n].x[2] = w[n].r * sin (PI / 4.);
    w[n].xcen[0] = w[n].xcen[2] = w[n].rcen * sin (PI / 4.);

    xfudge = (w[n].xcen[0] - w[n].x[0]);
    w[n].dfudge = XFUDGE * xfudge;
  }

  return (0);

}






/**********************************************************/
/**
 * @brief      Initialize some arrays in a domain for a spherical system that are used in
 * various interpolation routines
 *
 * @param [in] int  ndom   A domain which has a spherical grid defined
 * @param [in] WindPtr  w   the entire wind
 * @return     Alwasy returns 0
 *
 * @details
 *
 * This simple little routine just populates one dimensional
 * arrays that are used for interpolation.
 *
 * ### Notes ###
 *
 **********************************************************/

int
spherical_wind_complete (ndom, w)
     int ndom;
     WindPtr w;
{
  int i;
  int ndim, nstart;

  ndim = zdom[ndom].ndim;
  nstart = zdom[ndom].nstart;

  for (i = 0; i < ndim; i++)
    zdom[ndom].wind_x[i] = w[nstart + i].r;
  for (i = 0; i < ndim - 1; i++)
    zdom[ndom].wind_midx[i] = w[nstart + i].rcen;
  /* Add something plausible for the edges */
  zdom[ndom].wind_midx[ndim - 1] = 2. * zdom[ndom].wind_x[ndim - 1] - zdom[ndom].wind_midx[ndim - 2];
  zdom[ndom].wind_midz[0] = zdom[ndom].wind_z[0] / 2;

  return (0);
}



/*
 * What does RESOLUTION do?
 */

#define RESOLUTION   100

/**********************************************************/
/**
 * @brief      calculate the volume of a cell
 * 	allowing for the fact that a cell may be partially in the wind
 *
 * @param [in,out] WindPtr  w   a single wind cell to calculate the volume for
 * @return     Always returns 0
 *
 * @details
 *
 * The routine performs a simple 2d numerical
 * integration to find out what fraction of
 * a shell is in the wind of a specfic domain
 * and uses this to calculate the volume that
 * is in the wind.
 *
 * ### Notes ###
 *
 * The volume of edge cells are all set to 0
 *
 *
 **********************************************************/

int
spherical_cell_volume (WindPtr w)
{
  int i;
  double fraction;
  double num, denom;
  double r, theta;
  double dr, dtheta, x[3];
  double rmin, rmax;
  int kk, jj;
  int ndim;
  int ndom;

  const double thetamin = 0.0;
  const double thetamax = 0.5 * PI;

  i = w->nwind;
  ndom = w->ndom;
  ndim = zdom[ndom].ndim;
  rmin = zdom[ndom].wind_x[i];
  rmax = zdom[ndom].wind_x[i + 1];

  w->vol = 4. / 3. * PI * (rmax * rmax * rmax - rmin * rmin * rmin);

  if (i == ndim - 1 || i == ndim - 2)
  {
    fraction = 0.0;             /* Force outside edge volues to zero */
    jj = 0;
    kk = RESOLUTION;
  }
  else
  {                             /* Determine whether the cell is in the wind */
    num = denom = 0;
    jj = kk = 0;
    dr = (rmax - rmin) / RESOLUTION;
    dtheta = (thetamax - thetamin) / RESOLUTION;
    for (r = rmin + dr / 2; r < rmax; r += dr)
    {
      for (theta = thetamin + dtheta / 2; theta < thetamax; theta += dtheta)
      {
        denom += r * r * sin (theta);
        kk++;
        x[0] = r * sin (theta);
        x[1] = 0;
        x[2] = r * cos (theta);
        if (where_in_wind (x, &ndom) == W_ALL_INWIND)
        {
          num += r * r * sin (theta);
          jj++;
        }
      }
    }
    fraction = num / denom;
  }

  if (jj == 0)
  {
    w->inwind = W_NOT_INWIND;   // The cell is not in the wind
    w->vol = 0.0;
  }
  else if (jj == kk)
  {
    w->inwind = W_ALL_INWIND;   // The cell is completely in the wind
  }
  else
  {
    w->inwind = W_PART_INWIND;  //The cell is partially in the wind
    w->vol *= fraction;
  }

  return (0);
}





/**********************************************************/
/**
 * @brief      locates the element in wmain corrosponding to a position
 * 	when one is using spherical coordinates.
 *
 * @param [in] int  ndom   The domain number
 * @param [in] double  x[]   The postion
 * @return     the number of wind element associated with
 *  		a position.
 *
 * If the position is in the grid this will be a positive number.  If
 * s is inside the grid -1 will be returned, if x is outsice the domain
 * -2 will be returned
 *
 * @details
 *
 * ### Notes ###
 *  speherical_where_in grid does not tell you whether the x is in the wind or not.
 *
 * 	What one means by inside or outside the grid may well be different
 * 	for different coordinate systems.
 *
 * 	This routine is not normally used directly.  Instead it is called via
 * 	where_in_grid which determines the specific where_in_grid
 * 	routine to call, depending on the coordinate system.
 *
 **********************************************************/

int
spherical_where_in_grid (ndom, x)
     int ndom;
     double x[];
{
  int n;
  double r;
  double f;
  int ndim;


  ndim = zdom[ndom].ndim;

  r = length (x);

  /* Check to see if x is outside the region of the calculation */
  if (r > zdom[ndom].wind_x[ndim - 1])
  {
    return (-2);                /* x is outside grid */
  }
  else if (r < zdom[ndom].wind_x[0])
  {
    return (-1);                /*x is inside grid */
  }

  fraction (r, zdom[ndom].wind_x, ndim, &n, &f, 0);

  /* n is the position with this domain, so zdom[ndom].nstart is added get
   * to wmain
   */

  return (n + zdom[ndom].nstart);
}



/**********************************************************/
/**
 * @brief
 *
 * @param [in] int  n   -- Cell in which random position is to be generated
 * @param [out] double  x[]   A random position in the cell
 * @return     The retine generally returns W_ALL_IN_WIND
 *
 * @details
 *
 * The routine generates random positions within a spherical
 * cell and returns when it finds a position that is in the
 * wind region.
 *
 * ### Notes ###
 *
 * Comment: This routine has no protection against the possibility
 * that the cell is not at least partially in the wind.
 *
 **********************************************************/

int
spherical_get_random_location (n, x)
     int n;                     // Cell in which to create position
     double x[];                // Returned position
{
  int i, j;
  int inwind;
  double r, rmin, rmax;
  double theta, phi;
  int ndom, ndomain;

  ndom = wmain[n].ndom;
  wind_n_to_ij (ndom, n, &i, &j);
  rmin = zdom[ndom].wind_x[i];
  rmax = zdom[ndom].wind_x[i + 1];


  /* Generate a position which is both in the cell and in the wind */
  inwind = W_NOT_INWIND;
  while (inwind != W_ALL_INWIND)
  {
    r = (rmin * rmin * rmin) + (rmax * rmax * rmax - rmin * rmin * rmin) * random_number (0.0, 1.0);

    r = pow (r, (1. / 3.));
    theta = acos (random_number (-1.0, 1.0));

    phi = 2. * PI * random_number (0.0, 1.0);

/* Project from r, theta phi to x y z  */
    x[0] = r * cos (phi) * sin (theta);
    x[1] = r * sin (phi) * sin (theta);
    x[2] = r * cos (theta);
    inwind = where_in_wind (x, &ndomain);       /* Some photons will not be in the wind */
  }

  return (inwind);
}



/**********************************************************/
/**
 * @brief      extends the density to
 * regions just outside the wind regiions so that
 * extrapolations of density can be made there
 *
 * @param [in] int  ndom   The domain of interest
 * @param [in] WindPtr  w   The enetire wind
 * @return   Always returns 0
 *
 * @details
 *
 * In order for density to be continuous within the grid,
 * and to enable linear interpolation everwhere, we need
 * to define a density in the guard cells inside and outside
 * the wind region.
 *
 * This routine accomplishes this by assignining an element
 * in the plasma structure to the guard cells.
 *
 * ### Notes ###
 *
 * In principle one could interpolate on any of the prameters
 * in the plasma stucture, but in reality the main thing this
 * is used for is to assure tha ion denisites are continuous.
 *
 **********************************************************/

int
spherical_extend_density (ndom, w)
     int ndom;
     WindPtr w;
{

  int j, n, m;
  int ndim, nstart;

  ndim = zdom[ndom].ndim;
  nstart = zdom[ndom].nstart;
  /*
     Now we need to updated the densities immediately outside the wind so that the density interpolation in resonate will work.
     In this case all we have done is to copy the densities from the cell which is just in the wind (as one goes outward) to the
     cell that is just inside (or outside) the wind.
   */

  for (j = 0; j < ndim - 1; j++)
  {
    n = nstart + j;
    if (w[n].vol == 0 || (modes.partial_cells == PC_EXTEND && w[n].inwind == W_PART_INWIND))

//DEBUG    if (w[n].vol == 0)          // Then the grid point is not in the wind

    {
      m = n + 1;
      if (w[m].vol > 0)         // Then grid point n is just inside the wind
      {                         // To extend we copy  copy the densities to the grid cell n
        w[n].nplasma = w[m].nplasma;

      }
      else if (j > 0)
      {
        m = n - 1;
        if (w[m].vol > 0)       // The grid point is just outside the wind
        {
          w[n].nplasma = w[m].nplasma;

        }
      }
    }
  }

  return (0);

}
