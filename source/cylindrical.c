
/***********************************************************/
/** @file  cylindrical.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  Generic routines for handling cylindrical coordinate
 * systems
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief      calculate the distance to the far
 * boundary of the cell in which the photon bundle when dealing
 * with a cylindrical domain
 *
 * @param [in] ndom   The number of the domain of interest
 * @param [in] p   Photon pointer
 * @return     distance to the far boundary of the cell
 *
 * Negative numbers (and zero) should be
 * regarded as errors.
 *
 * @details
 *
 * This routine solves the quadratic equations that allow one
 * to determine the distance a photon can travel within a cylindrical
 * cell before hitting the edge of the cell.
 *
 *
 * ### Notes ###
 *
 * The distance if the smallest positive root of the quadratic
 * equation.
 *
 *
 *
 **********************************************************/

double
cylind_ds_in_cell (ndom, p)
     int ndom;
     PhotPtr p;


{

  int n, ix, iz, iroot;
  double a, b, c, root[2];
  double z1, z2, q;
  double smax;



  /* The next lines should be unnecessary; they check/update
   * the cell number for the photon.
   */


  if ((p->grid = n = where_in_grid (ndom, p->x)) < 0)
  {
    return (n);                 /* Photon was not in wind */
  }

  wind_n_to_ij (ndom, n, &ix, &iz);     /*Convert the index n to two dimensions */

  smax = VERY_BIG;

  /* Set up the quadratic equations in the radial rho direction */

  a = (p->lmn[0] * p->lmn[0] + p->lmn[1] * p->lmn[1]);
  b = 2. * (p->lmn[0] * p->x[0] + p->lmn[1] * p->x[1]);
  c = p->x[0] * p->x[0] + p->x[1] * p->x[1];

  iroot = quadratic (a, b, c - zdom[ndom].wind_x[ix] * zdom[ndom].wind_x[ix], root);    /* iroot will be the smallest positive root
                                                                                           if one exists or negative otherwise */

  if (iroot >= 0 && root[iroot] < smax)
    smax = root[iroot];

  iroot = quadratic (a, b, c - zdom[ndom].wind_x[ix + 1] * zdom[ndom].wind_x[ix + 1], root);

  if (iroot >= 0 && root[iroot] < smax)
    smax = root[iroot];

  /* At this point we have found how far the photon can travel in rho in its
     current direction.  Now we must worry about motion in the z direction  */

  z1 = zdom[ndom].wind_z[iz];
  z2 = zdom[ndom].wind_z[iz + 1];
  if (p->x[2] < 0)
  {                             /* We need to worry about which side of the plane the photon is on! */
    z1 *= (-1.);
    z2 *= (-1.);
  }

  if (p->lmn[2] != 0.0)
  {
    q = (z1 - p->x[2]) / p->lmn[2];
    if (q > 0 && q < smax)
      smax = q;
    q = (z2 - p->x[2]) / p->lmn[2];
    if (q > 0 && q < smax)
      smax = q;

  }

  return (smax);
}





/**********************************************************/
/**
 * @brief      defines the cells in a cylindrical grid
 *
 * @param [in] int  ndom   The domain number of interest
 * @param [in] WindPtr  w   The structure which defines the wind in Python
 * @return   Always returns 0
 *
 * @details
 *
 * Creates the boundaries of cells for a cylindrical domain.
 *
 * Depending on the value of log_linear for the domain, the
 * cells are linearly or logarithmically spaced in x and z
 *
 * ### Notes ###
 *
 **********************************************************/

int
cylind_make_grid (int ndom, WindPtr w)
{
  double dr, dz, dlogr, dlogz, xfudge;
  int i, j, n;
  DomainPtr one_dom;


  one_dom = &zdom[ndom];

  if (zdom[ndom].zmax == 0)
  {
    /* Check if zmax has been set, and if not set it to rmax */
    zdom[ndom].zmax = zdom[ndom].rmax;
  }


  Log ("cylind_make_grid: Making cylindrical grid %d\n", ndom);
  Log ("cylind_make_grid: rmax %e for domain %d\n", one_dom->rmax, ndom);

  /* In order to interpolate the velocity (and other) vectors out to zdom[ndom].rmax, we need
     to define the wind at least one grid cell outside the region in which we want photons
     to propagate.  */


  /* First calculate parameters that are to be calculated at the edge of the grid cell.  This is
     mainly the positions and the velocity */
  for (i = 0; i < one_dom->ndim; i++)
  {
    for (j = 0; j < one_dom->mdim; j++)
    {
      wind_ij_to_n (ndom, i, j, &n);
      w[n].x[1] = w[n].xcen[1] = 0;     /*The cells are all defined in the xz plane */

      /*Define the grid points */
      if (one_dom->log_linear == 1)
      {                         // linear intervals

        dr = one_dom->rmax / (one_dom->ndim - 3);
        dz = one_dom->zmax / (one_dom->mdim - 3);
        w[n].x[0] = i * dr;     /* The first zone is at the inner radius of
                                   the wind */
        w[n].x[2] = j * dz;
        w[n].xcen[0] = w[n].x[0] + 0.5 * dr;
        w[n].xcen[2] = w[n].x[2] + 0.5 * dz;

      }
      else
      {                         //logarithmic intervals

        dlogr = (log10 (one_dom->rmax / one_dom->xlog_scale)) / (one_dom->ndim - 3);
        dlogz = (log10 (one_dom->zmax / one_dom->zlog_scale)) / (one_dom->mdim - 3);

        if (dlogr <= 0)
        {
          Error ("cylindrical: dlogr %g is less than 0.  This is certainly wrong! Aborting\n", dlogr);
          Exit (0);
        }

        if (dlogz <= 0)
        {
          Error ("cylindrical: dlogz %g is less than 0.  This is certainly wrong! Aborting\n", dlogz);
          Exit (0);
        }

        if (i == 0)
        {
          w[n].x[0] = 0.0;
          w[n].xcen[0] = 0.5 * one_dom->xlog_scale;
        }
        else
        {
          w[n].x[0] = one_dom->xlog_scale * pow (10., dlogr * (i - 1));
          w[n].xcen[0] = 0.5 * one_dom->xlog_scale * (pow (10., dlogr * (i - 1)) + pow (10., dlogr * (i)));
        }

        if (j == 0)
        {
          w[n].x[2] = 0.0;
          w[n].xcen[2] = 0.5 * one_dom->zlog_scale;
        }
        else
        {
          w[n].x[2] = one_dom->zlog_scale * pow (10, dlogz * (j - 1));
          w[n].xcen[2] = 0.5 * one_dom->zlog_scale * (pow (10., dlogz * (j - 1)) + pow (10., dlogz * (j)));
        }
      }

      xfudge = fmin ((w[n].xcen[0] - w[n].x[0]), (w[n].xcen[2] - w[n].x[2]));
      w[n].dfudge = XFUDGE * xfudge;

    }
  }

  return (0);
}




/**********************************************************/
/**
 * @brief      Create arrays that are used for interpolating
 * quantities, such as density and velocity
 *
 * @param [in] int  ndom   The domain number of interest
 * @param [in] WindPtr  w   The entire wind
 * @return     Always returns 0
 *
 *
 * @details
 * This simple little routine just populates two one dimensional arrays
 * that are used for interpolation.
 *
 * ### Notes ###
 *
 * There is an analogous routine for each different type of coordinate
 * system
 *
 **********************************************************/

int
cylind_wind_complete (ndom, w)
     int ndom;
     WindPtr w;
{
  int i, j;
  int nstart, mdim, ndim;
  DomainPtr one_dom;

  one_dom = &zdom[ndom];
  nstart = one_dom->nstart;
  ndim = one_dom->ndim;
  mdim = one_dom->mdim;

  /* Finally define some one-d vectors that make it easier to locate a photon in the wind given that we
     have adoped a "rectangular" grid of points.  Note that rectangular does not mean equally spaced. */

  for (i = 0; i < ndim; i++)
    one_dom->wind_x[i] = w[nstart + i * mdim].x[0];

  for (j = 0; j < one_dom->mdim; j++)
    one_dom->wind_z[j] = w[nstart + j].x[2];

  for (i = 0; i < one_dom->ndim - 1; i++)
    one_dom->wind_midx[i] = 0.5 * (w[nstart + i * mdim].x[0] + w[nstart + (i + 1) * mdim].x[0]);

  for (j = 0; j < one_dom->mdim - 1; j++)
    one_dom->wind_midz[j] = 0.5 * (w[nstart + j].x[2] + w[nstart + j + 1].x[2]);

  /* Add something plausible for the edges */
  one_dom->wind_midx[one_dom->ndim - 1] = 2. * one_dom->wind_x[nstart + ndim - 1] - one_dom->wind_midx[nstart + ndim - 2];
  one_dom->wind_midz[one_dom->mdim - 1] = 2. * one_dom->wind_z[nstart + mdim - 1] - one_dom->wind_midz[nstart + mdim - 2];

  return (0);
}



#define RESOLUTION   1000

/**********************************************************/
/**
 * @brief      Calculates the volume for a cell in a cylindrical domain
 * 	allowing for the fact that it may not necessarily be entirely in
 * 	the wind
 *
 * @param [in,out] WindPtr  w   a single wind cell to calculate the volume for
 *
 * @return     Always returns 0
 *
 * @details
 * This is a brute-froce integration of the volume
 *
 * 	ksl--04aug--some of the machinations regarding what to to at the
 * 	edge of the wind seem bizarre, like a substiture for figuring out
 * 	what one actually should be doing.  However, volume > 0 is used
 * 	for making certain choices in the existing code, and so one does
 * 	need to be careful.
 *
 * ### Notes ###
 * Where_in grid does not tell you whether the photon is in the wind or not.
 *
 **********************************************************/

int
cylind_cell_volume (WindPtr w)
{
  int i, j;
  int jj, kk;
  double fraction, cell_volume;
  double num, denom;
  double r, z;
  double rmax, rmin;
  double zmin, zmax;
  double dr, dz, x[3];
  int n_inwind;
  int ndom, dom_check;
  DomainPtr one_dom;

  ndom = w->ndom;
  one_dom = &zdom[ndom];
  wind_n_to_ij (ndom, w->nwind, &i, &j);

  rmin = one_dom->wind_x[i];
  rmax = one_dom->wind_x[i + 1];
  zmin = one_dom->wind_z[j];
  zmax = one_dom->wind_z[j + 1];

  /* this is the full cell volume, which is adjusted by the fraction
     of the cell that is in the wind later if necessary leading factor of 2
     added to allow for volume above and below plane (SSMay04) */
  cell_volume = 2 * PI * (rmax * rmax - rmin * rmin) * (zmax - zmin);

  if (w->inwind == W_NOT_ASSIGNED)
  {
    if (one_dom->wind_type == IMPORT)
    {
      Error ("cylind_volumes: Shouldn't be redefining inwind in cylind_volumes with imported model.\n");
      Exit (0);
    }

    n_inwind = cylind_is_cell_in_wind (w->nwind);

    if (n_inwind == W_NOT_INWIND)
    {
      fraction = 0.0;
      jj = 0;
      kk = RESOLUTION * RESOLUTION;
    }
    else if (n_inwind == W_ALL_INWIND)
    {
      fraction = 1.0;
      jj = kk = RESOLUTION * RESOLUTION;
    }
    else
    {                           /* Determine whether the grid cell is in the wind */
      num = denom = 0;
      jj = kk = 0;
      dr = (rmax - rmin) / RESOLUTION;
      dz = (zmax - zmin) / RESOLUTION;
      for (r = rmin + dr / 2; r < rmax; r += dr)
      {
        for (z = zmin + dz / 2; z < zmax; z += dz)
        {
          denom += r * r;
          kk++;
          x[0] = r;
          x[1] = 0;
          x[2] = z;

          if (where_in_wind (x, &dom_check) == W_ALL_INWIND && dom_check == ndom)
          {
            num += r * r;       /* 0 implies in wind */
            jj++;
          }
        }
      }
      fraction = num / denom;
    }

    /* OK now make the final assignement of nwind and fix the volumes */
    /* XXX JM - not clear why these additional if statements are necessary */
    if (jj == 0)
    {
      w->inwind = W_NOT_INWIND; // The cell is not in the wind
      w->vol = 0.0;
    }
    else if (jj == kk)
    {
      w->inwind = W_ALL_INWIND; // All of cell is inwind
      w->vol = cell_volume;
    }
    else
    {
      w->inwind = W_PART_INWIND;        // Some of cell is inwind
      w->vol = cell_volume * fraction;
    }
  }
  /* JM 1711 -- the following two if statements are for if the inwind values are
     already assigned, for example by an imported model */
  else if (w->inwind == W_NOT_INWIND)
  {
    w->vol = 0.0;               /* need to zero volumes for cells not in the wind */
  }
  else if (w->inwind == W_ALL_INWIND)
  {
    w->vol = cell_volume;
  }

  return (0);
}



/**********************************************************/
/**
 * @brief      locates the grid position of the vector,
 * 	when one is using cylindrical coordinates.
 *
 * @param [in] ndom   The domain number
 * @param [in] x[]   A position
 * @return     the element number in wmain associated with
 * a position.  If the position is in the grid this will be a positive
 * integer
 *
 * if the number is -1 then the position is inside (has rho less than
 * rhomin) the grid.  If -2, then it is outside the grid.
 *
 * @details
 *
 * ### Notes ###
 *
 * The routine does not tell you whether the x is in the wind or not,
 * just that it is in the region covered by the grid.
 *
 **********************************************************/

int
cylind_where_in_grid (ndom, x)
     int ndom;
     double x[];
{
  int i, j, n;
  double z;
  double rho;
  double f;
  DomainPtr one_dom;

  one_dom = &zdom[ndom];

  z = fabs (x[2]);              /* This is necessary to get correct answer above
                                   and below plane */
  if (z == 0)
    z = 1.e4;                   //Force z to be positive  02feb ksl
  rho = sqrt (x[0] * x[0] + x[1] * x[1]);       /* This is distance from z
                                                   axis */
  /* Check to see if x is outside the region of the calculation */
  if (rho > one_dom->wind_x[one_dom->ndim - 1] || z > one_dom->wind_z[one_dom->mdim - 1])
  {
    return (-2);                /* x is outside grid */
  }

  if (rho < one_dom->wind_x[0])
    return (-1);

  fraction (rho, one_dom->wind_x, one_dom->ndim, &i, &f, 0);
  fraction (z, one_dom->wind_z, one_dom->mdim, &j, &f, 0);

  /* At this point i,j are just outside the x position */

  wind_ij_to_n (ndom, i, j, &n);

  /* n is the array element in wmain */
  return (n);
}




/**********************************************************/
/**
 * @brief      Generate a position at a random location within a specific
 * cell of a domain with cylindrical coordinates
 *
 * @param [in] int  n  Cell in wmain in which random position is to be generated
 * @param [out] double  x[] The random position generated by the routine
 * @return     A status flag indicating whether the position is or is not in the wind.
 *
 * @details
 * The routine generates a position somewhere in the (wind portion of the) volumme defined
 * by a cell.
 *
 * ### Notes ###
 *
 * The routine generates a position within the entire volume, checks
 * that this position is in the wind region of the domain, and if
 * not, continues to generate new positons until it finds one
 * in the wind region.
 *
 * Comment: The routine does not check to see whehter there
 * if the cell contains any valid wind region, and if this
 * is not the case, would end up in an infinite loop
 *
 *
 *
 **********************************************************/

int
cylind_get_random_location (n, x)
     int n;
     double x[];
{
  int i, j;
  int inwind;
  double r, rmin, rmax, zmin, zmax;
  double zz;
  double phi;
  int ndom, ndomain;
  DomainPtr one_dom;

  ndom = wmain[n].ndom;
  ndomain = -1;
  one_dom = &zdom[ndom];

  wind_n_to_ij (ndom, n, &i, &j);

  rmin = one_dom->wind_x[i];
  rmax = one_dom->wind_x[i + 1];
  zmin = one_dom->wind_z[j];
  zmax = one_dom->wind_z[j + 1];

  /* Generate a position which is both in the cell and in the wind */

  inwind = W_NOT_INWIND;
  while (inwind != W_ALL_INWIND || ndomain != ndom)
  {
    r = sqrt (rmin * rmin + random_number (0.0, 1.0) * (rmax * rmax - rmin * rmin));

/* Generate the azimuthal location */
    phi = 2. * PI * random_number (0.0, 1.0);
    x[0] = r * cos (phi);
    x[1] = r * sin (phi);



    x[2] = zmin + (zmax - zmin) * random_number (0.0, 1.0);
    inwind = where_in_wind (x, &ndomain);       /* Some photons will not be in the wind
                                                   because the boundaries of the wind split the grid cell */
  }

  zz = random_number (-0.5, 0.5);       //positions above are all at +z distances
  if (zz < 0)
    x[2] *= -1;                 /* The photon is in the bottom half of the wind */

  return (inwind);
}




/**********************************************************/
/**
 * @brief      extends the density to
 * 	regions just outside the wind regiions so that
 * 	extrapolations of density can be made there
 *
 * @param [in] int  ndom   The domain of interest
 * @param [in] WindPtr  w   The entire wind
 * @return     Always returns 0
 *
 * @details
 *
 * Densities and other parameters for wind cells are actually
 * in the PlasmaPtrs and elements wmain with any valid wind region
 * have an associated element in the Plasma structure, which
 * contains denisties, etc, for that cell.  For the purpose
 * of extapolating quantities, like density to the very edge
 * of the wind, this routine just assigns an element of the plasma
 * to a wind cell that is just outside the domain.  This allows
 * one to use the same linear interpolation routines throughout the
 * wind.
 *
 * ### Notes ###
 *
 * In principle one can interpolate over any variable contained
 * in the plasma pointer so that it is continuous within the wind.
 * In practive, the main thing that is interpolated are densities
 * of ions.
 *
 * Comment: One must be careful never to do anything to the values
 * in the Plasma element as a reult of anything happening outside
 * the wind.
 *
 **********************************************************/

int
cylind_extend_density (ndom, w)
     int ndom;
     WindPtr w;
{

  int i, j, n, m;
  int ndim, mdim;

  ndim = zdom[ndom].ndim;
  mdim = zdom[ndom].mdim;


  for (i = 0; i < ndim - 1; i++)
  {
    for (j = 0; j < mdim - 1; j++)
    {
      wind_ij_to_n (ndom, i, j, &n);
//DEBUG    if (w[n].vol == 0)
      if (w[n].vol == 0 || (modes.partial_cells == PC_EXTEND && w[n].inwind == W_PART_INWIND))

      {                         /*Then this grid point is not in the wind  */

        wind_ij_to_n (ndom, i + 1, j, &m);
        if (w[m].vol > 0)
        {                       /*Then the windcell in the +x direction is in the wind and
                                   we associate that plasma element with this grid cell
                                 */
          w[n].nplasma = w[m].nplasma;

        }
        else if (i > 0)
        {
          wind_ij_to_n (ndom, i - 1, j, &m);
          if (w[m].vol > 0)
          {                     /*Then the grid cell in the -x direction is in the wind and
                                   we associate that plasma element with this grid cell
                                 */
            w[n].nplasma = w[m].nplasma;

          }
        }
      }
    }
  }

  return (0);

}

/*

cylind_is_cell_in_wind (n)


Note that it simply calls where_in_wind multiple times.

History:
  11Aug	ksl	70b - Modified to incoporate torus
  		See sirocco.h for more complete explanation
		of how PART and ALL are related
  15sep ksl	Modified to ask the more refined question of
  		whether this cell is in the wind of the
		specific domain assigned to the cell.


*/


/**********************************************************/
/**
 * @brief      Check whether this cell is in the wind or not
 *
 * @param [in] int  n   The number of an element in the wind structure
 * @return     Various intergers depending on whether the cell is W_ALL_INWIND,
 * W_NOT_INWIND, W_PART_INWIND, etc.
 *
 * @details
 *
 * This routine performs is a robust check of whether a cell is in the wind or not,
 * by looking at the four corners of the cell as defined in wmain.
 * It was created to speed up the evaluation of the volumes for the wind.  It
 * checks each of the four boundaries of the wind to see whether any portions
 * of these are in the wind
 *
 * ### Notes ###
 *
 * This rooutine gets the domain number from wmain[n].ndom, and uses the
 * wind_x and wind_y postions to get the coreners, and then
 * calls where in wind multiple times.
 *
 * Comment: As with the case of a number of other routines the logic
 * of this is quite convoluted, reflecting the fact that much
 * of the code was written prior to the introduction of domains
 *
 **********************************************************/

int
cylind_is_cell_in_wind (n)
     int n;                     // cell number
{
  int i, j;
  double r, z, dr, dz;
  double rmin, rmax, zmin, zmax;
  double x[3];
  int ndom, ndomain;
  DomainPtr one_dom;

  ndom = wmain[n].ndom;
  one_dom = &zdom[ndom];
  wind_n_to_ij (ndom, n, &i, &j);

  /* First check if the cell is in the boundary */
  if (i >= (one_dom->ndim - 2) || j >= (one_dom->mdim - 2))
  {
    return (W_NOT_INWIND);
  }

  /* Assume that if all four corners are in the wind that the
     entire cell is in the wind.  check_corners also now checks
     that the corners are in the wind of a specific domain */

  if (check_corners_inwind (n) == 4)
  {
    return (W_ALL_INWIND);
  }

  /* So at this point, we have dealt with the easy cases
     of being in the region of the wind which is used for the
     boundary and when all the edges of the cell are in the wind.
   */


  rmin = one_dom->wind_x[i];
  rmax = one_dom->wind_x[i + 1];
  zmin = one_dom->wind_z[j];
  zmax = one_dom->wind_z[j + 1];

  dr = (rmax - rmin) / RESOLUTION;
  dz = (zmax - zmin) / RESOLUTION;

  /* Check inner and outer boundary in the z direction */

  x[1] = 0;

  for (z = zmin + dz / 2; z < zmax; z += dz)
  {
    x[2] = z;

    x[0] = rmin;

    if (where_in_wind (x, &ndomain) == W_ALL_INWIND && ndom == ndomain)
    {
      return (W_PART_INWIND);
    }

    x[0] = rmax;

    if (where_in_wind (x, &ndomain) == W_ALL_INWIND && ndom == ndomain)
    {
      return (W_PART_INWIND);
    }
  }


  /* Check inner and outer boundary in the z direction */

  for (r = rmin + dr / 2; r < rmax; r += dr)
  {

    x[0] = r;

    x[2] = zmin;

    if (where_in_wind (x, &ndomain) == W_ALL_INWIND && ndom == ndomain)
    {
      return (W_PART_INWIND);
    }

    x[2] = zmax;

    if (where_in_wind (x, &ndomain) == W_ALL_INWIND && ndom == ndomain)
    {
      return (W_PART_INWIND);
    }
  }

  /* If one has reached this point, then this wind cell is not in the wind */
  return (W_NOT_INWIND);
}
