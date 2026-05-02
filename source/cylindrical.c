
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
 * The z boundaries of the cell are taken directly from wmain[n].x[2]
 * (inner z, signed) and wmain[n].xmax[2] (outer z, signed).  Lower-hemisphere
 * cells have both values negative, so no sign flip is needed regardless of
 * the hemisphere the photon is in.
 *
 * ### Notes ###
 *
 * The distance is the smallest positive root of the quadratic
 * equation.
 *
 *
 *
 **********************************************************/

double
cylind_ds_in_cell (int ndom, PhotPtr p)
{

  int n, iroot;
  double a, b, c, root[2];
  double z1, z2, q;
  double smax;
  double rho_min, rho_max;



  /* The next lines should be unnecessary; they check/update
   * the cell number for the photon.
   */


  if ((p->grid = n = where_in_grid (ndom, p->x)) < 0)
  {
    return (n);                 /* Photon was not in wind */
  }

  smax = VERY_BIG;

  /* Cell rho and z boundaries come directly from the wind cell,
   * avoiding the need to convert n back to (i,j) indices. */
  rho_min = wmain[n].x[0];
  rho_max = wmain[n].xmax[0];

  /* Set up the quadratic equations in the radial rho direction */

  a = (p->lmn[0] * p->lmn[0] + p->lmn[1] * p->lmn[1]);
  b = 2. * (p->lmn[0] * p->x[0] + p->lmn[1] * p->x[1]);
  c = p->x[0] * p->x[0] + p->x[1] * p->x[1];

  iroot = quadratic (a, b, c - rho_min * rho_min, root);        /* iroot will be the smallest positive root
                                                                   if one exists or negative otherwise */

  if (iroot >= 0 && root[iroot] < smax)
    smax = root[iroot];

  iroot = quadratic (a, b, c - rho_max * rho_max, root);

  if (iroot >= 0 && root[iroot] < smax)
    smax = root[iroot];

  /* At this point we have found how far the photon can travel in rho in its
     current direction.  Now we must worry about motion in the z direction  */

  z1 = wmain[n].x[2];
  z2 = wmain[n].xmax[2];
  /* Lower-hemisphere cells have x[2] < 0 and xmax[2] < 0 already; no sign flip needed. */

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
 * @param [in] WindPtr  w   The structure which defines the wind in Sirocco
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
  int i, j, j_half, n;
  int mdim_half;                /* number of z cells per hemisphere */
  double zsign;                 /* +1 for upper hemisphere, -1 for lower */
  double z_j, z_jp1, z_jcen;    /* z boundaries and centre for cell j_half */
  DomainPtr one_dom;


  one_dom = &zdom[ndom];
  mdim_half = one_dom->mdim / 2;        /* mdim was doubled in setup_domains */

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

  /* Grid layout (Layout B): j=0..mdim_half-1 is the lower hemisphere (negative z, sorted
     from outermost to innermost), j=mdim_half..mdim-1 is the upper hemisphere (positive z,
     innermost to outermost).  wind_z is therefore monotonically increasing across the full
     range, which allows a single binary search to locate any photon in step 2.  */

  /* First calculate parameters that are to be calculated at the edge of the grid cell.  This is
     mainly the positions and the velocity */
  for (i = 0; i < one_dom->ndim; i++)
  {
    for (j = 0; j < one_dom->mdim; j++)
    {
      wind_ij_to_n (ndom, i, j, &n);
      w[n].x[1] = w[n].xcen[1] = 0;     /*The cells are all defined in the xz plane */

      /* Map j to a within-hemisphere index j_half (0 = innermost, mdim_half-1 = outermost)
         and the sign of z for this hemisphere. */
      if (j >= mdim_half)
      {
        j_half = j - mdim_half; /* upper hemisphere */
        zsign = 1.0;
      }
      else
      {
        j_half = mdim_half - 1 - j;     /* lower hemisphere: j=0 is outermost */
        zsign = -1.0;
      }

      /*Define the grid points */
      if (one_dom->log_linear == 1)
      {                         // linear intervals

        dr = one_dom->rmax / (one_dom->ndim - 3);
        dz = one_dom->zmax / (mdim_half - 3);   /* spacing based on cells per hemisphere */
        w[n].x[0] = i * dr;
        w[n].xcen[0] = w[n].x[0] + 0.5 * dr;

        z_j = j_half * dz;
        z_jp1 = (j_half + 1) * dz;
        z_jcen = z_j + 0.5 * dz;

      }
      else
      {                         //logarithmic intervals

        dlogr = (log10 (one_dom->rmax / one_dom->xlog_scale)) / (one_dom->ndim - 3);
        dlogz = (log10 (one_dom->zmax / one_dom->zlog_scale)) / (mdim_half - 3);

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

        if (j_half == 0)
        {
          z_j = 0.0;
          z_jp1 = one_dom->zlog_scale;
          z_jcen = 0.5 * one_dom->zlog_scale;
        }
        else
        {
          z_j = one_dom->zlog_scale * pow (10., dlogz * (j_half - 1));
          z_jp1 = one_dom->zlog_scale * pow (10., dlogz * j_half);
          z_jcen = 0.5 * one_dom->zlog_scale * (pow (10., dlogz * (j_half - 1)) + pow (10., dlogz * (j_half)));
        }
      }

      /* Apply hemisphere sign.  For the upper hemisphere x[2] is the inner (lower) boundary
         and xmax[2] the outer (upper) boundary — both positive.  For the lower hemisphere
         x[2] is the inner (more negative) boundary and xmax[2] is the outer (less negative)
         boundary; z_jp1 gives the larger magnitude so it becomes the more-negative x[2]. */
      if (zsign > 0)
      {
        w[n].x[2] = z_j;
        w[n].xcen[2] = z_jcen;
      }
      else
      {
        w[n].x[2] = -z_jp1;
        w[n].xcen[2] = -z_jcen;
      }

      xfudge = fmin ((w[n].xcen[0] - w[n].x[0]), fabs (w[n].xcen[2] - w[n].x[2]));
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
cylind_wind_complete (int ndom, WindPtr w)
{
  int i, j, n;
  int nstart, mdim, ndim;
  double xmax_rho, xmax_z;
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

  /* Add something plausible for the edges.  With both hemispheres present wind_midz also
     needs a plausible value at the lower (j=0) boundary of the lower hemisphere. */
  one_dom->wind_midx[one_dom->ndim - 1] = 2. * one_dom->wind_x[nstart + ndim - 1] - one_dom->wind_midx[nstart + ndim - 2];
  one_dom->wind_midz[one_dom->mdim - 1] = 2. * one_dom->wind_z[nstart + mdim - 1] - one_dom->wind_midz[nstart + mdim - 2];
  one_dom->wind_midz[0] = 2. * one_dom->wind_z[0] - one_dom->wind_midz[1];

  /* Set xmax (the far corner of each cell) using the now-populated wind_x and wind_z
     arrays.  For guard cells at the outer boundary of the grid the edge is extrapolated
     one grid spacing beyond the last defined value. */
  for (i = 0; i < ndim; i++)
  {
    xmax_rho = (i < ndim - 1) ? one_dom->wind_x[i + 1] : 2.0 * one_dom->wind_x[ndim - 1] - one_dom->wind_x[ndim - 2];
    for (j = 0; j < mdim; j++)
    {
      wind_ij_to_n (ndom, i, j, &n);
      xmax_z = (j < mdim - 1) ? one_dom->wind_z[j + 1] : 2.0 * one_dom->wind_z[mdim - 1] - one_dom->wind_z[mdim - 2];
      w[n].xmax[0] = xmax_rho;
      w[n].xmax[1] = 0.0;
      w[n].xmax[2] = xmax_z;
    }
  }

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

  rmin = w->x[0];
  rmax = w->xmax[0];
  zmin = w->x[2];
  zmax = w->xmax[2];

  /* Full cell volume for this single hemisphere cell; adjusted later by the
     fraction of the cell actually in the wind.  The historical factor-of-2 that
     accounted for both hemispheres is removed now that each hemisphere has its
     own cells. */
  cell_volume = PI * (rmax * rmax - rmin * rmin) * (zmax - zmin);

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
 * For parameterized models the routine uses signed z = x[2] directly and
 * searches the full wind_z[] array, which covers both hemispheres: lower
 * hemisphere in indices [0..mdim/2-1] (wind_z[] < 0) and upper hemisphere
 * in [mdim/2..mdim-1] (wind_z[] >= 0).  No folding to |z| is performed.
 *
 * For imported models the grid covers the upper hemisphere only (z >= 0),
 * so lower-hemisphere positions are folded via fabs(x[2]) to preserve
 * bilateral symmetry, as in the pre-doubled-grid code.
 *
 * ### Notes ###
 *
 * The routine does not tell you whether x is in the wind or not,
 * just that it is in the region covered by the grid.
 *
 **********************************************************/

int
cylind_where_in_grid (int ndom, double x[])
{
  int i, j, n;
  double z;
  double rho;
  double f;
  DomainPtr one_dom;

  one_dom = &zdom[ndom];

  /* Data-based hemisphere detection: if wind_z[0] >= 0 the grid only covers
     the upper hemisphere (either an imported single-hemisphere file or a
     legacy grid), so fold lower-hemisphere positions via fabs.  If wind_z[0]
     < 0 the grid covers both hemispheres (parameterized doubled grid or a
     two-hemisphere import), so use signed z directly. */
  if (zdom[ndom].wind_z[0] >= 0.0)
    z = fabs (x[2]);
  else
    z = x[2];                   /* signed z — each hemisphere has its own cells */
  rho = sqrt (x[0] * x[0] + x[1] * x[1]);       /* This is distance from z axis */

  /* Check to see if x is outside the region of the calculation.
     For parameterized grids wind_z[0] is the lower boundary of the lower
     hemisphere; for imported grids wind_z[0] is near zero. */
  if (rho > one_dom->wind_x[one_dom->ndim - 1] || z > one_dom->wind_z[one_dom->mdim - 1] || z < one_dom->wind_z[0])
  {
    return (-2);                /* x is outside grid */
  }

  if (rho < one_dom->wind_x[0])
    return (-1);

  fraction (rho, one_dom->wind_x, one_dom->ndim, &i, &f, 0);

  /* Search the full wind_z array — lower hemisphere in [0..mdim_half-1],
     upper hemisphere in [mdim_half..mdim-1] (parameterized); or the
     single-hemisphere imported grid. */
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
cylind_get_random_location (int n, double x[])
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

  /* For single-hemisphere models wind_z values are all >= 0, so x[2] from the
     loop above is always positive.  Randomly flip half the photons to z < 0
     to populate both hemispheres symmetrically.
     For two-hemisphere models wind_z[0] < 0 and x[2] already carries the
     correct sign from the cell boundaries — do not flip. */
  if (one_dom->wind_z[0] >= 0.0)
  {
    zz = random_number (-0.5, 0.5);
    if (zz < 0)
      x[2] *= -1;               /* The photon is in the bottom half of the wind */
  }

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
cylind_extend_density (int ndom, WindPtr w)
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
cylind_is_cell_in_wind (int n)
{
  int i, j, j_hemi;
  int mdim_half;
  double r, z, dr, dz;
  double rmin, rmax, zmin, zmax;
  double x[3];
  int ndom, ndomain;

  ndom = wmain[n].ndom;
  mdim_half = zdom[ndom].mdim / 2;
  wind_n_to_ij (ndom, n, &i, &j);

  /* Convert j to a within-hemisphere index for the guard-cell check.
     The outermost two cells of each hemisphere are guard cells. */
  j_hemi = (j >= mdim_half) ? (j - mdim_half) : (mdim_half - 1 - j);

  /* First check if the cell is in the boundary */
  if (i >= (zdom[ndom].ndim - 2) || j_hemi >= (mdim_half - 2))
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

  rmin = wmain[n].x[0];
  rmax = wmain[n].xmax[0];
  zmin = wmain[n].x[2];
  zmax = wmain[n].xmax[2];

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
