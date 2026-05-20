
/***********************************************************/
/** @file  cylindrical3d.c
 * @author ksl
 * @date   2025
 *
 * @brief  Routines for a full 3D cylindrical coordinate system (rho, phi, z).
 *
 * Cells are indexed as n = nstart + i*(mdim*pdim) + j*pdim + k where
 * i is the rho index, j is the z index (hemisphere-doubled, identical to
 * CYLIND), and k is the phi index (0..pdim-1 covering 0 to 2pi).
 *
 * Convention for per-cell position fields (matching CYLIND):
 *   x[0]    = rho_min  (cylindrical radius of inner boundary)
 *   x[1]    = 0        (unused; cells are not defined at a specific Cartesian y)
 *   x[2]    = z_min    (signed: negative for lower hemisphere)
 *   xcen[0] = rho_cen, xcen[1] = 0, xcen[2] = z_cen
 *   xmax[0] = rho_max, xmax[1] = 0, xmax[2] = z_max
 * Azimuthal extent is carried in the dedicated fields wmain[n].phi,
 * wmain[n].phicen, and wmain[n].phimax (all in radians, 0 to 2pi).
 *
 * Phase 1 provides: indexing helpers, make_grid, wind_complete,
 * where_in_grid, is_cell_in_wind, cell_volume, get_random_location,
 * extend_density.  ds_in_cell is a stub (Phase 2).
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

#define RESOLUTION 1000


/**********************************************************/
/**
 * @brief  Convert 3D grid indices (i, j, k) to a wmain array index.
 *
 * @param [in]  ndom   Domain number
 * @param [in]  i      Rho index  [0, ndim-1]
 * @param [in]  j      Z index    [0, mdim-1]  (mdim is hemisphere-doubled)
 * @param [in]  k      Phi index  [0, pdim-1]
 * @param [out] n      wmain array index
 * @return  0 on success, 1 if indices were clamped due to out-of-range input
 *
 **********************************************************/

int
wind_ijk_to_n (int ndom, int i, int j, int k, int *n)
{
  int ierror = 0;
  int ii = i, jj = j, kk = k;
  int ndim = zdom[ndom].ndim;
  int mdim = zdom[ndom].mdim;
  int pdim = zdom[ndom].pdim;

  if (ii >= ndim)
  {
    ii = ndim - 1;
    ierror = 1;
  }
  if (jj >= mdim)
  {
    jj = mdim - 1;
    ierror = 1;
  }
  if (kk >= pdim)
  {
    kk = pdim - 1;
    ierror = 1;
  }
  if (ierror)
  {
    Error ("wind_ijk_to_n: i %d j %d k %d out of range (%d %d %d) in domain %d\n", i, j, k, ndim, mdim, pdim, ndom);
  }

  *n = zdom[ndom].nstart + ii * (mdim * pdim) + jj * pdim + kk;
  return (ierror);
}


/**********************************************************/
/**
 * @brief  Convert a wmain index n to 3D grid indices (i, j, k).
 *
 * @param [in]  ndom   Domain number
 * @param [in]  n      wmain array index
 * @param [out] i      Rho index
 * @param [out] j      Z index
 * @param [out] k      Phi index
 *
 **********************************************************/

void
wind_n_to_ijk (int ndom, int n, int *i, int *j, int *k)
{
  int n_use = n - zdom[ndom].nstart;
  int pdim = zdom[ndom].pdim;
  int mdim = zdom[ndom].mdim;

  *i = n_use / (mdim * pdim);
  *j = (n_use / pdim) % mdim;
  *k = n_use % pdim;
}


/**********************************************************/
/**
 * @brief  Build the 3D cylindrical grid.
 *
 * @param [in] ndom   Domain number
 * @param [in] w      Pointer to the entire wmain wind array
 * @return  0
 *
 * @details
 * The rho and z grid layout is identical to cylind_make_grid:
 * logarithmic or linear spacing, hemisphere-doubled z axis.
 * A new outer loop over phi adds pdim azimuthal slices covering
 * 0 to 2pi uniformly.  All phi slices share the same (rho, z)
 * positions; the azimuthal extent is recorded in wmain[n].phi,
 * wmain[n].phicen, wmain[n].phimax.
 *
 * Cell position conventions follow CYLIND: x[0]=rho_min, x[1]=0,
 * x[2]=z_min (signed).  The phi fields carry the azimuthal bounds.
 *
 **********************************************************/

int
cylind3d_make_grid (int ndom, WindPtr w)
{
  double dr, dz, dlogr, dlogz, xfudge;
  double dphi;
  int i, j, k, j_half, n;
  int mdim_half;
  double zsign;
  double z_j, z_jp1, z_jcen;
  DomainPtr one_dom;

  one_dom = &zdom[ndom];
  mdim_half = one_dom->mdim / 2;

  if (one_dom->zmax == 0)
    one_dom->zmax = one_dom->rmax;

  Log ("cylind3d_make_grid: Making 3D cylindrical grid %d (pdim=%d)\n", ndom, one_dom->pdim);

  dphi = 2.0 * PI / one_dom->pdim;

  for (i = 0; i < one_dom->ndim; i++)
  {
    for (j = 0; j < one_dom->mdim; j++)
    {
      /* Hemisphere mapping — identical to cylind_make_grid */
      if (j >= mdim_half)
      {
        j_half = j - mdim_half;
        zsign = 1.0;
      }
      else
      {
        j_half = mdim_half - 1 - j;
        zsign = -1.0;
      }

      /* Rho and z positions */
      if (one_dom->log_linear == 1)
      {
        dr = one_dom->rmax / (one_dom->ndim - 3);
        dz = one_dom->zmax / (mdim_half - 3);

        w[one_dom->nstart + i * one_dom->mdim * one_dom->pdim].x[0] = i * dr;
        w[one_dom->nstart + i * one_dom->mdim * one_dom->pdim].xcen[0] = i * dr + 0.5 * dr;

        z_j = j_half * dz;
        z_jp1 = (j_half + 1) * dz;
        z_jcen = z_j + 0.5 * dz;
      }
      else
      {
        dlogr = (log10 (one_dom->rmax / one_dom->xlog_scale)) / (one_dom->ndim - 3);
        dlogz = (log10 (one_dom->zmax / one_dom->zlog_scale)) / (mdim_half - 3);

        if (dlogr <= 0)
        {
          Error ("cylind3d_make_grid: dlogr %g <= 0 for domain %d\n", dlogr, ndom);
          Exit (0);
        }
        if (dlogz <= 0)
        {
          Error ("cylind3d_make_grid: dlogz %g <= 0 for domain %d\n", dlogz, ndom);
          Exit (0);
        }

        double rho_i, rho_icen;
        if (i == 0)
        {
          rho_i = 0.0;
          rho_icen = 0.5 * one_dom->xlog_scale;
        }
        else
        {
          rho_i = one_dom->xlog_scale * pow (10., dlogr * (i - 1));
          rho_icen = 0.5 * one_dom->xlog_scale * (pow (10., dlogr * (i - 1)) + pow (10., dlogr * i));
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
          z_jcen = 0.5 * one_dom->zlog_scale * (pow (10., dlogz * (j_half - 1)) + pow (10., dlogz * j_half));
        }

        /* Store rho into the k=0 slot temporarily; copied to all k below */
        w[one_dom->nstart + i * one_dom->mdim * one_dom->pdim].x[0] = rho_i;
        w[one_dom->nstart + i * one_dom->mdim * one_dom->pdim].xcen[0] = rho_icen;
      }

      double rho_i = w[one_dom->nstart + i * one_dom->mdim * one_dom->pdim].x[0];
      double rho_icen = w[one_dom->nstart + i * one_dom->mdim * one_dom->pdim].xcen[0];

      double z_lo = (zsign > 0) ? z_j : -z_jp1;
      double z_cen = (zsign > 0) ? z_jcen : -z_jcen;

      /* Now fill all phi slices for this (i, j) */
      for (k = 0; k < one_dom->pdim; k++)
      {
        wind_ijk_to_n (ndom, i, j, k, &n);

        w[n].x[0] = rho_i;
        w[n].x[1] = 0.0;
        w[n].x[2] = z_lo;
        w[n].xcen[0] = rho_icen;
        w[n].xcen[1] = 0.0;
        w[n].xcen[2] = z_cen;

        w[n].phi = k * dphi;
        w[n].phicen = (k + 0.5) * dphi;
        w[n].phimax = (k + 1) * dphi;

        xfudge = fmin (fabs (rho_icen - rho_i), fabs (z_cen - z_lo));
        w[n].dfudge = XFUDGE * xfudge;
      }
    }
  }

  return (0);
}


/**********************************************************/
/**
 * @brief  Populate the 1D lookup arrays used for grid location.
 *
 * @param [in] ndom   Domain number
 * @param [in] w      Pointer to the entire wmain wind array
 * @return  0
 *
 * @details
 * Fills wind_x, wind_z, wind_midx, wind_midz (shared with CYLIND),
 * and wind_phi, wind_midphi (new for CYLIND3D).
 * Also sets xmax for each cell by stepping one grid position forward.
 *
 **********************************************************/

int
cylind3d_wind_complete (int ndom, WindPtr w)
{
  int i, j, k, n;
  int nstart, ndim, mdim, pdim;
  double dphi;
  DomainPtr one_dom;

  one_dom = &zdom[ndom];
  nstart = one_dom->nstart;
  ndim = one_dom->ndim;
  mdim = one_dom->mdim;
  pdim = one_dom->pdim;
  dphi = 2.0 * PI / pdim;

  /* Rho axis: read from k=0 column of each i-row */
  for (i = 0; i < ndim; i++)
    one_dom->wind_x[i] = w[nstart + i * mdim * pdim].x[0];

  /* Z axis: read from k=0, i=0 column */
  for (j = 0; j < mdim; j++)
    one_dom->wind_z[j] = w[nstart + j * pdim].x[2];

  /* Mid-rho */
  for (i = 0; i < ndim - 1; i++)
    one_dom->wind_midx[i] = 0.5 * (one_dom->wind_x[i] + one_dom->wind_x[i + 1]);
  one_dom->wind_midx[ndim - 1] = 2.0 * one_dom->wind_x[ndim - 1] - one_dom->wind_midx[ndim - 2];

  /* Mid-z */
  for (j = 1; j < mdim - 1; j++)
    one_dom->wind_midz[j] = 0.5 * (one_dom->wind_z[j] + one_dom->wind_z[j + 1]);
  one_dom->wind_midz[mdim - 1] = 2.0 * one_dom->wind_z[mdim - 1] - one_dom->wind_midz[mdim - 2];
  one_dom->wind_midz[0] = 2.0 * one_dom->wind_z[0] - one_dom->wind_midz[1];

  /* Phi axis */
  for (k = 0; k < pdim; k++)
    one_dom->wind_phi[k] = k * dphi;
  one_dom->wind_phi[pdim] = 2.0 * PI;

  for (k = 0; k < pdim; k++)
    one_dom->wind_midphi[k] = (k + 0.5) * dphi;

  /* xmax: step one grid position forward in each axis */
  for (i = 0; i < ndim; i++)
  {
    double xmax_rho = (i < ndim - 1) ? one_dom->wind_x[i + 1] : 2.0 * one_dom->wind_x[ndim - 1] - one_dom->wind_x[ndim - 2];
    for (j = 0; j < mdim; j++)
    {
      double xmax_z = (j < mdim - 1) ? one_dom->wind_z[j + 1] : 2.0 * one_dom->wind_z[mdim - 1] - one_dom->wind_z[mdim - 2];
      for (k = 0; k < pdim; k++)
      {
        wind_ijk_to_n (ndom, i, j, k, &n);
        w[n].xmax[0] = xmax_rho;
        w[n].xmax[1] = 0.0;
        w[n].xmax[2] = xmax_z;
      }
    }
  }

  return (0);
}


/**********************************************************/
/**
 * @brief  Locate a Cartesian position in the 3D cylindrical grid.
 *
 * @param [in] ndom   Domain number
 * @param [in] x[]    Cartesian position
 * @return  wmain index, or -1 (inside rho_min), or -2 (outside grid)
 *
 * @details
 * Converts x to (rho, phi, z) and binary-searches wind_x, wind_z,
 * wind_phi independently.  phi is mapped to [0, 2pi) via atan2.
 *
 **********************************************************/

int
cylind3d_where_in_grid (int ndom, double x[])
{
  int i, j, k, n;
  double rho, phi, z;
  double f;
  DomainPtr one_dom;

  one_dom = &zdom[ndom];

  rho = sqrt (x[0] * x[0] + x[1] * x[1]);
  z = x[2];                     /* signed z — hemisphere-doubled grid */

  phi = atan2 (x[1], x[0]);
  if (phi < 0.0)
    phi += 2.0 * PI;

  /* Bounds check */
  if (rho > one_dom->wind_x[one_dom->ndim - 1] || z > one_dom->wind_z[one_dom->mdim - 1] || z < one_dom->wind_z[0])
    return (-2);

  if (rho < one_dom->wind_x[0])
    return (-1);

  fraction (rho, one_dom->wind_x, one_dom->ndim, &i, &f, 0);
  fraction (z, one_dom->wind_z, one_dom->mdim, &j, &f, 0);
  fraction (phi, one_dom->wind_phi, one_dom->pdim + 1, &k, &f, 0);
  if (k >= one_dom->pdim)
    k = one_dom->pdim - 1;

  wind_ijk_to_n (ndom, i, j, k, &n);
  return (n);
}


/**********************************************************/
/**
 * @brief  Check whether cell n is in the wind.
 *
 * @param [in] n   wmain index
 * @return  W_ALL_INWIND, W_PART_INWIND, or W_NOT_INWIND
 *
 * @details
 * Wind occupancy depends only on (rho, z) for any axisymmetric model.
 * The phi index k is ignored; the check uses the same corner-sampling
 * logic as cylind_is_cell_in_wind, evaluated in the xz plane (phi=0).
 *
 **********************************************************/

int
cylind3d_is_cell_in_wind (int n)
{
  int i, j, k, j_hemi;
  int mdim_half;
  double r, z, dr, dz;
  double rmin, rmax, zmin, zmax;
  double x[3];
  int ndom, ndomain;

  ndom = wmain[n].ndom;
  mdim_half = zdom[ndom].mdim / 2;
  wind_n_to_ijk (ndom, n, &i, &j, &k);

  j_hemi = (j >= mdim_half) ? (j - mdim_half) : (mdim_half - 1 - j);

  if (i >= (zdom[ndom].ndim - 2) || j_hemi >= (mdim_half - 2))
    return (W_NOT_INWIND);

  rmin = wmain[n].x[0];
  rmax = wmain[n].xmax[0];
  zmin = wmain[n].x[2];
  zmax = wmain[n].xmax[2];

  /* Check 4 corners in the (rho, z) plane; if all in wind return W_ALL_INWIND early */
  int n_corners = 0;
  double xcorner[3];
  xcorner[1] = 0.0;
  xcorner[0] = rmin;
  xcorner[2] = zmin;
  if (where_in_wind (xcorner, &ndomain) == W_ALL_INWIND && ndomain == ndom)
    n_corners++;
  xcorner[0] = rmax;
  xcorner[2] = zmin;
  if (where_in_wind (xcorner, &ndomain) == W_ALL_INWIND && ndomain == ndom)
    n_corners++;
  xcorner[0] = rmin;
  xcorner[2] = zmax;
  if (where_in_wind (xcorner, &ndomain) == W_ALL_INWIND && ndomain == ndom)
    n_corners++;
  xcorner[0] = rmax;
  xcorner[2] = zmax;
  if (where_in_wind (xcorner, &ndomain) == W_ALL_INWIND && ndomain == ndom)
    n_corners++;

  if (n_corners == 4)
    return (W_ALL_INWIND);
  if (n_corners == 0)
    return (W_NOT_INWIND);

  dr = (rmax - rmin) / RESOLUTION;
  dz = (zmax - zmin) / RESOLUTION;

  x[1] = 0.0;

  for (z = zmin + dz / 2; z < zmax; z += dz)
  {
    x[2] = z;
    x[0] = rmin;
    if (where_in_wind (x, &ndomain) == W_ALL_INWIND && ndomain == ndom)
      return (W_PART_INWIND);
    x[0] = rmax;
    if (where_in_wind (x, &ndomain) == W_ALL_INWIND && ndomain == ndom)
      return (W_PART_INWIND);
  }

  for (r = rmin + dr / 2; r < rmax; r += dr)
  {
    x[0] = r;
    x[2] = zmin;
    if (where_in_wind (x, &ndomain) == W_ALL_INWIND && ndomain == ndom)
      return (W_PART_INWIND);
    x[2] = zmax;
    if (where_in_wind (x, &ndomain) == W_ALL_INWIND && ndomain == ndom)
      return (W_PART_INWIND);
  }

  return (W_NOT_INWIND);
}


/**********************************************************/
/**
 * @brief  Calculate the valid (in-wind) volume of a 3D cylindrical cell.
 *
 * @param [in,out] w   Single wind cell
 * @return  0
 *
 * @details
 * Full cell volume = 0.5*(rho_max^2 - rho_min^2) * |z_max - z_min| * delta_phi.
 * Multiplied by the in-wind fraction determined by random sampling in (rho, z).
 *
 **********************************************************/

int
cylind3d_cell_volume (WindPtr w)
{
  int jj, kk;
  double fraction_inwind, cell_volume;
  double num, denom;
  double r, z;
  double rmax, rmin, zmin, zmax, dphi;
  double dr, dz, x[3];
  int n_inwind;
  int ndom, dom_check;

  ndom = w->ndom;

  rmin = w->x[0];
  rmax = w->xmax[0];
  zmin = w->x[2];
  zmax = w->xmax[2];
  dphi = w->phimax - w->phi;

  cell_volume = 0.5 * (rmax * rmax - rmin * rmin) * fabs (zmax - zmin) * dphi;

  if (w->inwind == W_NOT_ASSIGNED)
  {
    if (zdom[ndom].wind_type == IMPORT)
    {
      Error ("cylind3d_cell_volume: should not reassign inwind for imported model\n");
      Exit (0);
    }

    n_inwind = cylind3d_is_cell_in_wind (w->nwind);

    if (n_inwind == W_NOT_INWIND)
    {
      fraction_inwind = 0.0;
      jj = 0;
      kk = RESOLUTION * RESOLUTION;
    }
    else if (n_inwind == W_ALL_INWIND)
    {
      fraction_inwind = 1.0;
      jj = kk = RESOLUTION * RESOLUTION;
    }
    else
    {
      num = denom = 0.0;
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
          x[1] = 0.0;
          x[2] = z;
          if (where_in_wind (x, &dom_check) == W_ALL_INWIND && dom_check == ndom)
          {
            num += r * r;
            jj++;
          }
        }
      }
      fraction_inwind = (denom > 0.0) ? num / denom : 0.0;
    }

    if (jj == 0)
    {
      w->inwind = W_NOT_INWIND;
      w->vol = 0.0;
    }
    else if (jj == kk)
    {
      w->inwind = W_ALL_INWIND;
      w->vol = cell_volume;
    }
    else
    {
      w->inwind = W_PART_INWIND;
      w->vol = cell_volume * fraction_inwind;
    }
  }
  else if (w->inwind == W_NOT_INWIND)
  {
    w->vol = 0.0;
  }
  else if (w->inwind == W_ALL_INWIND)
  {
    w->vol = cell_volume;
  }

  return (0);
}


/**********************************************************/
/**
 * @brief  Generate a random position within a CYLIND3D cell.
 *
 * @param [in]  n    wmain index
 * @param [out] x[]  Cartesian position
 * @return  inwind status
 *
 * @details
 * Samples rho uniformly in rho^2 (uniform area element), phi
 * uniformly in [phi_min, phi_max], and z uniformly in [z_min, z_max].
 * Rejects positions outside the wind cone.
 *
 **********************************************************/

int
cylind3d_get_random_location (int n, double x[])
{
  double r, rmin, rmax, zmin, zmax, phi_lo, phi_hi, phi;
  int inwind, ndomain, ndom;

  ndom = wmain[n].ndom;
  rmin = wmain[n].x[0];
  rmax = wmain[n].xmax[0];
  zmin = wmain[n].x[2];
  zmax = wmain[n].xmax[2];
  phi_lo = wmain[n].phi;
  phi_hi = wmain[n].phimax;

  ndomain = -1;
  inwind = W_NOT_INWIND;

  while (inwind != W_ALL_INWIND || ndomain != ndom)
  {
    r = sqrt (rmin * rmin + random_number (0.0, 1.0) * (rmax * rmax - rmin * rmin));
    phi = phi_lo + random_number (0.0, 1.0) * (phi_hi - phi_lo);
    x[0] = r * cos (phi);
    x[1] = r * sin (phi);
    x[2] = zmin + random_number (0.0, 1.0) * (zmax - zmin);
    inwind = where_in_wind (x, &ndomain);
  }

  return (inwind);
}


/**********************************************************/
/**
 * @brief  Extend density to cells just outside the wind.
 *
 * @param [in] ndom   Domain number
 * @param [in] w      Pointer to the entire wmain wind array
 * @return  0
 *
 * @details
 * For each out-of-wind cell, associates the nearest in-wind neighbour's
 * plasma index, scanning in the rho direction.  Phi slices at the same
 * (i, j) share the same in-wind status, so the loop is over (i, j, k).
 *
 **********************************************************/

int
cylind3d_extend_density (int ndom, WindPtr w)
{
  int i, j, k, n, m;
  int ndim, mdim, pdim;

  ndim = zdom[ndom].ndim;
  mdim = zdom[ndom].mdim;
  pdim = zdom[ndom].pdim;

  for (i = 0; i < ndim - 1; i++)
  {
    for (j = 0; j < mdim - 1; j++)
    {
      for (k = 0; k < pdim; k++)
      {
        wind_ijk_to_n (ndom, i, j, k, &n);
        if (w[n].vol == 0 || (modes.partial_cells == PC_EXTEND && w[n].inwind == W_PART_INWIND))
        {
          wind_ijk_to_n (ndom, i + 1, j, k, &m);
          if (w[m].vol > 0)
          {
            w[n].nplasma = w[m].nplasma;
          }
          else if (i > 0)
          {
            wind_ijk_to_n (ndom, i - 1, j, k, &m);
            if (w[m].vol > 0)
              w[n].nplasma = w[m].nplasma;
          }
        }
      }
    }
  }

  return (0);
}


/**********************************************************/
/**
 * @brief  Distance to the far boundary of a CYLIND3D cell.
 *
 * @param [in] ndom   Domain number
 * @param [in] p      Photon pointer (position p->x, direction p->lmn)
 * @return  Distance smax to the nearest cell boundary, or negative
 *          if the photon is not in the grid.
 *
 * @details
 * Checks five boundaries: inner rho cylinder, outer rho cylinder,
 * lower z plane, upper z plane, and (for pdim > 1) the two phi
 * half-planes bounding the azimuthal sector.
 *
 * For pdim == 1 the cell covers the full 2π azimuth and there are
 * no phi boundaries, so the function reduces to cylind_ds_in_cell.
 *
 * A phi half-plane at angle phi0 is the set of points satisfying
 *   y·cos(phi0) − x·sin(phi0) = 0   with   x·cos(phi0) + y·sin(phi0) > 0.
 * The intersection time with the ray (x0+lx·t, y0+ly·t) is
 *   t = (y0·cos(phi0) − x0·sin(phi0)) / (lx·sin(phi0) − ly·cos(phi0))
 * and the intersection is valid only when the rho at that point is > 0.
 *
 **********************************************************/

double
ds_phi_boundary (double phi0, PhotPtr p)
{
  double sphi = sin (phi0), cphi = cos (phi0);
  double denom = p->lmn[0] * sphi - p->lmn[1] * cphi;

  if (denom == 0.0)
    return (VERY_BIG);

  double t = (p->x[1] * cphi - p->x[0] * sphi) / denom;

  if (t <= 0.0)
    return (VERY_BIG);

  /* Intersection must be on the rho > 0 half (not through the axis) */
  double rho_at_t = (p->x[0] + p->lmn[0] * t) * cphi + (p->x[1] + p->lmn[1] * t) * sphi;
  if (rho_at_t <= 0.0)
    return (VERY_BIG);

  return (t);
}


double
cylind3d_ds_in_cell (int ndom, PhotPtr p)
{
  int n, iroot;
  double a, b, c, root[2];
  double z1, z2, q, t;
  double smax;
  double rho_min, rho_max;

  if ((p->grid = n = where_in_grid (ndom, p->x)) < 0)
    return (n);

  smax = VERY_BIG;

  rho_min = wmain[n].x[0];
  rho_max = wmain[n].xmax[0];

  /* Rho cylinder boundaries — same quadratic as cylind_ds_in_cell */
  a = p->lmn[0] * p->lmn[0] + p->lmn[1] * p->lmn[1];
  b = 2. * (p->lmn[0] * p->x[0] + p->lmn[1] * p->x[1]);
  c = p->x[0] * p->x[0] + p->x[1] * p->x[1];

  iroot = quadratic (a, b, c - rho_min * rho_min, root);
  if (iroot >= 0 && root[iroot] < smax)
    smax = root[iroot];

  iroot = quadratic (a, b, c - rho_max * rho_max, root);
  if (iroot >= 0 && root[iroot] < smax)
    smax = root[iroot];

  /* Z plane boundaries */
  z1 = wmain[n].x[2];
  z2 = wmain[n].xmax[2];

  if (p->lmn[2] != 0.0)
  {
    q = (z1 - p->x[2]) / p->lmn[2];
    if (q > 0 && q < smax)
      smax = q;
    q = (z2 - p->x[2]) / p->lmn[2];
    if (q > 0 && q < smax)
      smax = q;
  }

  /* Phi half-plane boundaries (skipped when pdim == 1: full azimuth, no phi walls) */
  if (zdom[ndom].pdim > 1)
  {
    t = ds_phi_boundary (wmain[n].phi, p);
    if (t < smax)
      smax = t;

    t = ds_phi_boundary (wmain[n].phimax, p);
    if (t < smax)
      smax = t;
  }

  if (smax <= 0.0)
  {
    Error ("cylind3d_ds_in_cell: smax %g <= 0 for cell %d\n", smax, n);
    smax = wmain[n].dfudge;
  }

  return (smax);
}
