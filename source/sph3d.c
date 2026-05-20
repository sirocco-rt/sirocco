
/***********************************************************/
/** @file  sph3d.c
 * @author ksl
 * @date   2026
 *
 * @brief  Routines for a full 3D spherical polar coordinate system (r, theta, phi).
 *
 * Cells are indexed as n = nstart + i*(mdim*pdim) + j*pdim + k where
 * i is the r index, j is the theta index (hemisphere-doubled, identical to
 * RTHETA), and k is the phi index (0..pdim-1 covering 0 to 2pi).
 *
 * Convention for per-cell position fields (matching RTHETA):
 *   r, rcen        inner radius and cell-centre radius
 *   theta, thetacen inner polar angle (degrees) and cell-centre angle
 *   x[0]  = r_min * sin(theta_min)   (xz-plane Cartesian x)
 *   x[1]  = 0
 *   x[2]  = r_min * cos(theta_min)   (xz-plane Cartesian z)
 *   xcen[0] = rcen * sin(thetacen), xcen[2] = rcen * cos(thetacen)
 *   xmax[0] = r_max * sin(theta_max), xmax[2] = r_max * cos(theta_max)
 * Azimuthal extent is carried in wmain[n].phi, wmain[n].phicen, wmain[n].phimax
 * (radians, 0 to 2pi).  The wcone / wcone_max structures hold the inner and
 * outer theta-boundary cones exactly as for RTHETA.
 *
 * Phase 1b provides: make_grid, wind_complete, where_in_grid,
 * is_cell_in_wind, cell_volume, get_random_location, extend_density.
 * sph3d_ds_in_cell is a Phase-2 stub that delegates to rtheta_ds_in_cell
 * (ignores phi boundaries).
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

#define RESOLUTION 1000

/* Macro for cone dzdr with polar-axis special cases */
#define SPH3D_CONE_DZDR(theta_rad) \
  ((theta_rad) < 1e-10 ? VERY_BIG : \
   (theta_rad) > M_PI - 1e-10 ? -VERY_BIG : \
   1.0 / tan (theta_rad))


/**********************************************************/
/**
 * @brief  Build the 3D spherical polar grid.
 *
 * @param [in] ndom   Domain number
 * @param [in] w      Pointer to the entire wmain wind array
 * @return  0
 *
 * @details
 * The r and theta grid layout is identical to rtheta_make_grid:
 * logarithmic or linear r spacing, linear theta spacing covering
 * both hemispheres (mdim doubled by setup_domains).  A new outer
 * loop over phi adds pdim azimuthal slices covering 0 to 2pi.
 * All phi slices at the same (i, j) share identical (r, theta, x, xcen);
 * phi/phicen/phimax carry the azimuthal extent.
 *
 **********************************************************/

int
sph3d_make_grid (int ndom, WindPtr w)
{
  double dr, dlogr, dtheta, dphi;
  double theta, thetacen, theta_rad, thetacen_rad;
  double r_ij, rcen_ij;
  double x0, x2, xcen0, xcen2, xfudge;
  int i, j, k, n;
  int ndim, mdim, pdim, mdim_half;

  if (zdom[ndom].zmax == 0)
    zdom[ndom].zmax = zdom[ndom].rmax;

  ndim = zdom[ndom].ndim;
  mdim = zdom[ndom].mdim;
  pdim = zdom[ndom].pdim;
  mdim_half = mdim / 2;

  dtheta = 90.0 / mdim_half;
  dphi = 2.0 * PI / pdim;

  for (i = 0; i < ndim; i++)
  {
    /* Radial position — same formula as rtheta_make_grid */
    if (zdom[ndom].log_linear == 1)
    {
      dr = (zdom[ndom].rmax - geo.rstar) / (ndim - 3);
      r_ij = geo.rstar + i * dr;
      rcen_ij = r_ij + 0.5 * dr;
    }
    else
    {
      dlogr = (log10 (zdom[ndom].rmax / geo.rstar)) / (ndim - 3);
      r_ij = geo.rstar * pow (10., dlogr * (i - 1));
      rcen_ij = 0.5 * geo.rstar * (pow (10., dlogr * i) + pow (10., dlogr * (i - 1)));
    }

    for (j = 0; j < mdim; j++)
    {
      theta = dtheta * j;
      thetacen = theta + 0.5 * dtheta;

      theta_rad = theta / RADIAN;
      thetacen_rad = thetacen / RADIAN;

      x0 = r_ij * sin (theta_rad);
      x2 = r_ij * cos (theta_rad);
      xcen0 = rcen_ij * sin (thetacen_rad);
      xcen2 = rcen_ij * cos (thetacen_rad);

      xfudge = fabs (xcen0 - x0);
      if (xfudge == 0.0)
        xfudge = fabs (xcen2 - x2);

      for (k = 0; k < pdim; k++)
      {
        wind_ijk_to_n (ndom, i, j, k, &n);

        w[n].r = r_ij;
        w[n].rcen = rcen_ij;
        w[n].theta = theta;
        w[n].thetacen = thetacen;

        w[n].x[0] = x0;
        w[n].x[1] = 0.0;
        w[n].x[2] = x2;
        w[n].xcen[0] = xcen0;
        w[n].xcen[1] = 0.0;
        w[n].xcen[2] = xcen2;

        w[n].phi = k * dphi;
        w[n].phicen = (k + 0.5) * dphi;
        w[n].phimax = (k + 1) * dphi;

        w[n].dfudge = XFUDGE * xfudge;
      }
    }
  }

  return (0);
}


/**********************************************************/
/**
 * @brief  Populate 1D lookup arrays and per-cell boundary data.
 *
 * @param [in] ndom   Domain number
 * @param [in] w      Pointer to the entire wmain wind array
 * @return  0
 *
 * @details
 * Fills wind_x (r boundaries), wind_z (theta boundaries in degrees),
 * wind_phi (phi boundaries 0..2pi), and their midpoint arrays.
 * Sets per-cell rmax, thetamax, phimax, xmax, and theta-cone structures
 * wcone / wcone_max, following the same convention as rtheta_wind_complete.
 *
 **********************************************************/

int
sph3d_wind_complete (int ndom, WindPtr w)
{
  int i, j, k, n;
  int nstart, ndim, mdim, pdim;
  double dphi;
  double r_outer, theta_outer_deg, theta_outer_rad;
  double theta_inner_deg, theta_inner_rad;
  DomainPtr one_dom;

  one_dom = &zdom[ndom];
  nstart = one_dom->nstart;
  ndim = one_dom->ndim;
  mdim = one_dom->mdim;
  pdim = one_dom->pdim;
  dphi = 2.0 * PI / pdim;

  /* Radial lookup array: read r from cell (i, j=0, k=0) */
  for (i = 0; i < ndim; i++)
    one_dom->wind_x[i] = w[nstart + i * mdim * pdim].r;

  /* Theta lookup array: read theta from cell (i=0, j, k=0) */
  for (j = 0; j < mdim; j++)
    one_dom->wind_z[j] = w[nstart + j * pdim].theta;

  /* Mid-r */
  for (i = 0; i < ndim - 1; i++)
    one_dom->wind_midx[i] = w[nstart + i * mdim * pdim].rcen;
  one_dom->wind_midx[ndim - 1] = 2.0 * one_dom->wind_x[ndim - 1] - one_dom->wind_midx[ndim - 2];

  /* Mid-theta */
  for (j = 0; j < mdim - 1; j++)
    one_dom->wind_midz[j] = w[nstart + j * pdim].thetacen;
  one_dom->wind_midz[mdim - 1] = 2.0 * one_dom->wind_z[mdim - 1] - one_dom->wind_midz[mdim - 2];

  /* Phi lookup array: uniform 0..2pi */
  for (k = 0; k < pdim; k++)
    one_dom->wind_phi[k] = k * dphi;
  one_dom->wind_phi[pdim] = 2.0 * PI;

  /* Mid-phi */
  for (k = 0; k < pdim; k++)
    one_dom->wind_midphi[k] = (k + 0.5) * dphi;

  /* Per-cell boundary data: xmax, rmax, thetamax, phimax, wcone, wcone_max */
  for (i = 0; i < ndim; i++)
  {
    r_outer = (i < ndim - 1) ? one_dom->wind_x[i + 1] : 2.0 * one_dom->wind_x[ndim - 1] - one_dom->wind_x[ndim - 2];

    for (j = 0; j < mdim; j++)
    {
      theta_outer_deg = (j < mdim - 1) ? one_dom->wind_z[j + 1] : 2.0 * one_dom->wind_z[mdim - 1] - one_dom->wind_z[mdim - 2];
      theta_outer_rad = theta_outer_deg / RADIAN;

      theta_inner_deg = one_dom->wind_z[j];
      theta_inner_rad = theta_inner_deg / RADIAN;

      for (k = 0; k < pdim; k++)
      {
        wind_ijk_to_n (ndom, i, j, k, &n);

        w[n].rmax = r_outer;
        w[n].thetamax = theta_outer_deg;
        w[n].phimax = (k + 1) * dphi;

        w[n].xmax[0] = r_outer * sin (theta_outer_rad);
        w[n].xmax[1] = 0.0;
        w[n].xmax[2] = r_outer * cos (theta_outer_rad);

        w[n].wcone.z = 0.0;
        w[n].wcone.dzdr = SPH3D_CONE_DZDR (theta_inner_rad);

        w[n].wcone_max.z = 0.0;
        w[n].wcone_max.dzdr = SPH3D_CONE_DZDR (theta_outer_rad);
      }
    }
  }

  return (0);
}


/**********************************************************/
/**
 * @brief  Locate a Cartesian position in the 3D spherical polar grid.
 *
 * @param [in] ndom   Domain number
 * @param [in] x[]    Cartesian position
 * @return  wmain index, or -1 (inside r_min), or -2 (outside grid)
 *
 * @details
 * Converts x to (r, theta, phi) then binary-searches wind_x, wind_z,
 * and wind_phi independently.  Uses signed theta (0-180 deg) so
 * lower-hemisphere photons (theta > 90 deg) are routed correctly,
 * unless the grid only covers the upper hemisphere, in which case
 * lower positions are folded via fabs(x[2]).
 *
 **********************************************************/

int
sph3d_where_in_grid (int ndom, double x[])
{
  int i, j, k, n;
  double r, theta, phi;
  double f;
  DomainPtr one_dom;

  one_dom = &zdom[ndom];

  r = length (x);

  /* Hemisphere detection: same logic as rtheta_where_in_grid */
  if (one_dom->wind_z[one_dom->mdim - 1] <= 90.0)
    theta = acos (fabs (x[2]) / r) * RADIAN;
  else
    theta = acos (x[2] / r) * RADIAN;

  phi = atan2 (x[1], x[0]);
  if (phi < 0.0)
    phi += 2.0 * PI;

  if (r > one_dom->wind_x[one_dom->ndim - 1])
    return (-2);
  if (r < one_dom->wind_x[0])
    return (-1);

  fraction (r, one_dom->wind_x, one_dom->ndim, &i, &f, 0);
  fraction (theta, one_dom->wind_z, one_dom->mdim, &j, &f, 0);
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
 * Wind occupancy for azimuthally symmetric models depends only on
 * (r, theta).  The phi index k is ignored; the check uses the same
 * 4-corner and boundary-scan logic as rtheta_is_cell_in_wind,
 * evaluated in the xz plane (phi = 0).
 *
 **********************************************************/

int
sph3d_is_cell_in_wind (int n)
{
  int i, j, k;
  int mdim_half, is_lower, j_hemi;
  double r, theta;
  double rmin, rmax, thetamin, thetamax, dr, dtheta;
  double x[3];
  int ndom, ndim, mdim, ndomain;

  ndom = wmain[n].ndom;
  wind_n_to_ijk (ndom, n, &i, &j, &k);
  ndim = zdom[ndom].ndim;
  mdim = zdom[ndom].mdim;
  mdim_half = mdim / 2;

  /* Outer radial guard */
  if (i >= (ndim - 2))
    return (W_NOT_INWIND);

  /* Theta-direction guards (same asymmetric logic as rtheta_is_cell_in_wind) */
  is_lower = (j >= mdim_half);
  j_hemi = is_lower ? (j - mdim_half) : j;
  if (is_lower && j_hemi >= (mdim_half - 2))
    return (W_NOT_INWIND);

  /* Check 4 corners in the (r, theta) plane (phi = 0) */
  if (check_corners_inwind (n) == 4)
    return (W_ALL_INWIND);

  rmin = zdom[ndom].wind_x[i];
  rmax = zdom[ndom].wind_x[i + 1];
  thetamin = zdom[ndom].wind_z[j] / RADIAN;
  thetamax = zdom[ndom].wind_z[j + 1] / RADIAN;

  dr = (rmax - rmin) / RESOLUTION;
  dtheta = (thetamax - thetamin) / RESOLUTION;

  x[1] = 0.0;

  for (theta = thetamin + dtheta / 2.; theta < thetamax; theta += dtheta)
  {
    x[0] = rmin * sin (theta);
    x[2] = rmin * cos (theta);
    if (where_in_wind (x, &ndomain) == W_ALL_INWIND)
      return (W_PART_INWIND);

    x[0] = rmax * sin (theta);
    x[2] = rmax * cos (theta);
    if (where_in_wind (x, &ndomain) == W_ALL_INWIND)
      return (W_PART_INWIND);
  }

  for (r = rmin + dr / 2.; r < rmax; r += dr)
  {
    x[0] = r * sin (thetamin);
    x[2] = r * cos (thetamin);
    if (where_in_wind (x, &ndomain) == W_ALL_INWIND)
      return (W_PART_INWIND);

    x[0] = r * sin (thetamax);
    x[2] = r * cos (thetamax);
    if (where_in_wind (x, &ndomain) == W_ALL_INWIND)
      return (W_PART_INWIND);
  }

  return (W_NOT_INWIND);
}


/**********************************************************/
/**
 * @brief  Calculate the valid (in-wind) volume of a SPH3D cell.
 *
 * @param [in,out] w   Single wind cell
 * @return  0
 *
 * @details
 * Full cell volume = (delta_phi / 3) * (r_max^3 - r_min^3) *
 *                   (cos(theta_min) - cos(theta_max)).
 * This is the RTHETA 2pi-integrated formula scaled by delta_phi / (2*pi).
 * Multiplied by the in-wind fraction from rtheta-style boundary sampling.
 *
 **********************************************************/

int
sph3d_cell_volume (WindPtr w)
{
  int i, j, k, jj, kk;
  double fraction_inwind, cell_volume;
  double num, denom;
  double r, theta;
  double dr, dtheta;
  double rmin, rmax, thetamin, thetamax, dphi;
  double x[3];
  int n_inwind, ndom, ndomain;

  ndom = w->ndom;
  n_inwind = w->nwind;
  wind_n_to_ijk (ndom, n_inwind, &i, &j, &k);

  rmin = zdom[ndom].wind_x[i];
  rmax = zdom[ndom].wind_x[i + 1];
  thetamin = zdom[ndom].wind_z[j] / RADIAN;
  thetamax = zdom[ndom].wind_z[j + 1] / RADIAN;
  dphi = w->phimax - w->phi;

  cell_volume = (dphi / 3.0) * (rmax * rmax * rmax - rmin * rmin * rmin) * (cos (thetamin) - cos (thetamax));

  if (w->inwind == W_NOT_INWIND || w->inwind == W_IGNORE)
  {
    w->vol = 0.0;
  }
  else if (w->inwind == W_ALL_INWIND)
  {
    w->vol = cell_volume;
  }
  else if (w->inwind == W_NOT_ASSIGNED)
  {
    if (zdom[ndom].wind_type == IMPORT)
    {
      Error ("sph3d_cell_volume: should not reassign inwind for imported model\n");
      Exit (EXIT_FAILURE);
    }

    n_inwind = sph3d_is_cell_in_wind (w->nwind);

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
      dtheta = (thetamax - thetamin) / RESOLUTION;
      for (r = rmin + dr / 2; r < rmax; r += dr)
      {
        for (theta = thetamin + dtheta / 2; theta < thetamax; theta += dtheta)
        {
          denom += r * r * sin (theta);
          kk++;
          x[0] = r * sin (theta);
          x[1] = 0.0;
          x[2] = r * cos (theta);
          if (where_in_wind (x, &ndomain) == W_ALL_INWIND)
          {
            num += r * r * sin (theta);
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

  return (0);
}


/**********************************************************/
/**
 * @brief  Generate a random position within a SPH3D cell that is in the wind.
 *
 * @param [in]  n    wmain index
 * @param [out] x[]  Cartesian position
 * @return  W_ALL_INWIND once a valid in-wind position has been found
 *
 * @details
 * Samples r uniformly in r^3 (uniform volume element), theta via
 * sin(theta) sampling (uniform solid angle), and phi uniformly in
 * [phi_min, phi_max].  Lower-hemisphere cells are handled by the sign
 * of cos(theta) — same as rtheta_get_random_location.
 *
 **********************************************************/

int
sph3d_get_random_location (int n, double x[])
{
  int i, j, k;
  int inwind, ndomain;
  double r, rmin, rmax;
  double s, cos_theta, phi, phi_lo, phi_hi;
  double sthetamin, sthetamax;
  int ndom, mdim_half, is_lower_hemi;

  ndom = wmain[n].ndom;
  wind_n_to_ijk (ndom, n, &i, &j, &k);

  mdim_half = zdom[ndom].mdim / 2;
  is_lower_hemi = (j >= mdim_half);

  rmin = zdom[ndom].wind_x[i];
  rmax = zdom[ndom].wind_x[i + 1];

  sthetamin = sin (zdom[ndom].wind_z[j] / RADIAN);
  sthetamax = sin (zdom[ndom].wind_z[j + 1] / RADIAN);

  phi_lo = wmain[n].phi;
  phi_hi = wmain[n].phimax;

  inwind = W_NOT_INWIND;
  while (inwind != W_ALL_INWIND)
  {
    r = pow (rmin * rmin * rmin + random_number (0.0, 1.0) * (rmax * rmax * rmax - rmin * rmin * rmin), 1.0 / 3.0);
    phi = phi_lo + random_number (0.0, 1.0) * (phi_hi - phi_lo);

    s = sthetamin + random_number (0.0, 1.0) * (sthetamax - sthetamin);
    cos_theta = sqrt (1.0 - s * s);
    if (is_lower_hemi)
      cos_theta = -cos_theta;

    x[0] = r * cos (phi) * s;
    x[1] = r * sin (phi) * s;
    x[2] = r * cos_theta;

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
 * plasma index, scanning in the r direction.  Phi slices at the same
 * (i, j) share the same in-wind status so the loop is over (i, j, k).
 *
 **********************************************************/

int
sph3d_extend_density (int ndom, WindPtr w)
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
 * @brief  Distance to the far boundary of a SPH3D cell.
 *
 * @param [in] ndom   Domain number
 * @param [in] p      Photon pointer
 * @return  Distance to nearest cell boundary (r, theta, or phi).
 *
 * @details
 * Checks six boundaries:
 *   - inner/outer sphere at wmain[n].r and wmain[n].rmax
 *   - inner/outer theta cone at wmain[n].wcone and wmain[n].wcone_max
 *   - inner/outer phi half-plane at wmain[n].phi and wmain[n].phimax
 *     (skipped when pdim==1: full-azimuth axisymmetric model)
 *
 **********************************************************/

double
sph3d_ds_in_cell (int ndom, PhotPtr p)
{
  int n;
  double s, smax;

  if ((p->grid = n = where_in_grid (ndom, p->x)) < 0)
  {
    Error ("sph3d_ds_in_cell: Photon not in grid when routine entered\n");
    return (n);
  }

  /* Radial sphere boundaries */
  smax = ds_to_sphere (wmain[n].r, p);
  s = ds_to_sphere (wmain[n].rmax, p);
  if (s < smax)
    smax = s;

  /* Theta cone boundaries */
  s = ds_to_cone (&wmain[n].wcone, p);
  if (s < smax)
    smax = s;
  s = ds_to_cone (&wmain[n].wcone_max, p);
  if (s < smax)
    smax = s;

  /* Phi half-plane boundaries (skipped when pdim==1: full azimuth) */
  if (zdom[ndom].pdim > 1)
  {
    s = ds_phi_boundary (wmain[n].phi, p);
    if (s < smax)
      smax = s;
    s = ds_phi_boundary (wmain[n].phimax, p);
    if (s < smax)
      smax = s;
  }

  if (smax <= 0)
  {
    int i, j, k;
    wind_n_to_ijk (ndom, p->grid, &i, &j, &k);
    Error ("sph3d_ds_in_cell: smax %g <= 0 wind cell %i (%d,%d,%d) inwind %d\n", smax, p->grid, i, j, k, wmain[p->grid].inwind);
  }

  return (smax);
}
