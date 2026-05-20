/***********************************************************/
/** @file  import_sph3d.c
 * @author ksl
 * @date   May, 2025
 *
 * @brief  Routines to read an arbitrary wind model in 3D
 * spherical polar coordinates (SPH3D).
 *
 * The import file format is:
 *   i  j  k  inwind  r  theta  phi  v_x  v_y  v_z  rho  [t_e  [t_r]]
 *
 * where:
 *   r     is the inner radial boundary of the cell (cm)
 *   theta is the inner polar-angle boundary (degrees, 0-180)
 *   phi   is the inner azimuthal boundary (radians, 0-2pi)
 *   v_x, v_y, v_z are Cartesian velocity components (cm/s)
 *   rho   is the mass density (g/cm^3)
 *
 * Cells must be written in row-major order (k varies fastest, then j, then i).
 * The phi grid must be uniformly spaced; sph3d_wind_complete will recompute
 * phimax as (k+1)*2pi/pdim.
 *
 * Lines where sscanf returns fewer than READ_NO_TEMP_3D columns are skipped
 * (handles comment and header lines).
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief      Read a 3D spherical polar import file.
 *
 * @param [in] ndom       Domain number
 * @param [in] filename   Path to the import file
 * @return     Always returns 0
 *
 * @details
 * Two-pass reader.  Pass 1 counts valid data lines and records the
 * maximum i, j, k indices.  Arrays are then allocated to the exact
 * required size.  Pass 2 fills them.
 *
 * After reading, the 1-D boundary arrays (wind_x for r, wind_z for
 * theta in degrees, wind_phi for phi in radians) and their midpoint
 * arrays are built from cells at fixed index slices.
 *
 **********************************************************/

int
import_sph3d (int ndom, char *filename)
{
  FILE *fptr;
  char line[LINELENGTH];
  int icell, jcell, kcell, inwind, n, ncell;
  double r, theta, phi, v_x, v_y, v_z, rho, t_e, t_r;
  int imax, jmax, kmax;
  struct Import *m;

  Log ("Reading a model in 3D spherical polar coordinates from %s\n", filename);

  /* Pass 1: count lines and track maximum indices */

  if ((fptr = fopen (filename, "r")) == NULL)
  {
    Error ("import_sph3d: cannot open %s\n", filename);
    Exit (1);
  }

  ncell = 0;
  imax = jmax = kmax = 0;
  while (fgets (line, LINELENGTH, fptr) != NULL)
  {
    n = sscanf (line, " %d %d %d %d %le %le %le %le %le %le %le %le %le",
                &icell, &jcell, &kcell, &inwind, &r, &theta, &phi, &v_x, &v_y, &v_z, &rho, &t_e, &t_r);
    if (n < READ_NO_TEMP_3D)
      continue;
    if (icell > imax)
      imax = icell;
    if (jcell > jmax)
      jmax = jcell;
    if (kcell > kmax)
      kmax = kcell;
    ncell++;
  }
  fclose (fptr);

  if (ncell == 0)
  {
    Error ("import_sph3d: no valid data lines in %s\n", filename);
    Exit (1);
  }

  /* Allocate arrays now that ncell is known */

  m = &imported_model[ndom];
  m->ncell = ncell;
  m->ndim = imax + 1;
  m->mdim = jmax + 1;
  m->pdim = kmax + 1;

  m->i = calloc (ncell, sizeof (int));
  m->j = calloc (ncell, sizeof (int));
  m->k = calloc (ncell, sizeof (int));
  m->inwind = calloc (ncell, sizeof (int));
  m->r = calloc (ncell, sizeof (double));
  m->theta = calloc (ncell, sizeof (double));
  m->phi = calloc (ncell, sizeof (double));
  m->v_x = calloc (ncell, sizeof (double));
  m->v_y = calloc (ncell, sizeof (double));
  m->v_z = calloc (ncell, sizeof (double));
  m->mass_rho = calloc (ncell, sizeof (double));
  m->t_e = calloc (ncell, sizeof (double));
  m->t_r = calloc (ncell, sizeof (double));

  /* wind_x: r boundaries (ndim values); wind_z: theta boundaries in degrees (mdim values);
   * wind_phi: phi boundaries (pdim+1 values);  midpoint arrays similarly. */
  m->wind_x = calloc (m->ndim + 1, sizeof (double));
  m->wind_z = calloc (m->mdim + 1, sizeof (double));
  m->wind_midx = calloc (m->ndim, sizeof (double));
  m->wind_midz = calloc (m->mdim, sizeof (double));
  m->wind_phi = calloc (m->pdim + 1, sizeof (double));
  m->wind_midphi = calloc (m->pdim, sizeof (double));

  if (!m->i || !m->j || !m->k || !m->inwind || !m->r || !m->theta || !m->phi
      || !m->v_x || !m->v_y || !m->v_z || !m->mass_rho || !m->t_e || !m->t_r
      || !m->wind_x || !m->wind_z || !m->wind_midx || !m->wind_midz || !m->wind_phi || !m->wind_midphi)
  {
    Error ("import_sph3d: memory allocation failed\n");
    Exit (1);
  }

  /* Pass 2: fill arrays */

  if ((fptr = fopen (filename, "r")) == NULL)
  {
    Error ("import_sph3d: cannot reopen %s\n", filename);
    Exit (1);
  }

  n = 0;
  m->init_temperature = TRUE;
  while (fgets (line, LINELENGTH, fptr) != NULL)
  {
    int nf = sscanf (line, " %d %d %d %d %le %le %le %le %le %le %le %le %le",
                     &icell, &jcell, &kcell, &inwind, &r, &theta, &phi, &v_x, &v_y, &v_z, &rho, &t_e, &t_r);
    if (nf < READ_NO_TEMP_3D)
      continue;

    m->i[n] = icell;
    m->j[n] = jcell;
    m->k[n] = kcell;
    m->inwind[n] = inwind;
    m->r[n] = r;
    m->theta[n] = theta;
    m->phi[n] = phi;
    m->v_x[n] = v_x;
    m->v_y[n] = v_y;
    m->v_z[n] = v_z;
    m->mass_rho[n] = rho;

    if (nf == READ_ELECTRON_TEMP_3D)
    {
      m->init_temperature = FALSE;
      m->t_e[n] = t_e;
      m->t_r[n] = 1.1 * t_e;
    }
    else if (nf >= READ_BOTH_TEMP_3D)
    {
      m->init_temperature = FALSE;
      m->t_e[n] = t_e;
      m->t_r[n] = t_r;
    }

    n++;
  }
  fclose (fptr);

  if (n != ncell)
  {
    Error ("import_sph3d: expected %d cells but read %d on second pass\n", ncell, n);
    Exit (1);
  }
  if (ncell != m->ndim * m->mdim * m->pdim)
  {
    Error ("import_sph3d: incomplete grid: %d x %d x %d = %d but ncell = %d\n",
           m->ndim, m->mdim, m->pdim, m->ndim * m->mdim * m->pdim, ncell);
    Exit (1);
  }

  /* Build 1-D coordinate arrays from cells at fixed index slices */

  for (n = 0; n < ncell; n++)
  {
    /* wind_x[i]: r corner from j=0, k=0 cells */
    if (m->j[n] == 0 && m->k[n] == 0)
      m->wind_x[m->i[n]] = m->r[n];

    /* wind_z[j]: theta corner (degrees) from i=0, k=0 cells */
    if (m->i[n] == 0 && m->k[n] == 0)
      m->wind_z[m->j[n]] = m->theta[n];

    /* wind_phi[k]: phi corner from i=0, j=0 cells */
    if (m->i[n] == 0 && m->j[n] == 0)
      m->wind_phi[m->k[n]] = m->phi[n];
  }

  /* Closing boundaries */
  {
    double delta;
    int nd = m->ndim, md = m->mdim, pd = m->pdim;

    delta = m->wind_x[nd - 1] - m->wind_x[nd - 2];
    m->wind_x[nd] = m->wind_x[nd - 1] + delta;

    /* Theta closing boundary: for a full (0-180) grid this extrapolates to 180;
     * for an upper-hemisphere-only grid it extrapolates to 90+delta. */
    delta = m->wind_z[md - 1] - m->wind_z[md - 2];
    m->wind_z[md] = m->wind_z[md - 1] + delta;

    /* Phi must close at 2*pi */
    if (pd > 1)
      m->wind_phi[pd] = 2.0 * M_PI;
    else
      m->wind_phi[1] = 2.0 * M_PI;
  }

  /* Cell-centre arrays */
  for (n = 0; n < m->ndim; n++)
    m->wind_midx[n] = 0.5 * (m->wind_x[n] + m->wind_x[n + 1]);

  for (n = 0; n < m->mdim; n++)
    m->wind_midz[n] = 0.5 * (m->wind_z[n] + m->wind_z[n + 1]);

  for (n = 0; n < m->pdim; n++)
    m->wind_midphi[n] = 0.5 * (m->wind_phi[n] + m->wind_phi[n + 1]);

  /* Push dimensions into zdom */
  zdom[ndom].ndim = m->ndim;
  zdom[ndom].mdim = m->mdim;
  zdom[ndom].pdim = m->pdim;
  zdom[ndom].ndim2 = m->ndim * m->mdim * m->pdim;

  Log ("import_sph3d: ndim=%d mdim=%d pdim=%d ncell=%d\n", m->ndim, m->mdim, m->pdim, m->ncell);

  return (0);
}


/**********************************************************/
/**
 * @brief      Set domain boundary parameters from a SPH3D imported model.
 *
 * @param [in] ndom   Domain number
 * @return     Always returns 0
 *
 * @details
 * Copies 1-D coordinate arrays into zdom and determines the spatial
 * extent of cells marked as in-wind.
 *
 **********************************************************/

int
import_sph3d_setup_boundaries (int ndom)
{
  int n;
  struct Import *m = &imported_model[ndom];
  double rmin, rmax, rho_min, rho_max, zmax;
  double rho, z_abs;

  /* Copy 1-D grids into zdom */
  for (n = 0; n < m->ndim; n++)
    zdom[ndom].wind_x[n] = m->wind_x[n];

  for (n = 0; n < m->mdim; n++)
    zdom[ndom].wind_z[n] = m->wind_z[n];

  for (n = 0; n <= m->pdim; n++)
    zdom[ndom].wind_phi[n] = m->wind_phi[n];

  for (n = 0; n < m->pdim; n++)
    zdom[ndom].wind_midphi[n] = m->wind_midphi[n];

  /* Determine spatial extent from in-wind cells */
  rmax = rho_max = zmax = 0.0;
  rmin = rho_min = VERY_BIG;

  for (n = 0; n < m->ncell; n++)
  {
    if (m->inwind[n] < 0)
      continue;

    double r_lo = m->r[n];
    double r_hi = m->wind_x[m->i[n] + 1];
    double theta_lo_rad = m->theta[n] / RADIAN;
    double theta_hi_rad = m->wind_z[m->j[n] + 1] / RADIAN;

    rho = r_hi * sin (theta_hi_rad);
    z_abs = fabs (r_hi * cos (theta_lo_rad));

    if (r_lo < rmin)
      rmin = r_lo;
    if (r_hi > rmax)
      rmax = r_hi;
    if (rho > rho_max)
      rho_max = rho;
    if (z_abs > zmax)
      zmax = z_abs;
  }

  if (rmin == VERY_BIG)
    rmin = m->wind_x[0];

  rho_min = rmin * sin (m->wind_z[0] / RADIAN);
  if (rho_min < 0.0)
    rho_min = 0.0;

  zdom[ndom].rmin = rmin;
  zdom[ndom].rmax = rmax;
  zdom[ndom].wind_rhomin_at_disk = rho_min;
  zdom[ndom].wind_rhomax_at_disk = rho_max;
  zdom[ndom].zmax = zmax;
  zdom[ndom].wind_thetamin = zdom[ndom].wind_thetamax = 0.0;

  return (0);
}


/**********************************************************/
/**
 * @brief      Populate wmain from a SPH3D imported model.
 *
 * @param [in] w      The wind array
 * @param [in] ndom   Domain number
 * @return     Always returns 0
 *
 * @details
 * Maps each imported cell (i, j, k) to its flat wind index via
 * wind_ijk_to_n and fills positions, velocities, and inwind flag.
 * Cell centres are taken from the pre-computed wind_mid* arrays.
 * The rmax, thetamax, wcone, and wcone_max fields are filled later
 * by sph3d_wind_complete.
 *
 **********************************************************/

int
sph3d_make_grid_import (WindPtr w, int ndom)
{
  int n, nn;
  double theta_rad, thetacen_rad;
  struct Import *m = &imported_model[ndom];

  for (n = 0; n < m->ncell; n++)
  {
    wind_ijk_to_n (ndom, m->i[n], m->j[n], m->k[n], &nn);

    w[nn].r = m->r[n];
    w[nn].theta = m->theta[n];
    w[nn].phi = m->phi[n];

    theta_rad = m->theta[n] / RADIAN;
    w[nn].x[0] = m->r[n] * sin (theta_rad);
    w[nn].x[1] = 0.0;
    w[nn].x[2] = m->r[n] * cos (theta_rad);

    w[nn].rcen = m->wind_midx[m->i[n]];
    w[nn].thetacen = m->wind_midz[m->j[n]];
    w[nn].phicen = m->wind_midphi[m->k[n]];

    thetacen_rad = w[nn].thetacen / RADIAN;
    w[nn].xcen[0] = w[nn].rcen * sin (thetacen_rad);
    w[nn].xcen[1] = 0.0;
    w[nn].xcen[2] = w[nn].rcen * cos (thetacen_rad);

    w[nn].v[0] = m->v_x[n];
    w[nn].v[1] = m->v_y[n];
    w[nn].v[2] = m->v_z[n];

    w[nn].inwind = m->inwind[n];
    if (w[nn].inwind == W_NOT_INWIND || w[nn].inwind == W_PART_INWIND)
      w[nn].inwind = W_IGNORE;
  }

  /* Copy 1-D grids into zdom (sph3d_wind_complete will also do this,
     but import_make_grid is called before wind_complete). */
  for (n = 0; n < zdom[ndom].ndim; n++)
    zdom[ndom].wind_x[n] = m->wind_x[n];

  for (n = 0; n < zdom[ndom].mdim; n++)
    zdom[ndom].wind_z[n] = m->wind_z[n];

  for (n = 0; n <= zdom[ndom].pdim; n++)
    zdom[ndom].wind_phi[n] = m->wind_phi[n];

  for (n = 0; n < zdom[ndom].pdim; n++)
    zdom[ndom].wind_midphi[n] = m->wind_midphi[n];

  return (0);
}


/**********************************************************/
/**
 * @brief      Velocity at an arbitrary position in a SPH3D imported model.
 *
 * @param [in] ndom   Domain number
 * @param [in] x      Position (3-vector, Cartesian)
 * @param [out] v     Velocity (3-vector, Cartesian cm/s)
 * @return     Speed at x
 *
 * @details
 * Uses coord_fraction for trilinear interpolation weights across the
 * 8 corners of the SPH3D cell, then combines corner velocities from wmain.
 *
 **********************************************************/

double
velocity_sph3d (int ndom, double *x, double *v)
{
  int j, nn;
  int nnn[8], nelem;
  double frac[8];
  double vv[3];
  double speed;

  coord_fraction (ndom, 0, x, nnn, frac, &nelem);

  for (j = 0; j < 3; j++)
  {
    vv[j] = 0;
    for (nn = 0; nn < nelem; nn++)
      vv[j] += wmain[zdom[ndom].nstart + nnn[nn]].v[j] * frac[nn];
  }

  speed = length (vv);
  v[0] = vv[0];
  v[1] = vv[1];
  v[2] = vv[2];

  return (speed);
}


/**********************************************************/
/**
 * @brief      Density at a position in a SPH3D imported model.
 *
 * @param [in] ndom   Domain number
 * @param [in] x      Position (3-vector, Cartesian)
 * @return     Mass density (g/cm^3) at x
 *
 * @details
 * Converts x to (r, theta, phi) and binary-searches the imported
 * 1-D grids to find the cell, returning its stored density.
 * No interpolation: density is defined at cell centres.
 *
 **********************************************************/

double
rho_sph3d (int ndom, double *x)
{
  struct Import *m = &imported_model[ndom];
  double r_pos, theta_pos, phi_pos;
  int i, j, k, n;

  r_pos = length (x);
  theta_pos = (r_pos > 0.0) ? acos (x[2] / r_pos) * RADIAN : 0.0;
  phi_pos = atan2 (x[1], x[0]);
  if (phi_pos < 0.0)
    phi_pos += 2.0 * M_PI;

  i = 0;
  while (i < m->ndim - 1 && r_pos >= m->wind_x[i + 1])
    i++;

  j = 0;
  while (j < m->mdim - 1 && theta_pos >= m->wind_z[j + 1])
    j++;

  k = 0;
  while (k < m->pdim - 1 && phi_pos >= m->wind_phi[k + 1])
    k++;

  n = (i * m->mdim + j) * m->pdim + k;

  return (m->mass_rho[n]);
}


/**********************************************************/
/**
 * @brief      Temperature at a position in a SPH3D imported model.
 *
 * @param [in] ndom        Domain number
 * @param [in] x           Position (3-vector, Cartesian)
 * @param [in] return_t_e  If TRUE return electron temperature, else radiation
 * @return     Temperature in Kelvin
 *
 **********************************************************/

double
temperature_sph3d (int ndom, double *x, int return_t_e)
{
  struct Import *m = &imported_model[ndom];
  double r_pos, theta_pos, phi_pos;
  int i, j, k, n;

  if (m->init_temperature)
    return (return_t_e ? 1.1 * zdom[ndom].twind : zdom[ndom].twind);

  r_pos = length (x);
  theta_pos = (r_pos > 0.0) ? acos (x[2] / r_pos) * RADIAN : 0.0;
  phi_pos = atan2 (x[1], x[0]);
  if (phi_pos < 0.0)
    phi_pos += 2.0 * M_PI;

  i = 0;
  while (i < m->ndim - 1 && r_pos >= m->wind_x[i + 1])
    i++;

  j = 0;
  while (j < m->mdim - 1 && theta_pos >= m->wind_z[j + 1])
    j++;

  k = 0;
  while (k < m->pdim - 1 && phi_pos >= m->wind_phi[k + 1])
    k++;

  n = (i * m->mdim + j) * m->pdim + k;

  return (return_t_e ? m->t_e[n] : m->t_r[n]);
}
