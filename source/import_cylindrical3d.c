/***********************************************************/
/** @file  import_cylindrical3d.c
 * @author ksl
 * @date   May, 2025
 *
 * @brief  Routines to read in an arbitrary wind model in
 * 3D cylindrical coordinates (CYLIND3D).
 *
 * The import file format is:
 *   i  j  k  inwind  x  z  phi  v_x  v_y  v_z  rho  [t_e  [t_r]]
 *
 * where x, z, phi are the cell corner coordinates (lower rho/z/phi
 * corner of the cell), and rho is the observer-frame mass density.
 * Cells must be in row-major order (k varies fastest, then j, then i).
 * The phi grid may be non-uniformly spaced but must be strictly increasing
 * in k order and span 0 to 2pi.
 * Lines where sscanf returns fewer than READ_NO_TEMP_3D values are
 * silently skipped (handles comment/header lines).
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"



/**********************************************************/
/**
 * @brief      Read a 3D cylindrical import file.
 *
 * @param [in] ndom       Domain number
 * @param [in] filename   Path to the import file
 * @return     Always returns 0
 *
 * @details
 * Two-pass reader.  Pass 1 counts valid lines and determines the
 * maximum i, j, k indices seen.  Arrays are then allocated to the
 * exact required size.  Pass 2 fills the arrays.
 *
 * After reading, the phi boundary array wind_phi and cell-centre
 * array wind_midphi are built from the unique phi corner values
 * found at i=0, j=0 (one entry per phi slice).
 *
 **********************************************************/

int
import_cylindrical3d (int ndom, char *filename)
{
  FILE *fptr;
  char line[LINELENGTH];
  int icell, jcell, kcell, inwind, n, ncell;
  double x, z, phi, v_x, v_y, v_z, rho, t_e, t_r;
  int imax, jmax, kmax;
  struct Import *m;

  Log ("Reading a model in 3D cylindrical coordinates from %s\n", filename);

  /* Pass 1: count lines and track max indices */

  if ((fptr = fopen (filename, "r")) == NULL)
  {
    Error ("import_cylindrical3d: cannot open %s\n", filename);
    Exit (1);
  }

  ncell = 0;
  imax = jmax = kmax = 0;
  while (fgets (line, LINELENGTH, fptr) != NULL)
  {
    n = sscanf (line, " %d %d %d %d %le %le %le %le %le %le %le %le %le",
                &icell, &jcell, &kcell, &inwind, &x, &z, &phi, &v_x, &v_y, &v_z, &rho, &t_e, &t_r);
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
    Error ("import_cylindrical3d: no valid data lines in %s\n", filename);
    Exit (1);
  }

  /* Allocate arrays now that we know ncell */

  m = &imported_model[ndom];
  m->ncell = ncell;
  m->ndim = imax + 1;
  m->mdim = jmax + 1;
  m->pdim = kmax + 1;

  m->i = calloc (ncell, sizeof (int));
  m->j = calloc (ncell, sizeof (int));
  m->k = calloc (ncell, sizeof (int));
  m->inwind = calloc (ncell, sizeof (int));
  m->x = calloc (ncell, sizeof (double));
  m->z = calloc (ncell, sizeof (double));
  m->phi = calloc (ncell, sizeof (double));
  m->v_x = calloc (ncell, sizeof (double));
  m->v_y = calloc (ncell, sizeof (double));
  m->v_z = calloc (ncell, sizeof (double));
  m->mass_rho = calloc (ncell, sizeof (double));
  m->t_e = calloc (ncell, sizeof (double));
  m->t_r = calloc (ncell, sizeof (double));

  /* wind_x: ndim entries; wind_z: mdim entries; wind_phi: pdim+1 boundaries;
   * wind_midx/midz/midphi: ndim/mdim/pdim centres */
  m->wind_x = calloc (m->ndim + 1, sizeof (double));
  m->wind_z = calloc (m->mdim + 1, sizeof (double));
  m->wind_midx = calloc (m->ndim, sizeof (double));
  m->wind_midz = calloc (m->mdim, sizeof (double));
  m->wind_phi = calloc (m->pdim + 1, sizeof (double));
  m->wind_midphi = calloc (m->pdim, sizeof (double));

  if (!m->i || !m->j || !m->k || !m->inwind || !m->x || !m->z || !m->phi
      || !m->v_x || !m->v_y || !m->v_z || !m->mass_rho || !m->t_e || !m->t_r
      || !m->wind_x || !m->wind_z || !m->wind_midx || !m->wind_midz || !m->wind_phi || !m->wind_midphi)
  {
    Error ("import_cylindrical3d: memory allocation failed\n");
    Exit (1);
  }

  /* Pass 2: fill arrays */

  if ((fptr = fopen (filename, "r")) == NULL)
  {
    Error ("import_cylindrical3d: cannot reopen %s\n", filename);
    Exit (1);
  }

  n = 0;
  m->init_temperature = TRUE;
  while (fgets (line, LINELENGTH, fptr) != NULL)
  {
    int nf = sscanf (line, " %d %d %d %d %le %le %le %le %le %le %le %le %le",
                     &icell, &jcell, &kcell, &inwind, &x, &z, &phi, &v_x, &v_y, &v_z, &rho, &t_e, &t_r);
    if (nf < READ_NO_TEMP_3D)
      continue;

    m->i[n] = icell;
    m->j[n] = jcell;
    m->k[n] = kcell;
    m->inwind[n] = inwind;
    m->x[n] = x;
    m->z[n] = z;
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

  /* Validate completeness */
  if (n != ncell)
  {
    Error ("import_cylindrical3d: expected %d cells but read %d on second pass\n", ncell, n);
    Exit (1);
  }
  if (ncell != m->ndim * m->mdim * m->pdim)
  {
    Error ("import_cylindrical3d: grid is incomplete: %d x %d x %d = %d but ncell = %d\n",
           m->ndim, m->mdim, m->pdim, m->ndim * m->mdim * m->pdim, ncell);
    Exit (1);
  }

  /* Build 1-D coordinate arrays from cells with the appropriate fixed indices */

  for (n = 0; n < ncell; n++)
  {
    /* wind_x[i]: rho corner, constant j=0 k=0 */
    if (m->j[n] == 0 && m->k[n] == 0)
      m->wind_x[m->i[n]] = m->x[n];

    /* wind_z[j]: z corner, constant i=0 k=0 */
    if (m->i[n] == 0 && m->k[n] == 0)
      m->wind_z[m->j[n]] = m->z[n];

    /* wind_phi[k]: phi corner, constant i=0 j=0 */
    if (m->i[n] == 0 && m->j[n] == 0)
      m->wind_phi[m->k[n]] = m->phi[n];
  }

  /* Closing boundary: extrapolate one cell beyond last corner */
  {
    double delta;
    int nd = m->ndim, md = m->mdim, pd = m->pdim;

    delta = m->wind_x[nd - 1] - m->wind_x[nd - 2];
    m->wind_x[nd] = m->wind_x[nd - 1] + delta;

    delta = m->wind_z[md - 1] - m->wind_z[md - 2];
    m->wind_z[md] = m->wind_z[md - 1] + delta;

    /* For phi the grid must close at 2*pi; the import file stores only the
     * lower boundary of each cell, so add 2*pi as the upper edge of the last
     * cell.  If pdim==1 (full torus collapsed to one phi cell) set 2*pi
     * directly. */
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

  /* Set ndim2 = total grid cells for domain */
  zdom[ndom].ndim = m->ndim;
  zdom[ndom].mdim = m->mdim;
  zdom[ndom].pdim = m->pdim;
  zdom[ndom].ndim2 = m->ndim * m->mdim * m->pdim;

  Log ("import_cylindrical3d: ndim=%d mdim=%d pdim=%d ncell=%d\n",
       m->ndim, m->mdim, m->pdim, m->ncell);

  return (0);
}


/**********************************************************/
/**
 * @brief      Set domain boundary parameters from a CYLIND3D imported model.
 *
 * @param [in] ndom   Domain number
 * @return     Always returns 0
 *
 * @details
 * Copies coordinate arrays into zdom and determines the spatial
 * extent of the wind-cells that are actually in the wind.
 *
 **********************************************************/

int
import_cylindrical3d_setup_boundaries (int ndom)
{
  int n;
  struct Import *m = &imported_model[ndom];
  double rmin, rmax, rho_min, rho_max, zmin, zmax, abs_zmax;

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
  rmax = rho_max = zmax = abs_zmax = 0;
  rmin = rho_min = zmin = VERY_BIG;

  for (n = 0; n < m->ncell; n++)
  {
    if (m->inwind[n] < 0)
      continue;

    double xn = m->x[n];
    double zn = m->z[n];
    double xnext = m->wind_x[m->i[n] + 1];
    double znext = m->wind_z[m->j[n] + 1];

    double r_inner = sqrt (xn * xn + zn * zn);
    double r_outer = sqrt (xnext * xnext + znext * znext);

    if (xnext > rho_max)
      rho_max = xnext;
    if (xn < rho_min)
      rho_min = xn;
    if (znext > zmax)
      zmax = znext;
    if (zn < zmin)
      zmin = zn;
    if (fabs (znext) > abs_zmax)
      abs_zmax = fabs (znext);
    if (fabs (zn) > abs_zmax)
      abs_zmax = fabs (zn);
    if (r_outer > rmax)
      rmax = r_outer;
    if (r_inner < rmin)
      rmin = r_inner;
  }

  zdom[ndom].wind_rhomin_at_disk = rho_min;
  zdom[ndom].wind_rhomax_at_disk = rho_max;
  zdom[ndom].zmax = abs_zmax;
  zdom[ndom].zmin = zmin;
  zdom[ndom].rmax = rmax;
  zdom[ndom].rmin = rmin;
  zdom[ndom].wind_thetamin = zdom[ndom].wind_thetamax = 0.0;

  zdom[ndom].windplane[0].x[0] = zdom[ndom].windplane[0].x[1] = 0;
  zdom[ndom].windplane[0].x[2] = zmin;
  zdom[ndom].windplane[0].lmn[0] = zdom[ndom].windplane[0].lmn[1] = 0;
  zdom[ndom].windplane[0].lmn[2] = 1;

  zdom[ndom].windplane[1].x[0] = zdom[ndom].windplane[1].x[1] = 0;
  zdom[ndom].windplane[1].x[2] = zmax;
  zdom[ndom].windplane[1].lmn[0] = zdom[ndom].windplane[1].lmn[1] = 0;
  zdom[ndom].windplane[1].lmn[2] = 1;

  return (0);
}


/**********************************************************/
/**
 * @brief      Populate wmain from a CYLIND3D imported model.
 *
 * @param [in] w      The wind array
 * @param [in] ndom   Domain number
 * @return     Always returns 0
 *
 * @details
 * Maps each imported cell (i, j, k) to its flat wind index via
 * wind_ijk_to_n and fills positions, velocities, and inwind flag.
 * Cell centres are taken from the pre-computed wind_mid* arrays.
 *
 **********************************************************/

int
cylindrical3d_make_grid_import (WindPtr w, int ndom)
{
  int n, nn;
  struct Import *m = &imported_model[ndom];

  for (n = 0; n < m->ncell; n++)
  {
    wind_ijk_to_n (ndom, m->i[n], m->j[n], m->k[n], &nn);

    w[nn].x[0] = m->x[n];
    w[nn].x[1] = 0.0;
    w[nn].x[2] = m->z[n];
    w[nn].phi = m->phi[n];
    w[nn].phimax = m->wind_phi[m->k[n] + 1];

    w[nn].xcen[0] = m->wind_midx[m->i[n]];
    w[nn].xcen[1] = 0.0;
    w[nn].xcen[2] = m->wind_midz[m->j[n]];
    w[nn].phicen = m->wind_midphi[m->k[n]];

    w[nn].v[0] = m->v_x[n];
    w[nn].v[1] = m->v_y[n];
    w[nn].v[2] = m->v_z[n];

    w[nn].inwind = m->inwind[n];
    if (w[nn].inwind == W_NOT_INWIND || w[nn].inwind == W_PART_INWIND)
      w[nn].inwind = W_IGNORE;
  }

  return (0);
}


/**********************************************************/
/**
 * @brief      Velocity at an arbitrary position in a CYLIND3D imported model.
 *
 * @param [in] ndom   Domain number
 * @param [in] x      Position (3-vector)
 * @param [out] v     Velocity (3-vector) in Cartesian coordinates
 * @return     Speed at x
 *
 * @details
 * Delegates to coord_fraction for bilinear/trilinear interpolation
 * weights, then combines the corner velocities from wmain.
 *
 **********************************************************/

double
velocity_cylindrical3d (int ndom, double *x, double *v)
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
 * @brief      Density at a position in a CYLIND3D imported model.
 *
 * @param [in] ndom   Domain number
 * @param [in] x      Position (3-vector)
 * @return     Mass density (cgs) at x
 *
 * @details
 * Finds the (i, j, k) cell containing x by binary search in the
 * rho, z, and phi axes and returns the stored density for that cell.
 * No interpolation: density is defined at cell centres.
 *
 **********************************************************/

double
rho_cylindrical3d (int ndom, double *x)
{
  struct Import *m = &imported_model[ndom];
  double rho_pos, z, phi_pos;
  int i, j, k, n;

  rho_pos = sqrt (x[0] * x[0] + x[1] * x[1]);
  z = x[2];
  phi_pos = atan2 (x[1], x[0]);
  if (phi_pos < 0.0)
    phi_pos += 2.0 * M_PI;

  /* Locate rho index */
  i = 0;
  while (i < m->ndim - 1 && rho_pos >= m->wind_x[i + 1])
    i++;

  /* Locate z index; for two-hemisphere grids z can be negative */
  j = 0;
  while (j < m->mdim - 1 && z >= m->wind_z[j + 1])
    j++;

  /* Locate phi index */
  k = 0;
  while (k < m->pdim - 1 && phi_pos >= m->wind_phi[k + 1])
    k++;

  n = (i * m->mdim + j) * m->pdim + k;

  return (m->mass_rho[n]);
}


/**********************************************************/
/**
 * @brief      Temperature at a position in a CYLIND3D imported model.
 *
 * @param [in] ndom        Domain number
 * @param [in] x           Position (3-vector)
 * @param [in] return_t_e  If TRUE return electron temperature, else radiation
 * @return     Temperature in Kelvin
 *
 **********************************************************/

double
temperature_cylindrical3d (int ndom, double *x, int return_t_e)
{
  struct Import *m = &imported_model[ndom];
  double rho_pos, z, phi_pos;
  int i, j, k, n;

  if (m->init_temperature)
  {
    return (return_t_e ? 1.1 * zdom[ndom].twind : zdom[ndom].twind);
  }

  rho_pos = sqrt (x[0] * x[0] + x[1] * x[1]);
  z = x[2];
  phi_pos = atan2 (x[1], x[0]);
  if (phi_pos < 0.0)
    phi_pos += 2.0 * M_PI;

  i = 0;
  while (i < m->ndim - 1 && rho_pos >= m->wind_x[i + 1])
    i++;

  j = 0;
  while (j < m->mdim - 1 && z >= m->wind_z[j + 1])
    j++;

  k = 0;
  while (k < m->pdim - 1 && phi_pos >= m->wind_phi[k + 1])
    k++;

  n = (i * m->mdim + j) * m->pdim + k;

  return (return_t_e ? m->t_e[n] : m->t_r[n]);
}
