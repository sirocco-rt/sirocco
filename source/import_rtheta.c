/***********************************************************/
/** @file  import_rtheta.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief
 * These are general routines to read in a model that
 * is in polar or rtheta coordinates.
 * ###Notes###
 *
 * There a various possibilities for how the
 * velocities could be entered.  One possibility
 * which is the way the zeus_sirocco models work
 * is for the velocity to be given in spherical
 * polar coordinates.

 * However, internally, sirocco uses xyz coordinates
 * for velocities (as measured in the xz plane),
 * and that is the model followed, here.  This also
 * makes these routines similar to those used
 * in imported cylindrical models.

 * This means that if the user provides a model
 * where velocities are in spherical polar coordinates
 * then one must translate them to the convention
 * here before the model is read in
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief      Read the an arbitrary wind model in polar coordinates
 *
 * @param [in] int  ndom   The domain for the imported model
 * @param [in] char *  filename   The file containing the model to import
 * @return   Always returns 0
 *
 * @details
 * This routine reads the data into a set of arrays.  It's
 * only purpose is to read in the data
 *
 * ### Notes ###
 * The basic data we need to read in are
 *
 * icell, jcell, r theta inwind v_x v_y v_z  rho (and optionally T)
 *
 * where
 *
 * * r is the radial coordinates
 * * theta is the angular coordinate measured from the z axis
 * * v_x, v_y, and v_z is the velocity in cartesian coordinates
 *      as measured in the x,z plane
 * * rho is the density in cgs units
 * * inwind defines whether or not a particular cell is actually
 * in the wind
 *
 * Guard cells are required at the outer boundaries. 
 *
 * This routine assumes the same conventions as used elsewhere
 * in Python, that is that the positions and velocities are given
 * at the edges of a cell, but that rho is given at the center.
 *
 **********************************************************/

int
import_rtheta (ndom, filename)
     int ndom;
     char *filename;
{
  FILE *fptr;
  char line[LINELENGTH];
  int n, icell, jcell, ncell, inwind;
  int jz, jx;
  double delta;
  double r, theta, v_x, v_y, v_z, rho, t_r, t_e;
  double mid_z;


  Log ("Reading a model %s in polar (r,theta) coordinates \n", filename);


  if ((fptr = fopen (filename, "r")) == NULL)
  {
    Error ("import_rtheta: No such file\n");
    Exit (0);
  }


  ncell = 0;
  while (fgets (line, LINELENGTH, fptr) != NULL)
  {
    n = sscanf (line, " %d %d %d %le %le %le %le %le %le %le %le", &icell, &jcell, &inwind, &r, &theta, &v_x, &v_y, &v_z, &rho, &t_e, &t_r);

    if (n < READ_NO_TEMP_2D)
    {
      continue;
    }
    else
    {
      imported_model[ndom].i[ncell] = icell;
      imported_model[ndom].j[ncell] = jcell;
      imported_model[ndom].inwind[ncell] = inwind;
      imported_model[ndom].r[ncell] = r;
      imported_model[ndom].theta[ncell] = theta;
      imported_model[ndom].v_x[ncell] = v_x;
      imported_model[ndom].v_y[ncell] = v_y;
      imported_model[ndom].v_z[ncell] = v_z;
      imported_model[ndom].mass_rho[ncell] = rho;

      if (n == READ_ELECTRON_TEMP_2D)
      {
        imported_model[ndom].init_temperature = FALSE;
        imported_model[ndom].t_e[ncell] = t_e;
        imported_model[ndom].t_r[ncell] = 1.1 * t_e;
      }
      else if (n == READ_BOTH_TEMP_2D)
      {
        imported_model[ndom].init_temperature = FALSE;
        imported_model[ndom].t_e[ncell] = t_e;
        imported_model[ndom].t_r[ncell] = t_r;
      }
      else
      {
        imported_model[ndom].init_temperature = TRUE;
      }

      ncell++;

      if (ncell > NDIM_MAX2D)
      {
        Error ("%s : %i : trying to read in more grid points than allowed (%i). Try changing NDIM_MAX and recompiling.\n", __FILE__,
               __LINE__, NDIM_MAX2D);
        Exit (1);
      }

    }
  }

  /* Set and check the dimensions of the grids to be set up.
   * Note that some assumptions are built into the way the grid
   * is read in, most notably that the last cell to be read in
   * defines the dimensions of the entire grid.
   */

  zdom[ndom].ndim = imported_model[ndom].ndim = icell + 1;
  zdom[ndom].mdim = imported_model[ndom].mdim = jcell + 1;
  imported_model[ndom].ncell = ncell;
  zdom[ndom].ndim2 = zdom[ndom].ndim * zdom[ndom].mdim;

  /* Check that the number of cells read in matches the number that was expected */

  if (ncell != imported_model[ndom].ndim * imported_model[ndom].mdim)
  {
    Error ("import_rtheta: The dimensions of the imported grid seem wrong % d x %d != %d\n", imported_model[ndom].ndim,
           imported_model[ndom].mdim, imported_model[ndom].ncell);
    exit (1);
  }


  jz = jx = 0;
  for (n = 0; n < imported_model[ndom].ncell; n++)
  {
    if (imported_model[ndom].i[n] == 0)
    {
      imported_model[ndom].wind_z[jz] = imported_model[ndom].theta[n];
      jz++;
    }
    if (imported_model[ndom].j[n] == 0)
    {
      imported_model[ndom].wind_x[jx] = imported_model[ndom].r[n];
      jx++;
    }
  }

  for (n = 0; n < jz - 1; n++)
  {
    mid_z = 0.5 * (imported_model[ndom].wind_z[n] + imported_model[ndom].wind_z[n + 1]);
    if (mid_z > 90)
      mid_z = 90;
    imported_model[ndom].wind_midz[n] = mid_z;
  }

  delta = (imported_model[ndom].wind_z[n - 1] - imported_model[ndom].wind_z[n - 2]);
  imported_model[ndom].wind_midz[n] = imported_model[ndom].wind_z[n - 1] + 0.5 * delta;

  for (n = 0; n < jx - 1; n++)
  {
    imported_model[ndom].wind_midx[n] = 0.5 * (imported_model[ndom].wind_x[n] + imported_model[ndom].wind_x[n + 1]);
  }


  delta = (imported_model[ndom].wind_x[n - 1] - imported_model[ndom].wind_x[n - 2]);
  imported_model[ndom].wind_midx[n] = imported_model[ndom].wind_x[n - 1] + 0.5 * delta;

  return (0);
}




/* ************************************************************************** */
/**
 * @brief   Set up the various domain boundaries for a rtheta coordinate
 *          system
 *
 * @param[in] int ndom         The domain of interest
 *
 * @return    Always returns 0
 *
 * @details
 *
 * This used to be contained within rtheta_make_grid_import, however, it
 * does not reply on any of the variables in that function and only relies on
 * the imported_model struct. Therefore, the boundary setup was moved into a
 * separate function so it could be done else where in the program flow.
 *
 * ************************************************************************** */

int
import_rtheta_setup_boundaries (int ndom)
{
  int n_inner, n_outer;
  double rho_max, rho_min, rho_inner, rho_outer, rmin, rmax;
  double zmin, zmax;
  double cell_inner[3], cell_outer[3];
  double x_inner, z_inner_next;
  double z_outer, x_outer_next;


  /* Now set up wind boundaries so they are harmless.
   * Note that the grid goes from near the pole
   * to the equator
   */

  rho_inner = rho_outer = 0;
  rmax = rho_max = zmax = 0;
  rmin = rho_min = zmin = VERY_BIG;

  for (n_inner = 0; n_inner < imported_model[ndom].ncell; n_inner++)
  {
    // Only update wind boundaries taking into account cells which are inwind
    if (imported_model[ndom].inwind[n_inner] >= 0)
    {
      n_outer = n_inner + imported_model[ndom].mdim;

      if (n_outer + 1 >= zdom[ndom].ndim2)
      {
        Error ("import_rtheta_setup_boundaries: attempting to access cell outside of grid");
        Exit (1);               // Something has gone wrong programmatically, best thing to do imo is exit
      }

      /*
       * Calculate the x, y, and z coordinates for the inner cell and the
       * next inner cell
       */

      x_inner = cell_inner[0] = imported_model[ndom].r[n_inner] * sin (imported_model[ndom].theta[n_inner] / RADIAN);
      cell_inner[1] = 0;
      cell_inner[2] = imported_model[ndom].r[n_inner] * cos (imported_model[ndom].theta[n_inner] / RADIAN);
      rho_inner = length (cell_inner);

      z_inner_next = imported_model[ndom].r[n_inner + 1] * cos (imported_model[ndom].theta[n_inner + 1] / RADIAN);

      /*
       * Calculate the x, y and z coordinates for the outer side of the cell
       * and for the next inner/outer cell
       */

      cell_outer[0] = imported_model[ndom].r[n_outer] * sin (imported_model[ndom].theta[n_outer] / RADIAN);
      cell_outer[1] = 0;
      z_outer = cell_outer[2] = imported_model[ndom].r[n_outer] * cos (imported_model[ndom].theta[n_outer] / RADIAN);

      if (n_outer < zdom[ndom].ndim2)
        rho_outer = length (cell_outer);

      x_outer_next = imported_model[ndom].r[n_outer + 1] * sin (imported_model[ndom].theta[n_outer + 1] / RADIAN);

      /*
       * Now we need to check if this cell is within the currently set boundaries.
       * If it is not, then we must update the boundaries to accommodate these
       * cells as they have been marked inwind, and hence should be inwind when
       * read into Python.
       */

      if (x_outer_next > rho_max)
      {
        rho_max = x_outer_next;
      }

      if (z_outer > zmax)
      {
        zmax = z_outer;
      }

      if (z_inner_next < zmin && z_inner_next > 0)
      {
        zmin = z_inner_next;
      }

      if (rho_outer > rmax)
      {
        rmax = rho_outer;
      }

      if (rho_min > x_inner)
      {
        rho_min = x_inner;
      }

      if (rmin > rho_inner)
      {
        rmin = rho_inner;
      }
    }
  }

  zdom[ndom].wind_rhomin_at_disk = rho_min;
  zdom[ndom].wind_rhomax_at_disk = rho_max;
  zdom[ndom].zmax = zmax;

  zdom[ndom].rmax = rmax;
  zdom[ndom].rmin = rmin;
  zdom[ndom].wind_thetamin = zdom[ndom].wind_thetamax = 0.;

  return 0;
}




/**********************************************************/
/**
 * @brief      Use the imported data to initialize various
 * portions of the Wind and Domain structures
 *
 * @param [out] WindPtr  w   The wind structure
 * @param [in] int  ndom   The domain for the imported model
 * @return     Always returns 0
 *
 * @details
 * This routine initializes the portions of the wind structure
 * using the imported model, specifically those portions having
 * to do with positions, and velocities.
 *
 * ### Notes ###
 * The routine also initializes wind_x and wind_z in the domain
 * structure, ans sets up wind_cones and other boundaries
 * intended to bound the wind.
 *
 **********************************************************/

int
rtheta_make_grid_import (w, ndom)
     WindPtr w;
     int ndom;
{
  int n, nn;
  double theta;

  /* As in the case of other models we assume that the grid has been
   * read in correctly and so now that the WindPtrs have been generated
   * we can put a lot of information in directly
   */

  for (n = 0; n < imported_model[ndom].ncell; n++)
  {
    wind_ij_to_n (ndom, imported_model[ndom].i[n], imported_model[ndom].j[n], &nn);
    w[nn].r = imported_model[ndom].r[n];
    w[nn].theta = theta = imported_model[ndom].theta[n];

    theta /= RADIAN;

    w[nn].x[0] = w[nn].r * sin (theta);
    w[nn].x[1] = 0;
    w[nn].x[2] = w[nn].r * cos (theta);
    w[nn].v[0] = imported_model[ndom].v_x[n];
    w[nn].v[1] = imported_model[ndom].v_y[n];
    w[nn].v[2] = imported_model[ndom].v_z[n];
    w[nn].inwind = imported_model[ndom].inwind[n];

    w[nn].thetacen = imported_model[ndom].wind_midz[imported_model[ndom].j[n]];
    if (w[nn].thetacen > 90)
      w[nn].thetacen = 90;
    theta = w[nn].thetacen / RADIAN;

    w[nn].rcen = imported_model[ndom].wind_midx[imported_model[ndom].i[n]];
    w[nn].xcen[0] = w[nn].rcen * sin (theta);
    w[nn].xcen[1] = 0;
    w[nn].xcen[2] = w[nn].rcen * cos (theta);

    /* Copy across the inwind variable to the wind pointer */
    w[nn].inwind = imported_model[ndom].inwind[n];


    /* 1812 - ksl - For imported models, one is either in the wind or not. But we need
     * to make sure the rest of the code knows that this cell is to be ignored in
     * this case. Adapted from the code in import_cylindrical */
    if (w[nn].inwind == W_NOT_INWIND || w[nn].inwind == W_PART_INWIND)
      w[nn].inwind = W_IGNORE;
  }

  /* Now add information used in zdom */

  for (n = 0; n < zdom[ndom].ndim; n++)
  {
    zdom[ndom].wind_x[n] = imported_model[ndom].wind_x[n];
  }



  for (n = 0; n < zdom[ndom].mdim; n++)
  {
    zdom[ndom].wind_z[n] = imported_model[ndom].wind_z[n];
  }

  /* The next line is necessary for calculating distances in a cell in rthota coordiatnes */

  rtheta_make_cones (ndom, w);

  return (0);
}


/* The next section calculates velocites.  We follow the hydro approach of
 * getting those velocities from the original grid.  This is really only
 * used for setting up the grid
 *
 * The code here is identical to that in velocity_cylindrical, which suggests
 * that it could be used for any regular 2d grid
 */


/**********************************************************/
/**
 * @brief      The velocity at any position in an imported
 * rtheat model
 *
 * @param [in] int  ndom   The domain for the imported model
 * @param [in] double *  x   A position (3 vector)
 * @param [out] double *  v   The calcuated velocity
 * @return     The speed at x
 *
 * @details
 * This routine interpolates on the values read in for the
 * imported model to give one a velocity
 *
 *
 * ### Notes ###
 *  In practice this routine is only used to initialize v in
 *  wind structure, and interpolation is only a convenient 
 *  way to do this.  
 *  
 *  This is consistent with the way velocities
 *  are treated throughout Python.
 *
 **********************************************************/

double
velocity_rtheta (ndom, x, v)
     int ndom;
     double *x, *v;
{
  int j;
  int nn;
  int nnn[4], nelem;
  double frac[4];
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

  /* Now copy the result into v, which is very necessary if refilling wmain.v */

  v[0] = vv[0];
  v[1] = vv[1];
  v[2] = vv[2];

  return (speed);
}



/**********************************************************/
/**
 * @brief      Get the density for an imported rtheta model at x
 *
 * @param [in] int  ndom   The domain for the imported model
 * @param [in] double *  x   A postion
 * @return     The density in cgs units is returned
 *
 * @details
 * This routine finds rho from the imported model
 * at a position x.  The routine does not interpolate rho, but
 * simply locates the cell associated with x
 *
 * ### Notes ###
 * This routine is really only used to initialize rho in the
 * Plasma structure.  In reality, once the Plasma structure is
 * initialized, we always interpolate within the plasma structure
 * and do not access the original data
 *
 * This routine is dependent on the assumption that x is in 
 * the center of a cell.
 *
 **********************************************************/

double
rho_rtheta (ndom, x)
     int ndom;
     double *x;
{
  double rho = 0;
  double r, z;
  int i, j, n;
  double ctheta, angle;

  r = length (x);

  z = fabs (x[2]);
  ctheta = z / r;
  angle = acos (ctheta) * RADIAN;

  i = 0;
  while (angle > imported_model[ndom].wind_z[i] && i < imported_model[ndom].mdim)
  {
    i++;
  }
  i--;

  j = 0;
  while (r > imported_model[ndom].wind_x[j] && j < imported_model[ndom].ndim)
  {
    j++;
  }
  j--;

  n = j * imported_model[ndom].mdim + i;

  rho = imported_model[ndom].mass_rho[n];

  return (rho);
}




/* ************************************************************************** */
/**
 * @brief      Get the temperature at a position x
 *
 * @param[in] int    ndom        The domain for the imported model
 * @param[in] double *x          A position (3d)
 * @param[in] int    return_t_e  If TRUE, the electron temperature is returned
 *
 * @return     The temperature in K
 *
 * @details
 *
 * ************************************************************************** */

double
temperature_rtheta (int ndom, double *x, int return_t_e)
{
  int i, j, n;
  double r, z;
  double temperature = 0;
  double ctheta, angle;

  if (imported_model[ndom].init_temperature)
  {
    if (return_t_e)
      temperature = 1.1 * zdom[ndom].twind;
    else
      temperature = zdom[ndom].twind;
  }
  else
  {
    r = length (x);
    z = fabs (x[2]);

    ctheta = z / r;
    angle = acos (ctheta) * RADIAN;

    i = 0;
    while (angle > imported_model[ndom].wind_z[i] && i < imported_model[ndom].mdim)
    {
      i++;
    }
    i--;

    j = 0;
    while (r > imported_model[ndom].wind_x[j] && j < imported_model[ndom].ndim)
    {
      j++;
    }
    j--;

    n = j * imported_model[ndom].mdim + i;

    if (return_t_e)
      temperature = imported_model[ndom].t_e[n];
    else
      temperature = imported_model[ndom].t_r[n];
  }

  return (temperature);
}
