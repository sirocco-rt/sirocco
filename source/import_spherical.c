
/***********************************************************/
/** @file  import_spherical.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief
 * General purpose routines for reading in an arbitrary wind model
 * in spherical coordinates.
 *
 * The basic data we need to read in are,
 *
 *     i r v_r mass_rho (and optionally T_r)
 *
 *  where,
 *
 *  * i are the element numbers (increasing outwards)
 *  * r is the radial coordinates
 *  * v_r is the velocity in the radial direction
 *  * mass_rho is the density in cgs units
 *  * T_r is the radiation temperature in Kelvin
 *
 * We assume that all of the physical quantities are centered, that is
 * we are assuming that we are giving mass_rho/T_r at the center of
 * a cell. However, r and v_r should be given at the edges of a cell.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief      Read an arbitrary wind model in spherical symmetry.
 *
 * @param [in out] int     ndom       The domain number for the imported model
 * @param [in out] char *  filename   The file containing the model to import
 * @return                            Always returns 0
 *
 * @details
 *
 * This routine just reads in the data and stores it in arrays
 *
 * ### Notes ###
 *
 * The basic data we need to read in are,
 *
 *     i r v_r mass_rho (and optionally T_r)
 *
 *  where,
 *
 *  * i are the element numbers (increasing outwards)
 *  * r is the radial coordinates
 *  * v_r is the velocity in the radial direction
 *  * mass_rho is the density in cgs units
 *  * T_r is the radiation temperature in Kelvin
 *
 * We assume that all of the physical quantities are centered, that is
 * we are assuming that we are giving mass_rho/T_r at the center of
 * a cell. However, r and v_r should be given at the edges of a cell.
 *
 **********************************************************/

int
import_1d (ndom, filename)
     int ndom;
     char *filename;
{
  FILE *fptr;
  char line[LINELENGTH];
  int n, icell, ncell;
  double r, v_r, mass_rho, t_r, t_e;

  Log ("Reading a 1d model %s\n", filename);

  if ((fptr = fopen (filename, "r")) == NULL)
  {
    Error ("import_1d: No such file\n");
    Exit (0);
  }

  ncell = 0;
  while (fgets (line, LINELENGTH, fptr) != NULL)
  {
    n = sscanf (line, " %d %le %le %le %le %le", &icell, &r, &v_r, &mass_rho, &t_e, &t_r);
    if (n < READ_NO_TEMP_1D)
    {
      continue;
    }
    else
    {
      imported_model[ndom].i[ncell] = icell;
      imported_model[ndom].r[ncell] = r;
      imported_model[ndom].v_r[ncell] = v_r;
      imported_model[ndom].mass_rho[ncell] = mass_rho;

      if (n == READ_ELECTRON_TEMP_1D)
      {
        imported_model[ndom].init_temperature = FALSE;
        imported_model[ndom].t_e[ncell] = t_e;
        imported_model[ndom].t_r[ncell] = 1.1 * t_e;
      }
      else if (n == READ_BOTH_TEMP_1D)
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

      if (ncell > NDIM_MAX)
      {
        Error ("import_1d: trying to read in more grid points than allowed (%i). Try changing NDIM_MAX and recompiling.\n", NDIM_MAX);
        Exit (1);
      }

    }
  }

  imported_model[ndom].ncell = ncell;

  /*
   * Check that each cell has its own unique radius and that the radius is
   * constantly increasing with grid cell
   */

  for (n = 1; n < imported_model[ndom].ncell; ++n)
  {
    if (imported_model[ndom].r[n] <= imported_model[ndom].r[n - 1])
    {
      Error ("import_1d: cell %i r %e < cell %i r %e. The grid radii must be constantly increasing in size. Exiting!\n", n,
             imported_model[ndom].r[n], n - 1, imported_model[ndom].r[n - 1]);
      Exit (1);
    }
  }

  /* Although much of the initialization of zdom can be postponed
   * one has to define mdim and ndim of zdom here, so that the correct
   * number of wind cells will be allocated */

  imported_model[ndom].ndim = zdom[ndom].ndim2 = zdom[ndom].ndim = imported_model[ndom].ncell;
  imported_model[ndom].mdim = zdom[ndom].mdim = 1;

  return (0);
}




/* ************************************************************************** */
/**
 * @brief   Set up the various domain boundaries for a spherical coordinate
 *          system
 *
 * @param[in] int ndom         The domain of interest
 *
 * @return    Always returns 0
 *
 * @details
 *
 * This used to be contained within spherical_make_grid_import, however, it
 * does not reply on any of the variables in that function and only relies on
 * the imported_model struct. Therefore, the boundary setup was moved into a
 * separate function so it could be done else where in the program flow.
 *
 * ************************************************************************** */

int
import_spherical_setup_boundaries (int ndom)
{
  zdom[ndom].wind_rhomin_at_disk = 0;
  zdom[ndom].rmin = imported_model[ndom].r[1];  // <- this assumes the 1st cell is a ghost cell
  zdom[ndom].wind_rhomax_at_disk = zdom[ndom].zmax = zdom[ndom].rmax = imported_model[ndom].r[imported_model[ndom].ncell - 2];  // <- this assumes the last 2 cells are ghost cells
  zdom[ndom].wind_thetamin = zdom[ndom].wind_thetamax = 0;

  return 0;
}




/* The next section contains routines to make the grids for imported models.

   We make some assumptions here.  We assume that every cell is in the wind.
   and we assume that r refers to the inside edge of the cell.

 * */


/**********************************************************/
/**
 * @brief      Use the imported data to initialize various
 * portions of the Wind and Domain structures
 *
 *
 * @param [out] WindPtr  w   The entire wind
 * @param [in] int  ndom   The domain for the imported model
 * @return     Always returns 0
 *
 * @details
 * This routine initializes the portions of the wind structure
 * using the imported model, specifically those portions having
 * to do with positions.
 *
 * ### Notes ###
 *
 * The velocity v_r is stored in v[0] of the wind array.
 * See velocity_1d for how this is interpreted to generate
 * a 3-d velocity
 *
 **********************************************************/

int
spherical_make_grid_import (w, ndom)
     WindPtr w;
     int ndom;
{

  int j, n;

  for (j = 0; j < imported_model[ndom].ncell; j++)
  {
    n = j + zdom[ndom].nstart;
    w[n].r = imported_model[ndom].r[j];
    /* Put the radial velocity in v[0] */
    w[n].v[0] = imported_model[ndom].v_r[j];
  }

  /* Need to define the midpoints of the grid */

  for (j = 0; j < zdom[ndom].ndim; j++)
  {
    n = j + zdom[ndom].nstart;
    if (j < zdom[ndom].ndim - 1)
    {
      w[n].rcen = 0.5 * (w[n].r + w[n + 1].r);
    }
    else
    {
      w[n].rcen = w[n].r * 1.005;
    }
    w[n].inwind = W_ALL_INWIND;
    w[n].x[1] = w[n].xcen[1] = 0.0;
    w[n].x[0] = w[n].x[2] = w[n].r * sin (PI / 4.);
    w[n].xcen[0] = w[n].xcen[2] = w[n].rcen * sin (PI / 4.);
  }

  /*
   * The first and last two cells are "guard" cells, as such they are not
   * considered to be in the wind.
   */

  w[0].inwind = w[imported_model[ndom].ncell - 1].inwind = w[imported_model[ndom].ncell - 2].inwind = W_NOT_INWIND;

  /* Since we assume all of the cells are in the wind in a spherical wind
   * we can use the standard routine to finish everything off
   */

  spherical_wind_complete (ndom, w);

  return (0);
}


/* The next section calculates velocities.
 *
 * One could follow the zeus_hydro approach of getting those velocities from the original grid.
 * but for consistency with the 2d case we get it by interpolating on values in the cells
 *
 *
 */


/**********************************************************/
/**
 * @brief      The velocity at any position in an imported spherical
 * model
 *
 *
 * @param [in] int  ndom   The domain of the imported model
 * @param [in] double *  x   A position (3d)
 * @param [out] double *  v   The velocity at x
 * @return     The speed at x
 *
 * @details
 * This routine interpolates on the values read in for the
 * imported model to give one a velocity
 *
 *
 * ### Notes ###

 * Note that v_r has been  stored (see sperical_make_grid_import)
 * in v_0 and this explains
 * the way the speed is calculated, and then translated 
 * to a velocity in a 3d space
 *
 **********************************************************/

double
velocity_1d (ndom, x, v)
     int ndom;
     double *x, *v;
{
  double speed;
  double r;
  int nelem, nn, nnn[4];
  double frac[4];

  r = length (x);
  coord_fraction (ndom, 0, x, nnn, frac, &nelem);

  /* For imported spherical/1d models the total velocity is
     stored in v[0], which explains the code below.  See
     spherical_make_grid_import. See issue #787
   */
  speed = 0;
  for (nn = 0; nn < nelem; nn++)
  {
    speed += wmain[zdom[ndom].nstart + nnn[nn]].v[0] * frac[nn];
  }

  v[0] = x[0] / r * speed;
  v[1] = x[1] / r * speed;
  v[2] = x[2] / r * speed;

  return (speed);
}




/**********************************************************/
/**
 * @brief      Get the density for an imported spherical model at x
 *
 * @param [in] int  ndom   The domain for the imported model
 * @param [in] double *  x   A position (3d)
 * @return     The density in cgs units is returned
 *
 * @details
 * This routine finds rho from the imported model
 * at a position x.  The routine does not interpolate rho, but
 * simply locates the cell associated with x
 *
 *
 * ### Notes ###
 * This routine is really only used to intialize rho in the
 * Plasma structure.  In reality, once the Plasma structure is
 * initialized we always interpolate within the plasma structure
 * and do not access the original data.
 *
 * This routine is peculiar because it depends on the assumption
 * that x is the center of the cell.
 *
 **********************************************************/

double
rho_1d (ndom, x)
     int ndom;
     double *x;
{
  double rho = 0;
  double r;
  int n;

  r = length (x);

  n = 0;
  while (r >= imported_model[ndom].r[n] && n < imported_model[ndom].ncell)
  {
    n++;
  }
  n--;

  if (n < imported_model[ndom].ncell)
  {
    rho = imported_model[ndom].mass_rho[n];
  }
  else
  {
    rho = imported_model[ndom].mass_rho[imported_model[ndom].ncell - 1];
  }

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
temperature_1d (int ndom, double *x, int return_t_e)
{
  int n;
  double r, temperature = 0;

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

    n = 0;
    while (r >= imported_model[ndom].r[n] && n < imported_model[ndom].ncell)
    {
      n++;
    }
    n--;

    if (n < imported_model[ndom].ncell)
    {
      if (return_t_e)
        temperature = imported_model[ndom].t_e[n];
      else
        temperature = imported_model[ndom].t_r[n];
    }
    else
    {
      if (return_t_e)
        temperature = imported_model[ndom].t_e[imported_model[ndom].ncell - 1];
      else
        temperature = imported_model[ndom].t_r[imported_model[ndom].ncell - 1];
    }
  }

  return temperature;
}
