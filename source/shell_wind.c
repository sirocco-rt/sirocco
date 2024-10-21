
/***********************************************************/
/** @file  shell_wind.c
 * @author ksl,nsh
 * @date   May, 2018
 *
 * @brief
 * Routines needed for a single shell wind model
 *
 * These routines were developed (by NSH) primarily for diagnostic purposes, e.g for calculating ionization models for
 * comparison with Cloudy.  The routines are in many ways similar to stellar wind model routines, but after a discussion
 * in 2015 (See issue #172), we decided to keep it.
 ***********************************************************/

/*
   This file was created in Feb 2011.
   The purpose is to have a model where we have a single shell of material.
   This is different from the stellar wind because we do not want the inner surface of the wind to touch the star.
   This requires tight control of the grid and makes for a very prescriptive model.  We also need a special grid,
   which is also stored in this file.

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

/**********************************************************/
/**
 * @brief Allocate memory for the wind coordinate arrays and for a shell wind
 *
 * @param [in] ndom The domain number to allocate for
 *
 * @details
 *
 * Memory is allocated for the `wind_x`, `wind_z`, `wind_midx` and `wind_midz`
 * arrays. The shell wind is a special wind case and the inputs and domain setup
 * are handled in a different order/way. This is why a shell wind has its own
 * domain allocation function in here.
 *
***********************************************************/

static int
allocate_shell_domain (int ndom)
{
  zdom[ndom].wind_x = calloc (zdom[ndom].ndim, sizeof (double));
  zdom[ndom].wind_midx = calloc (zdom[ndom].ndim, sizeof (double));
  if (zdom[ndom].wind_x == NULL || zdom[ndom].wind_midx == NULL)
  {
    Error ("allocate_domain_wind_coords: Problem allocating memory for x-dim for domain %d\n", ndom);
    Exit (EXIT_FAILURE);
  }

  /* Allocate z dimensions */
  zdom[ndom].wind_z = calloc (zdom[ndom].mdim, sizeof (double));
  zdom[ndom].wind_midz = calloc (zdom[ndom].mdim, sizeof (double));
  if (zdom[ndom].wind_z == NULL || zdom[ndom].wind_midz == NULL)
  {
    Error ("allocate_domain_wind_coords: Problem allocating memory for z-dim for domain %d\n", ndom);
    Exit (EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}

/**********************************************************/
/**
 * @brief      Get the inputs needed to specigy the shell_wind model.
 *
 * @param [in] int  ndom   The domain for the shell
 * @return     Alwyas returns 0
 *
 * @details
 * This routine gets the inputs needed for the shell_wind model and
 * translates them into those needed for the equvalent stellar wind
 * model. Several arrays specific to the shell_wind model are also
 * initialized.
 *
 * ### Notes ###
 *
 **********************************************************/

int
get_shell_wind_params (ndom)
     int ndom;
{
  double vtemp[3];
  double rhotemp[200];
  double rmin, dr, r[200];
  double postemp[3];
  double speedtemp;
  double cdensity;
  double shell_vmin, shell_vmax;        //Local variables to define shellwind
  double shell_rmin, shell_rmax;        //Local variables to define shellwind
  double factor;
  int i;

  Log ("Creating wind with a single shell for testing purposes\n");


/* In order to make life simple, we are first going to check that we are in spherical coordinates, if not change!! */

  if (zdom[ndom].coord_type != SPHERICAL)
  {
    Error ("For the shell type wind, we should be in spherical coordinates - changing....\n");
    zdom[ndom].coord_type = SPHERICAL;
  }

/* This may not be important, but we sould make sure that NDIM is 4... */

  if (zdom[ndom].ndim != 4)
  {
    Error ("For the shell type wind, we take control of the grid, and need NDIM to be the minimum - 4 - changing\n");
    zdom[ndom].ndim = 4;
    zdom[ndom].mdim = 1;
    zdom[ndom].ndim2 = 4;
  }

  allocate_shell_domain (ndom);
  zdom[ndom].stellar_wind_mdot = 1.e-6;
  zdom[ndom].rmin = geo.rstar;
  zdom[ndom].cl_beta = 1.0;
  shell_rmin = zdom[ndom].cl_rmin = 2.8e9;
  shell_vmin = zdom[ndom].cl_v_zero = 200e5;
  shell_vmax = zdom[ndom].cl_v_infinity = 3000e5;

  rddoub ("Shell.wind_mdot(msol/yr)", &zdom[ndom].stellar_wind_mdot);
  zdom[ndom].stellar_wind_mdot *= MSOL / YR;

  rddoub ("Shell.wind.radmin(cm)", &zdom[ndom].rmin);   /*Radius where wind begins */
  if (zdom[ndom].rmin < geo.rstar)
  {
    Error ("get_shell_wind_params: It is unreasonable to have the wind start inside the star!\n");
    Log ("Setting zdom[ndom].rmin to geo.rstar\n");
    zdom[ndom].rmin = geo.rstar;
  }
  zdom[ndom].cl_rmin = shell_rmin = zdom[ndom].rmin;

  zdom[ndom].rmax = 1.1 * zdom[ndom].rmin;
  rddoub ("Shell.wind.radmax(cm)", &zdom[ndom].rmax);   /*Radius where wind begins */

/*120130 NSH the next two lines have been modified to mean that the wind will end up as a CL wind,
 * but the v_0 and v_infinity will be calulated here from these two variables, which are now local */

  rddoub ("Shell.wind_v_at_rmin(cm)", &shell_vmin);     /* Velocity at base of the wind */
  rddoub ("Shell.wind.v_at_rmax(cm)", &shell_vmax);     /* Final speed of wind in units of escape velocity */


  rddoub ("Shell.wind.acceleration_exponent", &zdom[ndom].cl_beta);     /* Accleration scale exponent for a CL wind */
  Log ("Geo rmax = %f\n", zdom[ndom].rmax);

  shell_rmax = zdom[ndom].rmax;

  /* 120130 NSH These next lines invert the cl velocity equation to get the cl factors from the local shell factors */

  zdom[ndom].cl_v_zero = shell_vmin;
  factor = pow ((1 - (zdom[ndom].rmin / zdom[ndom].rmax)), zdom[ndom].cl_beta);
  zdom[ndom].cl_v_infinity = (shell_vmax - shell_vmin + shell_vmin * factor) / factor;


  /* Assign the generic parameters for the wind the generic parameters of the wind */

  zdom[ndom].wind_thetamin = 0.0;
  zdom[ndom].wind_thetamax = 90. / RADIAN;

  /* define the the variables that determine the gridding */
  zdom[ndom].wind_rhomin_at_disk = 0;
  zdom[ndom].wind_rhomax_at_disk = zdom[ndom].rmax;
  zdom[ndom].zmax = zdom[ndom].rmax;

  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
  {
    zdom[ndom].xlog_scale = 0.3 * geo.rstar;
    zdom[ndom].zlog_scale = 0.3 * geo.rstar;
  }

  /* Since this is a diagnostic routine, we write out some additional information to the log file.  Note that
   * much of this could be done as post-procssing step, and so it is possible that these messages are superflous. */

  Log ("shell rmin=%f shell rmax=%f\n", shell_rmin, shell_rmax);
  dr = (shell_rmax - shell_rmin) / 100.0000;
  Log ("dr= %e, root2= %10.30e\n", dr, pow (2.0, 0.5));
  rmin = shell_rmin - (dr);

  for (i = 0; i < 103; i++)
  {
    r[i] = rmin + i * dr;

    postemp[0] = postemp[2] = r[i] / pow (2.0, 0.5);
    postemp[1] = 0.0;
    speedtemp = stellar_velocity (ndom, postemp, vtemp);
    rhotemp[i] = stellar_rho (ndom, postemp) * rho2nh;
    Log_silent ("SHELL  ring=%i,x=%11.4e,r=%11.4e,speed=%11.4e,density=%11.4e\n", i, r[i] / pow (2.0, 0.5), r[i], speedtemp, rhotemp[i]);
  }

  cdensity = 0.0;
  for (i = 1; i < 100; i++)
  {
    cdensity += ((rhotemp[i] + rhotemp[i + 1]) / 2.) * dr;
  }
  Log ("Column density of hydrogen=%11.4e\n", cdensity);


  return (0);
}




/**********************************************************/
/**
 * @brief      defines the cells in a thin shell, a special case of a spherical wind.
 *
 * @param [in, out] WindPtr  w   The structure which defines the wind in Python
 * @param [in, out] int  ndom   The domain number
 * @return     Always retunrs 0
 *
 * @details
 *
 * The shell wind has 3 elements, one inside the shell, one outside, and
 * one outside and one shell exactly fitting the shell.
 *
 *
 * ### Notes ###
 *
 * The shell wind is intended for diagnostic purposes
 *
 *
 **********************************************************/

int
shell_make_grid (int ndom, WindPtr w)
{
  int n;
  int ndim;
  int nstart;

  ndim = zdom[ndom].ndim;
  nstart = zdom[ndom].nstart;


  w[nstart + 0].r = 0.999999 * zdom[ndom].rmin;
  w[nstart + 1].r = zdom[ndom].rmin;
  w[nstart + 2].r = zdom[ndom].rmax;
  w[nstart + 3].r = 1.000001 * zdom[ndom].rmax;



  w[nstart + 0].rcen = (w[nstart + 0].r + w[nstart + 1].r) / 2;
  w[nstart + 1].rcen = (w[nstart + 1].r + w[nstart + 2].r) / 2;
  w[nstart + 2].rcen = (w[nstart + 2].r + w[nstart + 3].r) / 2;
  w[nstart + 3].rcen = w[nstart + 2].rcen + (zdom[ndom].rmax - zdom[ndom].rmin);

  /* Now calculate the positions of these points in the xz plane.
     There is a choice about how one does this.   I have elected
     to assume that we want to calculate this at a 45 degree angle.
     in the hopes this will be a reasonable portion of the wind in
     a biconical flow.
   */
  for (n = 0; n < ndim; n++)
  {
    Log ("Shell_wind: cell %i:  inner edge = %2.20e, centre = %2.20e\n", n, w[nstart + n].r, w[nstart + n].rcen);
    w[nstart + n].x[1] = w[nstart + n].xcen[1] = 0.0;

    //NSH Slight change here, using 1/root2 give more accurate results than sin45.


    w[nstart + n].x[0] = w[nstart + n].x[2] = w[nstart + n].r / pow (2.0, 0.5);
    w[nstart + n].xcen[0] = w[nstart + n].xcen[2] = w[nstart + n].rcen / pow (2.0, 0.5);
  }


  return (0);

}
