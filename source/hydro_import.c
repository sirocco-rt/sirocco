/***********************************************************/
/** @file   hydro_import.c
 * @author KSL/NSH
 * @date   January 2018
 * @brief  Routines to import hydro-dynamic snapshots.
 *
 * These routines are designed to read the velocity and densities
 * of produced by Zeus  into the appropriate structures and
 * allow one to calculate the density at any point in the gridded space.
 *
 * Thesse routines were originally written by ksl to work with
 * models which Daniel Progra provided by they were extensively modified
 * by Nick for his work with Zeus
 *
 ***********************************************************/

#define MAXHYDRO 1000
#define IGHOST 0

double hydro_r_cent[MAXHYDRO];
double hydro_r_edge[MAXHYDRO];
double hydro_theta_cent[MAXHYDRO];
double hydro_theta_edge[MAXHYDRO];
int ihydro_r, ihydro_theta, j_hydro_thetamax, ihydro_mod;
double hydro_thetamax;          //The angle at which we want to truncate the theta grid
double v_r_input[MAXHYDRO * MAXHYDRO];
double v_theta_input[MAXHYDRO * MAXHYDRO];
double v_phi_input[MAXHYDRO * MAXHYDRO];
double rho_input[MAXHYDRO * MAXHYDRO];
double temp_input[MAXHYDRO * MAXHYDRO];


/*
typedef struct hydro_mod
{
  double rho;
  double v[3];
  double temp;
}
hydro_dummy, *HydroPtr;

HydroPtr hydro_ptr;
*/

#include  <stdio.h>
#include  <stdlib.h>
#include  <strings.h>
#include  <string.h>
#include  <math.h>

#include  "log.h"
#include  "atomic.h"
#include  "sirocco.h"

#define LINE 200






/**********************************************************/
/**
 * @brief	Sets up geometry ready to accept a hydro-dynamic model
 *
 * @param [in] ndom			The number of the domain which will be populated by the geometry
 * @return 					0 if successful
 *
 * Sets up the geometry - max and min - for the hydro-dynamic
 * snapshot. It uses limits read in from the file to fill in
 * parameters more usually read in from a .pf file
 *
 *
 * ###Notes###
***********************************************************/


int
get_hydro_wind_params (ndom)
     int ndom;
{
  Log ("Creating a wind model using a Hydro calculation = domain %i\n", ndom);

  get_hydro (ndom);

  /* updates_2d has special code related to zeus operations with zeus. The
   * next line preserves information related to the domain number
   */
  geo.hydro_domain_number = ndom;

/* Assign the generic parameters for the wind the generic parameters of the wind */
  zdom[ndom].rmin = hydro_r_edge[0];


/* ksl - One should not be defining geo.rmax on a domain basis.  It is calculated from the maximum valuse os zdom[].rmax for all of the domains */
  zdom[ndom].rmax = hydro_r_edge[ihydro_r] + 2.0 * (hydro_r_cent[ihydro_r] - hydro_r_edge[ihydro_r]);   //Set the outer edge of the wind to the outer edge of the final defined cell



  Log ("rmax=%e\n", zdom[ndom].rmax);
  zdom[ndom].wind_rhomin_at_disk = 0.0; //Set wind_rmin 0.0, otherwise wind cones dont work properly
  Log ("rho_min=%e\n", zdom[ndom].wind_rhomin_at_disk);
  zdom[ndom].wind_rhomax_at_disk = zdom[ndom].rmax;     //This is the outer edge of the
  Log ("rho_max=%e\n", zdom[ndom].wind_rhomax_at_disk);
  zdom[ndom].zmax = zdom[ndom].rmax;    //This is the outer edge of the
  Log ("zmax=%e\n", zdom[ndom].zmax);

  zdom[ndom].wind_thetamin = hydro_theta_edge[0];
  Log ("theta_min=%e\n", zdom[ndom].wind_thetamin);
  Log ("theta_max=%e\n", zdom[ndom].wind_thetamax);
  Log ("rmin=%e\n", zdom[ndom].rmin);
  Log ("rmax=%e\n", zdom[ndom].rmax);
  Log ("wind_rhomin_at_disk=%e\n", zdom[ndom].wind_rhomin_at_disk);
  Log ("wind_rhomax_at_disk=%e\n", zdom[ndom].wind_rhomax_at_disk);



  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
  {
    zdom[ndom].xlog_scale = 0.3 * geo.rstar;
    zdom[ndom].zlog_scale = 0.3 * geo.rstar;
  }

  return (0);
}


/**********************************************************/
/**
 * @brief	Reads in and decodes a hydro-dynamic model
 *
 * @param [in] ndom			The number of the domain which will be populated by the geometry
 * @return 					0 if successful
 *
 * This deals with importing the hydro-dynamic snapshot
 * It asks the user the filename and the maximum theta to be used
 * It also works out the maximum theta bin to be used, and
 * interpolates if this doesnt quite line up.
 *
 *
 * ###Notes###
***********************************************************/


int
get_hydro (ndom)
     int ndom;
{

  FILE *fopen (), *fptr;
  char datafile[LINE];
  char aline[LINE];
  char word[LINE];
  int i, j, k;
  double r, r_edge;
  double rho;
  double theta, theta_edge, temp;
  double vr, vtheta, vphi;
  int irmax, ithetamax, itest;
//OLD  int ndim, mdim;

/*Write something into the file name strings */

  strcpy (datafile, "hdf062.dat");


  for (k = 0; k < MAXHYDRO; k++)
  {
    hydro_r_cent[k] = 0;
    hydro_theta_cent[k] = 0;
  }

  rdstr ("Hydro.file", datafile);
  if ((fptr = fopen (datafile, "r")) == NULL)
  {
    Error ("Could not open %s\n", datafile);
    Exit (0);
  }


  hydro_thetamax = 89.9;


  rddoub ("Hydro.thetamax(degrees:negative_means_no_maximum)", &hydro_thetamax);

  //If we have set the maximum angle to a negaive value, we mean that we dont want to restrict.
  if (hydro_thetamax < 0.0)
    hydro_thetamax = VERY_BIG;
  hydro_thetamax = hydro_thetamax / RADIAN;



  ihydro_r = 0;
  j_hydro_thetamax = 0;         /* NSH 130605 to remove o3 compile error */
  irmax = 0;
  ithetamax = 0;                //Zero counters to check all cells actually have data

  while (fgets (aline, LINE, fptr) != NULL)
  {

    if (aline[0] != '#')
    {
      sscanf (aline, "%s", word);

      if (strncmp (word, "ir", 2) == 0)
      {
        Log ("We have a hydro title line, we might do something with this in the future \n");
      }
      else
      {
        itest =
          sscanf (aline, "%d %lf %lf %d %lf %lf %lf %lf %lf %lf %lf",
                  &i, &r, &r_edge, &j, &theta, &theta_edge, &vr, &vtheta, &vphi, &rho, &temp);
        if (itest != 11)        //We have an line which does not match what we expect, so quit
        {
          Error ("hydro.c data file improperly formatted\n");
          Exit (0);
        }
        // read read the r and theta coordinates into arrays
        hydro_r_edge[i] = r_edge;
        hydro_r_cent[i] = r;
        hydro_theta_cent[j] = theta;
        hydro_theta_edge[j] = theta_edge;
        //keep track of how many r and theta cells there are
        if (j > ithetamax)
          ithetamax = j;
        if (i > irmax)
          irmax = i;
        //If the value of theta in this cell, the edge, is greater than out theta_max, we want to make a note.
        if (hydro_theta_edge[j] > hydro_thetamax && hydro_theta_edge[j - 1] <= hydro_thetamax)
        {
          j_hydro_thetamax = j - 1;
          Log
            ("current theta  (%f) > theta_max  (%f) so setting j_hydro_thetamax=%i\n",
             theta * RADIAN, hydro_thetamax * RADIAN, j_hydro_thetamax);
        }
        //If theta is greater than thetamax, then we will replace rho with the last density above the disk
        else if (hydro_theta_edge[j] > hydro_thetamax)
        {
          /* NSH 130327 - for the time being, if theta is in the disk, replace with the last
             density above the disk */
          rho = rho_input[i * MAXHYDRO + j_hydro_thetamax];
        }
        rho_input[i * MAXHYDRO + j] = rho;
        temp_input[i * MAXHYDRO + j] = temp;
        v_r_input[i * MAXHYDRO + j] = vr;
        v_theta_input[i * MAXHYDRO + j] = vtheta;
        v_phi_input[i * MAXHYDRO + j] = vphi;
      }
    }
  }



  if (j_hydro_thetamax == 0 || j_hydro_thetamax == i - 1)
  {
    Log ("HYDRO j_hydro_thetamax never bracketed, using all data\n");
    ihydro_theta = ithetamax;
    zdom[ndom].wind_thetamax = 90. / RADIAN;
    hydro_thetamax = 90.0 / RADIAN;
//OLD    mdim = zdom[ndom].mdim = ihydro_theta + 2;
    zdom[ndom].mdim = ihydro_theta + 2;
  }
  else
  {
    Log
      ("HYDRO j_hydro_thetamax=%i, bracketing cells have theta = %f and %f\n",
       j_hydro_thetamax, hydro_theta_cent[j_hydro_thetamax] * RADIAN, hydro_theta_cent[j_hydro_thetamax + 1] * RADIAN);
    ihydro_theta = j_hydro_thetamax;
    zdom[ndom].wind_thetamax = hydro_thetamax;
//OLD    mdim = zdom[ndom].mdim = ihydro_theta + 2;
    zdom[ndom].mdim = ihydro_theta + 2;
  }




  if (hydro_r_edge[0] < geo.rstar)
  {
    Error ("Major problem, innermost edge of hydro radial grid begins inside geo.rstar\n");
    Exit (0);
  }

  ihydro_r = irmax;

  Log ("Read %d r values\n", ihydro_r);
  Log ("Read %d theta values\n", ihydro_theta);
  fclose (fptr);

  /* Set a couple of last tags */

  zdom[ndom].coord_type = RTHETA;       //At the moment we only deal with RTHETA - in the future we might want to do some clever stuff
//OLD  ndim = zdom[ndom].ndim = ihydro_r + 3;        //We need an inner radial cell to bridge the star and the inside of the wind, and an outer cell
  zdom[ndom].ndim = ihydro_r + 3;       //We need an inner radial cell to bridge the star and the inside of the wind, and an outer cell
  zdom[ndom].ndim2 = zdom[ndom].ndim * zdom[ndom].mdim; // Make ndim2 consistent with the individual dimensions


  return (0);
}






/**********************************************************/
/**
 * @brief	Interpolates on supplied data file to get velocity at all vertices
 *
 * @param [in] ndom			The number of the domain which will be populated by the geometry
 * @param [in] x			The x location we need a velocity for
 * @param [in] v			The velocity we are returning
 * @return 					The speed - relating the velocity we have computed
 *
 * The velocities in a zeus hydro-grid are defined on the middle of cell
 * edges, wheras in sirocco we require them at the verticies. This routine
 * interpolates on rhe supplied grid to make the sirocco grid.
 *
 *
 * ###Notes###
 *
***********************************************************/



double
hydro_velocity (ndom, x, v)
     int ndom;
     double x[];
     double v[];
{
  double length ();
  int ii, jj;
  int im, jm;
  double f1, f2;
  double r, theta;
  double v_r_interp, v_theta_interp, v_phi_interp;
  double speed;
  double xxx;

  if ((r = length (x)) == 0.0)
  {
    v[0] = v[1] = v[2] = 0.0;
    return (0.);
  };



  xxx = (sqrt (x[0] * x[0] + x[1] * x[1]) / r);
  if (xxx >= 1.)
  {
    theta = PI / 2.;
  }
  else
    theta = asin (xxx);


  /* We compute the radial velocity - in zeus this is defined at the
     cell centre in the theta dimension and the edge in the r dimension */


  im = jm = ii = jj = 0;
  f1 = f2 = 0.0;
  hydro_frac (r, hydro_r_edge, ihydro_r, &im, &ii, &f1);
  hydro_frac (theta, hydro_theta_cent, ihydro_theta, &jm, &jj, &f2);
  v_r_interp = hydro_interp_value (v_r_input, im, ii, jm, jj, f1, f2);


  /* We compute the theta velocity - in zeus this is defined at the
     cell centre in the r dimension and the edge in the theta dimension */

  im = jm = ii = jj = 0;
  f1 = f2 = 0.0;
  hydro_frac (r, hydro_r_cent, ihydro_r, &im, &ii, &f1);
  hydro_frac (theta, hydro_theta_edge, ihydro_theta, &jm, &jj, &f2);
  v_theta_interp = hydro_interp_value (v_theta_input, im, ii, jm, jj, f1, f2);


  /* We compute the phi velocity - in zeus this is defined at the
     cell centre in both the r and theta dimension */


  im = jm = ii = jj = 0;
  f1 = f2 = 0.0;
  hydro_frac (r, hydro_r_cent, ihydro_r, &im, &ii, &f1);
  hydro_frac (theta, hydro_theta_cent, ihydro_theta, &jm, &jj, &f2);
  v_phi_interp = hydro_interp_value (v_phi_input, im, ii, jm, jj, f1, f2);


  /* Convert to carteisian velocity */


  v[0] = v_r_interp * sin (theta) + v_theta_interp * cos (theta);
  v[1] = v_phi_interp;
  v[2] = v_r_interp * cos (theta) - v_theta_interp * sin (theta);


  /*Compute the speed */

  speed = sqrt (v[0] * v[0] + v[1] * v[1] * v[2] * v[2]);



  if (sane_check (speed))
  {
    Error ("hydro_velocity:sane_check v %e %e %e\n", v[0], v[1], v[2]);
  }
  return (speed);

}


/**********************************************************/
/**
 * @brief	Interpolates on supplied data file to get cell densities
 *
 * @param [in] x			The x location we need a density for
 * @return 	rrho			The density at the requested loc ation
 *
 * This is a routine which allows the density in a supplied
 * hydro-dynamic snapshot to be used to populate a sirocco
 * model. Often the grids overlap, but this should work whatever
 * the case.
 *
 *
 * ###Notes###
 * 12/99	-	Written by KSL
 * 12/04	-	Mod to prevent divide by zero when r=0.0, and NaN when xxx>1.0;
 * 04/13	-	Heavily modified by NSH
***********************************************************/



double
hydro_rho (x)
     double x[];
{
  double length ();
  int ii, jj;
  int im, jm;
  double r, theta;
  double rrho;
  double f1, f2;
  r = length (x);
  theta = asin (sqrt (x[0] * x[0] + x[1] * x[1]) / r);

  if ((hydro_r_cent[ihydro_r] - r) / r < -1e-6)
  {
    Log (" r outside hydro grid in hydro_rho %e > %e %e\n", r, hydro_r_cent[ihydro_r], (hydro_r_cent[ihydro_r] - r) / r);
    rrho = 1.e-23;
    return (rrho);
  }


  /* The density is a cell centred quantity in zeus */

  im = jm = ii = jj = 0;
  f1 = f2 = 0.0;
  hydro_frac (r, hydro_r_cent, ihydro_r, &im, &ii, &f1);
  hydro_frac (theta, hydro_theta_cent, ihydro_theta, &jm, &jj, &f2);

  rrho = hydro_interp_value (rho_input, im, ii, jm, jj, f1, f2);

  if (rrho < 1e-23)
    rrho = 1e-23;

  return (rrho);
}



/**********************************************************/
/**
 * @brief	calculates the wind temperature at any position x in cartesian
 * coordiantes
 *
 * @param [in] x			The x location we need a temperature for
 * @return 	temp			The temperature at the requested loc ation
 *
 * This is a routine which allows the temperature in a supplied
 * hydro-dynamic snapshot to be used to populate a sirocco
 * model. Often the grids overlap, but this should work whatever
 * the case.
 *
 *
 * ###Notes###
 * 04/13	-	NSH - began work
 * 01/15	-	NSH - temperature is now computed by a helper script, also now common code removed to subroutines
***********************************************************/



double
hydro_temp (x)
     double x[];
{
  double length ();
  int ii, jj;
  int im, jm;
  double r, theta, temp;
  double f1, f2;
  r = length (x);
  theta = asin (sqrt (x[0] * x[0] + x[1] * x[1]) / r);


  im = jm = ii = jj = 0;
  f1 = f2 = 0.0;


  if ((hydro_r_cent[ihydro_r] - r) / r < -1e-6)
  {
    Log (" r outside hydro grid in hydro_temp %e > %e %e\n", r, hydro_r_cent[ihydro_r], (hydro_r_cent[ihydro_r] - r) / r);
    temp = 1e4;
    return (temp);
  }


  /* The density is a cell centred quantity in zeus */



  hydro_frac (r, hydro_r_cent, ihydro_r, &im, &ii, &f1);
  hydro_frac (theta, hydro_theta_cent, ihydro_theta, &jm, &jj, &f2);

  temp = hydro_interp_value (temp_input, im, ii, jm, jj, f1, f2);



  return (temp);
}







/**********************************************************/
/**
 * @brief	Defines an r-theta grid to accept a Zeus hydro model
 *
 * @param [in] w			The wind pointer
 * @param [in] ndom			The number of the domain which will be populated by the geometry
 * @return 					0 if successful
 *
 * rtheta_make_zeus_grid defines the cells in a rtheta grid based
 * upon the coordinates that can be read in from a zeus
 * (the hydrocode used by Proga for some simulations)
 *
 * ###Notes###
***********************************************************/



int
rtheta_make_hydro_grid (int ndom, WindPtr w)
{
  double theta, thetacen, dtheta;
  int i, j, n;
  int ndim, mdim;

  ndim = zdom[ndom].ndim;
  mdim = zdom[ndom].mdim;


  for (i = 0; i < ndim; i++)
  {
    for (j = 0; j < mdim; j++)
    {
      wind_ij_to_n (ndom, i, j, &n);
      w[n].inwind = W_ALL_INWIND;
      if (i == 0)               // The inner edge of the grid should be geo.rstar
      {
        w[n].r = geo.rstar;     //So we set the inner edge to be the stellar (or QSO) radius
        w[n].rcen = (geo.rstar + hydro_r_edge[0]) / 2.0;        //It will be a big cell
        w[n].inwind = W_NOT_INWIND;
      }
      else if (i - 1 > ihydro_r)        // We are at the radial limit of the data, this last cell will be a ghost cell of our own
      {
        w[n].r = zdom[ndom].rmax;       // so we set the outer cell to be the edge of the wind
        w[n].rcen = zdom[ndom].rmax;    // And it has zero volume
        w[n].inwind = W_NOT_INWIND;
      }
      else
      {
        w[n].r = hydro_r_edge[i - 1];   //The -1 is to take account of the fact that i=0 cell is the inner part of the geometry inside the wind
        w[n].rcen = hydro_r_cent[i - 1];
      }
      dtheta = 2.0 * (hydro_theta_cent[j] - hydro_theta_edge[j]);



      if (hydro_theta_cent[j] + (dtheta / 2.0) > hydro_thetamax && hydro_theta_cent[j] - (dtheta / 2.0) < hydro_thetamax)       //This cell bridges the boundary - we reset it so the lower edge is at the disk boundary
      {
        theta = hydro_theta_edge[j];
        thetacen = ((hydro_thetamax + theta) / 2.0);
      }

      else if (j > ihydro_theta)        //We are setting up a cell past where there is any data
      {
        thetacen = hydro_thetamax;      //Set the center and the edge to the maximum extent of the data/interest
        theta = hydro_thetamax;
        w[n].inwind = W_NOT_INWIND;
      }
      else
      {
        thetacen = hydro_theta_cent[j];
        theta = hydro_theta_edge[j];
      }
      w[n].theta = theta * RADIAN;
      w[n].thetacen = thetacen * RADIAN;
      w[n].x[1] = w[n].xcen[1] = 0.0;
      w[n].x[0] = w[n].r * sin (theta);
      w[n].x[2] = w[n].r * cos (theta);
      w[n].xcen[0] = w[n].rcen * sin (thetacen);
      w[n].xcen[2] = w[n].rcen * cos (thetacen);

    }
  }

  /* Now set up the wind cones that are needed for calculating ds in a cell */

  rtheta_make_cones (ndom, w);


  /* OK finished successfuly */
  return (0);

}


/**********************************************************/
/**
 * @brief	Computes the volume for a cell in an r-theta hydro model
 *
 * @param [in,out] WindPtr  w   a single wind cell to calculate the volume for
 *
 * @return 					0 if successful
 *
 *   rtheta_zeus_volumes replaces rtheta_volumes for a zeus model.
 * We know wether cells are in the wimnd
 * so all we need to do is work out the volumes.
 *
 *
 * ###Notes###
***********************************************************/


int
rtheta_hydro_cell_volume (WindPtr w)
{
  int i, j;
  int ndom;
  double rmin, rmax, thetamin, thetamax;
  DomainPtr one_dom;

  ndom = w->ndom;
  one_dom = &zdom[ndom];
  wind_n_to_ij (ndom, w->nwind, &i, &j);

  if (w->inwind == W_ALL_INWIND)
  {

    rmin = one_dom->wind_x[i];
    rmax = one_dom->wind_x[i + 1];
    thetamin = one_dom->wind_z[j] / RADIAN;
    thetamax = one_dom->wind_z[j + 1] / RADIAN;

    //leading factor of 2 added to allow for volume above and below plane (SSMay04)
    w->vol = 2. * 2. / 3. * PI * (rmax * rmax * rmax - rmin * rmin * rmin) * (cos (thetamin) - cos (thetamax));

    if (w->vol == 0.0)
    {
      Log ("Found wind cell (%i) with no volume (%e) in wind, resetting\n", w->nwind, w->vol);
      w->inwind = W_NOT_INWIND;
    }
  }
  else
  {
    w->vol = 0.0;
  }

  return (0);
}




/**
 * @brief	Helper routine to interpolate on a grid
 *
 * @param [in] coord			The coordinate we want a value for
 * @param [in] coord_array		The array of grid points
 * @param [in] imax				The maximum extent of the array
 * @param [in] *cell1			One of the two bracketing cells (computed)
 * @param [in] *cell2			The other bracketing cell (computed)
 * @param [in] *frac			The fraction along the coord between the two cells (computed)
 * @return 					0 if successful

 *
 * hydro_frac replaces many lines of identical code, all of which
 * interpolate on the input grid to get value for the sirocco grid.
 * This code find out the fraction along a coord array where the centre
 * or edge of the sirocco grid cell lies on the hydro grid
 *
 *
 * ###Notes###
 *
 * ### Programming Comment ###
 * This is rather similar to the routine coord_frac, but that
 * loops over all the values in a domain. Here, because of the
 * way wind parameters are set up on a cell by cell basis, in
 * order to keep things the same between hydro and other wind
 * models we have routines that are called for a single cell.
***********************************************************/


int
hydro_frac (coord, coord_array, imax, cell1, cell2, frac)
     double coord;
     double coord_array[];
     int imax;
     int *cell1, *cell2;
     double *frac;
{
  int ii;
  ii = 0;
  *cell1 = 0;
  *cell2 = 0;

  while (coord_array[ii] < coord && ii < imax)  //Search through array, until array value is greater than your coordinate
    ii++;


  if (ii > imax)
  {                             // r is greater than anything in Proga's model
    *frac = 1;                  //If we are outside the model set fractional position to 1
    *cell1 = *cell2 = imax;     //And the bin to the outermost
    return (0);
  }
  else if (ii == 0)
  {
    *frac = 1;                  //Otherwise, we must be inside the innermost bin, so again, set fration to 1
    *cell1 = *cell2 = 0;        // And the bin to the innermost. Lines below to the same for theta.
    return (0);
  }
  else if (ii > 0)
  {                             //r is in the normal range

    *frac = (coord - coord_array[ii - 1]) / (coord_array[ii] - coord_array[ii - 1]);    //Work out fractional position in the ii-1th radial bin where you want to be
    *cell1 = ii - 1;            //This is the radial bin below your value of r
    *cell2 = ii;
    return (0);
  }

  return (0);
}



/**
 * @brief	Helper routine to interpolate on a grid
 *
 * @param [in] array			the array of input data
 * @param [in] im,ii			the two cells surrounding the cell in the first dim (r)
 * @param [in] jm,jj			the two cells surrounding the cell in the second dim (theta)
 * @param [in] f1,f2			the fraction between the two values in first and second dim
 * @return 		value			the interplated number

 *
 * hydro_interp_value replaces many lines of identical code, all of which
 * interpolate on the input grid to get value for the sirocco grid.
 * This returns the actuak value of an interpolated array
 *
 *
 * ###Notes###
 *
 * ### Programming Comment ###
 * This called after hydro_interp_frac and the programming
 * comment for that routine applies here too.
***********************************************************/



double
hydro_interp_value (array, im, ii, jm, jj, f1, f2)
     double array[];
     int im, ii;                //the two cells surrounding the cell in the first dim (r)
     int jm, jj;                //the two cells surrounding the cell in the second dim (theta)
     double f1, f2;             //the fraction between the two values in first and second dim
{
  double value;
  double d1, d2;


  d1 = array[im * MAXHYDRO + jm] + f1 * (array[ii * MAXHYDRO + jm] - array[im * MAXHYDRO + jm]);
  d2 = array[im * MAXHYDRO + jj] + f1 * (array[ii * MAXHYDRO + jj] - array[im * MAXHYDRO + jj]);
  value = d1 + f2 * (d2 - d1);


  return (value);
}





/**
 * @brief	Allows a previous windsave to be used in a new hydro simulation
 *
 * @param [in] ndom			The number of the domain which will be populated by the geometry
 * @return 					0 if successful

 *
 * hydro_restart is a subroutine which permits a previous wind save
 * file to be used in a hydro simulation. The density and temperature
 * for each cell are those from the hydro simulation. The ion fractions
 * are taken from the windsave file. The routine therefore changes *most*
 * things in the wind - only the ion fractions remain
 *
 *
 * ###Notes###
***********************************************************/





int
hydro_restart (ndom)
     int ndom;
{
  int n, nion;
  int nwind;
  double x[3];
  double old_density;
//OLD  int nstart, nstop, ndim2;
  int nstart, nstop;

  zdom[ndom].wind_type = 3;     //Temporarily set the wind type to hydro, so we can use the normal routines
  /* note that we will have passed the hydro domain number as default */
  nstart = zdom[ndom].nstart;
  nstop = zdom[ndom].nstop;
//OLD  ndim2 = zdom[ndom].ndim2;

  for (n = nstart; n < nstop; n++)
  {
    model_velocity (ndom, wmain[n].x, wmain[n].v);
    model_vgrad (ndom, wmain[n].x, wmain[n].v_grad);
  }


  for (nwind = zdom[ndom].nstart; nwind < zdom[ndom].nstop; nwind++)
  {
    if (wmain[nwind].vol > 0.0)
    {
      n = wmain[nwind].nplasma;
      stuff_v (wmain[nwind].xcen, x);
      old_density = plasmamain[n].rho;
      plasmamain[n].rho = model_rho (ndom, x) / zdom[ndom].fill;
      plasmamain[n].t_r = plasmamain[n].t_e = hydro_temp (x);

      for (nion = 0; nion < nions; nion++)      //Change the absolute number densities, fractions remain the same
      {
        plasmamain[n].density[nion] = plasmamain[n].density[nion] * (plasmamain[n].rho / old_density);
      }

      plasmamain[n].ne = get_ne (plasmamain[n].density);        //get the new electron density
      partition_functions (&plasmamain[n], NEBULARMODE_LTE_GROUND);     //set the level populations to ground state - this is because at the moment we dont know how to work out levels for cases that dont have a dilute BB radiation field. We need to set them to something however. Could do better in the future.

    }
  }
  /* Recreate the wind cones because these are not part of the windsave file */
  rtheta_make_cones (ndom, wmain);
  return (0);
}

/**********************************************************/
/**
 * @brief
 *
 * @param
 *
 *
 * ###Notes###
 *
***********************************************************/

void
create_hydro_output_files (void)
{
  FILE *fptr;
  FILE *fptr2;
  FILE *fptr3;
  FILE *fptr4;
  FILE *fptr5;

  int i;
  int ii;
  int j;
  int nwind;
  int nplasma;

  double v_th;
  double kappa_es;
  double t_opt;
  struct photon ptest;
  double t_UV;
  double t_Xray;
  double fhat[4];

  WindPtr w = wmain;
  double vol;

  Log ("Outputting heatcool file for connecting to zeus\n");
  fptr = fopen ("py_heatcool.dat", "w");
  fptr2 = fopen ("py_flux.dat", "w");
  fptr3 = fopen ("py_ion_data.dat", "w");
  fptr4 = fopen ("py_spec_data.dat", "w");
  fptr5 = fopen ("py_pcon_data.dat", "w");

  fprintf (fptr,
           "i j rcen thetacen vol temp xi ne heat_xray heat_comp heat_lines heat_ff cool_comp cool_lines cool_ff rho n_h rad_f_w rad_f_phi rad_f_z bf_f_w bf_f_phi bf_f_z\n");
  fprintf (fptr2, "i j F_vis_x F_vis_y F_vis_z F_UV_theta F_UV_phi F_UV_r F_Xray_x F_Xray_y F_Xray_z\n");       //directional flux by band

  fprintf (fptr3, "nions %i\n", nions);
  for (i = 0; i < nions; i++)
  {
    fprintf (fptr3, "ion %i %s %i %i\n", i, ele[ion[i].nelem].name, ion[i].z, ion[i].istate);
  }
  fprintf (fptr3, "nplasma %i\n", NPLASMA);

  fprintf (fptr4, "model %i\n", geo.ioniz_mode);
  fprintf (fptr4, "nbands %i\n", geo.nxfreq);
  fprintf (fptr4, "nplasma %i\n", NPLASMA);
  for (i = 0; i < geo.nxfreq + 1; i++)
    fprintf (fptr4, "%e ", geo.xfreq[i]);       //hard wired band edges
  fprintf (fptr4, "\n ");

  fprintf (fptr5, "nplasma %i\n", NPLASMA);

  for (nwind = zdom[geo.hydro_domain_number].nstart; nwind < zdom[geo.hydro_domain_number].nstop; nwind++)
  {
    if (wmain[nwind].inwind >= 0)
    {
      nplasma = wmain[nwind].nplasma;
      wind_n_to_ij (geo.hydro_domain_number, plasmamain[nplasma].nwind, &i, &j);
      i = i - 1;                //There is a radial 'ghost zone' in sirocco, we need to make our i,j agree with zeus
      vol = w[plasmamain[nplasma].nwind].vol;
      fprintf (fptr, "%d %d %e %e %e ", i, j, w[plasmamain[nplasma].nwind].rcen, w[plasmamain[nplasma].nwind].thetacen / RADIAN, vol);  //output geometric things
      fprintf (fptr, "%e %e %e ", plasmamain[nplasma].t_e, plasmamain[nplasma].xi, plasmamain[nplasma].ne);     //output temp, xi and ne to ease plotting of heating rates
      fprintf (fptr, "%e ", (plasmamain[nplasma].heat_photo + plasmamain[nplasma].heat_auger) / vol);   //Xray heating - or photoionization
      fprintf (fptr, "%e ", (plasmamain[nplasma].heat_comp) / vol);     //Compton heating
      fprintf (fptr, "%e ", (plasmamain[nplasma].heat_lines) / vol);    //Line heating 28/10/15 - not currently used in zeus
      fprintf (fptr, "%e ", (plasmamain[nplasma].heat_ff) / vol);       //FF heating 28/10/15 - not currently used in zeus
      fprintf (fptr, "%e ", (plasmamain[nplasma].cool_comp) / vol);     //Compton cooling
      fprintf (fptr, "%e ", (plasmamain[nplasma].lum_lines + plasmamain[nplasma].cool_rr + plasmamain[nplasma].cool_dr) / vol); //Line cooling must include all recombination cooling
      fprintf (fptr, "%e ", (plasmamain[nplasma].lum_ff) / vol);        //ff cooling
      fprintf (fptr, "%e ", plasmamain[nplasma].rho);   //density
      fprintf (fptr, "%e ", plasmamain[nplasma].rho * rho2nh);  //hydrogen number density
      fprintf (fptr, "%e ", plasmamain[nplasma].rad_force_es[0]);       //electron scattering radiation force in the w(x) direction
      fprintf (fptr, "%e ", plasmamain[nplasma].rad_force_es[1]);       //electron scattering radiation force in the phi(rotational) directionz direction
      fprintf (fptr, "%e ", plasmamain[nplasma].rad_force_es[2]);       //electron scattering radiation force in the z direction
      fprintf (fptr, "%e ", plasmamain[nplasma].rad_force_bf[0]);       //bound free scattering radiation force in the w(x) direction
      fprintf (fptr, "%e ", plasmamain[nplasma].rad_force_bf[1]);       //bound free scattering radiation force in the phi(rotational) direction
      fprintf (fptr, "%e \n", plasmamain[nplasma].rad_force_bf[2]);     //bound free scattering radiation force in the z direction
      fprintf (fptr2, "%d %d ", i, j);  //output geometric things
      fprintf (fptr2, "%e %e %e ", plasmamain[nplasma].F_vis[0], plasmamain[nplasma].F_vis[1], plasmamain[nplasma].F_vis[2]);   //directional flux by band
      fprintf (fptr2, "%e %e %e ", plasmamain[nplasma].F_UV[0], plasmamain[nplasma].F_UV[1], plasmamain[nplasma].F_UV[2]);      //directional flux by band
      fprintf (fptr2, "%e %e %e ", plasmamain[nplasma].F_Xray[0], plasmamain[nplasma].F_Xray[1], plasmamain[nplasma].F_Xray[2]);        //directional flux by band

      fprintf (fptr2, "\n");
      fprintf (fptr3, "%d %d ", i, j);  //output geometric things
      for (ii = 0; ii < nions; ii++)
        fprintf (fptr3, "%e ", plasmamain[nplasma].density[ii]);
      fprintf (fptr3, "\n");

      fprintf (fptr4, "%d %d ", i, j);  //output geometric things
      for (ii = 0; ii < geo.nxfreq; ii++)
        fprintf (fptr4, "%e %e %i %e %e %e %e ",
                 plasmamain[nplasma].fmin_mod[ii], plasmamain[nplasma].fmax_mod[ii], plasmamain[nplasma].spec_mod_type[ii],
                 plasmamain[nplasma].pl_log_w[ii], plasmamain[nplasma].pl_alpha[ii], plasmamain[nplasma].exp_w[ii],
                 plasmamain[nplasma].exp_temp[ii]);
      fprintf (fptr4, "\n ");


      //We need to compute the g factor for this cell and output it.

      v_th = pow ((2. * BOLTZMANN * plasmamain[nplasma].t_e / MPROT), 0.5);     //We need the thermal velocity for hydrogen
      stuff_v (w[plasmamain[nplasma].nwind].xcen, ptest.x);     //place our test photon at the centre of the cell
      ptest.grid = nwind;       //We need our test photon to know where it is
      kappa_es = THOMPSON * plasmamain[nplasma].ne / plasmamain[nplasma].rho;

      //First for the optical band (up to 4000AA)
      if (length (plasmamain[nplasma].F_vis) > 0.0)     //Only makes sense if flux in this band is non-zero
      {
        stuff_v (plasmamain[nplasma].F_vis, fhat);
        renorm (fhat, 1.);      //A unit vector in the direction of the flux - this can be treated as the lmn vector of a pretend photon
        stuff_v (fhat, ptest.lmn);      //place our test photon at the centre of the cell
        t_opt = kappa_es * plasmamain[nplasma].rho * v_th / fabs (dvwind_ds_cmf (&ptest));
      }
      else
        t_opt = 0.0;            //Essentually a flag that there is no way of computing t (and hence M) in this cell.

      //Now for the UV band (up to 4000AA->100AA)
      if (length (plasmamain[nplasma].F_UV) > 0.0)      //Only makes sense if flux in this band is non-zero
      {
        stuff_v (plasmamain[nplasma].F_UV, fhat);
        renorm (fhat, 1.);      //A unit vector in the direction of the flux - this can be treated as the lmn vector of a pretend photon
        stuff_v (fhat, ptest.lmn);      //place our test photon at the centre of the cell
        t_UV = kappa_es * plasmamain[nplasma].rho * v_th / fabs (dvwind_ds_cmf (&ptest));
      }
      else
        t_UV = 0.0;             //Essentually a flag that there is no way of computing t (and hence M) in this cell.


      //And finally for the Xray band (up to 100AA and up)
      if (length (plasmamain[nplasma].F_Xray) > 0.0)    //Only makes sense if flux in this band is non-zero
      {
        stuff_v (plasmamain[nplasma].F_Xray, fhat);
        renorm (fhat, 1.);      //A unit vector in the direction of the flux - this can be treated as the lmn vector of a pretend photon
        stuff_v (fhat, ptest.lmn);      //place our test photon at the centre of the cell
        t_Xray = kappa_es * plasmamain[nplasma].rho * v_th / fabs (dvwind_ds_cmf (&ptest));
      }
      else
        t_Xray = 0.0;           //Essentually a flag that there is no way of computing t (and hence M) in this cell.

      fprintf (fptr5, "%i %i %e %e %e %e %e %e %e\n", i, j, plasmamain[nplasma].t_e, plasmamain[nplasma].rho,
               plasmamain[nplasma].rho * rho2nh, plasmamain[nplasma].ne, t_opt, t_UV, t_Xray);
    }
  }
  fclose (fptr);
  fclose (fptr2);
  fclose (fptr3);
  fclose (fptr4);
  fclose (fptr5);
}
