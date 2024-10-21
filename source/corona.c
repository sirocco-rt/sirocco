
/***********************************************************/
/** @file  corona.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  Routines specific to setting up disk corona defined as
 * a region that rotates witht he underlying disk and has an 
 * expoential density distribution
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/** 
 * @brief      gets input data which is necessary to describe a corona above
 * 	the surface of the disk.
 *
 * @param [in] int  ndom   The domain number for the corona
 * @return  Always returns 0     
 *
 * @details
 * The routine initializes variables that are needed to define a corona (to
 * values that are appropriate for a WD system) and then queries the 
 * .pf or the user for the same variables, which are then stored in one
 * of the Domain elements
 *
 * The corona is assumed to exist between rmin and rmax, and between 0
 * and zmax.  Within this region it has a density specified by the density
 * at the base and an exponential scale height.
 *
 * Azimuthally the velocity will be set to the velocity of the disk, but one
 * can give it a outward velocity which is specifed as a fraction of the
 * azimuthal velocity
 *
 * ### Notes ###
 * The parameters obtained here are only used in the routines in corona.c
 * Initially, we define it as a gaussian ring with and exponetial density distribution.
 *
 * The terminology here is a bit confusing because rmin and rmax as read in to the routine
 * refer actually to the cylindrical rho, which is easy to get confused with other parameters like 
 * zdom[ndom].rmax
 * which really are radii for the central source.
 *
 **********************************************************/

int
get_corona_params (ndom)
     int ndom;
{
  Log ("Creating a corona above a disk\n");

  /* Start with reasonable values for a WD system */

  zdom[ndom].wind_thetamin = 0.0;
  zdom[ndom].wind_thetamax = 0.0;
  zdom[ndom].rmin = geo.rstar;

  zdom[ndom].corona_rmin = 1.e10;
  zdom[ndom].corona_rmax = 2.e10;
  zdom[ndom].corona_zmax = 2.e09;
  zdom[ndom].corona_base_density = 1.e13;
  zdom[ndom].corona_scale_height = 1.e9;

  rddoub ("Corona.radmin(cm)", &zdom[ndom].corona_rmin);        /*Radius where corona begins */
  if (zdom[ndom].corona_rmin < geo.rstar)
  {
    Error ("get_corona_params: It is unreasonable to have the corona start inside the star!\n");
    Log ("Setting geo.corona_rmin to geo.rstar\n");
    zdom[ndom].corona_rmin = geo.rstar;
  }
  rddoub ("Corona.radmax(cm)", &zdom[ndom].corona_rmax);        /*Radius where corona ends */
  rddoub ("Corona.zmax(cm)", &zdom[ndom].corona_zmax);  /*Veritical heighe where corona ends */
  rddoub ("Corona.base_den(cgs)", &zdom[ndom].corona_base_density);     /*Density at the base of the corona */
  rddoub ("Corona.scale_height(cm)", &zdom[ndom].corona_scale_height);  /*Scale height of corona */
  rddoub ("Corona.vel_frac", &zdom[ndom].corona_vel_frac);      /*fractional radial velocity of corona */

  zdom[ndom].rmin = zdom[ndom].corona_rmin;
  /* rmax here is the defines a radius beyond which this region does not exist, if the vertical height is large
   * compared to the horizontal size then one needs to include both */
  zdom[ndom].rmax = sqrt (zdom[ndom].corona_rmax * zdom[ndom].corona_rmax + zdom[ndom].corona_zmax * zdom[ndom].corona_zmax);
  zdom[ndom].zmax = zdom[ndom].corona_zmax;
  zdom[ndom].wind_rhomin_at_disk = zdom[ndom].corona_rmin;
  zdom[ndom].wind_rhomax_at_disk = zdom[ndom].corona_rmax;
  zdom[ndom].wind_thetamin = 0.0;
  zdom[ndom].wind_thetamax = 0.0;

  /* Set up wind planes for a layer with a specific height */

  zdom[ndom].windplane[0].x[0] = zdom[ndom].windplane[0].x[1] = zdom[ndom].windplane[0].x[2] = 0;
  zdom[ndom].windplane[0].lmn[0] = zdom[ndom].windplane[0].lmn[1] = 0;
  zdom[ndom].windplane[0].lmn[2] = 1;

  zdom[ndom].windplane[1].x[0] = zdom[ndom].windplane[0].x[1] = 0;
  zdom[ndom].windplane[1].x[2] = zdom[ndom].corona_zmax;
  zdom[ndom].windplane[1].lmn[0] = zdom[ndom].windplane[0].lmn[1] = 0;
  zdom[ndom].windplane[1].lmn[2] = 1;




  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
  {
    zdom[ndom].xlog_scale = 0.3 * zdom[ndom].corona_rmin;
    if (zdom[ndom].corona_scale_height < zdom[ndom].corona_zmax)
      zdom[ndom].zlog_scale = 0.3 * zdom[ndom].corona_scale_height;
    else if (zdom[ndom].mdim > 0)
    {
      zdom[ndom].zlog_scale = zdom[ndom].corona_zmax / zdom[ndom].mdim;
    }
    else
    {
      Error ("corona: Cannot define z coordinates unless zdom[ndom].mdim is defined. Aborting\n");
      Exit (0);
    }
  }


  zdom[ndom].wind_rhomin_at_disk = zdom[ndom].corona_rmin;
  zdom[ndom].wind_rhomax_at_disk = zdom[ndom].corona_rmax;
  zdom[ndom].wind_thetamin = 0;
  zdom[ndom].wind_thetamax = 0;

  return (0);
}



/**********************************************************/
/** 
 * @brief      double (x,v) calulates the v the wind at a position r
 * 	x
 *
 * @param [in] int  ndom   The domain number in which the corona is descibed
 * @param [in] double  x[]   the position where for the which one desires the velocity
 * @param [out] double  v[]   The calculated velocity in cartesian coordineates
 * @return     double v[]	The amplitude of the velocity is returned
 *
 * @details
 * Given a position, the azimuthal (v_phi) velocity is set to the velocity
 * of the underlying disk.  The outward (v_rho) velocity is set to a 
 * fraction of the disk velocity as specifed by corona_vel_frac
 *
 * ### Notes ###
 * The routine is valid even if x is not in the xz plane.
 *
 **********************************************************/

double
corona_velocity (ndom, x, v)
     int ndom;
     double x[], v[];
{
  double rho, speed;
  double xtest[3];

  rho = sqrt (x[0] * x[0] + x[1] * x[1]);
  if (rho > 0.0)
    speed = sqrt (GRAV * geo.mstar / rho);
  else
    speed = sqrt (GRAV * geo.mstar / geo.rstar);

  v[0] = -zdom[ndom].corona_vel_frac * speed;
  v[2] = 0.0;
  v[1] = speed;

  /* At this point we have calculated the velocity in the xz plane, which
   * is identical to the statement that we have calculated it in
   * cylindrical coordinates.  The next bit projects back to xyz 
   * coordinates if x was not originally in the xz plane.
   */
  if (x[1] != 0.0)
  {
    project_from_cyl_xyz (x, v, xtest);
    stuff_v (xtest, v);
  }


  return (speed);

}


/**********************************************************/
/** 
 * @brief      double (x) calculates the density of a corona at a position x
 *
 * @param [in] int  ndom   The domain where the corona is described
 * @param [in] double  x[]   the position at which one desires the denisty
 * @return     The density at x is returned in gram/cm**3
 *
 * @details
 * The function is simply and exponential distribution with a constant
 *         scale height
 *
 * ### Notes ###
 * The routine does not check whether the x is within the region
 * described by the corona
 *
 * @bug This and the other routines that describe a coronal model
 * do not handle the case of the vertically extended disk
 *
 **********************************************************/

double
corona_rho (ndom, x)
     int ndom;
     double x[];
{
  double rho;

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    Error ("corona_rho: Quitting. Need to think more about coronal model more with vertically extended disk\n");
    Exit (0);
  }

  rho = zdom[ndom].corona_base_density * exp (-(fabs (x[2])) / zdom[ndom].corona_scale_height);

  if (rho < 1.e-10)
    rho = 1.e-10;               // A floor to the density appears to be needed for some of the
  // ionization calculations

  rho = rho / rho2nh;

  return (rho);
}
