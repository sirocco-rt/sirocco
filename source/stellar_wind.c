
/***********************************************************/
/** @file  stellar_wind.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  Routines needed to define a spherical wind with a Castor and 
 * Larmers (1979) velocity law
 *
 * The Castor and Lamers (1979)  velocity law has
 *
 * v(r)=V_o + (V_infinity-V_o) (1-R/r)**beta
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/** 
 * @brief      gets input data which is necessary for a Castor and Larmors'
 * 	(1979) description of an O star wind, 
 *
 * @param [in] int  ndom   The domain number
 * @return     Always returns zero
 *
 *
 * ###Notes###
 *
 * Although it seems natural to assume that the wind starts at the photosphere, this
 * is not required (and in Mauche and Raymond's spherical model for a CV wind they set
 * the inner radius of the wind at several WD radii.  That possibility is allowed for
 *
 * 	
 *
 **********************************************************/

int
get_stellar_wind_params (ndom)
     int ndom;
{
  Log ("Creating a wind model for a Star\n");

  /* Initialize some of the parameters */
  zdom[ndom].stellar_wind_mdot = 1.e-6;
  zdom[ndom].rmin = geo.rstar;
  zdom[ndom].cl_beta = 1.0;
  zdom[ndom].cl_rmin = 2.8e9;
  zdom[ndom].cl_v_zero = 200e5;
  zdom[ndom].cl_v_infinity = 3000e5;

  /* Get the parameters needed to define the wind */

  rddoub ("Stellar_wind.mdot(msol/yr)", &zdom[ndom].stellar_wind_mdot);
  zdom[ndom].stellar_wind_mdot *= MSOL / YR;

  rddoub ("Stellar_wind.radmin(cm)", &zdom[ndom].rmin); /*Radius where wind begins */
  if (zdom[ndom].rmin < geo.rstar)
  {
    Error ("get_stellar_wind_params: It is unreasonable to have the wind start inside the star!\n");
    Log ("Setting zdom[ndom].rmin to geo.rstar\n");
    zdom[ndom].rmin = geo.rstar;
  }
  zdom[ndom].rmax = 10 * zdom[ndom].rmin;
  rddoub ("Stellar_wind.radmax(cm)", &zdom[ndom].rmax); /*Radius where wind ends */
  zdom[ndom].cl_rmin = zdom[ndom].rmin;
  rddoub ("Stellar_wind.vbase(cm)", &zdom[ndom].cl_v_zero);     /* Velocity at base of the wind */
  rddoub ("Stellar_wind.v_infinity(cm)", &zdom[ndom].cl_v_infinity);    /* Final speed of wind in units of escape velocity */

  rddoub ("Stellar_wind.acceleration_exponent", &zdom[ndom].cl_beta);   /* Accleration scale exponent */

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

  return (0);
}





/**********************************************************/
/** 
 * @brief      Calbulate the wind velocity at a position
 *
 * @param [in] int  ndom   The domain number
 * @param [in] double  x[]   the postion (in cartesian coordinates)
 * @param [out] double  v[]  the resulting velocity (in cartesain coordiantes
 * @return   the speed at the position               
 * 	
 *
 * The model is that of Castor and Lamars 1979  where
 * 
 * v(r)=V_o + (V_infinity-V_o) (1-R/r)**beta
 * 	
 * The values of the individiual constants should all be part of the structure zdom[ndom].
 * 
 * - V_o:  			zdom[ndom].cl_v_zero;		velocity at base of wind 
 * - V_infinity:		zdom[ndom].cl_v_infinity;	the velocity at infinity
 * - R:			zdom[ndom].cl_rmin;	       	the inner radius of the wind
 * - beta:			zdom[ndom].cl_beta;		power law exponent for the velocity law
 *
 * ###Notes###
 *
 * v is set to V_o inside of rstar or w_rscale, even though the should not really exist there!
 * 
 * Since this is a radial wind, it is not affected by the disk at all. In particular there
 * are no modifications to this routine if the disk is vertically extended.
 *
 **********************************************************/

double
stellar_velocity (ndom, x, v)
     int ndom;
     double x[], v[];
{
  double r, speed, zzz;
  double length ();

  if ((r = length (x)) == 0.0)
  {
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    return (0.0);
  }


  if (r <= geo.rstar || r <= zdom[ndom].cl_rmin)
    speed = zdom[ndom].cl_v_zero;
  else
  {
    zzz = pow (1. - zdom[ndom].cl_rmin / r, zdom[ndom].cl_beta);
    speed = zdom[ndom].cl_v_zero + (zdom[ndom].cl_v_infinity - zdom[ndom].cl_v_zero) * zzz;
  }
  v[0] = speed * x[0] / r;
  v[1] = speed * x[1] / r;
  v[2] = speed * x[2] / r;

  return (speed);

}


/**********************************************************/
/** 
 * @brief      Calculate the density of an stellar_wind at a position x
 *
 * @param [in] int  ndom   The domain number
 * @param [in] double  x[]   the position where for the which one desires the density
 * @return     The density at x is returned in gram/cm**3
 *
 * rho=mdot/(4PI*r*r*v);
 *
 * ###Notes###
 *
 *
 **********************************************************/

double
stellar_rho (ndom, x)
     int ndom;
     double x[];
{
  double r, rho, v[3];

  r = length (x);
  if (r < zdom[ndom].cl_rmin)
  {
    rho = zdom[ndom].stellar_wind_mdot / (4. * PI * r * r * zdom[ndom].cl_v_zero);
  }
  else
  {
    rho = zdom[ndom].stellar_wind_mdot / (4. * PI * r * r * stellar_velocity (ndom, x, v));
  }

  return (rho);
}
