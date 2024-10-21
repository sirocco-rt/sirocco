
/***********************************************************/
/** @file  homologous.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  Routines describing a homologous flow in Python.
 *
 * A homlogous flow is a flow in which the velocity is proportional
 * to the distance from the central source.  In this case, the
 * density is taken to have the form of a power law.
 * 
 * ###Notes###
 * Homlogous flows were added to allow comparisons with SN codes,
 * such as Tardis
 * 
 * These routines were largely adapted from those associated with stellar
 * winds, and some of the variables use those that are associated with stellar
 * winds.  With domains, it would be clearer to give them there own variable
 * names.  Note that the maximum radius is not, as it should be defined here, but
 * rather relies on zdmo[ndom].rmax.  This should be fixed.
 * 
 ***********************************************************/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"




/**********************************************************/
/** 
 * @brief      gets input data which is necessary for a homologous expansion law
 *    		v(r)= r / t0
 *
 * @param [in] int  ndom   The domain number
 * @return     Always returns 0
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
get_homologous_params (ndom)
     int ndom;
{
  DomainPtr one_dom;
  one_dom = &zdom[ndom];

  Log ("Creating a homolgous wind model in domain %d\n", ndom);


  /* Initialize some of the relevant parameters */


  one_dom->stellar_wind_mdot = 100.;
  one_dom->rmin = geo.rstar;
  one_dom->rmax = 2.2464e15;
  one_dom->cl_v_zero = 200e5;
  one_dom->cl_beta = 7.0;

  one_dom->stellar_wind_mdot /= MSOL / YR;
  rddoub ("Homologous.boundary_mdot(msol/yr)", &one_dom->stellar_wind_mdot);
  one_dom->stellar_wind_mdot *= MSOL / YR;


  rddoub ("Homologous.radmin(cm)", &one_dom->rmin);     /*Radius where wind begins */
  rddoub ("Homologous.radmax(cm)", &one_dom->rmax);     /*Radius where wind ends */

  if (one_dom->rmin < geo.rstar)
  {
    Error ("get_homologous_params: It is unreasonable to have the wind start inside the star!\n");
    Log ("Setting one_dom->rmin to geo.rstar\n");
    one_dom->rmin = geo.rstar;
  }
  one_dom->cl_rmin = one_dom->rmin;

  rddoub ("Homologous.vbase(cm)", &one_dom->cl_v_zero); /* Velocity at base of the wind */
  rddoub ("Homologous.density_exponent", &one_dom->cl_beta);    /* Density law exponent */


  /* Assign the generic parameters for the wind the generic parameters of the wind */
  one_dom->wind_thetamin = 0.0;
  one_dom->wind_thetamax = 90. / RADIAN;

  /* define the the variables that determine the gridding */
  one_dom->wind_rhomin_at_disk = 0;
  /* JM 1710 -- fixed these to use domain value rather than geo value see #305 */
  one_dom->wind_rhomax_at_disk = one_dom->rmax;
  one_dom->zmax = one_dom->rmax;

  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
  {
    one_dom->xlog_scale = 0.3 * geo.rstar;
    one_dom->zlog_scale = 0.3 * geo.rstar;
  }

  return (0);
}




/**********************************************************/
/** 
 * @brief      Calulate the velocity the wind at a position 
 *
 * @param [in, out] int  ndom   The domain number
 * @param [in, out] double  x[]   the position (in cartesian coordinates)                   
 * @param [in, out] double  v[]   the calculated velocity (in cartesian coordinates)
 * @return     The amplitude of the velocity is returned
 * 	
 *
 * The model is homologous expansion, namely
 * 
 * v(r)=vmin r / R
 * 	
 * The values of the individiual constants should all be part of the structure geo.
 * 
 * - Vmin:  			cl_v_zero;		velocity at base of wind 
 * - R:			cl_rmin;	       	the inner radius of the wind
 *
 * ###Notes###
 *
 *
 **********************************************************/

double
homologous_velocity (ndom, x, v)
     int ndom;
     double x[], v[];
{
  double r, speed;

  DomainPtr one_dom;
  one_dom = &zdom[ndom];

  if ((r = length (x)) == 0.0)
  {
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    return (0.0);
  }


  if (r <= geo.rstar || r <= one_dom->cl_rmin)
    speed = one_dom->cl_v_zero;
  else
  {
    speed = r * one_dom->cl_v_zero / one_dom->cl_rmin;
  }
  v[0] = speed * x[0] / r;
  v[1] = speed * x[1] / r;
  v[2] = speed * x[2] / r;

  return (speed);

}



/**********************************************************/
/** 
 * @brief      Calculate the density of a homologous flow at a position
 *
 * @param [in] int  ndom   The domain number
 * @param [in] double  x[]   the position where for the which one desires the density
 * @return     The density at x is returned in gram/cm**3
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

double
homologous_rho (ndom, x)
     int ndom;
     double x[];
{
  double r, rho, length ();
  DomainPtr one_dom;
  one_dom = &zdom[ndom];


  r = length (x);

  rho =
    one_dom->stellar_wind_mdot / (4. * PI * one_dom->cl_rmin * one_dom->cl_rmin *
                                  one_dom->cl_v_zero) * pow ((r / one_dom->cl_rmin), -1. * one_dom->cl_beta);

  return (rho);
}
