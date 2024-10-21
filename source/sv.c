
/***********************************************************/
/** @file  sv.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  The collection of routines needed to define a Shlossman-Vitello wind 
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "sirocco.h"

int sdom;
int sv_zero_r_ndom;


/**********************************************************/
/** 
 * @brief      gets input data which is necessary for a Shlossman & Vitello 
 * 	description of the wind
 *
 * @param [in] int  ndom   The domain number
 * @return     Always returns 0
 *
 *  Gets mdot for this wind component, initialize the parameters needed to define a SV wind, and then request those
 *  parameters as inputs.   Calculate the normalization factor needed to convert the global mass loss rate to the
 *  mass loss rate per unit area.
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
get_sv_wind_params (ndom)
     int ndom;
{
  double windmin, windmax, theta_min, theta_max;
  char answer[LINELENGTH];

  Log ("Creating an SV wind in domain %d\n", ndom);

  zdom[ndom].wind_mdot = 0.1 * geo.disk_mdot / (MSOL / YR);     // Convert to MSOL/YR for easy of data entry
  rddoub ("Wind.mdot(msol/yr)", &zdom[ndom].wind_mdot);
  zdom[ndom].wind_mdot *= MSOL / YR;


  zdom[ndom].sv_thetamin = 20. / RADIAN;
  zdom[ndom].sv_thetamax = 65. / RADIAN;
  zdom[ndom].sv_gamma = 1.;
  /* JM -- this can be altered in advanced mode below */
  zdom[ndom].sv_v_zero_mode = FIXED;
  zdom[ndom].sv_v_zero = 6e5;
  zdom[ndom].sv_r_scale = 7e10; /*Accleration length scale for wind */
  zdom[ndom].sv_alpha = 1.5;    /* Accleration scale exponent */
  zdom[ndom].sv_v_infinity = 3; /* Final speed of wind in units of escape velocity */
  zdom[ndom].sv_lambda = 0.0;   /* Mass loss rate exponent */

  windmin = 4;
  windmax = 12;
  rddoub ("SV.diskmin(units_of_rstar)", &windmin);
  rddoub ("SV.diskmax(units_of_rstar)", &windmax);


  zdom[ndom].sv_rmin = windmin * geo.rstar;
  zdom[ndom].sv_rmax = windmax * geo.rstar;

  theta_min = zdom[ndom].sv_thetamin * RADIAN;
  theta_max = zdom[ndom].sv_thetamax * RADIAN;
  rddoub ("SV.thetamin(deg)", &theta_min);
  rddoub ("SV.thetamax(deg)", &theta_max);
  zdom[ndom].sv_thetamin = theta_min / RADIAN;
  zdom[ndom].sv_thetamax = theta_max / RADIAN;

  rddoub ("SV.mdot_r_exponent", &zdom[ndom].sv_lambda); /* Mass loss rate exponent */
  rddoub ("SV.v_infinity(in_units_of_vescape", &zdom[ndom].sv_v_infinity);      /* Final speed of wind in units of escape velocity */

  rddoub ("SV.acceleration_length(cm)", &zdom[ndom].sv_r_scale);        /*Acceleration length scale for wind */
  rddoub ("SV.acceleration_exponent", &zdom[ndom].sv_alpha);    /* Acceleration scale exponent */
  rddoub ("SV.gamma(streamline_skew;1=usually)", &zdom[ndom].sv_gamma); /* Parameter controlling how concentrated the streamlines
                                                                           are toward theta_min.  Large values concentrate toward
                                                                           the polls, must be greater than 0 */
  if (zdom[ndom].sv_gamma <= 0)
  {
    Error ("SV.gamma must be greater than 0. Large values skew streamlines to poles, small values to equator\n");
    exit (0);
  }

  /* allow the user to pick whether they set v0 by a fixed value or by the sound speed */
  strcpy (answer, "fixed");
  zdom[ndom].sv_v_zero_mode = rdchoice ("SV.v_zero_mode(fixed,sound_speed)", "0,1", answer);
  if (zdom[ndom].sv_v_zero_mode == FIXED)
    rddoub ("SV.v_zero(cm/s)", &zdom[ndom].sv_v_zero);
  else if (zdom[ndom].sv_v_zero_mode == SOUND_SPEED)
    rddoub ("SV.v_zero(multiple_of_sound_speed)", &zdom[ndom].sv_v_zero);

/* Assign the generic parameters for the wind the generic parameters of the wind */

  zdom[ndom].rmin = geo.rstar;
  zdom[ndom].wind_rhomin_at_disk = zdom[ndom].sv_rmin;
  zdom[ndom].wind_rhomax_at_disk = zdom[ndom].sv_rmax;
  zdom[ndom].wind_thetamin = zdom[ndom].sv_thetamin;
  zdom[ndom].wind_thetamax = zdom[ndom].sv_thetamax;


  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
  {
    zdom[ndom].xlog_scale = zdom[ndom].sv_rmin;
    zdom[ndom].zlog_scale = geo.rstar;
  }


  /*Now calculate the normalization factor for the wind */

  sdom = ndom;

  zdom[ndom].mdot_norm = num_int (sv_wind_mdot_integral, zdom[ndom].sv_rmin, zdom[ndom].sv_rmax, 1e-6);

  return (0);
}



/**********************************************************/
/** 
 * @brief      Calulate the v of a Schlossman Vitello wind at a position
 * 	x
 *
 * @param [in] double  x[]   the position for which one desires the velocity
 * @param [out] double  v[]  the calcualted velocity at that postion
 * @param [in] int  ndom the domain number 
 * @return     The amplitude of the velocity is returned
 * 	
 *
 * ###Notes###
 *
 * In addition to the standard SV law where the base velocity
 * is fixed, this supports a base velocity law which follows
 * the KWD prescription for the base velocity.
 *
 * For a vertically extended disk, we want the base of the
 * wind to have same base velocity as that for a flat disk
 * in the poloidal direction and proper rotational velocity
 * for that particular radial velocity.  
 *
 **********************************************************/

double
sv_velocity (x, v, ndom)
     double x[], v[];
     int ndom;
{
  double r, rzero, theta, speed;
  double ldist, zzz, v_escape, vl = 0.0;
  struct photon ptest;
  double xtest[3];
  double s;
  double vzero;
//OLD  double ldist_orig, rzero_orig;
  int hit_disk;

  zzz = v_escape = vzero = -99.;


//OLD  rzero_orig = rzero = sv_find_wind_rzero (ndom, x);
  rzero = sv_find_wind_rzero (ndom, x);
  theta = sv_theta_wind (ndom, rzero);

  r = sqrt (x[0] * x[0] + x[1] * x[1]);
//OLD  ldist_orig = ldist = sqrt ((r - rzero) * (r - rzero) + x[2] * x[2]);
  ldist = sqrt ((r - rzero) * (r - rzero) + x[2] * x[2]);

  /* Calculate the poloidal distance for a vertically extended disk ksl 111124 */
  if (geo.disk_type == DISK_VERTICALLY_EXTENDED && rzero < zdom[ndom].sv_rmax)
  {
    init_dummy_phot (&ptest);
    ptest.x[0] = rzero;         // Define ptest to be the footpoint extended to xy plane
    ptest.x[1] = 0.0;
    ptest.x[2] = EPSILON;
    ptest.lmn[0] = sin (theta); // ptest direction is along the stream line
    ptest.lmn[1] = 0.0;
    ptest.lmn[2] = cos (theta);
    s = ds_to_disk (&ptest, 1, &hit_disk);
    if (hit_disk)
    {
      move_phot (&ptest, s);    // Now move the test photon to  disk surface
      xtest[0] = r;             // Define xtest in the +z plane
      xtest[1] = 0;
      xtest[2] = fabs (x[2]);
      vsub (ptest.x, xtest, xtest);     // Poloidal distance is just the distance between these two points.
      ldist = length (xtest);
      rzero = sqrt (ptest.x[0] * ptest.x[0] + ptest.x[1] * ptest.x[1]);
    }
  }

  /* Depending on the mode, we can set the initial velocity as fixed or by the sound speed. See #482. */
  if (zdom[ndom].sv_v_zero_mode == FIXED)
    vzero = zdom[ndom].sv_v_zero;
  else
    vzero = kn_vzero (rzero) * zdom[ndom].sv_v_zero;

  if (ldist > 0)
  {
    zzz = pow (ldist / zdom[ndom].sv_r_scale, zdom[ndom].sv_alpha);

    if (rzero < geo.rstar)
      v_escape = sqrt (2. * GRAV * geo.mstar / geo.rstar);
    else
      v_escape = sqrt (2. * GRAV * geo.mstar / rzero);

    vl = vzero + (zdom[ndom].sv_v_infinity * v_escape - vzero) * zzz / (1. + zzz);
  }

  v[0] = vl * sin (theta);

  if (r > 0)
    v[1] = sqrt (GRAV * geo.mstar * rzero) / r;
  else
    v[1] = 0;

  v[2] = vl * cos (theta);

  if (x[2] < 0)
    v[2] *= (-1);

  /* At this point we have calculated the velocity in the xz plane, which
   * is identical to the statement that we have calculated it in
   * cylindrical coordinates.  Now project the velocity back to  
   * coordinates if x was not originally in the xz plane.
   */

  if (x[1] != 0.0)
  {
    project_from_cyl_xyz (x, v, xtest);
    stuff_v (xtest, v);
  }


  speed = (sqrt (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
  if (sane_check (speed))
  {
    Error ("sv_velocity: x %g %g %g v %f %f %f\n", x[0], x[1], x[2], v[0], v[1], v[2]);
    Error ("sv_velocity: rzero %g theta %f ldist %f zzz %f v_escape %f vl %f\n", rzero, theta, ldist, zzz, v_escape, vl);
  }


  return (speed);

}



/**********************************************************/
/** 
 * @brief      Calculate the density of an sv_wind at a position x
 *
 * @param [in] int  ndom   The domain number
 * @param [in] double  x[] The position where for the which one desires the denisty
 * @return     The density at x in gram/cm**3
 *
 *
 * ###Notes###
 *
 * 52 -- If the disk is thick, our intention is that the stream lines
 * 	trace back to the disk, and have the same properties at that radius
 * 	as they would have had except for a vertical offset.
 *
 **********************************************************/

double
sv_rho (ndom, x)
     double x[];
     int ndom;
{
  double r, rzero, theta;
  double ldist;
  double dmdot_da;
  double dtheta_drzero, dr_drzero;

  double v[3], rho;
  struct photon ptest;
  double xtest[3];
  double s;
  int hit_disk;

  sv_velocity (x, v, ndom);

  rzero = sv_find_wind_rzero (ndom, x);
  theta = sv_theta_wind (ndom, rzero);

  r = sqrt (x[0] * x[0] + x[1] * x[1]);
  ldist = sqrt ((r - rzero) * (r - rzero) + x[2] * x[2]);



  /* Calculate the poloidal distance for a vertically extended disk */
  if (geo.disk_type == DISK_VERTICALLY_EXTENDED && rzero < zdom[ndom].sv_rmax)
  {
    init_dummy_phot (&ptest);
    ptest.x[0] = rzero;         // Define ptest to be the footpoint extended to xy plane
    ptest.x[1] = 0.0;
    ptest.x[2] = EPSILON;
    ptest.lmn[0] = sin (theta); // ptest direction is along the stream line
    ptest.lmn[1] = 0.0;
    ptest.lmn[2] = cos (theta);
    s = ds_to_disk (&ptest, 1, &hit_disk);
    if (hit_disk)
    {
      move_phot (&ptest, s);    // Now move the test photon to  disk surface
      xtest[0] = r;             // Define xtest in the +z plane
      xtest[1] = 0;
      xtest[2] = fabs (x[2]);
      vsub (ptest.x, xtest, xtest);     // Poloidal distance is just the distance between these two points.
      ldist = length (xtest);
      rzero = sqrt (ptest.x[0] * ptest.x[0] + ptest.x[1] * ptest.x[1]);
    }
  }


/* Reduced by a factor of 2 to account for radiation on both sides of the disk */
  dmdot_da = zdom[ndom].wind_mdot * pow (rzero, zdom[ndom].sv_lambda) * cos (theta) / zdom[ndom].mdot_norm / 2.;

/* Although the current definition of sv_theta_wind is continuous, the derivative is not continuous accross the
   outer boundary of the wind and thus dtheta_drzero would appear to change discontinuously.   This created
   large apparent density jumps at the outside edge of the wind. We can't allow that and so we force the derivative to equal
   that at the edge of the wind if you try to calculate the density rho.  ksl 97apr23 */

  if (rzero > zdom[ndom].sv_rmax)
    rzero = zdom[ndom].sv_rmax;
  dtheta_drzero = (sv_theta_wind (ndom, rzero) - sv_theta_wind (ndom, (1. - EPSILON) * rzero)) / (EPSILON * rzero);

  dr_drzero = 1. + ldist * dtheta_drzero / cos (theta);
  rho = rzero * dmdot_da / (dr_drzero * r * v[2]);

  return (rho);
}





/**********************************************************/
/** 
 * @brief      Locates the radial position in the disk from which the 
 * 	stream line arises.
 *
 * @param [in] int  ndom   The domain number
 * @param [in] double  p[]   A position in cartesian coordinates
 *
 * @return     the radius of the footpoint of the stream line
 *
 * ###Notes###
 *
 * In our models, which follows Shlossman and Vitello the wind is expected to arise in a cone.  
 * 	However, we have also defined r_o in a sensible way for any position regardless of whether it
 * 	is in the cone or not.   The routine should be accurate even if p[2] is negative, i.e if we
 * 	are trying to calculate r in the -z hemisphere.	
 * 
 * 	More specifically, if the point is a distance drho outside of the windcone then rzero
 * 	will be rmax+drho. Alternatively if the point is inside the wind cone then, rzero
 * 	will rmin*rho/rho_min where rho_min is the minimum value of rho to be in the windcone
 * 	at that height.
 *
 *
 *
 **********************************************************/

double
sv_find_wind_rzero (ndom, p)
     int ndom;
     double p[];                /* Note that p is a 3 vector and not a photon structure */
{
  double x, z;
  double rho_min, rho_max, rho;
  int ierr;




  /* thetamin and theta max are defined w.r.t  z axis */

  z = fabs (p[2]);              /* This is necessary to get correct answer above
                                   and below plane */

  if (z == 0)
  {
    x = (sqrt (p[0] * p[0] + p[1] * p[1]));     // If p is in the xy plane, there is no need to think further
    return (x);
  }


  sv_zero_init (p);             /* This initializes the routine sv_zero_r.  It is not
                                   actually needed unless zbrent is called, but it
                                   does allow you to check your answer otherwize
                                 */
  /* The next lines provide a graceful answer in the case where the
   * position is actually outside the wind so that rzero returned is
   * continuous
   */

  rho_min = zdom[ndom].sv_rmin + z * tan (zdom[ndom].sv_thetamin);
  rho_max = zdom[ndom].sv_rmax + z * tan (zdom[ndom].sv_thetamax);
  rho = sqrt (p[0] * p[0] + p[1] * p[1]);

  if (rho <= rho_min)
  {
    x = zdom[ndom].sv_rmin * rho / rho_min;
    return (x);
  }
  if (rho >= rho_max)
  {
    x = zdom[ndom].sv_rmax + rho - rho_max;
    return (x);
  }

  /* 100 here means that zbrent will end if one has a guess of rzero which is
     correct at 100 cm */

  /* change the global variable sv_zero_r_ndom before we call zbrent */
  sv_zero_r_ndom = ndom;
//  x = zbrent (sv_zero_r, zdom[ndom].sv_rmin, zdom[ndom].sv_rmax, 100.);
  x = zero_find (sv_zero_r, zdom[ndom].sv_rmin, zdom[ndom].sv_rmax, 100., &ierr);

  if (ierr)
  {
    Error ("sv_find_rzero: zero_find error at rho %.3e and z %.3e\n", rho, z);

  }


  return (x);

}


double zero_p[3];


/**********************************************************/
/** 
 * @brief      A setup routine for finding the footpoint of a stream line in a SV model
 *
 *
 * @param [in] double  p[]   A position in cartesian coordinates
 * @return     The routine simply returns 0
 *
 * One of two routines used to find the position on the disk from which a steamline
 * arises. This routine is just a setup routine, to make the position an external
 * variable so it can be accessed by sv_zero_r 
 *
 * ###Notes###
 *
 * The routine stores the position in zero_p
 *
 **********************************************************/

int
sv_zero_init (p)
     double p[];
{
  stuff_v (p, zero_p);
  zero_p[2] = fabs (zero_p[2]); /* Required to get correct 
                                   answer below (in -z ) the disk */
  return (0);
}




/**********************************************************/
/** 
 * @brief      Routine used by zbrent to find a footpoint of a streamline in an SV model
 *
 * @param [in] double  r   A position along the surface of the disk
 * @return     0 if r is the footpoint
 *
 * This routine is used to test whether a guess of r_zero is correct.  If
 * you have the answer exactly then sv_zero_r will return 0
 *
 * sv_zero_r is the routine called by
 * the numerical recipes routine zbrent to walk down on the actual value
 * of r in the disk.   
 *
 * sv_zero_r returns the difference
 * between rho actual (sqrt x*x + y *y ) and rho_guessed
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

double
sv_zero_r (double r, void *params)
{
  double theta;
  double rho, rho_guess;

  /* sv_zero_r_ndom should be set to the domain number before calling sv_zero_r */

  theta = sv_theta_wind (sv_zero_r_ndom, r);

  rho = sqrt (zero_p[0] * zero_p[0] + zero_p[1] * zero_p[1]);
  rho_guess = r + tan (theta) * zero_p[2];
  return (rho_guess - rho);

}



/**********************************************************/
/** 
 * @brief      find the angle at which the wind emerges from at a specific
 * 	radius on surface of the disk
 *
 * @param [in] int  ndom   The domain number
 * @param [in] double  r   a radial distance on the surface of the disk
 * @return   An angle (in radians) 
 *
 * As long as r is between sv_rmin and sv_rmax, sv_theta_wind returns the
 * angle given by the SV prescription.
 * 	
 * Inside sv_rmin, it returns a value which smoothly goes from thetamin to 0
 * as the radius goes to 0.
 * 	
 * Outside sv_rmax, it returns sv_thetamax
 *
 *
 * ###Notes###
 *
 * Theta is measuered from the z axis.
 *
 **********************************************************/

double
sv_theta_wind (ndom, r)
     int ndom;
     double r;
{
  double theta;


  if (r <= zdom[ndom].sv_rmin)
    return (atan (tan (zdom[ndom].sv_thetamin * r / zdom[ndom].sv_rmin)));
  if (r >= zdom[ndom].sv_rmax)
    return (zdom[ndom].sv_thetamax);

  theta = zdom[ndom].sv_thetamin +
    (zdom[ndom].sv_thetamax - zdom[ndom].sv_thetamin) *
    pow ((r - zdom[ndom].sv_rmin) / (zdom[ndom].sv_rmax - zdom[ndom].sv_rmin), zdom[ndom].sv_gamma);

  return (theta);

}




/**********************************************************/
/** 
 * @brief      The integrand required to calculate the normalization 
 *             factor between the global mass loss rate and the mass 
 *              loss per unit area at a particular place on the disk
 *
 * @param [in] double  r        A position (radius) on the disk
 * @param [in] void    params   An extra (unused) variable to make it paletable for the gsl integrator
 * @return     The value of the integrand 
 *
 * The SV model is defined in terms of a mass loss rate per unit area.  
 * mdot is the total mass loss rate from the disk.  In order to connect 
 * the local and global rates one
 * must integrate the local mass loss rate and renormalize this.  
 *
 * ###Notes###
 *
 * The integration is carried out in get get_sv_wind_params
 *
 **********************************************************/

double
sv_wind_mdot_integral (double r, void *params)
{
  double x;

  x = 2 * PI * pow (r, zdom[sdom].sv_lambda + 1.) * cos (sv_theta_wind (sdom, r));
  return (x);

}
