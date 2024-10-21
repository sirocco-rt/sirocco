
/***********************************************************/
/** @file  roche.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief      These routines are used to describe the Roche geometry and are intended to allow one
 *    to calculate whether a photon intersects the surface of the star
 *
 *
 * The important routines are
 *
 *    - binary_basics  which uses the masses of the primary and the secondary to calculate
 *    the position of various Lagrange points, etc, and the size of a pillbox
 *    which encloses the Roche lobe of the secondary.  It must be called before
 *    the other routines can be used!
 *
 *    - hit_secondary   determines whether a photon defined by p_roche will hit the secondary
 *    or not.  If it misses, hit_secondary returns 0, if it hits, the return is
 *    P_SEC (or hit secondary as contained in sirocco.h).
 *
 * ###Notes###
 *
 * The secondary is located along the x, i.e. 0, axis
 *
 *    To utilize the routines you must stuff the photon of interest into the external photon "p_roche".
 *    This is necessary because the program makes use of the gsl-based routines find_rzero    
 *    and find_function_miniumum to find zeroes and minima respectively.  
 *
 *    The routines adopt the same approach as in Keith Horne's routines dealing with the Roche
 *    geometry, that is a pillbox is defined around the secondary.  The radius of the pillbox is
 *    determined by the maximum size of the secondary perpendicular to the line of sight between
 *    the two stars.  One end of the pillbox is defined by the L1 point.  The other end of the pillbox
 *    is defined by the backside of the secondary.  One first tests to find out whether a photon
 *    intersects the pillbox, and then one checks more carefully to see whether the photon actually
 *    hits the secondary (using the intersections with the pillbox to define the limits of the careful
 *    search.
 *

 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

/* Given the mass of the primary and secondary and the period calculate some other basic parameters
 for the binary system, including the L1 point. */


/**********************************************************/
/** p_roche provides a way to pass a photon structure to roche() and roche_derive which
 * are in turn used  to find the zero either to dphi_ds
 * or to phi along the path lenght of the photon. p_roche is intialized in binary_basics
 **********************************************************/
struct photon p_roche;


/**********************************************************/
/**
 * @brief      From the masses and period, setup the information one needs to calculate the position of the Roche lobe of the secondary
 *    to calculate whether a photon intersects the surface of the star
 *
 * @return     Always returns 0
 *
 *
 * Notes
 * Initailizes p_roche
 *
 *
 **********************************************************/

int
binary_basics ()
{
  double x;
  void *dummy_par = NULL;
  int ierr;

  /* Calculate a the separation of the primary and the secondary */

  x = GRAV * geo.period * geo.period * (geo.mstar + geo.m_sec) / (16. * PI * PI);
  geo.a = pow (x, 1. / 3.);

  /*Calculate the mass ratio.  */

  geo.q = geo.m_sec / geo.mstar;

  /* Find the position of the L1 point with respect to the primary */

  geo.l1 = 0;

  /* p_roche is used here to provide an external way to pass a photon structure to find_l1 */
  /* p_roche is (near) origin of primary and pointing at secondary */

  p_roche.x[0] = 1000;
  p_roche.x[1] = p_roche.x[2] = 0;
  p_roche.lmn[0] = 1;
  p_roche.lmn[1] = p_roche.lmn[2] = 0;
  p_roche.frame = F_OBSERVER;

  geo.l1 = x = zero_find (dphi_ds, 0.01 * geo.a, 0.99 * geo.a, geo.a / 1000., &ierr);

  if (ierr)
  {
    Error ("binary_basics: Did not find L1\n");
  }

  geo.l1_from_m2 = geo.a - geo.l1;

  plane_m2_near.x[0] = geo.l1 + 10000.;
  plane_m2_near.x[1] = plane_m2_near.x[2] = 0;
  plane_m2_near.lmn[0] = 1;
  plane_m2_near.lmn[1] = plane_m2_near.lmn[2] = 0;

  /* Similarly, find l2, the Lagrangian point behind the secondary */

  geo.l2 = x = zero_find (dphi_ds, 1.01 * geo.a, 10.0 * geo.a, geo.a / 1000., &ierr);

  if (ierr)
  {
    Error ("binary_basics: Did not find L2\n");
  }


  geo.phi = 0;
  geo.phi = phi (geo.l1, dummy_par);    /* Set geo.phi so phi will be zero on Roche lobes */


  /* Set geo.r2_far to be the radius of the secondary on the backside of the secondary */

  geo.r2_far = x = zero_find (phi, 1.01 * geo.a, geo.l2, geo.a / 1000., &ierr) - geo.a;


  if (ierr)
  {
    Error ("binary_basics: Did not find r2_far\n");
  }


  /* Define a plane on the backside of the secondary */

  plane_m2_far.x[0] = geo.r2_far + geo.a;
  plane_m2_far.x[1] = plane_m2_far.x[2] = 0;
  plane_m2_far.lmn[0] = 1.0;
  plane_m2_far.lmn[1] = plane_m2_far.lmn[2] = 0;


  Log ("binary_basics: a=%8.2e l1=%8.2e l2=%8.2e l1_from_m2=%8.2e r2_far %8.2e\n", geo.a, geo.l1, geo.l2, geo.l1_from_m2, geo.r2_far);
  Log ("binary_basics: m2_near=%8.2e m2_far=%8.2e \n", plane_m2_near.x[0], plane_m2_far.x[0]);

  /* Calculate the maximum half width of the secondary in the plane of the orbit */

  geo.r2_width = roche2_half_width ();
  Log ("binary_basics: r2_width=%6.2e\n", geo.r2_width);
  return (0);
}






/**********************************************************/
/**
 * @brief      Find out whether a photon will hit the secondary
 *
 * @param [in] PhotPtr  p   A photon, for which we have a possiton and adirection of travel
 * @return     0 if the photon would miss the secondary, P_SEC otherwise
 *
 * The routine chacks whether a photon hits the Roche lobe of the secondary star.
 *
 * ###Notes###
 *
 * The routine first checks if the photon encounters a pill_box that surrouds the secondary.  Photons
 * that do not encounter this pill box do not hit the secondary.  If  the photon enters the pill box, then
 * the routine determines whether the photon hits the roche surface of the secondary, by finding the minimum
 * value of the Rooche potential. If this is greater than 0, then the photon has missed the secondary. If less
 * thn zero then it has hit the secondary.
 *
 * The routine does not in most cases calculate where the photon hits the secondary, only if it hits
 * the secondary.
 *
 **********************************************************/

int
hit_secondary (p)
     PhotPtr p;
{
  double smin, smax, smid, s;
  double potential;
  double idelt = 1e-2;
  void *dummy_par = NULL;       //A variable required (but not set) for calls to phi

  /* Check to see if the photon encountered the pill box */
  if (pillbox (p, &smin, &smax) == VERY_BIG)
    return (0);
  stuff_phot (p, &p_roche);


  smid = 0.5 * (smin + smax);
  if (phi (smid, dummy_par) > phi (smin, dummy_par) || phi (smid, dummy_par) > phi (smax, dummy_par))   //value of the function at the midpoint is greater than the end(s) problem
  {
    if (phi (smin + idelt * (smax - smin), dummy_par) > phi (smin, dummy_par) && phi (smax - idelt * (smax - smin), dummy_par) < phi (smax, dummy_par)) //monotonically rising? - no minimum
    {
      potential = fmin (phi (smin, dummy_par), phi (smax, dummy_par));
    }

    else if (phi (smin + idelt * (smax - smin), dummy_par) < phi (smin, dummy_par) && phi (smax - idelt * (smax - smin), dummy_par) > phi (smax, dummy_par))    //monotonically rising? - no minimum
    {
      potential = fmin (phi (smin, dummy_par), phi (smax, dummy_par));
    }
    else if (phi (smin + idelt * (smax - smin), dummy_par) > phi (smin, dummy_par) && phi (smax - idelt * (smax - smin), dummy_par) > phi (smax, dummy_par))    //hump - probably no minimum
    {
      potential = fmin (phi (smin, dummy_par), phi (smax, dummy_par));
    }
    else if (phi (smin + idelt * (smax - smin), dummy_par) < phi (smin, dummy_par) && phi (smax - idelt * (smax - smin), dummy_par) < phi (smax, dummy_par))    //there most likely is a minimum - need to find a sensibly midpoint 
    {
      if (phi (smin, dummy_par) < phi (smax, dummy_par))
      {
        smid = smin + idelt * (smax - smin);    //the function at this point, just inside smin should be lower than the function at smin
      }
      else
      {
        smid = smax - idelt * (smax - smin);    //the function at this point, just inside smax should be lower than the function at smax
      }
      potential = find_function_minimum (smin, smid, smax, phi, 0.0001, &s);    //call the minimuser            
    }
    else
    {
      potential = 1.;           //goodness only knows, lets say the photon missed
    }
  }
  else
  {
    potential = find_function_minimum (smin, smid, smax, phi, 0.0001, &s);
  }


  if (potential > 0.0)
    return (0);


  return (P_SEC);
}




/**********************************************************/
/**
 * @brief       Find out whether the photon hits the pillbox which surrounds the Roche surface of the secondary.
 *
 * @param [in] PhotPtr  p   The photon whose trajectory we wish to examine
 * @param [out] double *  smin   If the photon, hits the pillbox, this is how far the photon must travel to enter the pillbox
 * @param [out] double *  smax   If the photon, hits the pillbox, this is how far the photon must travel to exit the pillbox
 * @return
 * pillbox returns VERY_BIG if the photon p has not intersected the pillbox.  It return the distance
 * to near edge the pillbox if the photon p hits the pillbox, or 0 if the photon is inside the pillbox.  
 *
 * The purpose of this routine is to separate photons which may hit the Roche surface and those which certainly do not
 * by seeing if the photon avoids a pillbox surrounding the secondary.  If the photons does not enter the pillbox then
 * one dos not have to explore in detail whether the photon hits the Roche surrvace.
 *
 * If the photon does pass through the secondary, then the routine returns the distances of entry and exit which
 * are used to limit the range of root finding that then must take place.
 *
 * ###Notes###
 *
 * The pillbox is defined by a plane at L1, a plane on the backside of the secondary star and a cylinder whose radius
 * is the maximum (half) width of the star in the plane of the orbit.
 *
 *
 **********************************************************/

double
pillbox (p, smin, smax)
     PhotPtr p;
     double *smin, *smax;
{
  double x1, x2;
  double a, b, c;
  double root[2], ss[2];
  int i, n;
  struct photon pp;

  n = 0;

/* If the photon is along the x-axis, then it will not intersect the cylinder anywhere, and
the only possibility is that it hit the endcaps of the cylinder. But normally, it will
not be along the x axis and one must find if there are any intersections with the cylinder
and determine where they are.  This calculation is carried out in the if section below. */

  if ((a = (p->lmn[0] * p->lmn[0])) < 1.0)
  {
    a = 1. - a;
    b = 2. * (p->x[1] * p->lmn[1] + p->x[2] * p->lmn[2]);
    c = p->x[1] * p->x[1] + p->x[2] * p->x[2] - geo.r2_width * geo.r2_width;

/* After "quadratic", root contains the distances the photon pp would need to travel to
hit the edges of the cylinder. i is an index to the smallest positive root, if one exists.
If both of the roots are negative or imaginary then i will be negative, and we know
that the ray did not hit the cylinder while travelling in the positive direction.  In
that case return  VERY_BIG

Note that the fact that the ray does hit the cylinder going in the positive direction
does not necessarily mean that it hits the pillbox.
*/

    if ((i = quadratic (a, b, c, root)) < 0)
      return (VERY_BIG);

/* If the intersection is in the part of the cylinder between the caps it is probably legitimate.
So tranport pp to the intersection of the first root and check where that lies. And then
do the same thing for the second root.  The "probably refers to the fact that a negative
root is really only possible if the photon is already in the pillbox  */

    stuff_phot (p, &pp);
    move_phot (&pp, root[0]);

    if (geo.l1 <= pp.x[0] && pp.x[0] <= plane_m2_far.x[0])
    {
      ss[n] = root[0];
      n++;
    }


    stuff_phot (p, &pp);
    move_phot (&pp, root[1]);   /* So pp is now located at intersection of the second root */
    if (geo.l1 <= pp.x[0] && pp.x[0] <= plane_m2_far.x[0])
    {
      ss[n] = root[1];
      n++;
    }

  }

/* Now find out if the photon hits the caps. This is done by translating to the
planes of the two end caps and then checking the cylindrical radius. Negative
distances are valid only if the photon is in the pillbox already */

  x1 = ds_to_plane (&plane_m2_near, p, FALSE);  /* Calculate the distance to the L1 plane */
  stuff_phot (p, &pp);
  move_phot (&pp, x1);          /* So pp is now located at the l1 plane */
  if ((pp.x[1] * pp.x[1] + pp.x[2] * pp.x[2]) <= geo.r2_width * geo.r2_width)
  {
    ss[n] = x1;
    n++;
  }

  x2 = ds_to_plane (&plane_m2_far, p, FALSE);   /* Calculate the distance to the plane behind the the star.  */
  stuff_phot (p, &pp);
  move_phot (&pp, x2);          /* So pp is now located at the r2_far plane */
  if ((pp.x[1] * pp.x[1] + pp.x[2] * pp.x[2]) <= geo.r2_width * geo.r2_width)
  {
    ss[n] = x2;
    n++;
  }

/* Finally sort it all out.  If n == 0, or if both s[0] and s[2] are negative
then the photon did not hit the pillbox ksl 02jan */

  if ((n == 0) || (ss[0] < 0 && ss[1] < 0))
    return (VERY_BIG);

  if (n == 2)
  {
    if (ss[0] < ss[1])
    {
      *smin = ss[0];
      *smax = ss[1];
    }
    else
    {
      *smin = ss[1];
      *smax = ss[0];
    }

    /* If the photon is inside the pillbox, we set the minimum distance to 0,
       we only want to consider travel in the forward direction.
     */

    if (*smin < 0)
    {
      *smin = 0;
    }

    return (*smin);
  }

  Error ("pillbox %d interfaces to pillbox is impossible\n", n);
  return (VERY_BIG);

}


int phi_init = 0;
double phi_gm1, phi_gm2, phi_3, phi_4;


/**********************************************************/
/**
 * @brief      Calculate the Roche potential at position indicated by the photon p_roche and the length s
 *
 * @param [in] double  s   The distance from the current position of the photon to the point where one wishes to calculate the Roche potential
 *              void * params  An unused variable required to make the function compatible with the gsl routine used to minimise it
 * @return     The Roche potential at the position given by the photon moved by a distance s
 *
 * Given a photon stored in p_roche and a distance s to move the photone, phi returns the roche potentail at that point
 *
 * ###Notes###
 * This routine is used in a number of locations for calculating either where the
 * Roche surface is, or whether a photon has or will cross the Roche
 * lobe of the secondary star
 * 
 * The external variable p_roche gives the initial position and direction
 * of a photon, and s the distance along the path length at 
 * which the potential is calculated.
 *
 *
 **********************************************************/

double
phi (double s, void *params)
{
  struct photon pp;
  double x1, x2, z, z1, z2, z3;
  double length ();

  if (phi_init == 0)
  {
    phi_gm1 = GRAV * geo.mstar;
    phi_gm2 = GRAV * geo.m_sec;
    phi_3 = (phi_gm1 + phi_gm2) / (geo.a * geo.a * geo.a);
    phi_3 = 0.5 * phi_3;
    phi_4 = geo.a * geo.m_sec / (geo.mstar + geo.m_sec);
    phi_init++;
  }

  stuff_phot (&p_roche, &pp);
  move_phot (&pp, s);

  if ((x1 = length (&pp.x[0])) == 0)
    return (-VERY_BIG);

  z1 = -phi_gm1 / x1;
  z3 = -phi_3 * ((pp.x[0] - phi_4) * (pp.x[0] - phi_4) + pp.x[1] * pp.x[1]);

  pp.x[0] -= geo.a;             /* Here we make pp refer to the positions w.r.t. the secondary */
  if ((x2 = length (&pp.x[0])) == 0)
    return (-VERY_BIG);

  z2 = -phi_gm2 / x2;

  z = z1 + z2 + z3 - geo.phi;

  //printf ("foo %e %e %e %e \n", z1, z2, z3, geo.phi);

  if (z < -1.e20)
    z = -1.e20;

  return (z);

}

#define EPS 10000.


/**********************************************************/
/**
 * @brief      Calculate the gradient of the Roche potential as the photon travels
 *
 * @param [in] double  s   The distance the photon has traveled from its initial position
 *             void * params  An unused variable required to make the function compatible with the gsl routine used to minimise it

 * @return     the first derivative of the potential
 *
 * This is a brute force calculation of the first derivative of the Roche potential
 *
 * ###Notes###
 *
 * XXXX This seems awfully coarse and it is not a two-sided derivative
 *
 *
 **********************************************************/

double
dphi_ds (double s, void *params)
{
  double phi (), x1, x2;
  void *dummy_par = NULL;
  double dx, z;
  if ((dx = 0.001 * geo.a) < EPS)
    dx = EPS;

  x2 = phi (s + dx, dummy_par);
  x1 = phi (s, dummy_par);
  z = (x2 - x1) / dx;
  return (z);

}


/**********************************************************/
/**
 * @brief      Calculate the half width of the Roche potential at a position x along the x axis
 *
 * @param [in] double  x   A position (distance from the primary) along the positive x axis                            
 *             void params  unused variable required to present the correct function to gsl
 *
 * @return     the half width (or actually the negative of the half width)
 *
 * Calculate the half width of the Roche lobe at the position x along the
 *   x axis in the plane of the orbit.
 *
 * ###Notes###
 * There is no real guarantee that this would
 * work if you were outside L2 or L3.
 *
 * The routine return -rho because it is used as part of a search for
 * the maximum value of the half width, and that search is based on
 * finding a minimum and not a maximum value.
 *
 **********************************************************/


double
roche_width (double x, void *params)
{
  double rho, smax;
  int ierr;
  void *dummy_par = NULL;

  p_roche.x[1] = p_roche.x[2] = 0.0;
  p_roche.lmn[0] = 0;
  p_roche.lmn[1] = 1;
  p_roche.lmn[2] = 0;

//OLD  geo.phi = 0;
//OLD  p_roche.x[0] = geo.l1;
//OLD  geo.phi = phi (0, dummy_par);


  if (x < geo.l1)
    smax = geo.l1;
  else
    smax = geo.l1_from_m2;

  p_roche.x[0] = x;

  rho = zero_find (phi, 1000., smax, geo.a / 1000., &ierr);
  if (rho < 0 || ierr)
  {
    p_roche.x[0] = geo.l1;
    Error ("roche_width : zero_find failure x=%6.2e rho %.2e phi(l2) %.2e\n", x, rho, phi (0, dummy_par));
  }
  return (-rho);
}


/**********************************************************/
/**
 * @brief      Find the half width of the secondary Roche lobe
 *
 * @return     The half width is returned
 *
 * Find the maximum half width of the Roche lobe of the secondary.   
 *
 * ###Notes###
 *
 *
 *
 **********************************************************/

double
roche2_half_width ()
{
  double xmin, xmax, xmid, xbest;
  double rmin;
  double roche_width ();

  xmin = geo.l1 + 1.e5;
  xmax = geo.r2_far + geo.a - 1.e5;
  xmid = geo.a;

  Log_silent ("roche2_half_width: Search from %6.2e %6.2e %6.2e\n", xmin, xmid, xmax);


  rmin = find_function_minimum (xmin, xmid, xmax, roche_width, 0.0001, &xbest);

  Log_silent ("roche2_half_width: Max half width of %6.2e is located at %6.2e\n", -rmin, xbest);

  return (-rmin);
}
