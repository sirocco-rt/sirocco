
/***********************************************************/
/** @file  phot_util.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief
 * This file contains utilities which can be used to work with photon structues.  It also
 * contains routines to calculate the distances to intersections of spheres and cones.  The
 * routines included here should not be tightly tied to a particular wind model.   The one
 * exception at present is ds_to_cone which calculates the distance to a bifurcated cone.
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "atomic.h"
#include "sirocco.h"

int init_stuff_phot = 0;
size_t sizeofphot;

/**********************************************************/
/**
 * @brief      Simply initialize some of the variables in a 
 *             single photon used for some special purpose 
 *             calculation
 *
 *
 * @param [in] PhotPtr  p   The photon bundle to be intialized
 * @return     Always returns 0
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/


int
init_dummy_phot (p)
     PhotPtr p;
{
  p->x[0] = p->x[1] = p->x[2] = 0.0;
  p->lmn[0] = p->lmn[1] = 0;
  p->lmn[2] = 1;
  p->freq = p->freq_orig = 0.0;
  p->frame = F_OBSERVER;
  p->origin = p->origin_orig = PTYPE_DUMMY;
  p->np = -1;
  p->tau = p->ds = 0;

  return (0);


}


/**********************************************************/
/**
 * @brief      Simply copy one photon bundle into another
 *
 * @param [in] PhotPtr  pin   The photon bundle to be copied
 * @param [out] PhotPtr  pout   The PhotPtr that contains the output
 * @return     Always returns 0
 *
 * @details
 * This just a simple utility to copy a photon bundle
 *
 * ### Notes ###
 * All variables are (or should be) copied.
 *
 **********************************************************/

int
stuff_phot (pin, pout)
     PhotPtr pin, pout;
{
  pout->x[0] = pin->x[0];
  pout->x[1] = pin->x[1];
  pout->x[2] = pin->x[2];

  pout->lmn[0] = pin->lmn[0];
  pout->lmn[1] = pin->lmn[1];
  pout->lmn[2] = pin->lmn[2];

  pout->w = pin->w;
  pout->w_orig = pin->w_orig;
  pout->freq = pin->freq;
  pout->freq_orig = pin->freq_orig;
  pout->tau = pin->tau;

  pout->istat = pin->istat;
  pout->frame = pin->frame;
  pout->nres = pin->nres;
  pout->line_res = pin->line_res;
  pout->nmacro = pin->nmacro;
  pout->nrscat = pin->nrscat;
  pout->nscat = pin->nscat;
  pout->nnscat = pin->nnscat;
  pout->grid = pin->grid;
  pout->origin = pin->origin;
  pout->origin_orig = pin->origin_orig;
  pout->nnscat = pin->nnscat;
  pout->np = pin->np;

  pout->path = pin->path;

  pout->ds = pin->ds;

  return (0);
}



/**********************************************************/
/**
 * @brief      move a photon by a distance ds
 *
 * @param [in, out] PhotPtr  pp   A photon bundle
 * @param [in] double  ds   The distance to move the photon
 * @return     Always returns 0
 *
 * @details
 * The routine simply updates the position x of the
 * photon using the current direction and the distance.
 *
 * The path length for the photon is also updated
 * for use in reverberation calculations
 *
 * ### Notes ###
 *
 **********************************************************/

int
move_phot (pp, ds)
     PhotPtr pp;
     double ds;
{
  int ierr;

  ierr = check_frame (pp, F_OBSERVER, "move_phot");

  pp->x[0] += pp->lmn[0] * ds;
  pp->x[1] += pp->lmn[1] * ds;
  pp->x[2] += pp->lmn[2] * ds;

  pp->ds += ds;
  pp->path += fabs (ds);
  return (ierr);
}



/**********************************************************/
/**
 * @brief      Compares two photons to determine whether they have
 *  	the identical position and direction.
 *
 * @param [in] PhotPtr  p1   The first photon bundle
 * @param [in] PhotPtr  p2   The second phton bundle
 * @return     Returns 0 if the position and direction are identical,
 * 1 otherwise
 *
 * @details
 *
 * ### Notes ###
 * The routine is used in calculate_ds
 *
 **********************************************************/

int
comp_phot (p1, p2)
     PhotPtr p1, p2;
{
  if (p1->x[0] != p2->x[0] || p1->lmn[0] != p2->lmn[0])
    return (1);
  if (p1->x[1] != p2->x[1] || p1->lmn[1] != p2->lmn[1])
    return (1);
  if (p1->x[2] != p2->x[2] || p1->lmn[2] != p2->lmn[2])
    return (1);
  return (0);
}





/**********************************************************/
/**
 * @brief      A routine to calculate the distance a photon must travel
 *   to hit a cone.
 *
 * @param [in] ConePtr  cc   A structure which describes a cone
 * @param [in] struct photon *  p   A photon
 * @return     The distance to the cone
 *
 * @details
 *
 * ### Notes ###
 *
 *
 * The cone is defined by  z=z_o + dz/drho *drho
 * The photon is defined by X=X_o + LMN s
 *
 * Unless dz/drho = 0, this turns into a simple quadratic equation that must
 * be solved.  We work in the northen hemisphere only.
 *
 * There are issues concerning what to do so one crosses the disk plane.  For
 * finding how far a photon can travel in a cell, we are *not* interested
 * in the possibility that the photon hits the cone on the other side of
 * the disk plane, but for a photon travelling out of the wind we are
 * very interested in this possibility.
 *
 * This version is intended to allow cylvar coordinates, in old versions (56d)
 * prior to cylvar coordiantes, this routine was called ds_to_windcone
 **********************************************************/

double
ds_to_cone (cc, p)
     ConePtr cc;
     struct photon *p;
{
  double dz, dzdr2;
  double a, b, c, root[2];
  double s_to_zero;             /* The path length to the xy plane */
  int i;
  struct photon pp;

  /* First of all let's work only in the "northern" hemisphere */
  stuff_phot (p, &pp);

  if (pp.x[2] < 0.0)
  {                             /*move the photon to the northen hemisphere */
    pp.x[2] = -pp.x[2];
    pp.lmn[2] = -pp.lmn[2];
  }

  /* Set up and solve the quadratic equation that gives the cone intercept */

  dzdr2 = cc->dzdr * cc->dzdr;
  dz = pp.x[2] - cc->z;

  a = dzdr2 * (pp.lmn[0] * pp.lmn[0] + pp.lmn[1] * pp.lmn[1]) - (pp.lmn[2] * pp.lmn[2]);
  b = 2. * (dzdr2 * (pp.lmn[0] * pp.x[0] + pp.lmn[1] * pp.x[1]) - pp.lmn[2] * dz);
  c = dzdr2 * (pp.x[0] * pp.x[0] + pp.x[1] * pp.x[1]) - dz * dz;

  i = quadratic (a, b, c, root);        /* root[i] is the smallest positive root unless i is
                                           negative in which case either both roots were negative or both roots were imaginary */

  /* Calculate the positive path length to the xy plane.  If either the
     photon is travelling in the xy plane or if the intercept is in the negative
     direction set the path length to infinity */

  if (pp.x[2] * pp.lmn[2] >= 0)
    s_to_zero = VERY_BIG;
  else
    s_to_zero = (-pp.x[2] / pp.lmn[2]);


  if (i >= 0 && root[i] < s_to_zero)
    return (root[i]);           /*Because that implies
                                   the ray hit the cone before hitting the xy plane */
  return (s_to_zero);


}



/**********************************************************/
/**
 * @brief      Calculate the path length along a line of sight defined by
 * 	a photon p to a sphere centered on the origin.
 *
 * @param [in] double  r   The radius of the sphere
 * @param [in] struct photon *  p   The photon
 * @return     The distance to the sphere.
 * 	If
 * 	the photon does not hit the sphere return a large number VERY_BIG
 *
 * @details
 * The routine simply solves the quadratic equation for calculating
 * the distance to the sphere.  It returns the smalles positive root
 *
 * ### Notes ###
 *
 *
 **********************************************************/

double
ds_to_sphere (r, p)
     double r;
     struct photon *p;
{
  int i;
  double a, b, c, root[2];

  a = 1.;
  b = 2. * dot (p->x, p->lmn);
  c = dot (p->x, p->x) - r * r;

  i = quadratic (a, b, c, root);
/* root[i] is the smallest positive root unless i is
negative in which case either both roots were negative or
both roots were imaginary */


  if (i >= 0)
    return (root[i]);           /*Because that implies
                                   the ray hit the sphere */

  return (VERY_BIG);
}




/**********************************************************/
/**
 * @brief      Find the positive distance to
 *    	a sphere centered at x with radius r
 *
 * @param [in] double  x[]   The center of the sphere
 * @param [in] double  r   The radius of the sphere
 * @param [in] struct photon *  p   A photon
 * @return     The distance to the sphere.
 * 	If
 * 	the photon does not hit the sphere return a large number VERY_BIG
 *
 *
 * @details
 * The routine simply solves the quadratic equation for calculating
 * the distance to the sphere.  It returns the smalles positive root
 *
 *
 * ### Notes ###
 * It is not clear we need two routines like this, one which
 * assumes that the sphere is centered on 0 and one which does
 * not.
 *
 **********************************************************/

double
ds_to_sphere2 (x, r, p)
     double x[], r;
     struct photon *p;
{
  double a, b, c, root[2], delta[3];
  int i;
  double dot ();

  vsub (p->x, x, delta);

  a = 1.;
  b = 2. * dot (delta, p->lmn);
  c = dot (delta, delta) - r * r;

  i = quadratic (a, b, c, root);        /* root[i] is the smallest positive root unless i is
                                           negative in which case either both roots were negative or both roots were imaginary */


  if (i >= 0)
    return (root[i]);           /*Because that implies
                                   the ray hit the sphere */

  return (VERY_BIG);
}




/**********************************************************/
/**
 * @brief      This solves a simple  quadratic (or if a is zero linear) equation.
 *
 * @param [in] double  a   The constant which multiples x**2 in the quadratic equation
 * @param [in] double  b   The constant which multiples x
 * @param [in] double  c   The constant
 * @param [out] double  r[]   The roots of the quadratic equation
 * @return   The returns depend on the a status tha depends on the nature of the solutions
 * The return is set up
 *    to make it easy to identify the smallest positive root if one exists.  The routine returns
 *    a negative number if both roots are negative or imaginary.  More specifically the routine returns
 *  *  -1 -> both roots are imaginary
 *  * -2 -> both roots are negative or 0
 *  * 0 -> the first root is the smallest positive  root
 *  * 1 -> the second root is the smallest positive root
 *
 * @details
 * The equation that is solved is ax**2+bx+c=0.
 *
 * ### Notes ###
 *
 **********************************************************/

int
quadratic (a, b, c, r)
     double a, b, c, r[];
{
  double q, z;

  if (a == 0.0)
  {                             /* Then it's not really a quadratic but we can solve it
                                   anyway */
    if (b == 0.0)
    {
      r[0] = r[1] = -99.;
      return (-1);              /* The roots are extremely imaginary, since both a a b were 0 */
    }

    r[0] = r[1] = (-c / b);     // Then it was a linear equation. Setting both roots to the same thing could be a problem ksl
    if (r[0] < 0.0)
      return (-2);              /* Generally speaking we are not interested in
                                   negative distances */
    else
      return (0);
  }

  if ((q = b * b - 4. * a * c) < 0.0)
  {
    r[0] = r[1] = -99.;
    return (-1);                /* both roots are imaginary */
  }
  q = sqrt (q);
  z = 0.5 / a;

  r[0] = (-b - q) * z;
  r[1] = (-b + q) * z;


  if (r[0] > 0.0 && (r[0] < r[1] || r[1] <= 0.0))
    return (0);                 /* r[0] is smallest positive root */
  if (r[1] > 0.0 && (r[1] < r[0] || r[0] <= 0.0))
    return (1);                 /* r[1] is smallest positive root */
  return (-2);                  /* both roots are negative */

  /* x1 should be the smallest positive root for most applications */
}





/**********************************************************/
/**
 * @brief      calculates the distance of a photon must travel to hit the plane.
 *
 * @param [in] struct plane *  pl   A structure which descibes a plane in terms of a position and vector normal to the surface
 * @param [in] struct photon *  p   A photon
 * @param [in] int force_postive_z  Whether to force the photon to be moving towards the positve z quadrant (appropriate for planes that 
               are symmetric about z=0, but only defined once, such as the windplanes.)
 * @return     The distance to the plane, or VERY_BIG if the photon does not hit the plane (in the positive direction)
 *
 * @details
 * Calculate the distance a photon must travel to hit a plane.  The plane, just for simplicity,
 * is defined as a photon structure since you need both a point and a direction to define it.
 * A plane can be defined by a position x_p and a normal lmn_p.  If the photon ray is then
 * defined by x=x_v+s lmn_v and the allowed values of s are determined by the equation
 *
 * (x_v+s lmn_v - x_p)  .  lmn_p=0 where . implies the dot-product.   The routine returns
 * VERY_BIG if the photon ray does not intersect the plane .
 *
 * ### Notes ###
 *
 **********************************************************/

double
ds_to_plane (pl, p, force_positive_z)
     struct plane *pl;
     struct photon *p;
     int force_positive_z;
{
  double denom, diff[3], numer;
  struct photon ptest;

  stuff_phot (p, &ptest);
  if (ptest.lmn[2] < 0 && (force_positive_z == TRUE))
  {
    ptest.x[2] = -ptest.x[2];   /* force the photon to be in the positive x,z quadrant */
    ptest.lmn[2] = -ptest.lmn[2];       /* force the photon to moving in the positive z direction */
  }

  if ((denom = dot (ptest.lmn, pl->lmn)) == 0)
    return (VERY_BIG);

  vsub (pl->x, ptest.x, diff);

  numer = dot (diff, pl->lmn);

  return (numer / denom);
}




/**********************************************************/
/**
 * @brief      This routine calculates the distance a photon has to be moved to
 * 	reach the point of closest approach to a point described by x and also
 * 	calculates the distance of closest approach (i.e. the impact parameter).
 *
 * @param [in] double  x[]   A position in space
 * @param [in] struct photon *  p   A photon bundle
 * @param [out] double *  impact_parameter   The distance of closest approach of the photon
 * @return     The distance the photon would have to be moved to reach the point of closest approach,
 * which can be negative
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

double
ds_to_closest_approach (x, p, impact_parameter)
     double x[];                /* point for which impact parameter is calculated */
     struct photon *p;          /* Photon ptr of interest */
     double *impact_parameter;  /* distance of ray to point a closest approach */
{
  double diff[3], s, result[3];
  double length (), dot ();

  vsub (p->x, x, diff);
  s = -dot (diff, p->lmn);

  vmove (diff, p->lmn, s, result);
  *impact_parameter = length (result);

  return (s);
}




/**********************************************************/
/**
 * @brief      calculate the distance a photon must travel to encouter a cylinder
 *   that is centered on the z axis.
 *
 * @param [in] double  rho   the size of the cylinder
 * @param [in] struct photon *  p   a photon ptr
 * @return     The smallest positive distance to the cylinder. If the
 *   ray does not hit the cylinder a VERY_BIG is
 *   returned.
 *
 * @details
 * Given a position and a direction the distance to a cylinder,
 * solve the simple quadratic equation that gives the distance
 * to the cylinder.
 *
 * ### Notes ###
 *
 **********************************************************/

double
ds_to_cylinder (rho, p)
     double rho;
     struct photon *p;
{
  double a, b, c, root[2];
  int i;

  a = 1. - p->lmn[2] * p->lmn[2];
  b = 2. * (p->lmn[0] * p->x[0] + p->lmn[1] * p->x[1]);
  c = (p->x[0] * p->x[0] + p->x[1] * p->x[1]) - rho * rho;

  i = quadratic (a, b, c, root);

/* root[i] is the smallest positive root unless i is
negative in which case either both roots were negative or
both roots were imaginary */

  if (i >= 0)
    return (root[i]);
  return (VERY_BIG);
}
