
/***********************************************************/
/** @file   cv.c
 * @author  ksl
 * @date   2019 February  
 * @brief  Calculate parameter to use a default inputs of
 * models of catalclysmic variables
 *
 * These routines are largely independent of Sirocco as a whole
 * and are only intended to intialize variables to be used in 
 * Sirocco
***********************************************************/

#include	<math.h>
#include	<stdio.h>
#include	<strings.h>
#define  PI  3.141592
#define MSOL    1.989e33
#define GRAV    6.67e-8
#define VLIGHT       2.997e10
#define YEAR    3.1556925e7
#define SIGMA   5.66956e-5


#define LINELENGTH 132

/* Simple program to estimate the WD radius from the mass in msol */


/**********************************************************/
/**
 * @brief  Calculate the radius of a WD from its mass using a standard
 * mass radius relation 
 * @param  [in] double m  The mass of the wd in gm
 *
 * @return The predicted radius of the WD (in cm)
 * 
 *
 * @details     
 *
 * This is the mass radius relationship given by 
 * Anderson, N 1988 (ApJ 225, 266)
 *
**********************************************************/
double
wdrad (m)
     double m;
{
  double r;

  m /= MSOL;                   

  r = pow ((m / 1.458), 4. / 3.);
  r = pow (1 - r, 0.47);
  r = 8.91e8 * pow (m, -1. / 3.) * r;

  return (r);
}


/**********************************************************/
/**
 * @brief   Calculate a plausible value for the outer radius of the disk
 * @param  [in] double m1  mass of the WD in grams
 * @param  [in] double m2  mass of the secondary in grams
 * @param  [in] double period  period of the system in seconds 
 *
 * @return An estimate of the size of the disk in cm
 * 
 *
 * @details     
 * 
 * The routine calculates the size of the Roche lobe of the primary
 * and returns 0.9 of this value as a plausible estimate of the 
 * size of the disk
 *
 *
 * ### Notes ###
 *
 * Various other choices are possible, e.g the circularization radius
 * as discussed in Frank, King and Raine (Accretion Power in Astophysics)
 * 
**********************************************************/
double
diskrad (m1, m2, period)
     double m1, m2, period;
{

  double t2p, x, a;
  double rlobe1, q;
  double roche2 ();

  q = m2 / m1;

  t2p = period / (2 * PI);
  x = GRAV * t2p * t2p * (m1 + m2);

  a = pow (x, 1. / 3.);

  rlobe1 = roche2 (1. / q, a);
  return (rlobe1 * 0.9);


}


/**********************************************************/
/**
 * @brief   The approximate size of the Roche lobe of the secondary
 * @param  [in] double q  The mass ratio
 * @param  [in] double a  The binary separation
 *
 * @return The approximate size of the Roche lobe of the secondary in cm
 * 
 *
 * @details     
 * See Frank, King & Raine, equation 4.6a 
 *
**********************************************************/
double
roche2 (q, a)
     double q, a;

{
  double rouche;
  double x, y;


  x = 0.49 * pow (q, 2. / 3.);
  y = 0.6 * pow (q, 2. / 3.) + log (1 + pow (q, 1. / 3.));
  rouche = a * x / y;
  return (rouche);
}


/**********************************************************/
/**
 * @brief   Calculate the log of the gravity of a star given the mass and radius
 * @param  [in] double mass  The mass of the star in grams
 * @param  [in,out] double rwd  The radius of the star in cm
 *
 * @return The log of the gravity, a number typically between 4 and 8 (in cgs unists)
 * 
 *
 * @details     
 *
 * A very straightforward calculation of the gravity of a star
 *
 * ### Notes ###
 *
 * It is likely this code is duplicated elsewhere in Sirocco
 * 
**********************************************************/
double
logg (mass, rwd)
     double mass, rwd;
{
  double gravity;
  gravity = GRAV * mass / (rwd * rwd);
  gravity = log10 (gravity);

  return (gravity);
}
