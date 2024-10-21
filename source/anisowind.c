
/***********************************************************/
/** @file  anisowind.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  Routines to implement anisotropic scattering in
 * the wind
 *
 *
 * Python supports two ways to determine a new photon
 * direction when a photon scatters in thie wind.  These
 * include 
 *
 * * SCATTER_MODE_ISOTROPIC -isotropic scattering (randvec),
 * * SCATTER_MODE_THERMAL - thermally-broadened anisotropic scattering
 * (randwind_thermal_trapping).
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief the routine which chooses
 *   a new anisotropic direction in geo.scatter_mode = SCATTER_MODE_THERMAL
 *
 * @param [in,out] PhotPtr  p   The photon being scattered
 * @param [out] int *  nnscat   The number of times the photon
 * scattered internally before escaping the local scattering region
 *
 * @return     0 for success. Also modifies the photon ptr p
 *   to reflect new direction (p->lmn), and nnscat, which
 *   should be copied to the photon structure after calling
 *   this routine.
 *
 * @details
 * This routine uses a rejection method to choose a direction
 * so that the probability distribution of directions generated
 * reflects the probability of escape along each direction in accordance
 * with the sobolev optical depth.
 *
 * ### Notes ###
 * The resonance that caused the scatter must be stored in the
 * photon bundle.
 *
 * The name of the routine is something of a misnomer. Pure
 * sobolev optical depths are used.  The temperature in the
 * cell does not come into the calculation.
 *
 * 1904 - ksl - In certain situations, see issue #505, the photon can
 * get trapped in the while loop.  To prevent his from hanging
 * the routine bails after NSCAT_MAX trips through the while loop
 * and chooses the direction with the lowest tau it has found. 
 * From a physics perspective, this is not an ideal treatment of
 * the problem ans so if this occurs a lot of times, something
 * better needs to be done.
 **********************************************************/

#define NSCAT_MAX 10000


int
randwind_thermal_trapping (p, nnscat)
     PhotPtr p;
     int *nnscat;
{
  double tau_norm, p_norm;
  double tau, dvds, z, ztest;
  double z_prime[3], z_min[3];
  int nscat;
  double tau_min = VERY_BIG;
  double dvds_max, dvds_test = 0;
  WindPtr one;

  /* find the wind pointer for the photon */
  one = &wmain[p->grid];

  /* we want to normalise our rejection method by the escape
     probability along the vector of maximum velocity gradient.
     First find the sobolev optical depth along that vector 

   */


  dvds_max = get_dvds_max (p);
  tau_norm = sobolev (one, p->x, -1.0, lin_ptr[p->nres], dvds_max);

  /* then turn into a probability. p_norm is the probability a
     photon will escape along the line so sight with maximum
     velocity graident, which is the easiest way for aphoton
     to escape.
   */
  p_norm = p_escape_from_tau (tau_norm);

  /* Throw error if p_norm is 0 */
  if (p_norm <= 0)
    Error ("randwind_thermal_trapping: p_norm is %8.4e in cell %i", p_norm, one->nplasma);

  /* JM 1406 -- we increment nnscat here, and it is recorded in the photon
     structure. This is done because we actuall have to multiply the photon weight
     by 1/mean escape probability- which is nnscat. this is done in trans_phot.c
     before extract is called.
   */
  *nnscat = *nnscat - 1;
  nscat = 0;

  ztest = 1.0;
  z = 0.0;

  /* rejection method loop, which chooses direction and also calculates nnscat */
  while (ztest > z && nscat < NSCAT_MAX)
  {
    *nnscat = *nnscat + 1;
    randvec (z_prime, 1.0);
    stuff_v (z_prime, p->lmn);

    dvds = dvwind_ds_cmf (p);

    if (dvds > dvds_test)
    {
      dvds_test = dvds;
    }

    if (dvds > dvds_max)
    {
      Error ("randwind_thermal trapping: dvds (%e) > dvds_max (%e) ratio %e in grid %d at %e %e %e\n",
             dvds, dvds_max, dvds / dvds_max, p->grid, p->x[0], p->x[1], p->x[2]);
    }
    tau = sobolev (one, p->x, -1.0, lin_ptr[p->nres], dvds);

    if (tau < tau_min)
    {
      stuff_v (z_prime, z_min);
      tau_min = tau;
    }

    z = p_escape_from_tau (tau);        /* probability of escaping in a given (random) direction */
    /* generate random number, normalised by p_norm as dvds_max is worked out with a sample of directions) */
    ztest = random_number (0.0, 1.0) * p_norm;

    nscat++;
  }

  if (nscat == NSCAT_MAX)
  {
    stuff_v (z_min, p->lmn);    // copy to photon pointer
    Error
      ("randwind_thermal_trapping: photon %d needed > %5d directions in cell %4d at %9.2e %9.2e %9.2e with p_norm %9.2e zmin %9.2e tau_min %9.2e tau_norm %9.2e dvds_test %9.2e dvds_max %9.2e\n",
       p->np, nscat, p->grid, p->x[0], p->x[1], p->x[2], p_norm, z, tau_min, tau_norm, dvds_test, dvds_max);
  }

  return (0);
}
