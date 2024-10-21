/***********************************************************/
/** @file  random.c
 *  @author ksl
 * @date   May, 2018
 *
 * @brief  Various general purpose routines for generating random
 * numbers including varous routines for generating 
 * randomly oriented vectors
 *
 * These routines should be kept SEPARATE from routines that 
 * require the Python specific
 * structures in sirocco.h so that it is possible to test 
 * them more easily.
 *
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/stat.h>
#include <string.h>
#include <errno.h>

#include "constants.h"
#include "math_struc.h"
#include "math_proto.h"
#include "log.h"

/* A basis is defined such that if x is a 3 vector as expressed an unprimed cartesian coordinate
   frame, and if y is the same vector in some rotated frame, then
   x[i] = a[i][j] y[j]
   Basis is defined in sirocco.h
 */

gsl_rng *rng = NULL;            // pointer to a global random number generator
char rngsave_file[LINELENGTH];


/**********************************************************/
/** 
 *
 * @brief Get a 3 vector "a" whose direction will be random and whose
 * lenght will be r
 *
 * @param [out] double a[] The resulting 3 vector      
 * @param [in] double  r   desired length of the vector
 * @return     Always retursn 0                       
 *
 * @details
 *
 * The current routine creates a random vector
 * in spherical coordinates, noting that an element of solid angle i
 * s sin(theta) d(theta) d(phi).  But sin(theta)d(theta) is simply 
 * -d(cos(theta)), and so a random vector is created by generating
 * a random number between 0 and 2 PI representing phi and between
 * -1 and 1 representing cos(phi).   
 *
 * ### Notes ###
 *
 *
 **********************************************************/


int
randvec (a, r)
     double a[], r;
{

  double costheta, sintheta, phi, sinphi, cosphi;

  phi = 2. * PI * random_number (0.0, 1.0);
  sinphi = sin (phi);
  cosphi = cos (phi);
  costheta = random_number (-1.0, 1.0);
  sintheta = sqrt (1. - costheta * costheta);
  a[0] = r * cosphi * sintheta;
  a[1] = r * sinphi * sintheta;
  a[2] = r * costheta;

  return (0);

}


double zzz[] = { 0.0, 0.0, 1.0 };


int init_vcos = 0;

/**********************************************************/
/** 
 * @brief Create a photon direction "lmn" in the hemisphere with the 
 * vector "north pointing to the north pole based on  the Eddingon
 * approximation
 *
 * @param [out] double lmn[] The resulting 3 vector containing the correct 
 * direction cosines      
 * @param [in] double  north [] The direction  of local north
 * @return     Always retursn 0                       
 *
 * @details
 *
 * Create a photon direction "lmn" in the hemisphere with the vector "north" pointing to the "north
 * pole" pole of the hemisphere in the case where the photon is originating in a photosphere.
 * Another way of saying this is that north is the normal to surface at the point
 * at which the photon originates.  
 * 
 * The photon directions will be distributed according to the Eddington 
 * approximation
 *
 *
 * ### Notes ###
 * The routine calls vcos for generating the cdf that is used
 * vcos contains the Eddington 
 * approximation (aka linear limb darkening)
 *
 * Jumps were added to include more points near 90 degrees.
 *
 *
 **********************************************************/

int
randvcos (lmn, north)
     double lmn[], north[];
{
  double x[3];                  /* the photon direction in the rotated frame */
  double l, m, n;               /* the individual direction cosines in the rotated frame */
  double q, jumps[5];
  struct basis nbasis;
  int echeck;
  double phi;

  /*** 
   * ### Programming Comment ###
   * pdf_gen_from_func still uses jumps, so this is OK, but it may not be
   * necessary as PDFSTEPS has been increased to 10000 in cdf.c  180715 ksl.
   */

  if (init_vcos == 0)
  {
    jumps[0] = 0.01745;
    jumps[1] = 0.03490;
    jumps[2] = 0.05230;
    jumps[3] = 0.06976;
    jumps[4] = 0.08716;

    if ((echeck = cdf_gen_from_func (&cdf_vcos, &vcos, 0., 1., 5, jumps)) != 0)
    {
      Error ("Randvcos: return from cdf_gen_from_func %d\n", echeck);;
    }
    init_vcos = 1;
  }


  n = cdf_get_rand (&cdf_vcos);
  q = sqrt (1. - n * n);


  phi = 2. * PI * random_number (0.0, 1.0);
  l = q * cos (phi);
  m = q * sin (phi);

/* So at this point we have the direction cosines in a frame in which
the z axis is along the normal to the surface.  We must now put the
direction in the cartesian frame.  If north is in the +-z direction
this is simple. Otherwise one must do a coordinate rotation. */

  if (north[0] == 0 && north[1] == 0)
  {
    lmn[0] = l;
    lmn[1] = m;
    if (north[2] > 0)
      lmn[2] = n;
    else
      lmn[2] = -n;
  }
  else
  {
    create_basis (north, zzz, &nbasis); /* Create a basis with the first axis in 
                                           direction of "north" and the second axis in the 
                                           yz plane */
    x[0] = n;
    x[1] = l;
    x[2] = m;

    project_from (&nbasis, x, lmn);     /* Project the vector back to the standard
                                           frame */
  }
  return (0);

}



/**********************************************************/
/** 
 *
 * @brief get the probability density associated with Eddington limb darkening 
 * for cos theta
 *
 * @param [in] double x       cos theta
 * @param [in] void * params  Unused parameters for passing to the GSL integrator
 * @return     The probability density of the Eddington approximation at cos theta
 *
 * @details
 *
 *
 * ### Notes ###
 *
 * See Hubeny & Mihalas Equation 17.17  
 *
 * The extra factor of x arises from the fact that we also need to
 * account for the geometric factor at differnt angles.  If
 * there were no limb darkening, the second term, the whole
 * term in parenthesis could simply be dropped.
 * 
 * The extra factor of x arises from the fact that we want to
 * account for  the probability density for
 * all azimuthal angles
 *
 *
 **********************************************************/


double
vcos (double x, void *params)
{
  double a, b;
  double z;

  (void) params;

  a = 0.5;
  b = 1.5;
  z = x * (a * (1. + b * x));
  return (z);
}


int init_vdipole = 0;

/**********************************************************/
/** 
 *
 * @brief Create a photon direction "lmn" in with the 
 * vector "north" pointing in the direction of the the photon
 * before scttering
 *
 * @param [out] double lmn[] The resulting 3 vector containg the correct 
 * direction cosines      
 * @param [in] double  north [] The direction of the photon before scatterin
 * @return     Always returns 0                       
 *
 * @details
 *
 * 
 * The photon directions will be distributed according to the dipole    
 * approximation
 *
 *
 * ### Notes ###
 * The routine calls vdipole for generating the cdf that is used
 *
 * Jumps were added to include more points near 90 degrees.
 *
 * This routine was adapted from randvcos and retains some of 
 * that routines terminology for clarity.
 *
 *
 **********************************************************/

int
randvdipole (lmn, north)
     double lmn[], north[];
{
  double x[3];                  /* the photon direction in the rotated frame */
  double l, m, n;               /* the individual direction cosines in the rotated frame */
  double q, jumps[10];
  struct basis nbasis;
  int echeck;
  double phi;

  /*** 
   * ### Programming Comment ###
   * pdf_gen_from_func still uses jumps, so this is OK, but it may not be
   * necessary as PDFSTEPS has been increased to 10000 in cdf.c  180715 ksl.
   */

  if (init_vdipole == 0)
  {

    jumps[0] = 0.00010;
    jumps[1] = 0.00030;
    jumps[2] = 0.00100;
    jumps[3] = 0.00300;
    jumps[4] = 0.00500;

    jumps[5] = 1. - 0.00500;
    jumps[6] = 1. - 0.00300;
    jumps[7] = 1. - 0.00100;
    jumps[8] = 1. - 0.00030;
    jumps[9] = 1. - 0.00010;


    if ((echeck = cdf_gen_from_func (&cdf_vdipole, &vdipole, -1., 1., 10, jumps)) != 0)
    {
      Error ("Randvcos: return from cdf_gen_from_func %d\n", echeck);;
    }
    init_vdipole = 1;
    cdf_to_file (&cdf_vdipole, "Dipole");
  }


  n = cdf_get_rand (&cdf_vdipole);


  q = sqrt (1. - n * n);


  phi = 2. * PI * random_number (0.0, 1.0);
  l = q * cos (phi);
  m = q * sin (phi);

/* So at this point we have the direction cosines in a frame in which
the z axis we want to be in the direction of the initial travel    
direction in the cartesian frame.  If north is in the +-z direction
this is simple. Otherwise one must do a coordinate rotation. */

  if (north[0] == 0 && north[1] == 0)
  {
    lmn[0] = l;
    lmn[1] = m;
    if (north[2] > 0)
      lmn[2] = n;
    else
      lmn[2] = -n;

    Log ("ZZZZ  %e \n", n);
  }
  else
  {
    create_basis (north, zzz, &nbasis); /* Create a basis with the first axis in 
                                           direction of "north" and the second axis in the 
                                           yz plane */
    x[0] = n;                   /* Because the photon is travelling in this direction */
    x[1] = l;
    x[2] = m;

    project_from (&nbasis, x, lmn);     /* Project the vector back to the standard
                                           frame */

    renorm (lmn, 1.0);          /* Eliminate roundoff errors */

  }


  return (0);

}




/**********************************************************/
/** 
 *
 * @brief get the probablity density associated with dipole radiation
 * for cos theta
 *
 * @param [in] double x       cos theta
 * @param [in] void * params  Unused parameters for passing to the GSL integrator
 * @return     The probability density of the Eddingtog approximation at cos theta
 *
 * @details
 *
 *
 * ### Notes ###
 *
 * The probability densits is not  normalized, necessarily.
 *
 *
 **********************************************************/


double
vdipole (double cos_theta, void *params)
{


  double pdf;

  pdf = 0.5 * (1 + cos_theta * cos_theta);


  return (pdf);
}



/**********************************************************/
/** 
 * @brief	Sets up a random number generator 
 *
 * @param [in] seed  The seed to set up the generator
 * @return 	     0
 *
 * Sets up a random number generator. The resulting generator
 * is addressed by the pointer rng, which is set up as a local
 * variable at the top of the file. The type of generator is
 * set in the call to gsl_rng_alloc - currently a meursenne
 * twister
 *
 * ###Notes###
 * 2/18	-	Written by NSH
 * The number generator uses  GSL Meursenne twirster
***********************************************************/


int
init_rand (seed)
     int seed;
{
  rng = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_set (rng, seed);
  return (0);
}

/**********************************************************/
/**
 * @brief  Initialise the RNG directory structure.
 *
 * @details
 *
 * Creates a hidden folder named .rng_root where root.rng_save files are
 * stored.
 *
 **********************************************************/

void
init_rng_directory (root, rank)
     char *root;
     int rank;
{
  int err;
  char dir_name[LINELENGTH];
  char file_name[LINELENGTH];

  sprintf (dir_name, ".rng_%.100s/", root);

  err = mkdir (dir_name, 0777);
  if (err)
  {
    if (errno != EEXIST)
    {
      perror ("init_rng_directory");
      Exit (1);
    }
  }

  sprintf (file_name, "%.50s%.50s_%d.rng_save", dir_name, root, rank);
  strcpy (rngsave_file, file_name);
}

/**********************************************************/
/**
 * @brief  Dump the current GSL RNG state to a binary file.
 *
 * @details
 *
 * Saves the RNG state to the file root.gsl_save. The point of this
 * is mostly for debugging purposes, being able to keep the same RNG if
 * restarting the model.
 *
 **********************************************************/

void
save_gsl_rng_state ()
{
  FILE *file;

  if ((file = fopen (rngsave_file, "w")) == NULL)
  {
    Error ("save_gsl_rng_state: unable to open %s to write RNG to file\n", rngsave_file);
    return;
  }

  if (gsl_rng_fwrite (file, rng))
  {
    Error ("save_gsl_rng_state: gsl_rng_fwrite failed to write RNG state to file\n");
  }
  {
    Log ("GSL RNG state saved to %s\n", rngsave_file);
  }

  if (fclose (file))
  {
    Error ("save_gsl_rng_state: there was a problem when closing %s\n", rngsave_file);
  }
}

/**********************************************************/
/**
 * @brief  Load a dump GSL RNG save into the current RNG.
 *
 * @details
 *
 * Retreive a dumped RNG state from the file root.gsl_save. The point of this
 * is mostly for debugging purposes, being able to keep the same RNG if
 * restarting the model.
 *
 * gsl_rng_fread will modify rng, the global rng variable, so we do not need
 * to re-set the rng system. However note that we still need to allocate
 * memory for it. The previous rng memory is deallocated if, for some reason,
 * it is not a NULL pointer.
 *
 **********************************************************/

void
reload_gsl_rng_state ()
{
  FILE *file;

  if (rng != NULL)
    gsl_rng_free (rng);

  rng = gsl_rng_alloc (gsl_rng_mt19937);

  if ((file = fopen (rngsave_file, "r")) == NULL)
  {
    Error ("reload_gsl_rng_state: unable to open %s so using a new seed\n", rngsave_file);
    return;
  }

  if (gsl_rng_fread (file, rng))
  {
    Error ("reload_gsl_rng_state: gsl_rng_fread failed to read the RNG state from %s so using a new seed\n", rngsave_file);
  }
  else
  {
    Log ("GSL RNG state loaded from %s\n", rngsave_file);
  }

  if (fclose (file))
  {
    Error ("reload_gsl_rng_state: there was a problem when closing %s\n", rngsave_file);
  }
}

/**********************************************************/
/** 
 * @brief	Gets a random number from the generator set up in init_rand
 *
 * @param [in] min			The minimum value to be generated 
 * @param [in] max			The maximum value to be generated 
 * @return [out] x 			The generated number
 *
 * Produces a number from min to max (exclusive).
 *
 * ###Notes###
 * 2/18	-	Written by NSH
***********************************************************/


double
random_number (double min, double max)
{
  double num = gsl_rng_uniform_pos (rng);
  double x = min + ((max - min) * num);
  return (x);
}
