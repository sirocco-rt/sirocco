
/***********************************************************/
/** @file  photon_gen.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  Primary routines for creating photons for use in
 * the radiative transfer calculation.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


#define PRINT_OFF 0
#define PRINT_ON  1


/* THE NEXT FEW ROUTINES PERTAIN ONLY TO THE DISK */




/**********************************************************/
/**
 * @brief      Generate disk photons
 *
 * @param [out] PhotPtr  p   The entire photon structure
 * @param [in] double  weight   The weight of photons to generate
 * @param [in] double  f1   The minimum frequency
 * @param [in] double  f2   The maxnimum frequency
 * @param [in] int  spectype   The spectrum type to generate: 0 is bb
 * @param [in] int  istart   The starting point in p for generating photons
 * @param [in] int  nphot   The number of photons to generate
 * @return     Always returns 0
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

int
photo_gen_disk (p, weight, f1, f2, spectype, istart, nphot)
     PhotPtr p;
     double weight;
     double f1, f2;
     int spectype;
     int istart, nphot;
{

  double freqmin, freqmax;
  int i, iend;
  double planck ();
  double r, z, theta, phi;
  int nring = 0;
  double north[3];
  double fcol;

  if ((iend = istart + nphot) > NPHOT)
  {
    Error ("photo_gen_disk: iend %d > NPHOT %d\n", iend, NPHOT);
    Exit (0);
  }
  if (f2 < f1)
  {
    Error ("photo_gen_disk: Can't do anything if f2 %g < f1 %g\n", f2, f1);
    Exit (0);
  }
  Log_silent ("photo_gen_disk creates nphot %5d photons from %5d to %5d \n", nphot, istart, iend);
  freqmin = f1;
  freqmax = f2;
  for (i = istart; i < iend; i++)
  {
    p[i].origin = PTYPE_DISK;   // identify this as a disk photon
    p[i].frame = F_LOCAL;
    p[i].w = weight;
    p[i].istat = p[i].nscat = p[i].nrscat = p[i].nmacro = 0;
    p[i].tau = 0;
    p[i].nres = p[i].line_res = -1;     // It's a continuum photon
    p[i].nnscat = 1;
    p[i].np = istart + i;
    if (geo.reverb_disk == REV_DISK_UNCORRELATED)
      p[i].path = 0;            //If we're assuming disk photons are uncorrelated, leave them at 0


    /* First set the photon postion for the  special case of searchlight mode */
    if (modes.searchlight && geo.ioniz_or_extract == CYCLE_EXTRACT)
    {
      if (i == istart)
      {
        // we need to find the ring associated with the searchlight x
        nring = 0;
        while (disk.r[nring] < geo.searchlight_x[0] && nring < NRINGS - 1)
          nring++;

        if (nring == NRINGS)
        {
          Error ("disk_photon_gen: Trying to generate searchlight photons beyond the disk\n");
          Exit (1);
        }
      }
      stuff_v (geo.searchlight_lmn, p[i].lmn);
      stuff_v (geo.searchlight_x, p[i].x);
    }
    else
    {

/* The ring boundaries are defined so that an equal number of photons are
 * generated in each ring.  Howver, there is a possibility that the number
 * of photons to be generated is small, and therefore we, we still randomly
 * generate photon.  04march -- ksl
 */
      nring = random_number (0.0, 1.0) * (NRINGS - 1);


      if ((nring < 0) || (nring > NRINGS - 2))
      {
        Error ("photon_gen: photon launch out of bounds. nring = %d\n", nring);
        Exit (0);
      }

      disk.nphot[nring]++;

/* The next line is really valid only if dr is small.  Otherwise one
 * should account for the area.  But haven't fixed this yet ?? 04Dec
 */

      r = disk.r[nring] + (disk.r[nring + 1] - disk.r[nring]) * random_number (0.0, 1.0);

      /* Generate a photon in the plane of the disk a distance r */


      phi = 2. * PI * random_number (0.0, 1.0);

      p[i].x[0] = r * cos (phi);
      p[i].x[1] = r * sin (phi);


      z = 0.0;
      north[0] = 0;
      north[1] = 0;
      north[2] = 1;

      /* Correct photon direction of a vertically extended disk
       */

      if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
      {
        if (r == 0)
          theta = 0;
        else
        {
          z = zdisk (r);
          theta = atan ((zdisk (r * (1. + EPSILON)) - z) / (EPSILON * r));
        }
        north[0] = (-cos (phi) * sin (theta));
        north[1] = (-sin (phi) * sin (theta));
        north[2] = cos (theta);

      }

      if (random_number (-0.5, 0.5) > 0.0)
      {                         /* Then the photon emerges in the upper hemisphere */
        p[i].x[2] = (z + EPSILON);
      }
      else
      {
        p[i].x[2] = -(z + EPSILON);
        north[2] *= -1;
      }

      randvcos (p[i].lmn, north);
    }


    /* Note that the next bit of code is almost duplicated in photo_gen_star.  It's
     * possilbe this should be collected into a single routine   080518 -ksl
     */

    if (spectype == SPECTYPE_BB)
    {
      p[i].freq = planck (disk.t[nring], freqmin, freqmax);
    }
    else if (spectype == SPECTYPE_BB_FCOL)
    {
      fcol = disk_colour_correction (disk.t[nring]);
      p[i].freq = planck (fcol * disk.t[nring], freqmin, freqmax);
    }
    else if (spectype == SPECTYPE_UNIFORM)
    {
      p[i].freq = random_number (freqmin, freqmax);
    }

    else if (spectype == SPECTYPE_MONO)
    {
      p[i].w = 1.;
      p[i].freq = geo.mono_freq;
    }

    else
    {                           /* Then we will use a model which was read in */
      p[i].freq = one_continuum (spectype, disk.t[nring], log10 (disk.g[nring]), freqmin, freqmax);
    }

    if (p[i].freq < freqmin || freqmax < p[i].freq)
    {
      Error_silent ("photo_gen_disk: phot no. %d freq %g out of range %g %g\n", i, p[i].freq, freqmin, freqmax);
    }
    /* Now Doppler shift this. Use convention of dividing when going from rest
       to moving frame */


    // if (modes.save_photons && geo.ioniz_or_extract == CYCLE_EXTRACT)
    //   save_photons (&p[i], "gen_local");

    /* When given the same input photons the transformation is made in place */
    if (local_to_observer_frame_disk (&p[i], &p[i]))
    {
      Error ("photon_gen: Frame Error\n");
    }
    // if (modes.save_photons && geo.ioniz_or_extract == CYCLE_EXTRACT)
    //   save_photons (&p[i], "gen_obs");




  }


  return (0);
}




/**********************************************************/
/**
 * @brief      write information about the way photons are
 * created and absropted by the disk to a file
 *
 * @param [in out] char  filename[]   The name of the file
 * @param [in out] char  mode[]   A switch to determine whther to
 * write a new file or append to an existing file
 * @return     Always returns 0
 *
 * @details
 * The routine normal writes the disk.diag file.  It provides
 * information not only about the rings used, but also about
 * the number of photons that hit each ring,
 *
 * ### Notes ###
 *
 *
 **********************************************************/

int
disk_photon_summary (filename, mode)
     char filename[], mode[];
{
  FILE *fopen (), *ptr;
  int n;
  double x;
  if (mode[0] == 'a')
    ptr = fopen (filename, "a");
  else
    ptr = fopen (filename, "w");
  fprintf (ptr, "Ring     r      t      nphot   dN/dr\n");
  for (n = 0; n < NRINGS; n++)
  {
    if (n == 0)
      x = 0;
    else
      x = disk.nphot[n] / (disk.r[n] - disk.r[n - 1]);
    fprintf (ptr, "%d %8.2e %8.2e %8d %8.2e %8d\n", n, disk.r[n], disk.t[n], disk.nphot[n], x, disk.nhit[n]);
  }

  fclose (ptr);
  return (0);
}

/* THESE ROUTINES ARE FOR THE BOUNDARY LAYER */
