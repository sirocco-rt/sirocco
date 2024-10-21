
/***********************************************************/
/** @file  disk_init.c
 * @author ksl
 * @date   April, 2020
 *
 * @brief  Primary routines for initializing the disk
 * as on a luminoisity weighted basis
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


#define STEPS 100000



/**********************************************************/
/**
 * @brief      calculates the total luminosity and the luminosity between freqqmin and freqmax
 * 	of the disk.  More importantly  divides the disk into annulae such that each
 * 	annulus contributes and equal amount to the lumionosity of the disk (within the frequency
 * 	limits).  Thus  initializes the structure "disk".
 *
 * @param [in] double  rmin   The minimum radius of the disk
 * @param [in] double  rmax   The maximum radius of the disk
 * @param [in] double  m   mass of central object
 * @param [in] double  mdot   mass accretion rate
 * @param [in] double  freqmin   The minimum frequency
 * @param [in] double  freqmax   The maximum frequency
 * @param [in] int  ioniz_or_extract   A flag indicating whether this is an ionization or
 * a detailed spectral cycle (used to determine the spectral type to use)
 * @param [out] double *  ftot   The band limited luminosity in the freqency interval
 * @return     the total luminosity of the disk
 *
 * @details
 * This routine assumes the temperature distribution for the disk is
 * that of a simple Shakura-Sunyaev disk, and uses this to determine
 * the band limited luminosity of the disk.  Additionally, it divides
 * the disk in the rings of the same band-limited luminosity, so that
 * equal numbers of photons can be generated from each ring.  (The
 * reason the disk has to be initilized mulitple times is because
 * the rings are different for different freqency intervals.)
 *
 * ### Notes ###
 * The information needed to generate photons from the disk is stored
 * in the disk structure.
 * The positional parameters x and v are at the edge of the ring,
 * but many of the other parameters (like temperature) are at the mid point.
 *
 * The calculation of the disk rings (which depends on the area) does not make
 * any allowances for a vertical disk structure.
 *
 *
 **********************************************************/

double
disk_init (rmin, rmax, m, mdot, freqmin, freqmax, ioniz_or_extract, ftot)
     double rmin, rmax, m, mdot, freqmin, freqmax, *ftot;
     int ioniz_or_extract;
{
  double t;
  double log_g;
  double v, dr, r;
  double logdr, logrmin, logrmax, logr;
  double f, ltot;
  double q1;
  int nrings;
  int spectype;
  double emit;
  double factor, fcol;




  /*
   * Compute the apparent luminosity of the disk.  This is not
   * actually used to determine how annulae are set up.  It is just
   * used to populate geo.ltot. It can change if photons hitting the
   * disk are allowed to raise the temperature
   */

  logrmax = log (rmax);
  logrmin = log (rmin);
  logdr = (logrmax - logrmin) / STEPS;

  for (nrings = 0; nrings < NRINGS; nrings++)
  {
    disk.nphot[nrings] = 0;
    disk.nphot[nrings] = 0;
    disk.r[nrings] = 0;
    disk.t[nrings] = 0;
    disk.nhit[nrings] = 0;
    disk.heat[nrings] = 0;
    disk.ave_freq[nrings] = 0;
    disk.w[nrings] = 0;
    disk.t_hit[nrings] = 0;
  }




  ltot = 0;

  for (logr = logrmin; logr < logrmax; logr += logdr)
  {
    r = exp (logr);
    dr = exp (logr + logdr) - r;
    t = teff ((r + 0.5 * dr) / rmin);
    ltot += t * t * t * t * (2. * r + dr) * dr;
  }
  geo.lum_disk_init = ltot *= 2. * STEFAN_BOLTZMANN * PI;


  /* Now establish the type of spectrum to create */

  if (ioniz_or_extract == CYCLE_EXTRACT)
    spectype = geo.disk_spectype;
  else
    spectype = geo.disk_ion_spectype;

  /* Next compute the band limited luminosity ftot */

  /*
   * The area of an annulus is  PI*((r+dr)**2-r**2) = PI * (2. * r +dr) *
   * dr. The extra factor of two arises because the disk radiates from
   * both of its sides.
   */

  q1 = 2. * PI;

  (*ftot) = 0;


  for (logr = logrmin; logr < logrmax; logr += logdr)
  {
    r = exp (logr);
    dr = exp (logr + logdr) - r;
    t = teff ((r + 0.5 * dr) / rmin);
    log_g = (geff ((r + 0.5 * dr) / rmin));
    v = sqrt (GRAV * geo.mstar / r);
    v /= VLIGHT;
    if (rel_mode == REL_MODE_FULL)
    {
      factor = sqrt (1. - v * v);
    }
    else
    {
      factor = 1.0;
    }

    if (spectype > -1)
    {
      emit = emittance_continuum (spectype, freqmin, freqmax, t, log_g);
    }
    else if (spectype == SPECTYPE_BB_FCOL)
    {
      fcol = disk_colour_correction (t);
      emit = emittance_bb (freqmin, freqmax, fcol * t) / pow (fcol, 4.0);
    }
    else
    {
      emit = emittance_bb (freqmin, freqmax, t);
    }
    (*ftot) += emit * (2. * r + dr) * dr * factor;
  }

  (*ftot) *= q1;



  /*
   * If *ftot is 0 in this energy range then all the photons come
   * elsewhere, e. g. the star or BL
   */

  if ((*ftot) < EPSILON)
  {
    Log_silent ("disk_init: Warning! Disk does not radiate enough to matter in this wavelength range\n");
    return (ltot);
  }
  /*
   * Now find the boundaries of the each annulus, which depends on the
   * band limited flux. Note that disk.v is calculated at the
   * boundaries, because vdisk() interporlates on the actual radius.
   */

  disk.r[0] = rmin;
  disk.v[0] = sqrt (GRAV * geo.mstar / rmin);
  nrings = 1;
  f = 0;

  for (logr = logrmin; logr < logrmax; logr += logdr)
  {
    r = exp (logr);
    dr = exp (logr + logdr) - r;
    t = teff ((r + 0.5 * dr) / rmin);
    log_g = (geff ((r + 0.5 * dr) / rmin));
    v = sqrt (GRAV * geo.mstar / r);
    v /= VLIGHT;

    if (rel_mode == REL_MODE_FULL)
    {
      factor = sqrt (1. - v * v);
    }
    else
    {
      factor = 1.0;
    }

    if (spectype > -1)
    {
      emit = emittance_continuum (spectype, freqmin, freqmax, t, log_g);
    }
    else if (spectype == SPECTYPE_BB_FCOL)
    {
      fcol = disk_colour_correction (t);
      emit = emittance_bb (freqmin, freqmax, fcol * t) / pow (fcol, 4.0);
    }
    else
    {
      emit = emittance_bb (freqmin, freqmax, t);
    }

    f += q1 * emit * (2. * r + dr) * dr * factor;

    /*
     * EPSILON to assure that roundoffs don't affect result of if
     * statement
     */
    if (f / (*ftot) * (NRINGS - 1) >= nrings)
    {
      if (r <= disk.r[nrings - 1])
        //If the radius we have reached is smaller than or equal to the last assigned radius - we make a tiny annulus
      {
        r = disk.r[nrings - 1] * (1. + 1.e-10);
      }
      disk.r[nrings] = r;
      disk.v[nrings] = sqrt (GRAV * geo.mstar / r);
      nrings++;
      if (nrings >= NRINGS)
      {
        break;
      }
    }
  }
  if (nrings < NRINGS - 1)
  {
    Error ("error: disk_init: Integration on setting r boundaries got %d nrings instead of %d\n", nrings, NRINGS - 1);
    Exit (0);
  }
  disk.r[NRINGS - 1] = exp (logrmax);
  disk.v[NRINGS - 1] = sqrt (GRAV * geo.mstar / disk.r[NRINGS - 1]);


  /* Now calculate the temperature and gravity of the annulae */

  for (nrings = 0; nrings < NRINGS - 1; nrings++)
  {
    r = 0.5 * (disk.r[nrings + 1] + disk.r[nrings]);
    disk.t[nrings] = teff (r / rmin);
    disk.g[nrings] = geff (r / rmin);
  }

  /* Wrap up by zeroing other parameters */
  for (nrings = 0; nrings < NRINGS; nrings++)
  {
    disk.nphot[nrings] = 0;
    disk.nhit[nrings] = 0;
    disk.heat[nrings] = 0;
    disk.ave_freq[nrings] = 0;
    disk.w[nrings] = 0;
    disk.t_hit[nrings] = 0;
  }
  geo.lum_disk = ltot;
  return (ltot);
}




/**********************************************************/
/**
 * @brief      Initialize a structure (qdisk) for recording information about photons/energy impinging
 * 	on the disk, which is stored in a disk structure called qdisk.
 *
 * @return     Always return zero
 *
 *
 * ###Notes###
 *
 * The information stored in qdisk can be used to modify the effective temperature
 * of the disk
 *
 * The reason the qdisk is needed as well as disk structure is that whenever the
 * 	wavelength bands are changed the radii in the disk structure are recalibratied.
 * 	We want qdisk to have fixed boundaries when this happens.
 *
 *
 * The annular are defined differently than in disk_init.  Here they are simply
 * logarithmically spaced; there they are spaced so that equal amounts of emission
 * are emitted from each annulus.
 *
 **********************************************************/

int
qdisk_init (rmin, rmax, m, mdot)
     double rmin, rmax, m, mdot;
{
  int nrings;
  double log_rmin, log_rmax, dlog_r, log_r;
  double r;

  log_rmin = log10 (disk.r[0]);
  log_rmax = log10 (disk.r[NRINGS - 1]);
  dlog_r = (log_rmax - log_rmin) / (NRINGS - 1);


  for (nrings = 0; nrings < NRINGS; nrings++)
  {
    log_r = log_rmin + dlog_r * nrings;
    qdisk.r[nrings] = pow (10, log_r);
  }


  for (nrings = 0; nrings < NRINGS; nrings++)
  {
    if (nrings < NRINGS - 1)
    {
      r = 0.5 * (qdisk.r[nrings + 1] + qdisk.r[nrings]);
    }
    else
    {
      r = qdisk.r[nrings];
    }
    qdisk.t[nrings] = teff (r / rmin);
    qdisk.g[nrings] = geff (r / rmin);
    qdisk.v[nrings] = sqrt (GRAV * geo.mstar / r);
    qdisk.heat[nrings] = 0.0;
    qdisk.nphot[nrings] = 0;
    qdisk.nhit[nrings] = 0;
    qdisk.w[nrings] = 0;
    qdisk.ave_freq[nrings] = 0;
    qdisk.t_hit[nrings] = 0;
  }

  return (0);
}


/**********************************************************/
/**
 * @brief      Reinitialize portions a structure (qdisk) for recording information about 
 *      photons/energy impinging
 * 	on the disk.
 *
 * @return     Always return zero
 *
 *
 * ###Notes###
 *
 * The structure qdisk records information about how many photons hit the disk
 * and what the impact of these photons would be if these
 * photons are used to modify the temperatures structure
 * of the disk
 *
 * This routine records information about the 
 * photons which were emitted from the disk 
 * and zeros various arrays that will record information
 * about photons that will hit the disk
 * 
 *
 **********************************************************/

int
qdisk_reinit (p)
     PhotPtr p;
{
  int nphot, i, n;
  double rho;
  struct photon pp;
  for (n = 0; n < NRINGS; n++)
  {
    qdisk.emit[n] = qdisk.nhit[n] = qdisk.heat[n] = qdisk.nphot[n] = qdisk.w[n] = qdisk.ave_freq[n] = 0;
  }

  for (nphot = 0; nphot < NPHOT; nphot++)
  {
    stuff_phot (&p[nphot], &pp);
    if (pp.origin_orig == PTYPE_DISK)
    {
      rho = sqrt (pp.x[0] * pp.x[0] + pp.x[1] * pp.x[1]);

      i = 0;
      while (rho > qdisk.r[i] && i < NRINGS - 1)
        i++;
      i--;                      /* So that the heating refers to the heating between i and i+1 */

      qdisk.emit[i] += pp.w;
      qdisk.nphot[i] += 1;


    }
  }

  return (0);
}


/**********************************************************/
/**
 * @brief      Save the information in the qdisk structure to a file
 *
 * @param [in] char *  diskfile   Name of the file whihc is writteen
 * @param [in] int  ichoice   A switch that just says whether the weighted
 *       frequency has been converted to a real frequency (0 = No, 1=Yes)
 * @return     Always returns 0
 *
 * The routine reformats the data about disk heating which has
 * been accumulated and writes it to a file
 *
 * ###Notes###
 *
 *
 * The file that is produced should be readable as an astropy
 * table
 *
 * Here is a definition of the columnts
 *
 * * r            the radius of the wirng
 * * v            the velocity of the ring
 * * zdisk        the zheight ot the rign
 * * t_disk       the temperature derived from viscous heating or the input grid
 * * g            the log of the gravity
 * * emit         the luminosity of photons emitted by the ring in this cycle
 * * nemit        the number photons emitted by the ring in this cycle
 * * heat         the total luminosity of photons that hit the ring
 * * nhit         the number of photon bundles hitting a ring
 * * ehit/emit    the ratio of the energy of photons that hit the ring/to those
 *              to those that are emitted by the ring
 * * t_heat       the temperature calculated from the energy of the
 *              photons hitting a ring
 * * t_freq       The temperature calculated from the weighted frequency of 
 *              photons hitting the ring
 * * W_freq       Assuming t_freq for the temperature, this gives ther
 *              ratio of the energy flux to LTE
 * * t_rerad      The temperature required to radiate away both the
 *              the viscous heating and the heating from photons
 *              in the next cycle
* 
 **********************************************************/

int
qdisk_save (diskfile, ichoice)
     char *diskfile;
     int ichoice;
{
  FILE *qptr;
  int n;
  double area, theat, ttot;
  double ratio[NRINGS];
  qptr = fopen (diskfile, "w");
  fprintf (qptr,
           "         r          v      zdisk    t_disk         g       emit      nemit      heat       nhit ehit/emit    t_heat    t_freq     W_freq   t_rerad\n");

  for (n = 0; n < NRINGS; n++)
  {
    if (n < NRINGS - 1)
    {
      // The factor of 2 comes from the fact that the disk has two sides
      area = (2. * PI * (qdisk.r[n + 1] * qdisk.r[n + 1] - qdisk.r[n] * qdisk.r[n]));
    }
    else
    {
      area = 0;
    }
    if (qdisk.nhit[n] > 0)
    {

      /* During photon transport ave_freq is actually the w x frequency */
      if (ichoice == 0)
        qdisk.ave_freq[n] /= qdisk.heat[n];

      qdisk.t_hit[n] = PLANCK * qdisk.ave_freq[n] / (BOLTZMANN * 3.832);
      //Basic conversion from freq to T
      qdisk.w[n] = qdisk.heat[n] / (4. * PI * STEFAN_BOLTZMANN * area * qdisk.t_hit[n] * qdisk.t_hit[n] * qdisk.t_hit[n] * qdisk.t_hit[n]);
      ratio[n] = 99.;
      if (qdisk.emit[n] > 0.0)
        ratio[n] = qdisk.heat[n] / qdisk.emit[n];
      theat = qdisk.heat[n] / area;
      theat = pow (theat / STEFAN_BOLTZMANN, 0.25);
      //theat is temperature if no internal energy production

    }
    else
    {
      qdisk.ave_freq[n] = 0.0;
      qdisk.t_hit[n] = 0.0;
      qdisk.w[n] = 0.0;
      ratio[n] = 0.0;
      theat = 0;
    }
    ttot = pow (qdisk.t[n], 4) + pow (theat, 4);
    ttot = pow (ttot, 0.25);
    fprintf (qptr,
             "%9.4e %9.4e %0.4e %8.3e %8.3e  %8.3e %10d %8.3e %10d %8.3e %8.3e %8.3e %10.3e %8.3e\n",
             qdisk.r[n], qdisk.v[n], zdisk (qdisk.r[n]), qdisk.t[n], qdisk.g[n],
             qdisk.emit[n], qdisk.nphot[n], qdisk.heat[n], qdisk.nhit[n], ratio[n], theat, qdisk.t_hit[n], qdisk.w[n], ttot);
  }

  fclose (qptr);
  return (0);
}





/**********************************************************/
/**
 * @brief      Read the temperature profile from a file
 *
 * @param [in] char *  tprofile   Name of the input file
 * @return     0
 *
 * @details
 *
 * Each line of the input file
 * a radius and a temperature in the first two columns.
 * An optional 3rd column can contain a gravity
 * for use with stellar atmospheres models
 *
 * Comment lines (and other lines) that can not
 * be parsed are ignored, but will be printed out
 * to make sure that nothing is amiss.
 *
 * The radius values shoule be in cm
 * The temperature in degrees Kelvin
 *
 * The minimum and maximum radius of the disk is set to the minimum
 * and maximum radius of the disk.  e 
 *
 * ###Notes###
 *
 * 210226 - In the distant past, the units were different, that is
 * the r values were in units of 10**11 cm and the temperature
 * values were in 1000 of degrees.  At one time, the temperature
 * profile did not have to be for the entire disk, and one used
 * mdot and a standard model to describe the disk.  This has been
 * changed so that if one reads in a non standard temperature
 * profile, it must include the entire disk..
 *
 **********************************************************/

int
read_non_standard_disk_profile (tprofile)
     char *tprofile;
{

  FILE *fopen (), *fptr;
  int n;
  float r, t, g;
  int one, two;

  char *line;
  size_t buffsize = LINELENGTH;

  if ((fptr = fopen (tprofile, "r")) == NULL)
  {
    Error ("read_non_standard_disk_profile: Could not open filename %s\n", tprofile);
    Exit (1);
  }
  line = (char *) malloc (buffsize * sizeof (char));
  blmod.n_blpts = 0;
  one = 0;
  two = 0;


  while (getline (&line, &buffsize, fptr) > 0)
  {
    n = sscanf (line, "%g %g  %g", &r, &t, &g);
    if (n >= 2)
    {
      blmod.r[blmod.n_blpts] = r;
      blmod.t[blmod.n_blpts] = t;
      if (n == 3)
      {
        blmod.g[blmod.n_blpts] = g;
        two += 1;
      }
      else
      {
        blmod.g[blmod.n_blpts] = -9999.;
        one += 1;
      }
      blmod.n_blpts += 1;
    }
    else
    {
      Error ("read_non_standard_disk_file: Could not convert a line in %s, OK if comment\n", tprofile);
    }

    if (blmod.n_blpts == NBLMODEL)
    {
      Error ("read_non_standard_disk_file: More than %d points in %s; increase NBLMODEL\n", NBLMODEL, tprofile);
      Exit (1);

    }
  }
  if (one == blmod.n_blpts)
  {
    blmod.n_params = 1;
  }
  else if (two == blmod.n_blpts)
  {
    blmod.n_params = 2;
  }
  else
  {
    Error ("read_non_standard_disk_file: Inconsistent input lines: one %d, two %d pts %d\n", one, two, blmod.n_blpts);
  }

  for (n = 0; n < blmod.n_blpts; n++)
  {
    if (blmod.n_params == 1)
    {
      blmod.n_params = 1;
      Log ("Disk: r %.3e t %.3e \n", blmod.r[n], blmod.t[n]);
    }
    else
    {
      Log ("Disk: r %.3e t %.3e g %.3e\n", blmod.r[n], blmod.t[n], blmod.g[n]);
    }
  }


  fclose (fptr);

  geo.disk_rad_min = 0;
  if (blmod.r[0] > geo.rstar)
  {
    geo.disk_rad_min = blmod.r[0];
    geo.disk_type = DISK_WITH_HOLE;
    Log ("The temperature profile begins at %e, which is larger than central object %e\n", geo.disk_rad_min, geo.rstar);
    Log ("Treating the disk as having a hole through which photon can pass\n");
  }
  if (blmod.r[0] < geo.rstar)
  {
    Error ("The temperature profile which was read in begins at %e, which is less than central object %e\n", blmod.r[0], geo.rstar);
    double frac;
    frac = (geo.rstar - blmod.r[0]) / geo.rstar;
    if (frac > 0.01 || blmod.r[1] < geo.rstar)
    {
      Error ("This does not appear to be a rounding error so exiting\n");
      Exit (0);
    }
    else
    {
      Error ("This appears to be a rounding error, so continuining\n");
      blmod.r[0] = geo.rstar;
      geo.disk_rad_min = geo.rstar;;
    }
  }

  geo.disk_rad_max = blmod.r[blmod.n_blpts - 1];
  geo.disk_mdot = 0;


  return (0);
}
