/***********************************************************/
/** @file  polarization.c
 * @author ksl
 * @date   October, 2026
 *
 * @brief
 * Routines for following the linear polarization of photons
 * and for constructing Q and U spectra for the observers
 *
 * ###Notes###
 *
 * Polarization is followed only when modes.polarization is
 * TRUE (@Spectrum.polarization), and only matters in the
 * spectral cycles.
 *
 * Each photon carries q = Q/I and u = U/I, where I is the
 * photon weight w.  The reference axis is not stored.  It is
 * always the default axis for the photon's direction n,
 *
 *     e1 = normalize(z x n),  or normalize(x x n) near the poles,
 *     e2 = n x e1,
 *
 * so that (e1, e2, n) is right handed.  q > 0 means the electric
 * vector is along e1, and u > 0 means along (e1 + e2)/sqrt(2).
 * Away from the poles e2 is the projected +z axis, i.e. North
 * on the sky.
 *
 * Observer spectra are referred to North on the sky, the
 * projected +z axis, so Q > 0 is an electric vector parallel to
 * the projected symmetry axis (IAU convention).
 *
 * The conventions are described in more detail in
 * Sirocco-2026-Polarization-Conventions, and were prototyped in
 * py_progs/polarization.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

#define POL_EPS  1e-6           /* |z x n| below which the x axis is used to define e1 */


/**********************************************************/
/**
 * @brief      Construct the default polarization basis for a direction
 *
 * @param [in] double  lmn[]   A unit vector, the direction of a photon
 * @param [out] double  e1[]   The reference axis, q > 0 along it
 * @param [out] double  e2[]   n x e1
 * @return     Always returns 0
 *
 * @details
 * e1 = normalize(z x n), or normalize(x x n) if n is within
 * about 1e-6 of the z axis.
 *
 **********************************************************/

int
pol_basis (double lmn[], double e1[], double e2[])
{
  double zhat[3] = { 0, 0, 1 };
  double xhat[3] = { 1, 0, 0 };

  cross (zhat, lmn, e1);
  if (length (e1) < POL_EPS)
  {
    cross (xhat, lmn, e1);
  }
  renorm (e1, 1.0);
  cross (lmn, e1, e2);

  return (0);
}



/**********************************************************/
/**
 * @brief      Rotate the reference axis of a Stokes vector
 *
 * @param [in, out] double  *q   Q/I
 * @param [in, out] double  *u   U/I
 * @param [in] double  alpha   The angle, measured from e1 toward e2,
 * of the new reference axis
 * @return     Always returns 0
 *
 **********************************************************/

int
pol_rotate (double *q, double *u, double alpha)
{
  double c, s, qq;

  c = cos (2. * alpha);
  s = sin (2. * alpha);

  qq = (*q) * c + (*u) * s;
  *u = -(*q) * s + (*u) * c;
  *q = qq;

  return (0);
}



/**********************************************************/
/**
 * @brief      The angle of an axis t, perpendicular to lmn, measured from
 * the default axis e1 toward e2
 *
 * @param [in] double  lmn[]   The direction of the photon
 * @param [in] double  t[]   A unit vector perpendicular to lmn
 * @return     The angle in radians
 *
 **********************************************************/

double
pol_axis_angle (double lmn[], double t[])
{
  double e1[3], e2[3];

  pol_basis (lmn, e1, e2);
  return (atan2 (dot (t, e2), dot (t, e1)));
}



/**********************************************************/
/**
 * @brief      Express the polarization of a photon against another axis
 *
 * @param [in] PhotPtr  p   The photon, whose q and u refer to the default axis
 * @param [in] double  t[]   A unit vector perpendicular to p->lmn
 * @param [out] double  *qt   Q/I against t
 * @param [out] double  *ut   U/I against t
 * @return     Always returns 0
 *
 **********************************************************/

int
pol_to_axis (PhotPtr p, double t[], double *qt, double *ut)
{
  *qt = p->q;
  *ut = p->u;
  if (*qt != 0 || *ut != 0)
  {
    pol_rotate (qt, ut, pol_axis_angle (p->lmn, t));
  }

  return (0);
}



/**********************************************************/
/**
 * @brief      Set the polarization of a photon from values measured against
 * another axis
 *
 * @param [in, out] PhotPtr  p   The photon; p->q and p->u are set
 * against the default axis for p->lmn
 * @param [in] double  t[]   A unit vector perpendicular to p->lmn
 * @param [in] double  qt   Q/I against t
 * @param [in] double  ut   U/I against t
 * @return     Always returns 0
 *
 **********************************************************/

int
pol_from_axis (PhotPtr p, double t[], double qt, double ut)
{
  p->q = qt;
  p->u = ut;
  if (qt != 0 || ut != 0)
  {
    pol_rotate (&p->q, &p->u, -pol_axis_angle (p->lmn, t));
  }

  return (0);
}



/**********************************************************/
/**
 * @brief      North on the sky for a direction
 *
 * @param [in] double  lmn[]   The direction toward the observer
 * @param [out] double  north[]   The projected +z axis, a unit vector
 * perpendicular to lmn
 * @return     Always returns 0
 *
 * @details
 * Away from the poles this is e2 of the default basis.  At the
 * pole, where the projected z axis is undefined, e2 of the default
 * basis is used, so that the result is still well defined.  For
 * the observers, pol_spectra_init uses the inclination and phase
 * instead, which is well defined at the pole too.
 *
 **********************************************************/

int
pol_north (double lmn[], double north[])
{
  double e1[3];

  pol_basis (lmn, e1, north);

  return (0);
}



/**********************************************************/
/**
 * @brief      Allocate (the first time) and zero the Q and U spectra,
 * and set North on the sky for each observer
 *
 * @param [in] int  nangle   The number of observers
 * @param [in] double  angle[]   The inclinations of the observers (deg)
 * @param [in] double  phase[]   The orbital phases of the observers
 * @return     Always returns 0
 *
 * @details
 * Does nothing unless modes.polarization is TRUE.  It is
 * called from spectrum_init.
 *
 * North is e_N = -dn/di, the projected +z axis, which is defined
 * at i = 0 as well.  The azimuth follows the convention used
 * for xxspec[].lmn in spectrum_init.
 *
 **********************************************************/

int
pol_spectra_init (int nangle, double angle[], double phase[])
{
  int m, i;
  double inc, psi;

  if (modes.polarization == FALSE)
  {
    return (0);
  }

  if (xxpol == NULL)
  {
    npolspec = nangle;
    if ((xxpol = calloc (npolspec, sizeof (spec_pol_dummy))) == NULL)
    {
      Error ("pol_spectra_init: Could not allocate memory for %d polarization spectra\n", npolspec);
      Exit (EXIT_FAILURE);
    }
    pol_spectra_allocate ();
  }

  for (m = 0; m < npolspec; m++)
  {
    inc = angle[m] / RADIAN;
    psi = -phase[m] * 360. / RADIAN;
    xxpol[m].north[0] = -cos (inc) * cos (psi);
    xxpol[m].north[1] = -cos (inc) * sin (psi);
    xxpol[m].north[2] = sin (inc);

    for (i = 0; i < NWAVE_MAX; i++)
    {
      xxpol[m].q[i] = xxpol[m].u[i] = 0;
      xxpol[m].lq[i] = xxpol[m].lu[i] = 0;
      xxpol[m].q_wind[i] = xxpol[m].u_wind[i] = 0;
      xxpol[m].lq_wind[i] = xxpol[m].lu_wind[i] = 0;
    }
  }

  return (0);
}



/**********************************************************/
/**
 * @brief      Allocate the arrays of xxpol
 *
 * @return     Always returns 0
 *
 * @details
 * xxpol itself, with npolspec elements, must already exist.
 *
 **********************************************************/

int
pol_spectra_allocate (void)
{
  int m;

  for (m = 0; m < npolspec; m++)
  {
    xxpol[m].q = calloc (NWAVE_MAX, sizeof (double));
    xxpol[m].u = calloc (NWAVE_MAX, sizeof (double));
    xxpol[m].lq = calloc (NWAVE_MAX, sizeof (double));
    xxpol[m].lu = calloc (NWAVE_MAX, sizeof (double));
    xxpol[m].q_wind = calloc (NWAVE_MAX, sizeof (double));
    xxpol[m].u_wind = calloc (NWAVE_MAX, sizeof (double));
    xxpol[m].lq_wind = calloc (NWAVE_MAX, sizeof (double));
    xxpol[m].lu_wind = calloc (NWAVE_MAX, sizeof (double));
    if (xxpol[m].q == NULL || xxpol[m].u == NULL || xxpol[m].lq == NULL || xxpol[m].lu == NULL
        || xxpol[m].q_wind == NULL || xxpol[m].u_wind == NULL || xxpol[m].lq_wind == NULL || xxpol[m].lu_wind == NULL)
    {
      Error ("pol_spectra_allocate: Unable to allocate memory for xxpol[%d] with %d bins\n", m, NWAVE_MAX);
      Exit (EXIT_FAILURE);
    }
  }

  return (0);
}



/**********************************************************/
/**
 * @brief      Add the Q and U of a photon to the spectra of an observer
 *
 * @param [in] PhotPtr  p   The photon
 * @param [in] double  north[]   North on the sky, perpendicular to p->lmn
 * @param [in] int  nspec   The spectrum, xxspec[nspec], to which the photon
 * contributes; nspec >= MSPEC
 * @param [in] double  w   The weight with which the photon is added to I
 * @param [in] int  k   The bin in the linear spectrum
 * @param [in] int  k1   The bin in the logarithmic spectrum
 * @param [in] int  iwind   TRUE if the photon also goes into the wind spectra
 * @return     Always returns 0
 *
 * @details
 * Does nothing unless modes.polarization is TRUE.
 *
 **********************************************************/

int
pol_tally (PhotPtr p, double north[], int nspec, double w, int k, int k1, int iwind)
{
  int m;
  double qs, us;

  if (modes.polarization == FALSE || (p->q == 0 && p->u == 0))
  {
    return (0);
  }

  m = nspec - MSPEC;
  if (m < 0 || m >= npolspec)
  {
    Error ("pol_tally: spectrum %d has no polarization spectrum\n", nspec);
    return (0);
  }

  pol_to_axis (p, north, &qs, &us);

  xxpol[m].q[k] += w * qs;
  xxpol[m].u[k] += w * us;
  xxpol[m].lq[k1] += w * qs;
  xxpol[m].lu[k1] += w * us;

  if (iwind)
  {
    xxpol[m].q_wind[k] += w * qs;
    xxpol[m].u_wind[k] += w * us;
    xxpol[m].lq_wind[k1] += w * qs;
    xxpol[m].lu_wind[k1] += w * us;
  }

  return (0);
}



/**********************************************************/
/**
 * @brief      Multiply all of the Q and U spectra by a factor
 *
 * @param [in] double  factor   The factor
 * @return     Always returns 0
 *
 * @details
 * Used when a run with spectral cycles is restarted with a
 * different number of cycles (see spectrum_restart_renormalise).
 *
 * ### Notes ###
 * spectrum_restart_renormalise does not renormalize the wind
 * (f_wind, lf_wind) spectra, so the wind Q and U spectra are
 * not renormalized here either, to keep Q/I consistent.
 *
 **********************************************************/

int
pol_spectra_renormalise (double factor)
{
  int m, i;

  if (modes.polarization == FALSE)
  {
    return (0);
  }

  for (m = 0; m < npolspec; m++)
  {
    for (i = 0; i < NWAVE_EXTRACT; i++)
    {
      xxpol[m].q[i] *= factor;
      xxpol[m].u[i] *= factor;
      xxpol[m].lq[i] *= factor;
      xxpol[m].lu[i] *= factor;
    }
  }

  return (0);
}



/**********************************************************/
/**
 * @brief      Write the Q and U spectra to an open spec_save file
 *
 * @param [in] FILE  *fptr   The file
 * @return     The number of arrays written
 *
 **********************************************************/

int
pol_spectra_write (FILE *fptr)
{
  int m;
  int count = 0;

  for (m = 0; m < npolspec; m++)
  {
    count += (int) fwrite (xxpol[m].north, sizeof (double), 3, fptr);
    count += (int) fwrite (xxpol[m].q, sizeof (double), NWAVE_MAX, fptr);
    count += (int) fwrite (xxpol[m].u, sizeof (double), NWAVE_MAX, fptr);
    count += (int) fwrite (xxpol[m].lq, sizeof (double), NWAVE_MAX, fptr);
    count += (int) fwrite (xxpol[m].lu, sizeof (double), NWAVE_MAX, fptr);
    count += (int) fwrite (xxpol[m].q_wind, sizeof (double), NWAVE_MAX, fptr);
    count += (int) fwrite (xxpol[m].u_wind, sizeof (double), NWAVE_MAX, fptr);
    count += (int) fwrite (xxpol[m].lq_wind, sizeof (double), NWAVE_MAX, fptr);
    count += (int) fwrite (xxpol[m].lu_wind, sizeof (double), NWAVE_MAX, fptr);
  }

  return (count);
}



/**********************************************************/
/**
 * @brief      Allocate and read the Q and U spectra from an open spec_save file
 *
 * @param [in] FILE  *fptr   The file
 * @param [in] int  nangle   The number of observers
 * @return     The number of arrays read
 *
 * @details
 * Called by spec_read when the spec_save file has polarization
 * spectra.  It sets modes.polarization to TRUE, so that a
 * restarted run continues to follow the polarization.
 *
 **********************************************************/

int
pol_spectra_read (FILE *fptr, int nangle)
{
  int m;
  int count = 0;

  modes.polarization = TRUE;
  npolspec = nangle;
  if ((xxpol = calloc (npolspec, sizeof (spec_pol_dummy))) == NULL)
  {
    Error ("pol_spectra_read: Could not allocate memory for %d polarization spectra\n", npolspec);
    Exit (EXIT_FAILURE);
  }
  pol_spectra_allocate ();

  for (m = 0; m < npolspec; m++)
  {
    count += (int) fread (xxpol[m].north, sizeof (double), 3, fptr);
    count += (int) fread (xxpol[m].q, sizeof (double), NWAVE_MAX, fptr);
    count += (int) fread (xxpol[m].u, sizeof (double), NWAVE_MAX, fptr);
    count += (int) fread (xxpol[m].lq, sizeof (double), NWAVE_MAX, fptr);
    count += (int) fread (xxpol[m].lu, sizeof (double), NWAVE_MAX, fptr);
    count += (int) fread (xxpol[m].q_wind, sizeof (double), NWAVE_MAX, fptr);
    count += (int) fread (xxpol[m].u_wind, sizeof (double), NWAVE_MAX, fptr);
    count += (int) fread (xxpol[m].lq_wind, sizeof (double), NWAVE_MAX, fptr);
    count += (int) fread (xxpol[m].lu_wind, sizeof (double), NWAVE_MAX, fptr);
  }

  return (count);
}



/**********************************************************/
/**
 * @brief      Sum and normalize the Q and U spectra across MPI ranks
 *
 * @return     Always returns 0
 *
 * @details
 * The counterpart of normalize_spectra_across_ranks, which calls
 * it.  Does nothing unless modes.polarization is TRUE and this is
 * a spectral cycle.
 *
 **********************************************************/

int
pol_spectra_reduce (void)
{
#ifdef MPI_ON
  int i, m, j;
  int nbuf;
  double *buffer;
  double *arr[8];

  if (modes.polarization == FALSE || geo.ioniz_or_extract != CYCLE_EXTRACT)
  {
    return (0);
  }

  nbuf = 8 * npolspec * NWAVE_MAX;
  if ((buffer = calloc (nbuf, sizeof (double))) == NULL)
  {
    Error ("pol_spectra_reduce: Unable to allocate the communication buffer\n");
    Exit (EXIT_FAILURE);
  }

  for (m = 0; m < npolspec; m++)
  {
    arr[0] = xxpol[m].q;
    arr[1] = xxpol[m].u;
    arr[2] = xxpol[m].lq;
    arr[3] = xxpol[m].lu;
    arr[4] = xxpol[m].q_wind;
    arr[5] = xxpol[m].u_wind;
    arr[6] = xxpol[m].lq_wind;
    arr[7] = xxpol[m].lu_wind;
    for (j = 0; j < 8; j++)
    {
      for (i = 0; i < NWAVE_MAX; i++)
      {
        buffer[(m * 8 + j) * NWAVE_MAX + i] = arr[j][i] / np_mpi_global;
      }
    }
  }

  MPI_Allreduce (MPI_IN_PLACE, buffer, nbuf, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (m = 0; m < npolspec; m++)
  {
    arr[0] = xxpol[m].q;
    arr[1] = xxpol[m].u;
    arr[2] = xxpol[m].lq;
    arr[3] = xxpol[m].lu;
    arr[4] = xxpol[m].q_wind;
    arr[5] = xxpol[m].u_wind;
    arr[6] = xxpol[m].lq_wind;
    arr[7] = xxpol[m].lu_wind;
    for (j = 0; j < 8; j++)
    {
      for (i = 0; i < NWAVE_MAX; i++)
      {
        arr[j][i] = buffer[(m * 8 + j) * NWAVE_MAX + i];
      }
    }
  }

  free (buffer);
#endif
  return (0);
}
