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



/**********************************************************/
/**
 * @brief      Whether photons should be scattered with polarization
 *
 * @return     TRUE if modes.polarization is set and this is a spectral cycle
 *
 * @details
 * The ionization cycles are left exactly as they were, including
 * the sequence of random numbers.
 *
 **********************************************************/

int
pol_active (void)
{
  return (modes.polarization && geo.ioniz_or_extract == CYCLE_EXTRACT);
}



/**********************************************************/
/**
 * @brief      Draw cos(Theta) from the Thomson (dipole) distribution
 *
 * @return     cos(Theta), distributed as (3/8)(1 + mu^2)
 *
 * @details
 * The CDF is inverted exactly: mu^3 + 3 mu - (8 xi - 4) = 0 has the
 * single real root
 *
 *     mu = cbrt(A + sqrt(A^2 + 1)) + cbrt(A - sqrt(A^2 + 1)),  A = 4 xi - 2
 *
 * This is also the marginal distribution of Theta for polarized light,
 * because the terms in phi average to zero.
 *
 **********************************************************/

double
pol_thomson_cos (void)
{
  double a, r;

  a = 4. * random_number (0.0, 1.0) - 2.;
  r = sqrt (a * a + 1.);
  return (cbrt (a + r) + cbrt (a - r));
}



/**********************************************************/
/**
 * @brief      The electron-scattering phase matrix applied to a Stokes vector
 *
 * @param [in] double  c   cos(Theta), Theta the scattering angle
 * @param [in] double  eps   nu_out/nu_in; 1 for Thomson scattering
 * @param [in] double  q_perp   Q/I of the incoming photon, referred to e_perp,
 * the normal to the scattering plane
 * @param [in] double  u_perp   U/I, referred to e_perp
 * @param [out] double  *q_out   Q/I of the scattered photon, referred to e_perp
 * @param [out] double  *u_out   U/I of the scattered photon, referred to e_perp
 * @return     D/P11, the ratio of the scattered intensity in this direction to
 * that for unpolarized light
 *
 * @details
 * With e_perp as the reference axis the (Fano) matrix has
 *
 *     P11 = eps + 1/eps - s^2,  P12 = P21 = s^2,  P22 = 1 + c^2,  P33 = 2c
 *
 * which for eps = 1 is the Thomson matrix, P11 = 1 + c^2.  The
 * scattered intensity is proportional to D = P11 + s^2 q_perp.
 *
 **********************************************************/

double
pol_es_matrix (double c, double eps, double q_perp, double u_perp, double *q_out, double *u_out)
{
  double s2, p11, d;

  s2 = 1. - c * c;
  p11 = eps + 1. / eps - s2;
  d = p11 + s2 * q_perp;

  *q_out = (s2 + (1. + c * c) * q_perp) / d;
  *u_out = 2. * c * u_perp / d;

  return (d / p11);
}



/**********************************************************/
/**
 * @brief      Choose the azimuth of an electron scatter and set the new
 * direction and polarization of the photon
 *
 * @param [in, out] PhotPtr  p   The photon, with its incoming direction and
 * polarization; on return its new direction and polarization
 * @param [in] double  c   cos(Theta), already drawn from the unpolarized
 * distribution (which is the exact marginal distribution of Theta)
 * @param [in] double  eps   nu_out/nu_in; 1 for Thomson scattering
 * @return     Always returns 0
 *
 * @details
 * phi, the angle of e_perp measured from the photon's e1 toward e2,
 * is drawn from the conditional distribution given Theta,
 *
 *     p(phi | Theta) ~ 1 + a (q cos 2phi + u sin 2phi),  a = s^2 / P11,
 *
 * by rejection with the envelope 1 + a sqrt(q^2 + u^2).  The new
 * direction is c n + s (n x e_perp), so e_perp is perpendicular to
 * both directions.  The weight does not change, because the
 * direction is drawn from the exact phase function.
 *
 **********************************************************/

int
pol_scatter_es (PhotPtr p, double c, double eps)
{
  double e1[3], e2[3], e_perp[3], ne[3], n_new[3];
  double s, s2, a, pmax, phi, f;
  double q_perp, u_perp, q_out, u_out;
  int i;

  if (c > 1.)
    c = 1.;
  else if (c < -1.)
    c = -1.;

  s2 = 1. - c * c;
  s = sqrt (s2);
  a = s2 / (eps + 1. / eps - s2);
  pmax = 1. + a * sqrt (p->q * p->q + p->u * p->u);

  do
  {
    phi = random_number (0.0, 2. * M_PI);
    f = 1. + a * (p->q * cos (2. * phi) + p->u * sin (2. * phi));
  }
  while (random_number (0.0, 1.0) * pmax > f);

  pol_basis (p->lmn, e1, e2);
  for (i = 0; i < 3; i++)
  {
    e_perp[i] = cos (phi) * e1[i] + sin (phi) * e2[i];
  }
  cross (p->lmn, e_perp, ne);
  for (i = 0; i < 3; i++)
  {
    n_new[i] = c * p->lmn[i] + s * ne[i];
  }
  renorm (n_new, 1.0);

  q_perp = p->q;
  u_perp = p->u;
  pol_rotate (&q_perp, &u_perp, phi);
  pol_es_matrix (c, eps, q_perp, u_perp, &q_out, &u_out);

  stuff_v (n_new, p->lmn);
  pol_from_axis (p, e_perp, q_out, u_out);

  return (0);
}



/**********************************************************/
/**
 * @brief      Polarization of an electron-scattered photon peeled off toward
 * an observer
 *
 * @param [in] PhotPtr  p_in   The photon before the scatter
 * @param [in, out] PhotPtr  pp   The photon to be extracted, already pointed
 * toward the observer and reweighted for unpolarized scattering (and, for
 * Compton scattering, with its frequency changed)
 * @return     Always returns 0
 *
 * @details
 * The weight of pp is multiplied by D/P11, and its q and u are set
 * from the phase matrix.  For unpolarized incoming light D/P11 = 1, so
 * the weight is unchanged, but the scattered light is polarized.
 *
 * If the observer is (almost) exactly forward or backward, the
 * scattering plane is undefined, and e1 of the incoming photon is used
 * for e_perp.
 *
 **********************************************************/

int
pol_peel_es (PhotPtr p_in, PhotPtr pp)
{
  double e_perp[3], e2[3];
  double c, eps;
  double q_perp, u_perp, q_out, u_out;

  if (modes.polarization == FALSE)
  {
    return (0);
  }

  c = dot (p_in->lmn, pp->lmn);
  if (c > 1.)
    c = 1.;
  else if (c < -1.)
    c = -1.;
  eps = pp->freq / p_in->freq;

  cross (p_in->lmn, pp->lmn, e_perp);
  if (length (e_perp) < POL_EPS)
  {
    pol_basis (p_in->lmn, e_perp, e2);
  }
  renorm (e_perp, 1.0);

  q_perp = p_in->q;
  u_perp = p_in->u;
  pol_rotate (&q_perp, &u_perp, pol_axis_angle (p_in->lmn, e_perp));
  pp->w *= pol_es_matrix (c, eps, q_perp, u_perp, &q_out, &u_out);

  pol_from_axis (pp, e_perp, q_out, u_out);

  return (0);
}



/**********************************************************/
/**
 * @brief      Re-express the polarization of a photon whose direction has been
 * changed by a Lorentz transformation
 *
 * @param [in] double  n_in[]   The direction before the transformation
 * @param [in] double  q   Q/I before the transformation, against the default
 * axis for n_in
 * @param [in] double  u   U/I before the transformation
 * @param [in] double  v[]   The velocity (cm/s) of the frame being transformed into
 * @param [in, out] PhotPtr  p_out   The transformed photon, whose lmn has
 * already been changed; its q and u are set
 * @return     Always returns 0
 *
 * @details
 * Aberration is conformal, so q and u are unchanged against a reference
 * axis that is carried through the transformation like the electric field,
 *
 *     e' = gamma (e + beta x (n x e)) - gamma^2/(gamma+1) beta (beta . e),
 *
 * but that axis is not the default axis for the new direction, so q and
 * u are rotated onto it.  e' is made exactly perpendicular to the new
 * direction before it is used.
 *
 **********************************************************/

int
pol_lorentz (double n_in[], double q, double u, double v[], PhotPtr p_out)
{
  double e1[3], e2[3], b[3], bxb[3], beta[3], e[3];
  double gamma, bdote;
  int i;

  p_out->q = q;
  p_out->u = u;
  if (modes.polarization == FALSE || (q == 0 && u == 0))
  {
    return (0);
  }

  for (i = 0; i < 3; i++)
  {
    beta[i] = v[i] / VLIGHT;
  }
  gamma = 1. / sqrt (1. - dot (beta, beta));

  pol_basis (n_in, e1, e2);
  cross (n_in, e1, b);
  cross (beta, b, bxb);
  bdote = dot (beta, e1);
  for (i = 0; i < 3; i++)
  {
    e[i] = gamma * (e1[i] + bxb[i]) - gamma * gamma / (gamma + 1.) * beta[i] * bdote;
  }

  /* Make e exactly perpendicular to the new direction */
  bdote = dot (e, p_out->lmn);
  for (i = 0; i < 3; i++)
  {
    e[i] -= bdote * p_out->lmn[i];
  }
  renorm (e, 1.0);

  pol_from_axis (p_out, e, q, u);

  return (0);
}
