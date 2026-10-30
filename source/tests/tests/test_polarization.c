/** ********************************************************************************************************************
 *
 *  @file test_polarization.c
 *  @author ksl
 *  @date October 2026
 *
 *  @brief Unit tests for the geometry of the polarization routines in polarization.c
 *
 * ****************************************************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <CUnit/CUnit.h>

#include "../../atomic.h"
#include "../../sirocco.h"
#include "../assert.h"

#define POL_TOL 1e-12

static double test_dirs[][3] = {
  {0.0, 0.0, 1.0},
  {0.0, 0.0, -1.0},
  {1.0, 0.0, 0.0},
  {0.3, -0.4, 0.866025403784},
  {-0.5, 0.5, -0.707106781187},
};

#define NDIRS (int) (sizeof (test_dirs) / sizeof (test_dirs[0]))

/** *******************************************************************************************************************
 *
 * @brief The default basis (e1, e2, n) is orthonormal and right handed, including along the z axis
 *
 * ****************************************************************************************************************** */

void
test_pol_basis (void)
{
  int i;
  double n[3], e1[3], e2[3], c[3];

  for (i = 0; i < NDIRS; i++)
  {
    stuff_v (test_dirs[i], n);
    renorm (n, 1.0);
    pol_basis (n, e1, e2);

    CU_ASSERT_DOUBLE_EQUAL (length (e1), 1.0, POL_TOL);
    CU_ASSERT_DOUBLE_EQUAL (length (e2), 1.0, POL_TOL);
    CU_ASSERT_DOUBLE_EQUAL (dot (e1, n), 0.0, POL_TOL);
    CU_ASSERT_DOUBLE_EQUAL (dot (e2, n), 0.0, POL_TOL);
    CU_ASSERT_DOUBLE_EQUAL (dot (e1, e2), 0.0, POL_TOL);

    cross (e1, e2, c);
    CU_ASSERT_DOUBLE_EQUAL (dot (c, n), 1.0, POL_TOL);
  }
}

/** *******************************************************************************************************************
 *
 * @brief Rotations of the reference axis: by pi is the identity, by pi/4 takes Q to -U, and +alpha then -alpha
 * returns the starting state
 *
 * ****************************************************************************************************************** */

void
test_pol_rotate (void)
{
  double q, u;

  q = 0.3;
  u = -0.2;
  pol_rotate (&q, &u, M_PI);
  CU_ASSERT_DOUBLE_EQUAL (q, 0.3, POL_TOL);
  CU_ASSERT_DOUBLE_EQUAL (u, -0.2, POL_TOL);

  q = 1.0;
  u = 0.0;
  pol_rotate (&q, &u, M_PI / 4.);
  CU_ASSERT_DOUBLE_EQUAL (q, 0.0, POL_TOL);
  CU_ASSERT_DOUBLE_EQUAL (u, -1.0, POL_TOL);

  q = 0.3;
  u = -0.2;
  pol_rotate (&q, &u, 0.7);
  pol_rotate (&q, &u, -0.7);
  CU_ASSERT_DOUBLE_EQUAL (q, 0.3, POL_TOL);
  CU_ASSERT_DOUBLE_EQUAL (u, -0.2, POL_TOL);
}

/** *******************************************************************************************************************
 *
 * @brief pol_from_axis and pol_to_axis are inverses, and the degree of polarization does not change
 *
 * ****************************************************************************************************************** */

void
test_pol_axis_round_trip (void)
{
  int i;
  double e1[3], e2[3], t[3];
  double qt, ut;
  struct photon p;

  init_dummy_phot (&p);

  for (i = 0; i < NDIRS; i++)
  {
    stuff_v (test_dirs[i], p.lmn);
    renorm (p.lmn, 1.0);

    /* An axis at 40 degrees from e1 */
    pol_basis (p.lmn, e1, e2);
    t[0] = cos (0.698) * e1[0] + sin (0.698) * e2[0];
    t[1] = cos (0.698) * e1[1] + sin (0.698) * e2[1];
    t[2] = cos (0.698) * e1[2] + sin (0.698) * e2[2];

    pol_from_axis (&p, t, 0.25, 0.1);
    CU_ASSERT_DOUBLE_EQUAL (p.q * p.q + p.u * p.u, 0.25 * 0.25 + 0.1 * 0.1, POL_TOL);

    pol_to_axis (&p, t, &qt, &ut);
    CU_ASSERT_DOUBLE_EQUAL (qt, 0.25, POL_TOL);
    CU_ASSERT_DOUBLE_EQUAL (ut, 0.1, POL_TOL);

    /* Against e1 itself nothing changes */
    pol_to_axis (&p, e1, &qt, &ut);
    CU_ASSERT_DOUBLE_EQUAL (qt, p.q, POL_TOL);
    CU_ASSERT_DOUBLE_EQUAL (ut, p.u, POL_TOL);
  }
}

/** *******************************************************************************************************************
 *
 * @brief North on the sky for the observers
 *
 * @details
 *
 * North must be a unit vector perpendicular to the observer direction (constructed as in spectrum_init) with a
 * positive z component, and away from the pole it must equal pol_north for that direction.  At i = 0 it must still
 * be defined.  A photon whose electric vector is along the default e1 (in the plane of the disk) must give Q < 0
 * in the sky frame, and one with its electric vector along North must give Q > 0.
 *
 * ****************************************************************************************************************** */

void
test_pol_sky_north (void)
{
  int m;
  double angle[3] = { 0.0, 60.0, 85.0 };
  double phase[3] = { 0.0, 0.25, 0.6 };
  double n[3], north[3];
  double qs, us;
  struct photon p;

  modes.polarization = TRUE;
  NWAVE_MAX = 10;
  pol_spectra_init (3, angle, phase);

  for (m = 0; m < 3; m++)
  {
    n[0] = sin (angle[m] / RADIAN) * cos (-phase[m] * 360. / RADIAN);
    n[1] = sin (angle[m] / RADIAN) * sin (-phase[m] * 360. / RADIAN);
    n[2] = cos (angle[m] / RADIAN);

    CU_ASSERT_DOUBLE_EQUAL (length (xxpol[m].north), 1.0, POL_TOL);
    CU_ASSERT_DOUBLE_EQUAL (dot (xxpol[m].north, n), 0.0, POL_TOL);
    CU_ASSERT (xxpol[m].north[2] >= 0.0);

    if (angle[m] > 0)
    {
      pol_north (n, north);
      CU_ASSERT_DOUBLE_EQUAL (dot (north, xxpol[m].north), 1.0, POL_TOL);
    }

    init_dummy_phot (&p);
    stuff_v (n, p.lmn);
    p.q = 0.5;                  /* along the default e1 */
    p.u = 0.0;
    pol_to_axis (&p, xxpol[m].north, &qs, &us);
    CU_ASSERT (qs < 0);

    pol_from_axis (&p, xxpol[m].north, 0.5, 0.0);
    pol_to_axis (&p, xxpol[m].north, &qs, &us);
    CU_ASSERT_DOUBLE_EQUAL (qs, 0.5, POL_TOL);
    CU_ASSERT_DOUBLE_EQUAL (us, 0.0, POL_TOL);
  }

  modes.polarization = FALSE;
}

/** *******************************************************************************************************************
 *
 * @brief The electron-scattering matrix: unpolarized light scattered through 90 degrees is fully polarized
 * perpendicular to the scattering plane, forward scattering keeps the polarization, and Klein-Nishina reduces the
 * polarization of the scattered light
 *
 * ****************************************************************************************************************** */

void
test_pol_es_matrix (void)
{
  double q, u, r;

  r = pol_es_matrix (0.0, 1.0, 0.0, 0.0, &q, &u);
  CU_ASSERT_DOUBLE_EQUAL (r, 1.0, POL_TOL);
  CU_ASSERT_DOUBLE_EQUAL (q, 1.0, POL_TOL);
  CU_ASSERT_DOUBLE_EQUAL (u, 0.0, POL_TOL);

  r = pol_es_matrix (1.0, 1.0, 0.3, -0.2, &q, &u);
  CU_ASSERT_DOUBLE_EQUAL (r, 1.0, POL_TOL);
  CU_ASSERT_DOUBLE_EQUAL (q, 0.3, POL_TOL);
  CU_ASSERT_DOUBLE_EQUAL (u, -0.2, POL_TOL);

  /* Light polarized in the scattering plane does not scatter through 90 degrees */
  r = pol_es_matrix (0.0, 1.0, -1.0, 0.0, &q, &u);
  CU_ASSERT_DOUBLE_EQUAL (r, 0.0, POL_TOL);

  /* q = s^2 / (1 + c^2) for Thomson; smaller for eps < 1 */
  pol_es_matrix (0.5, 1.0, 0.0, 0.0, &q, &u);
  CU_ASSERT_DOUBLE_EQUAL (q, 0.75 / 1.25, POL_TOL);
  pol_es_matrix (0.5, 0.5, 0.0, 0.0, &q, &u);
  CU_ASSERT_DOUBLE_EQUAL (q, 0.75 / (0.5 + 2.0 - 0.75), POL_TOL);
}

/** *******************************************************************************************************************
 *
 * @brief Peel-off: an unpolarized photon peeled at 90 degrees keeps its weight and is fully polarized, with the
 * electric vector perpendicular to the scattering plane
 *
 * ****************************************************************************************************************** */

void
test_pol_peel_es (void)
{
  struct photon p_in, pp;
  double t[3] = { 0, 1, 0 };
  double qt, ut;

  modes.polarization = TRUE;
  init_dummy_phot (&p_in);
  p_in.lmn[0] = 1;
  p_in.lmn[1] = 0;
  p_in.lmn[2] = 0;
  p_in.freq = 1e15;
  p_in.w = 1.0;

  stuff_phot (&p_in, &pp);
  pp.lmn[0] = 0;
  pp.lmn[1] = 0;
  pp.lmn[2] = 1;

  pol_peel_es (&p_in, &pp);
  CU_ASSERT_DOUBLE_EQUAL (pp.w, 1.0, POL_TOL);

  /* The scattering plane is the xz plane, so the electric vector is along y */
  pol_to_axis (&pp, t, &qt, &ut);
  CU_ASSERT_DOUBLE_EQUAL (qt, 1.0, POL_TOL);
  CU_ASSERT_DOUBLE_EQUAL (ut, 0.0, POL_TOL);
  modes.polarization = FALSE;
}

/** *******************************************************************************************************************
 *
 * @brief Lorentz transformation: the polarization is carried through a boost and its inverse unchanged, and a
 * pure boost along the photon direction does not change it
 *
 * ****************************************************************************************************************** */

void
test_pol_lorentz (void)
{
  struct photon p;
  double v[3] = { 0.2 * VLIGHT, -0.3 * VLIGHT, 0.1 * VLIGHT };
  double vback[3];
  double q0, u0;
  double e1[3], e2[3], qa, ua, qb, ub;

  int rel_mode_save;

  rel_mode_save = rel_mode;
  rel_mode = REL_MODE_FULL;     /* REL_MODE_LINEAR does not change the direction */
  modes.polarization = TRUE;
  init_dummy_phot (&p);
  p.lmn[0] = 0.6;
  p.lmn[1] = 0.0;
  p.lmn[2] = 0.8;
  p.freq = 1e15;
  p.w = 1.0;
  p.q = q0 = 0.3;
  p.u = u0 = 0.15;

  lorentz_transform (&p, &p, v);
  CU_ASSERT (fabs (p.lmn[0] - 0.6) > 0.01);     /* the boost really aberrates the photon */
  CU_ASSERT (fabs (p.q - q0) > 1e-3);   /* and the default axis changes with it */
  CU_ASSERT_DOUBLE_EQUAL (p.q * p.q + p.u * p.u, q0 * q0 + u0 * u0, POL_TOL);

  rescale (v, -1, vback);
  lorentz_transform (&p, &p, vback);
  CU_ASSERT_DOUBLE_EQUAL (p.lmn[0], 0.6, 1e-10);
  CU_ASSERT_DOUBLE_EQUAL (p.lmn[2], 0.8, 1e-10);
  CU_ASSERT_DOUBLE_EQUAL (p.q, q0, 1e-10);
  CU_ASSERT_DOUBLE_EQUAL (p.u, u0, 1e-10);

  /* A boost along the photon direction changes neither the direction nor the axes */
  p.lmn[0] = 0.6;
  p.lmn[1] = 0.0;
  p.lmn[2] = 0.8;
  p.q = q0;
  p.u = u0;
  rescale (p.lmn, 0.3 * VLIGHT, v);
  pol_basis (p.lmn, e1, e2);
  pol_to_axis (&p, e1, &qa, &ua);
  lorentz_transform (&p, &p, v);
  pol_to_axis (&p, e1, &qb, &ub);
  CU_ASSERT_DOUBLE_EQUAL (qa, qb, 1e-10);
  CU_ASSERT_DOUBLE_EQUAL (ua, ub, 1e-10);

  modes.polarization = FALSE;
  rel_mode = rel_mode_save;
}

/** *******************************************************************************************************************
 *
 * @brief Lorentz transformation, analytic check: for a pure boost the polarization angle measured from the normal
 * m = n x v of the plane containing the photon direction and the velocity is unchanged
 *
 * @details
 * The aberrated direction stays in the (n, v) plane, so m is perpendicular to it before and after.  For E parallel
 * to m, E' = gamma (1 - beta.n) E, and for E in the plane, E' has the same magnitude and lies in the plane, so the
 * map from E to E' is a rotation plus a uniform scaling and the angle from m is preserved.  In REL_MODE_LINEAR the
 * direction is not aberrated and q, u must be unchanged.
 *
 * ****************************************************************************************************************** */

void
test_pol_lorentz_plane (void)
{
  struct photon p;
  double vels[][3] = {
    {0.3, 0.0, 0.0},
    {0.0, -0.5, 0.2},
    {-0.4, 0.3, -0.6},
    {0.05, 0.02, -0.01},
  };
  double pols[][2] = {
    {0.5, 0.0},
    {0.0, -0.7},
    {0.3, 0.4},
  };
  double n[3], v[3], m[3], qa, ua, qb, ub;
  int i, j, k, rel_mode_save;

  rel_mode_save = rel_mode;
  modes.polarization = TRUE;
  rel_mode = REL_MODE_FULL;

  for (i = 0; i < NDIRS; i++)
  {
    for (j = 0; j < (int) (sizeof (vels) / sizeof (vels[0])); j++)
    {
      stuff_v (test_dirs[i], n);
      renorm (n, 1.0);
      rescale (vels[j], VLIGHT, v);
      cross (n, v, m);
      if (length (m) < 1e-3 * length (v))
        continue;
      renorm (m, 1.0);

      for (k = 0; k < (int) (sizeof (pols) / sizeof (pols[0])); k++)
      {
        init_dummy_phot (&p);
        stuff_v (n, p.lmn);
        p.freq = 1e15;
        p.w = 1.0;
        p.q = pols[k][0];
        p.u = pols[k][1];

        pol_to_axis (&p, m, &qa, &ua);
        lorentz_transform (&p, &p, v);
        CU_ASSERT_DOUBLE_EQUAL (dot (p.lmn, m), 0.0, 1e-12);
        pol_to_axis (&p, m, &qb, &ub);
        CU_ASSERT_DOUBLE_EQUAL (qa, qb, 1e-10);
        CU_ASSERT_DOUBLE_EQUAL (ua, ub, 1e-10);
      }
    }
  }

  /* REL_MODE_LINEAR: no aberration, so neither the direction nor q, u change */
  rel_mode = REL_MODE_LINEAR;
  init_dummy_phot (&p);
  p.lmn[0] = 0.6;
  p.lmn[1] = 0.0;
  p.lmn[2] = 0.8;
  p.freq = 1e15;
  p.w = 1.0;
  p.q = 0.3;
  p.u = 0.15;
  rescale (vels[2], VLIGHT, v);
  lorentz_transform (&p, &p, v);
  CU_ASSERT_DOUBLE_EQUAL (p.lmn[0], 0.6, 1e-15);
  CU_ASSERT_DOUBLE_EQUAL (p.q, 0.3, 1e-15);
  CU_ASSERT_DOUBLE_EQUAL (p.u, 0.15, 1e-15);

  modes.polarization = FALSE;
  rel_mode = rel_mode_save;
}

/** *******************************************************************************************************************
 *
 * @brief Create the polarization test suite
 *
 * ****************************************************************************************************************** */

void
create_polarization_test_suite (void)
{
  CU_pSuite suite = CU_add_suite ("Polarization", NULL, NULL);

  if (suite == NULL)
  {
    fprintf (stderr, "Failed to create `Polarization` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }

  if ((CU_add_test (suite, "Default basis", test_pol_basis) == NULL) ||
      (CU_add_test (suite, "Rotation of the reference axis", test_pol_rotate) == NULL) ||
      (CU_add_test (suite, "Axis round trip", test_pol_axis_round_trip) == NULL) ||
      (CU_add_test (suite, "North on the sky", test_pol_sky_north) == NULL) ||
      (CU_add_test (suite, "Electron-scattering matrix", test_pol_es_matrix) == NULL) ||
      (CU_add_test (suite, "Electron-scattering peel-off", test_pol_peel_es) == NULL) ||
      (CU_add_test (suite, "Lorentz transformation", test_pol_lorentz) == NULL) ||
      (CU_add_test (suite, "Lorentz transformation, angle from the (n, v) plane", test_pol_lorentz_plane) == NULL))
  {
    fprintf (stderr, "Failed to add tests to `Polarization` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
