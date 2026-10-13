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
      (CU_add_test (suite, "North on the sky", test_pol_sky_north) == NULL))
  {
    fprintf (stderr, "Failed to add tests to `Polarization` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
