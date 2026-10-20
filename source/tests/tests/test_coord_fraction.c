/** ********************************************************************************************************************
 *
 *  @file test_coord_fraction.c
 *  @date October 2026
 *
 *  @brief Unit tests for coord_fraction and phi_fraction on 3D grids (CYLIND3D, SPH3D)
 *
 * ****************************************************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CUnit/CUnit.h"

#include "../assert.h"

#include "../../atomic.h"
#include "../../sirocco.h"

#define NDIM_T 5
#define MDIM_T 6
#define PDIM_T 8
#define TOL 1e-12

static DomainPtr zdom_save;
static domain_dummy test_dom;
static double wx[NDIM_T], wz[MDIM_T], wmx[NDIM_T], wmz[MDIM_T];
static double wphi[PDIM_T + 1], wmphi[PDIM_T];


/** *******************************************************************************************************************
 *
 * @brief Set up a single test domain with a regular grid in (r, z or theta) and phi
 *
 * @param [in] coord_type  CYLIND, CYLIND3D, RTHETA or SPH3D
 * @param [in] pdim        number of phi cells (1 for 2d grids)
 *
 * ****************************************************************************************************************** */

static void
setup_test_domain (int coord_type, int pdim)
{
  int i;
  double dphi;

  memset (&test_dom, 0, sizeof (test_dom));
  test_dom.coord_type = coord_type;
  test_dom.ndim = NDIM_T;
  test_dom.mdim = MDIM_T;
  test_dom.pdim = pdim;

  for (i = 0; i < NDIM_T; i++)
  {
    wx[i] = 1.0 + i;
    wmx[i] = 1.5 + i;
  }
  for (i = 0; i < MDIM_T; i++)
  {
    if (coord_type == RTHETA || coord_type == SPH3D)
    {
      wz[i] = 15. * i;          /* theta, degrees, upper hemisphere only */
      wmz[i] = 15. * i + 7.5;
    }
    else
    {
      wz[i] = 0.5 * i;
      wmz[i] = 0.5 * i + 0.25;
    }
  }
  dphi = 2. * PI / pdim;
  for (i = 0; i <= pdim; i++)
    wphi[i] = i * dphi;
  for (i = 0; i < pdim; i++)
    wmphi[i] = (i + 0.5) * dphi;

  test_dom.wind_x = wx;
  test_dom.wind_z = wz;
  test_dom.wind_midx = wmx;
  test_dom.wind_midz = wmz;
  test_dom.wind_phi = wphi;
  test_dom.wind_midphi = wmphi;

  zdom = &test_dom;
}


/* Position at cylindrical (rho, phi, z) */
static void
cyl_pos (double rho, double phi, double z, double x[])
{
  x[0] = rho * cos (phi);
  x[1] = rho * sin (phi);
  x[2] = z;
}


/* The phi slice of a 3D element index */
static int
slice_of (int n)
{
  return n % zdom[0].pdim;
}


/** *******************************************************************************************************************
 *
 * @brief The weights sum to 1 and point to two adjacent phi slices
 *
 * ****************************************************************************************************************** */

void
test_weights_sum_to_one (void)
{
  int ichoice, n, nn, nelem, ii[MAX_COORD_FRAC];
  double frac[MAX_COORD_FRAC], x[3], sum, phi;

  setup_test_domain (CYLIND3D, PDIM_T);
  for (ichoice = 0; ichoice < 2; ichoice++)
  {
    for (n = 0; n < 200; n++)
    {
      phi = 2. * PI * n / 200. + 0.001;
      cyl_pos (2.3 + 0.01 * n, phi, 0.7, x);
      coord_fraction (0, ichoice, x, ii, frac, &nelem);
      CU_ASSERT_EQUAL_FATAL (nelem, 8);
      sum = 0;
      for (nn = 0; nn < nelem; nn++)
      {
        CU_ASSERT_FATAL (frac[nn] >= -TOL);
        CU_ASSERT_FATAL (ii[nn] >= 0 && ii[nn] < NDIM_T * MDIM_T * PDIM_T);
        sum += frac[nn];
      }
      CU_ASSERT_DOUBLE_EQUAL_FATAL (sum, 1.0, TOL);
      /* elements 0-3 are in one slice, 4-7 in the next one round */
      for (nn = 1; nn < 4; nn++)
        CU_ASSERT_EQUAL_FATAL (slice_of (ii[nn]), slice_of (ii[0]));
      for (nn = 4; nn < 8; nn++)
        CU_ASSERT_EQUAL_FATAL (slice_of (ii[nn]), (slice_of (ii[0]) + 1) % PDIM_T);
    }
  }
}


/** *******************************************************************************************************************
 *
 * @brief A position on a phi slice (vertices) or at a slice centre (centres) gives weight to that slice only,
 * and a position half way between two slices splits the weight equally
 *
 * ****************************************************************************************************************** */

void
test_on_slice (void)
{
  int k, nn, nelem, ii[MAX_COORD_FRAC];
  double frac[MAX_COORD_FRAC], x[3], w_lo, w_hi;

  setup_test_domain (CYLIND3D, PDIM_T);
  for (k = 0; k < PDIM_T; k++)
  {
    /* vertex values sit at wind_phi[k] */
    cyl_pos (2.3, wphi[k] + 1e-9, 0.7, x);
    coord_fraction (0, 0, x, ii, frac, &nelem);
    w_lo = w_hi = 0;
    for (nn = 0; nn < nelem; nn++)
    {
      if (slice_of (ii[nn]) == k)
        w_lo += frac[nn];
    }
    CU_ASSERT_DOUBLE_EQUAL_FATAL (w_lo, 1.0, 1e-6);

    /* centre values sit at wind_midphi[k] */
    cyl_pos (2.3, wmphi[k] + 1e-9, 0.7, x);
    coord_fraction (0, 1, x, ii, frac, &nelem);
    w_lo = 0;
    for (nn = 0; nn < nelem; nn++)
    {
      if (slice_of (ii[nn]) == k)
        w_lo += frac[nn];
    }
    CU_ASSERT_DOUBLE_EQUAL_FATAL (w_lo, 1.0, 1e-6);

    /* half way between vertex slices k and k+1 */
    cyl_pos (2.3, wmphi[k], 0.7, x);
    coord_fraction (0, 0, x, ii, frac, &nelem);
    w_lo = w_hi = 0;
    for (nn = 0; nn < nelem; nn++)
    {
      if (slice_of (ii[nn]) == k)
        w_lo += frac[nn];
      else if (slice_of (ii[nn]) == (k + 1) % PDIM_T)
        w_hi += frac[nn];
    }
    CU_ASSERT_DOUBLE_EQUAL_FATAL (w_lo, 0.5, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL_FATAL (w_hi, 0.5, 1e-6);
  }
}


/** *******************************************************************************************************************
 *
 * @brief Interpolation wraps from the last phi slice back to slice 0
 *
 * ****************************************************************************************************************** */

void
test_phi_wrap (void)
{
  int k, kk;
  double x[3], dphi, step;

  setup_test_domain (CYLIND3D, PDIM_T);
  step = 2. * PI / PDIM_T;

  /* Vertices: just below 2 pi lies between slice pdim-1 and slice 0 */
  cyl_pos (2.3, 2. * PI - 0.25 * step, 0.7, x);
  phi_fraction (0, 0, x, &k, &kk, &dphi);
  CU_ASSERT_EQUAL_FATAL (k, PDIM_T - 1);
  CU_ASSERT_EQUAL_FATAL (kk, 0);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (dphi, 0.75, 1e-6);

  /* Negative angles from atan2 are handled the same way */
  cyl_pos (2.3, -0.25 * step, 0.7, x);
  phi_fraction (0, 0, x, &k, &kk, &dphi);
  CU_ASSERT_EQUAL_FATAL (k, PDIM_T - 1);
  CU_ASSERT_EQUAL_FATAL (kk, 0);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (dphi, 0.75, 1e-6);

  /* Centres: below the first centre lies between slice pdim-1 and slice 0 */
  cyl_pos (2.3, 0.25 * step, 0.7, x);
  phi_fraction (0, 1, x, &k, &kk, &dphi);
  CU_ASSERT_EQUAL_FATAL (k, PDIM_T - 1);
  CU_ASSERT_EQUAL_FATAL (kk, 0);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (dphi, 0.75, 1e-6);

  /* Centres: above the last centre also wraps */
  cyl_pos (2.3, 2. * PI - 0.25 * step, 0.7, x);
  phi_fraction (0, 1, x, &k, &kk, &dphi);
  CU_ASSERT_EQUAL_FATAL (k, PDIM_T - 1);
  CU_ASSERT_EQUAL_FATAL (kk, 0);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (dphi, 0.25, 1e-6);
}


/** *******************************************************************************************************************
 *
 * @brief With one phi cell a 3D grid gives exactly the 2D result; the (r, z or theta) weights are unchanged
 * by the phi interpolation
 *
 * ****************************************************************************************************************** */

void
test_pdim_one_matches_2d (void)
{
  int types3d[2] = { CYLIND3D, SPH3D };
  int types2d[2] = { CYLIND, RTHETA };
  int t, n, nn, ichoice, nelem2, nelem3, ii2[MAX_COORD_FRAC], ii3[MAX_COORD_FRAC];
  double frac2[MAX_COORD_FRAC], frac3[MAX_COORD_FRAC], x[3], s2[4], s3;

  for (t = 0; t < 2; t++)
  {
    for (ichoice = 0; ichoice < 2; ichoice++)
    {
      for (n = 0; n < 50; n++)
      {
        cyl_pos (1.2 + 0.07 * n, 0.3 + 0.11 * n, 0.05 + 0.04 * n, x);

        setup_test_domain (types2d[t], 1);
        coord_fraction (0, ichoice, x, ii2, frac2, &nelem2);

        /* pdim = 1: identical elements and weights */
        setup_test_domain (types3d[t], 1);
        coord_fraction (0, ichoice, x, ii3, frac3, &nelem3);
        CU_ASSERT_EQUAL_FATAL (nelem3, nelem2);
        for (nn = 0; nn < nelem2; nn++)
        {
          CU_ASSERT_EQUAL_FATAL (ii3[nn], ii2[nn]);
          CU_ASSERT_DOUBLE_EQUAL_FATAL (frac3[nn], frac2[nn], TOL);
        }

        /* pdim > 1: summing the two phi slices recovers the 2d weights */
        for (nn = 0; nn < 4; nn++)
          s2[nn] = frac2[nn];
        setup_test_domain (types3d[t], PDIM_T);
        coord_fraction (0, ichoice, x, ii3, frac3, &nelem3);
        CU_ASSERT_EQUAL_FATAL (nelem3, 8);
        for (nn = 0; nn < 4; nn++)
        {
          s3 = frac3[nn] + frac3[nn + 4];
          CU_ASSERT_DOUBLE_EQUAL_FATAL (s3, s2[nn], TOL);
          /* and the (i, j) of each element matches the 2d one */
          CU_ASSERT_EQUAL_FATAL (ii3[nn] / PDIM_T, ii2[nn]);
          CU_ASSERT_EQUAL_FATAL (ii3[nn + 4] / PDIM_T, ii2[nn]);
        }
      }
    }
  }
}


/** *******************************************************************************************************************
 *
 * @brief A quantity that is linear in phi between two slices is interpolated exactly
 *
 * ****************************************************************************************************************** */

void
test_linear_in_phi (void)
{
  int n, nn, nelem, ii[MAX_COORD_FRAC];
  double frac[MAX_COORD_FRAC], x[3], phi, val, expect, step;

  setup_test_domain (CYLIND3D, PDIM_T);
  step = 2. * PI / PDIM_T;
  for (n = 0; n < 100; n++)
  {
    /* stay within slices 2 and 3 so the slice index is linear in phi */
    phi = 2 * step + step * (n + 0.5) / 100.;
    cyl_pos (2.3, phi, 0.7, x);
    coord_fraction (0, 0, x, ii, frac, &nelem);
    val = 0;
    for (nn = 0; nn < nelem; nn++)
      val += slice_of (ii[nn]) * frac[nn];
    expect = phi / step;
    CU_ASSERT_DOUBLE_EQUAL_FATAL (val, expect, 1e-6);
  }
}


int
coord_fraction_suite_init (void)
{
  zdom_save = zdom;
  return EXIT_SUCCESS;
}

int
coord_fraction_suite_teardown (void)
{
  zdom = zdom_save;
  return EXIT_SUCCESS;
}

void
create_coord_fraction_test_suite (void)
{
  char *suite_name = "coord_fraction 3D";

  CU_pSuite suite = CU_add_suite (suite_name, coord_fraction_suite_init, coord_fraction_suite_teardown);

  if (suite == NULL)
  {
    fprintf (stderr, "Failed to create `%s' suite\n", suite_name);
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }

  if ((CU_add_test (suite, "Weights sum to one", test_weights_sum_to_one) == NULL) ||
      (CU_add_test (suite, "Positions on phi slices", test_on_slice) == NULL) ||
      (CU_add_test (suite, "Periodic wrap in phi", test_phi_wrap) == NULL) ||
      (CU_add_test (suite, "pdim = 1 matches 2d", test_pdim_one_matches_2d) == NULL) ||
      (CU_add_test (suite, "Linear in phi", test_linear_in_phi) == NULL))
  {
    fprintf (stderr, "Failed to add tests to `%s'\n", suite_name);
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
