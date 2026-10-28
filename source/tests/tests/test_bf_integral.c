/** ********************************************************************************************************************
 *
 *  @file test_bf_integral.c
 *  @date October 2026
 *
 *  @brief Unit tests for bf_cross_section_integral, the integral of a bound-free cross section times a weight
 *
 *  The cross sections are made up here, so no atomic data are needed. Power-law cross sections with weights for which
 *  the integral is analytic test the quadrature itself; a cross section with many resonances, like the TOPbase cross
 *  sections, is compared with a fine sum over the same interpolated cross section.
 *
 * ****************************************************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CUnit/CUnit.h"

#include "../../atomic.h"
#include "../../sirocco.h"

#define NPTS_MAX 2000

static double xfreq[NPTS_MAX], xlog_freq[NPTS_MAX], xx[NPTS_MAX], xlog_x[NPTS_MAX];
static Topbase_phot xcont;

/** *******************************************************************************************************************
 *
 * @brief Fill the test continuum from np frequencies and cross sections
 *
 * ****************************************************************************************************************** */

static void
set_cross_section (int np, const double *freq, const double *x)
{
  int i;

  memset (&xcont, 0, sizeof (xcont));
  xcont.np = np;
  xcont.freq = xfreq;
  xcont.log_freq = xlog_freq;
  xcont.x = xx;
  xcont.log_x = xlog_x;
  for (i = 0; i < np; i++)
  {
    xfreq[i] = freq[i];
    xx[i] = x[i];
    xlog_freq[i] = log (freq[i]);
    xlog_x[i] = (x[i] > 0) ? log (x[i]) : -999.;
  }
}

/** *******************************************************************************************************************
 *
 * @brief The cross section in segment i, interpolated as sigma_phot does (a power law, or linear if either point is
 * zero)
 *
 * ****************************************************************************************************************** */

static double
interpolated_sigma (int i, double nu)
{
  if (xx[i] > 0 && xx[i + 1] > 0)
    return exp (xlog_x[i] + (xlog_x[i + 1] - xlog_x[i]) * (log (nu) - xlog_freq[i]) / (xlog_freq[i + 1] - xlog_freq[i]));
  return xx[i] + (xx[i + 1] - xx[i]) * (nu - xfreq[i]) / (xfreq[i + 1] - xfreq[i]);
}

/* Weights: 1, and nu^2 exp (a (nu_0 - nu)) as in alpha_sp */

static double
unit_weight (double nu, void *params)
{
  return 1.0;
}

struct exp_weight_params
{
  double nu0, a;
};

static double
exp_weight (double nu, void *params)
{
  struct exp_weight_params *p = (struct exp_weight_params *) params;
  return nu * nu * exp (p->a * (p->nu0 - nu));
}

/** *******************************************************************************************************************
 *
 * @brief A power-law cross section, sigma = s0 (nu / nu0)^-3, with unit weight
 *
 * ****************************************************************************************************************** */

static void
test_power_law (void)
{
  double freq[6], x[6], nu0 = 1e15, s0 = 1e-18, expected, result;
  int i;

  for (i = 0; i < 6; i++)
  {
    freq[i] = nu0 * pow (2., i);
    x[i] = s0 * pow (freq[i] / nu0, -3.);
  }
  set_cross_section (6, freq, x);

  expected = s0 * nu0 / 2. * (1. - pow (freq[5] / nu0, -2.));
  result = bf_cross_section_integral (&xcont, freq[0], freq[5], 0.0, unit_weight, NULL);
  CU_ASSERT_DOUBLE_EQUAL (result / expected, 1.0, 1e-10);
}

/** *******************************************************************************************************************
 *
 * @brief sigma proportional to nu^-2, weight nu^2 exp (a (nu_0 - nu)), so the integrand is a pure exponential
 *
 * The test covers temperatures for which the exponential falls by very little, and by very much, across the
 * cross section, including the case where the integral stops early.
 *
 * ****************************************************************************************************************** */

static void
test_exponential (void)
{
  double freq[20], x[20], nu0 = 3e15, s0 = 6e-18, expected, result, t;
  struct exp_weight_params p;
  int i, k;

  for (i = 0; i < 20; i++)
  {
    freq[i] = nu0 * pow (1.3, i);
    x[i] = s0 * pow (freq[i] / nu0, -2.);
  }
  set_cross_section (20, freq, x);

  for (k = 0; k < 6; k++)
  {
    t = pow (10., 3. + k);
    p.nu0 = nu0;
    p.a = H_OVER_K / t;
    expected = s0 * nu0 * nu0 * (-expm1 (-p.a * (freq[19] - nu0))) / p.a;
    result = bf_cross_section_integral (&xcont, freq[0], freq[19], p.a, exp_weight, &p);
    CU_ASSERT_DOUBLE_EQUAL (result / expected, 1.0, 1e-8);
  }
}

/** *******************************************************************************************************************
 *
 * @brief A cross section with many narrow resonances, compared with a fine sum
 *
 * ****************************************************************************************************************** */

static void
test_resonances (void)
{
  double freq[1200], x[1200], nu0 = 2e15, s0 = 1e-18, f, expected, result, t, df, nu;
  struct exp_weight_params p;
  int i, k, n, nsub = 2000;

  for (i = 0; i < 1200; i++)
  {
    freq[i] = nu0 * pow (1.004, i);
    f = (i % 7 == 3) ? 300. : ((i % 13 == 5) ? 0.01 : 1.);
    x[i] = s0 * f * pow (freq[i] / nu0, -2.5);
  }
  set_cross_section (1200, freq, x);

  for (k = 0; k < 3; k++)
  {
    t = pow (10., 4. + k);
    p.nu0 = nu0;
    p.a = H_OVER_K / t;

    expected = 0;
    for (i = 0; i < 1199; i++)
    {
      df = (freq[i + 1] - freq[i]) / nsub;
      for (n = 0; n < nsub; n++)
      {
        nu = freq[i] + (n + 0.5) * df;
        expected += interpolated_sigma (i, nu) * exp_weight (nu, &p) * df;
      }
    }
    result = bf_cross_section_integral (&xcont, freq[0], freq[1199], p.a, exp_weight, &p);
    CU_ASSERT_DOUBLE_EQUAL (result / expected, 1.0, 1e-4);
  }
}

/** *******************************************************************************************************************
 *
 * @brief Integrals over adjacent bands add up to the integral over both, and limits outside the cross section
 * contribute nothing
 *
 * ****************************************************************************************************************** */

static void
test_bands (void)
{
  double freq[5] = { 1e15, 1.5e15, 2e15, 4e15, 8e15 };
  double x[5] = { 1e-18, 3e-18, 5e-19, 2e-19, 1e-19 };
  double whole, part1, part2, outside;
  struct exp_weight_params p;

  set_cross_section (5, freq, x);
  p.nu0 = freq[0];
  p.a = H_OVER_K / 3e4;

  whole = bf_cross_section_integral (&xcont, freq[0], freq[4], p.a, exp_weight, &p);
  part1 = bf_cross_section_integral (&xcont, freq[0], 1.77e15, p.a, exp_weight, &p);
  part2 = bf_cross_section_integral (&xcont, 1.77e15, freq[4], p.a, exp_weight, &p);
  CU_ASSERT_DOUBLE_EQUAL ((part1 + part2) / whole, 1.0, 1e-10);

  outside = bf_cross_section_integral (&xcont, 0.5 * freq[0], 2. * freq[4], p.a, exp_weight, &p);
  CU_ASSERT_DOUBLE_EQUAL (outside / whole, 1.0, 1e-12);

  CU_ASSERT_EQUAL (bf_cross_section_integral (&xcont, 1e14, 5e14, p.a, exp_weight, &p), 0.0);
}

/** *******************************************************************************************************************
 *
 * @brief A segment with a zero cross section is interpolated linearly, and an infinite rate gives zero
 *
 * ****************************************************************************************************************** */

static void
test_zero_and_infinite (void)
{
  double freq[3] = { 1e15, 2e15, 3e15 };
  double x[3] = { 0.0, 2e-18, 2e-18 };
  double result;

  set_cross_section (3, freq, x);

  /* triangle from 0 to 2e-18 over 1e15, then a rectangle */
  result = bf_cross_section_integral (&xcont, freq[0], freq[2], 0.0, unit_weight, NULL);
  CU_ASSERT_DOUBLE_EQUAL (result / (0.5 * 2e-18 * 1e15 + 2e-18 * 1e15), 1.0, 1e-12);

  result = bf_cross_section_integral (&xcont, freq[0], freq[2], INFINITY, unit_weight, NULL);
  CU_ASSERT_EQUAL (result, 0.0);
}

void
create_bf_integral_test_suite (void)
{
  char *suite_name = "bf cross-section integral";

  CU_pSuite suite = CU_add_suite (suite_name, NULL, NULL);

  if (suite == NULL)
  {
    fprintf (stderr, "Failed to create `%s' suite\n", suite_name);
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }

  if ((CU_add_test (suite, "Power law", test_power_law) == NULL) ||
      (CU_add_test (suite, "Exponential", test_exponential) == NULL) ||
      (CU_add_test (suite, "Resonances", test_resonances) == NULL) ||
      (CU_add_test (suite, "Bands", test_bands) == NULL) ||
      (CU_add_test (suite, "Zero and infinite", test_zero_and_infinite) == NULL))
  {
    fprintf (stderr, "Failed to add tests to `%s'\n", suite_name);
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
