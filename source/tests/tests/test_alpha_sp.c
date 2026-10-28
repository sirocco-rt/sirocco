/** ********************************************************************************************************************
 *
 *  @file test_alpha_sp.c
 *  @date October 2026
 *
 *  @brief Unit tests for the macro-atom bound-free rate integrals bf_rate_integrals, alpha_sp_all and alpha_sp
 *
 *  bf_rate_integrals integrates the log-log interpolated cross section of a continuum times nu^2 and the Maxwellian
 *  and Planck factors with Gauss-Legendre rules. These tests hold it to an independent adaptive quadrature (GSL cquad
 *  of sigma_phot itself, one call per table interval) for every photoionization continuum of the 20-level H + He
 *  macro-atom data set, over the temperatures a wind reaches, and check that no temperature, however extreme, gives
 *  an infinite, NaN or negative integral, and that the integrals alpha_sp_all remembers are the fresh ones. The
 *  free-bound emissivities and recombination rates of xinteg_fb are held to the same kind of reference of
 *  fb_topbase_partial.
 *
 * ****************************************************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <gsl/gsl_integration.h>
#include <CUnit/CUnit.h>

#include "../../atomic.h"
#include "../../sirocco.h"

#define N_TEST_TEMPERATURES 13

/* the accuracy asked of the reference, cquad at 1e-13, with margin; bf_rate_integrals itself is good to ~1e-13 */
#define ALPHA_SP_TOLERANCE 1e-10

static struct topbase_phot *ref_cont;
static double ref_beta, ref_beta_r;
static int ref_weight;

static double
reference_integrand (double freq, void *params)
{
  const double fthresh = ref_cont->freq[0];
  double f = sigma_phot (ref_cont, freq) * freq * freq * exp (ref_beta * (fthresh - freq));

  (void) params;
  if (ref_beta_r > 0)
    f /= expm1 (ref_beta_r * freq);
  if (ref_weight == 1)
    f *= freq / fthresh;
  else if (ref_weight == 2)
    f *= (freq - fthresh) / fthresh;
  return f;
}

/** *******************************************************************************************************************
 *
 * @brief The reference: adaptive quadrature of sigma_phot times the factors
 *
 * @details
 * One cquad call per -ln DBL_EPSILON e-folds of the exponential factor within each table interval: over a whole
 * interval many e-folds wide cquad can miss the narrow peak at its start by more than its tolerance. The integral
 * stops where the exponential factors have fallen by DBL_EPSILON^2.
 *
 * ****************************************************************************************************************** */

static double
reference_integral (struct topbase_phot *cont_ptr, double freq_max, double beta, double beta_r, int weight)
{
  gsl_integration_cquad_workspace *ws = gsl_integration_cquad_workspace_alloc (200);
  gsl_function F = { reference_integrand, NULL };
  double total = 0, result, abserr;
  size_t neval;
  int i;

  ref_cont = cont_ptr;
  ref_beta = beta;
  ref_beta_r = beta_r;
  ref_weight = weight;
  if (beta + beta_r > 0)
    freq_max = fmin (freq_max, cont_ptr->freq[0] - 2 * log (DBL_EPSILON) / (beta + beta_r));
  for (i = 0; i < cont_ptr->np - 1; i++)
  {
    const double a = cont_ptr->freq[i], b = fmin (cont_ptr->freq[i + 1], freq_max);
    const int pieces = b > a ? 1 + (int) ((beta + beta_r) * (b - a) / -log (DBL_EPSILON)) : 0;
    int k;
    for (k = 0; k < pieces; k++)
    {
      gsl_integration_cquad (&F, a + (b - a) * k / pieces, a + (b - a) * (k + 1) / pieces, 0, 1e-13, ws, &result, &abserr, &neval);
      total += result;
    }
  }
  gsl_integration_cquad_workspace_free (ws);
  return total;
}

/** *******************************************************************************************************************
 *
 * @brief alpha_sp_all and alpha_sp are the reference integral over the whole table times the constants
 *
 * ****************************************************************************************************************** */

static void
test_alpha_sp_matches_reference (void)
{
  struct plasma xplasma;
  double alpha[3];
  int n_cont, n_temp, ichoice;

  memset (&xplasma, 0, sizeof xplasma);
  for (n_cont = 0; n_cont < nphot_total; n_cont++)
  {
    struct topbase_phot *cont_ptr = &phot_top[n_cont];
    const double g_upper = cont_ptr->macro_info == TRUE && geo.macro_simple == FALSE ? xconfig[cont_ptr->uplev].g : ion[cont_ptr->nion + 1].g;
    for (n_temp = 0; n_temp < N_TEST_TEMPERATURES; n_temp++)
    {
      const double t_e = pow (10., 3.3 + 0.25 * n_temp + 0.0137);       /* 2e3 K to 2e6 K */
      xplasma.state.t_e = t_e;
      alpha_sp_all (cont_ptr, &xplasma, alpha);
      for (ichoice = 0; ichoice < 3; ichoice++)
      {
        const double reference = reference_integral (cont_ptr, cont_ptr->freq[cont_ptr->np - 1], H_OVER_K / t_e, 0, ichoice)
          * xconfig[cont_ptr->nlev].g / g_upper * pow (t_e, -1.5) * ALPHA_SP_CONSTANT;
        CU_ASSERT_TRUE (alpha[ichoice] > 0.0);
        CU_ASSERT_DOUBLE_EQUAL (alpha[ichoice], reference, ALPHA_SP_TOLERANCE * reference);
        CU_ASSERT_DOUBLE_EQUAL (alpha_sp (cont_ptr, &xplasma, ichoice), alpha[ichoice], 0.0);
      }
    }
  }
}

/** *******************************************************************************************************************
 *
 * @brief With a Planck factor (stimulated recombination, photoionization) the integrals are the reference
 *
 * @details
 * T_r is T_e / 2 and T_e, and the photoionization integrals have T_r = T_e and no Maxwellian. Only where h
 * nu_threshold / k T_r < 200: deeper in the Wien tail cquad itself is unreliable; extreme temperatures are
 * covered by test_extreme_temperatures.
 *
 * ****************************************************************************************************************** */

static void
test_stimulated_matches_reference (void)
{
  double sum[3];
  int n_cont, n_temp, n_rad, ichoice;

  for (n_cont = 0; n_cont < nphot_total; n_cont++)
  {
    struct topbase_phot *cont_ptr = &phot_top[n_cont];
    const double fthresh = cont_ptr->freq[0], flast = cont_ptr->freq[cont_ptr->np - 1];
    for (n_temp = 0; n_temp < N_TEST_TEMPERATURES; n_temp++)
    {
      const double t_e = pow (10., 3.3 + 0.25 * n_temp + 0.0137);
      for (n_rad = 0; n_rad < 3; n_rad++)
      {
        const double t_r = n_rad == 0 ? 0.5 * t_e : t_e, beta = n_rad < 2 ? H_OVER_K / t_e : 0;
        if (H_OVER_K * fthresh / t_r > 200)
          continue;
        bf_rate_integrals (cont_ptr, fthresh, flast, beta, H_OVER_K / t_r, sum);
        for (ichoice = 0; ichoice < 3; ichoice++)
        {
          const double reference = reference_integral (cont_ptr, flast, beta, H_OVER_K / t_r, ichoice);
          CU_ASSERT_DOUBLE_EQUAL (sum[ichoice], reference, ALPHA_SP_TOLERANCE * reference);
        }
      }
    }
  }
}

/** *******************************************************************************************************************
 *
 * @brief Band-limited integrals of two adjacent bands add up to the integral over both
 *
 * ****************************************************************************************************************** */

static void
test_band_limited_adds_up (void)
{
  struct plasma xplasma;
  int n_cont;

  memset (&xplasma, 0, sizeof xplasma);
  xplasma.state.t_e = 1.7e4;
  for (n_cont = 0; n_cont < nphot_total; n_cont++)
  {
    struct topbase_phot *cont_ptr = &phot_top[n_cont];
    const double fthresh = cont_ptr->freq[0], flast = cont_ptr->freq[cont_ptr->np - 1], fmid = fthresh * 1.3;
    const double full = scaled_alpha_sp_integral_band_limited (cont_ptr, &xplasma, 0, fthresh, flast);
    const double lo = scaled_alpha_sp_integral_band_limited (cont_ptr, &xplasma, 0, fthresh, fmid);
    const double hi = scaled_alpha_sp_integral_band_limited (cont_ptr, &xplasma, 0, fmid, flast);
    CU_ASSERT_DOUBLE_EQUAL (lo + hi, full, 1e-12 * full);
    CU_ASSERT_TRUE (lo > 0 && hi >= 0);
  }
}

/** *******************************************************************************************************************
 *
 * @brief Extreme temperatures give finite, non-negative integrals; zero temperature gives zero
 *
 * ****************************************************************************************************************** */

static void
test_extreme_temperatures (void)
{
  const double temps[] = { 1e-3, 1, 30, 1e8, 1e10, 1e12 };
  const int n_temps = sizeof temps / sizeof temps[0];
  double sum[3];
  int n_cont, i, j, ichoice;

  for (n_cont = 0; n_cont < nphot_total; n_cont++)
  {
    struct topbase_phot *cont_ptr = &phot_top[n_cont];
    const double fthresh = cont_ptr->freq[0], flast = cont_ptr->freq[cont_ptr->np - 1];
    for (i = 0; i < n_temps; i++)
    {
      for (j = 0; j <= n_temps; j++)
      {
        const double beta_r = j < n_temps ? H_OVER_K / temps[j] : 0;
        bf_rate_integrals (cont_ptr, fthresh, flast, H_OVER_K / temps[i], beta_r, sum);
        for (ichoice = 0; ichoice < 3; ichoice++)
          CU_ASSERT_TRUE (isfinite (sum[ichoice]) && sum[ichoice] >= 0);
      }
      bf_rate_integrals (cont_ptr, fthresh, flast, 0, H_OVER_K / temps[i], sum);
      for (ichoice = 0; ichoice < 3; ichoice++)
        CU_ASSERT_TRUE (isfinite (sum[ichoice]) && sum[ichoice] >= 0);
    }
    bf_rate_integrals (cont_ptr, fthresh, flast, INFINITY, 0, sum);
    for (ichoice = 0; ichoice < 3; ichoice++)
      CU_ASSERT_EQUAL (sum[ichoice], 0.0);
    bf_rate_integrals (cont_ptr, fthresh, flast, 0, INFINITY, sum);
    for (ichoice = 0; ichoice < 3; ichoice++)
      CU_ASSERT_EQUAL (sum[ichoice], 0.0);
  }
}

/** *******************************************************************************************************************
 *
 * @brief alpha_sp_all gives the same rates, bit for bit, whether it integrates afresh or remembers
 *
 * @details
 * Three temperatures, more than are remembered, in an order that returns to each of them after one and after two
 * others, so that every call is either a remembered entry or a replacement of one; the full-band
 * scaled_alpha_sp_integral_band_limited shares the remembered integrals.
 *
 * ****************************************************************************************************************** */

static void
test_memo_is_exact (void)
{
  const double temps[] = { 9.1e3, 1.3e4, 2.7e4 };
  const int order[] = { 0, 0, 1, 0, 1, 2, 0, 2, 2, 1, 0, 1 };
  struct plasma xplasma;
  double alpha[3], fresh[3][3], sum[3];
  int n_cont, i, ichoice;

  memset (&xplasma, 0, sizeof xplasma);
  for (n_cont = 0; n_cont < nphot_total; n_cont++)
  {
    struct topbase_phot *cont_ptr = &phot_top[n_cont];
    for (i = 0; i < 3; i++)
    {
      alpha_sp_forget ();
      xplasma.state.t_e = temps[i];
      alpha_sp_all (cont_ptr, &xplasma, fresh[i]);
    }
    alpha_sp_forget ();
    for (i = 0; i < (int) (sizeof order / sizeof order[0]); i++)
    {
      xplasma.state.t_e = temps[order[i]];
      alpha_sp_all (cont_ptr, &xplasma, alpha);
      for (ichoice = 0; ichoice < 3; ichoice++)
        CU_ASSERT_EQUAL (alpha[ichoice], fresh[order[i]][ichoice]);
      bf_rate_integrals (cont_ptr, cont_ptr->freq[0], cont_ptr->freq[cont_ptr->np - 1], H_OVER_K / xplasma.state.t_e, 0, sum);
      CU_ASSERT_EQUAL (scaled_alpha_sp_integral_band_limited (cont_ptr, &xplasma, 0, 0, VERY_BIG), sum[0]);
    }
  }
  alpha_sp_forget ();
}

/* fb_topbase_partial takes its record, temperature and choice from these (recomb.c) */
extern struct topbase_phot *fb_xtop;
extern double fbt, log_fbt;
extern int fbfr;

static double
fb_reference_integrand (double freq, void *params)
{
  (void) params;
  return fb_topbase_partial (freq);
}

/** @brief cquad of fb_topbase_partial over lo ... hi for the record cont_ptr, pieces as in reference_integral */

static double
fb_reference_integral (struct topbase_phot *cont_ptr, double t, double lo, double hi, int fb_choice)
{
  gsl_integration_cquad_workspace *ws = gsl_integration_cquad_workspace_alloc (200);
  gsl_function F = { fb_reference_integrand, NULL };
  const double big_b = H_OVER_K / t;
  double total = 0, result, abserr;
  size_t neval;
  int i, k;

  fb_xtop = cont_ptr;
  fbt = t;
  log_fbt = log (t);
  fbfr = fb_choice;
  hi = fmin (hi, cont_ptr->freq[0] - 2 * log (DBL_EPSILON) / big_b);
  for (i = 0; i < cont_ptr->np - 1; i++)
  {
    const double a = fmax (lo, cont_ptr->freq[i]), b = fmin (cont_ptr->freq[i + 1], hi);
    const int pieces = b > a ? 1 + (int) (big_b * (b - a) / -log (DBL_EPSILON)) : 0;
    for (k = 0; k < pieces; k++)
    {
      gsl_integration_cquad (&F, a + (b - a) * k / pieces, a + (b - a) * (k + 1) / pieces, 0, 1e-13, ws, &result, &abserr, &neval);
      total += result;
    }
  }
  gsl_integration_cquad_workspace_free (ws);
  return total;
}

/** *******************************************************************************************************************
 *
 * @brief xinteg_fb, for each free-bound choice, is the reference integral of fb_topbase_partial, and two halves of
 * a band add up to the whole
 *
 * @details
 * Every ion with photoionization records, 2e3 K to 6e8 K (below 2e3 K cquad itself misses 1e-10 on some records,
 * by up to 2e-9 at 100 K against a long-double Gauss-Legendre integral, while xinteg_fb stays within 2e-12 of it);
 * the macro-atom records count as simple ones here
 * (geo.macro_simple), as xinteg_fb otherwise leaves them to the macro-atom machinery.
 *
 * ****************************************************************************************************************** */

static void
test_free_bound_matches_reference (void)
{
  const int choices[] = { FB_FULL, FB_REDUCED, FB_RATE };
  const int macro_simple = geo.macro_simple;
  int nion, n_temp, c, n;

  geo.macro_simple = TRUE;
  for (nion = 0; nion < nions; nion++)
  {
    int nmin, nmax;
    if (ion[nion].phot_info > 0)
      nmin = ion[nion].ntop_first, nmax = nmin + ion[nion].ntop;
    else if (ion[nion].phot_info == 0)
      nmin = ion[nion].nxphot, nmax = nmin + 1;
    else
      continue;
    for (n_temp = 0; n_temp < 12; n_temp++)
    {
      const double t = pow (10., 3.3 + 0.5 * n_temp + 0.0137);
      for (c = 0; c < 3; c++)
      {
        const double value = xinteg_fb (t, 0, VERY_BIG, nion, choices[c]);
        const double fmid = 1.3 * phot_top[nmin].freq[0];
        double reference = 0;
        for (n = nmin; n < nmax; n++)
          reference += fb_reference_integral (&phot_top[n], t, 3e12, 3e18, choices[c]);
        CU_ASSERT_TRUE (value > 0);
        CU_ASSERT_DOUBLE_EQUAL (value, reference, ALPHA_SP_TOLERANCE * reference);
        CU_ASSERT_DOUBLE_EQUAL (xinteg_fb (t, 0, fmid, nion, choices[c]) + xinteg_fb (t, fmid, VERY_BIG, nion, choices[c]), value, 1e-12 * value);
      }
    }
  }
  geo.macro_simple = macro_simple;
}

/** *******************************************************************************************************************
 *
 * @brief Read the 20-level H + He macro-atom data set through a `data` link to $SIROCCO/xdata
 *
 * ****************************************************************************************************************** */

static int
suite_init (void)
{
  char target[LINELENGTH];
  struct stat sb;
  const char *sirocco_env = getenv ("SIROCCO");

  if (sirocco_env == NULL)
  {
    fprintf (stderr, "Failed to find SIROCCO environment variable\n");
    return EXIT_FAILURE;
  }
  snprintf (target, LINELENGTH, "%s/xdata", sirocco_env);
  if (!(stat (target, &sb) == EXIT_SUCCESS && S_ISDIR (sb.st_mode)))
  {
    perror ("Unable to find atomic data directory");
    return EXIT_FAILURE;
  }
  if (symlink (target, "data") != EXIT_SUCCESS && errno != EEXIST)
  {
    perror ("Unable to created symbolic link for atomic data for test case");
    return EXIT_FAILURE;
  }

  geo.macro_simple = FALSE;
  if (get_atomic_data ("data/h20_hetop_standard80.dat") != 0 || nphot_total <= 0)
  {
    fprintf (stderr, "Failed to read the H20 + He atomic data\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Create a CUnit test suite for alpha_sp
 *
 * ****************************************************************************************************************** */

void
create_alpha_sp_test_suite (void)
{
  CU_pSuite suite = CU_add_suite ("Macro-atom recombination integrals", suite_init, NULL);

  if (suite == NULL)
  {
    fprintf (stderr, "Failed to create `Macro-atom recombination integrals` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }

  if ((CU_add_test (suite, "alpha_sp matches adaptive quadrature", test_alpha_sp_matches_reference) == NULL) ||
      (CU_add_test (suite, "Stimulated and photoionization integrals match adaptive quadrature", test_stimulated_matches_reference) == NULL) ||
      (CU_add_test (suite, "Band-limited integrals add up", test_band_limited_adds_up) == NULL) ||
      (CU_add_test (suite, "Extreme temperatures", test_extreme_temperatures) == NULL) ||
      (CU_add_test (suite, "Remembered integrals are the fresh ones", test_memo_is_exact) == NULL) ||
      (CU_add_test (suite, "Free-bound integrals match adaptive quadrature", test_free_bound_matches_reference) == NULL))
  {
    fprintf (stderr, "Failed to add tests to `Macro-atom recombination integrals` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
