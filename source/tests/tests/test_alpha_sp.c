/** ********************************************************************************************************************
 *
 *  @file test_alpha_sp.c
 *  @date October 2026
 *
 *  @brief Unit tests for the macro-atom spontaneous recombination rate alpha_sp
 *
 *  alpha_sp is held to an independent adaptive quadrature (GSL cquad of sigma_phot itself, one call per table
 *  interval) for every photoionization continuum of the 20-level H + He macro-atom data set, over the electron
 *  temperatures a wind reaches.
 *
 * ****************************************************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <gsl/gsl_integration.h>
#include <CUnit/CUnit.h>

#include "../../atomic.h"
#include "../../sirocco.h"

#define ALPHA_SP_CONSTANT 5.79618e-36
#define N_TEST_TEMPERATURES 13

/* the 10000-point quadrature of alpha_sp is good to 2e-3 */
#define ALPHA_SP_TOLERANCE 2e-3

static struct topbase_phot *ref_cont;
static double ref_beta;
static int ref_weight;

static double
reference_integrand (double freq, void *params)
{
  const double fthresh = ref_cont->freq[0];
  double f = sigma_phot (ref_cont, freq) * freq * freq * exp (ref_beta * (fthresh - freq));

  (void) params;
  if (ref_weight == 1)
    f *= freq / fthresh;
  else if (ref_weight == 2)
    f *= (freq - fthresh) / fthresh;
  return f;
}

/** *******************************************************************************************************************
 *
 * @brief The reference: adaptive quadrature of sigma_phot times the factors, one call per table interval
 *
 * ****************************************************************************************************************** */

static double
reference_integral (struct topbase_phot *cont_ptr, double freq_max, double beta, int weight)
{
  gsl_integration_cquad_workspace *ws = gsl_integration_cquad_workspace_alloc (200);
  gsl_function F = { reference_integrand, NULL };
  double total = 0, result, abserr;
  size_t neval;
  int i;

  ref_cont = cont_ptr;
  ref_beta = beta;
  ref_weight = weight;
  for (i = 0; i < cont_ptr->np - 1; i++)
  {
    const double a = cont_ptr->freq[i], b = fmin (cont_ptr->freq[i + 1], freq_max);
    if (b > a)
    {
      gsl_integration_cquad (&F, a, b, 0, 1e-13, ws, &result, &abserr, &neval);
      total += result;
    }
  }
  gsl_integration_cquad_workspace_free (ws);
  return total;
}

/** *******************************************************************************************************************
 *
 * @brief alpha_sp is the reference integral, up to ALPHA_MATOM_NUMAX_LIMIT k T_e / h above the threshold, times
 * the constants
 *
 * ****************************************************************************************************************** */

static void
test_alpha_sp_matches_reference (void)
{
  struct plasma xplasma;
  int n_cont, n_temp, ichoice;

  memset (&xplasma, 0, sizeof xplasma);
  for (n_cont = 0; n_cont < nphot_total; n_cont++)
  {
    struct topbase_phot *cont_ptr = &phot_top[n_cont];
    const double g_upper = cont_ptr->macro_info == TRUE && geo.macro_simple == FALSE ? xconfig[cont_ptr->uplev].g : ion[cont_ptr->nion + 1].g;
    for (n_temp = 0; n_temp < N_TEST_TEMPERATURES; n_temp++)
    {
      const double t_e = pow (10., 3.3 + 0.25 * n_temp + 0.0137);       /* 2e3 K to 2e6 K */
      const double flast = fmin (cont_ptr->freq[cont_ptr->np - 1], cont_ptr->freq[0] + t_e * ALPHA_MATOM_NUMAX_LIMIT / H_OVER_K);
      xplasma.state.t_e = t_e;
      for (ichoice = 0; ichoice < 3; ichoice++)
      {
        const double reference = reference_integral (cont_ptr, flast, H_OVER_K / t_e, ichoice)
          * xconfig[cont_ptr->nlev].g / g_upper * pow (t_e, -1.5) * ALPHA_SP_CONSTANT;
        const double alpha = alpha_sp (cont_ptr, &xplasma, ichoice);
        CU_ASSERT_TRUE (alpha > 0.0);
        CU_ASSERT_DOUBLE_EQUAL (alpha, reference, ALPHA_SP_TOLERANCE * reference);
      }
    }
  }
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

  if (CU_add_test (suite, "alpha_sp matches adaptive quadrature", test_alpha_sp_matches_reference) == NULL)
  {
    fprintf (stderr, "Failed to add tests to `Macro-atom recombination integrals` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
