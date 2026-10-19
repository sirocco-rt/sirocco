/** ********************************************************************************************************************
 *
 *  @file test_macro_estimators.c
 *  @date September 2026
 *
 *  @brief Unit test for the per-cycle reset of the macro-atom bound-free estimators
 *
 *  On the SV AGN macro-atom model (test_data/define_wind/agn_macro), `init_macro_rad_properties` must zero gamma,
 *  gamma_e, alpha_st and alpha_st_e of every upward bound-free jump. These arrays are indexed by bfu_indx_first;
 *  zeroing alpha_st/alpha_st_e at bfd_indx_first left stale values. (The MPI broadcast of alpha_st_e is not unit
 *  tested: the unit-test driver refuses more than one rank.)
 *
 * ****************************************************************************************************************** */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <CUnit/CUnit.h>

#include "../../atomic.h"
#include "../../sirocco.h"
#include "../unit_test.h"

static char ATOMIC_DATA_DEST[LINELENGTH];
static int MODEL_READY = FALSE;

static void
test_init_resets_bf_estimators (void)
{
  int n, j, level, k;

  if (!MODEL_READY)
    CU_FAIL_FATAL ("AGN macro model not set up");

  for (n = 0; n < NPLASMA; ++n)
  {
    for (j = 0; j < size_gamma_est; ++j)
    {
      macromain[n].est.gamma[j] = macromain[n].est.gamma_e[j] = 1.0;
      macromain[n].est.alpha_st[j] = macromain[n].est.alpha_st_e[j] = 1.0;
    }
  }

  init_macro_rad_properties ();

  for (n = 0; n < NPLASMA; ++n)
  {
    for (level = 0; level < nlevels_macro; ++level)
    {
      for (k = 0; k < xconfig[level].n_bfu_jump; ++k)
      {
        j = xconfig[level].bfu_indx_first + k;
        CU_ASSERT_EQUAL_FATAL (macromain[n].est.gamma[j], 0.0);
        CU_ASSERT_EQUAL_FATAL (macromain[n].est.gamma_e[j], 0.0);
        CU_ASSERT_EQUAL_FATAL (macromain[n].est.alpha_st[j], 0.0);
        CU_ASSERT_EQUAL_FATAL (macromain[n].est.alpha_st_e[j], 0.0);
      }
    }
  }
}

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
  snprintf (ATOMIC_DATA_DEST, LINELENGTH, "data");
  if (symlink (target, ATOMIC_DATA_DEST) != EXIT_SUCCESS && errno != EEXIST)
  {
    perror ("Unable to created symbolic link for atomic data for test case");
    return EXIT_FAILURE;
  }

  rel_mode = REL_MODE_FULL;
  if (setup_model_grid ("agn_macro", ATOMIC_DATA_DEST))
  {
    fprintf (stderr, "Unable to initialise AGN Macro model\n");
    return EXIT_FAILURE;
  }
  define_wind ();
  MODEL_READY = nlevels_macro > 0 && size_gamma_est > 0 && NPLASMA > 0;

  return EXIT_SUCCESS;
}

static int
suite_teardown (void)
{
  if (MODEL_READY)
    cleanup_model ("agn_macro");
  MODEL_READY = FALSE;
  return EXIT_SUCCESS;
}

void
create_macro_estimators_test_suite (void)
{
  CU_pSuite suite = CU_add_suite ("Macro-atom estimators", suite_init, suite_teardown);

  if (suite == NULL)
  {
    fprintf (stderr, "Failed to create `Macro-atom estimators` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }

  if (CU_add_test (suite, "init_macro_rad_properties resets every bfu estimator", test_init_resets_bf_estimators) == NULL)
  {
    fprintf (stderr, "Failed to add tests to `Macro-atom estimators` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
