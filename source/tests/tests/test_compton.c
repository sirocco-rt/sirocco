/** ********************************************************************************************************************
 *
 *  @file test_compton.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief
 *
 * ****************************************************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <CUnit/CUnit.h>

#include "../../atomic.h"
#include "../../sirocco.h"
#include "../assert.h"

/** *******************************************************************************************************************
 *
 * @brief Test the Compton formula, which calculates the fractional energy change
 *
 * @details
 *
 * In this test, a (photon) frequency of 5e16 Hz is used. Testing this function is a bit trickier, as we need access
 * to some (static) global variables.
 *
 * ****************************************************************************************************************** */

void
test_compton_func (void)
{
  const double test_frequency = 5e18;
  const double cross_section_test = 0.5;
  const double energy_ratio = (PLANCK * test_frequency) / (MELEC * VLIGHT * VLIGHT);
  const double f_min = 1.0;
  const double f_max = 1 + (2 * energy_ratio);
  const double mid_f_max = 0.5 * (f_min + f_max);
  const double three_quarters_f_max = 0.75 * (f_min + f_max);
  const double cross_section_max = sigma_compton_partial (f_max, energy_ratio);
  set_comp_func_values (cross_section_test, cross_section_max, energy_ratio);

  CU_ASSERT_DOUBLE_EQUAL_FATAL (compton_func (f_min, NULL), -cross_section_test, EPSILON);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (compton_func (mid_f_max, NULL), 2.130244e-02, EPSILON);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (compton_func (three_quarters_f_max, NULL), 1.457703e+02, EPSILON);
}

/** *******************************************************************************************************************
 *
 * @brief Test the Compton cooling cross-section formula
 *
 * @details
 *
 * This contains two tests. One test where the frequency is below a threshold, where the cross section is equal to 1.0.
 * The other test has a cross-section < 1.
 *
 * ****************************************************************************************************************** */

void
test_compton_beta (void)
{
  CU_ASSERT_DOUBLE_EQUAL_FATAL (compton_beta (1e16), 1.0, EPSILON);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (compton_beta (1e17), 0.9991066689339109, EPSILON);
}

/** *******************************************************************************************************************
 *
 * @brief Test the Compton heating cross-section formula
 *
 * @details
 *
 * This contains two tests. One test where the frequency is below a threshold, where the cross section is equal to 1.0.
 * The other test has a cross-section < 1.
 *
 * ****************************************************************************************************************** */

void
test_compton_alpha (void)
{
  CU_ASSERT_DOUBLE_EQUAL_FATAL (compton_alpha (1e16), 1.0, EPSILON);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (compton_alpha (1e17), 0.9964273280931072, EPSILON);
}

/** *******************************************************************************************************************
 *
 * @brief Test the Klein-Nishina cross-section formula
 *
 * @details
 *
 * Four tests are performed. Two tests for frequencies below the threshold, which should return the Thompson scattering
 * cross-section. The remaining two tests will be corrected by the KN formula.
 *
 * ****************************************************************************************************************** */

void
test_klein_nishina (void)
{
  /* These frequencies should be below the KM threshold and be equal to the
     Thompson cross-section */
  CU_ASSERT_DOUBLE_EQUAL_FATAL (klein_nishina (1e14), THOMPSON, EPSILON);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (klein_nishina (1e15), THOMPSON, EPSILON);

  /* These frequency values will need the KN correction  --
     I think Epsilon needs to be a bit smaller in this case... */
  CU_ASSERT_DOUBLE_EQUAL_FATAL (klein_nishina (2e16), 6.650136664443343e-25, 1e-30);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (klein_nishina (1e18), 6.546940139261951e-25, 1e-30);
}

/** *******************************************************************************************************************
 *
 * @brief Create a CUnit test suite for Compton processes
 *
 * @details
 *
 * This function will create a test suite for functions related to Compton processes, e.g. scattering or thermal
 * velocity distribution.
 *
 * ****************************************************************************************************************** */

void
create_compton_test_suite (void)
{
  CU_pSuite suite = CU_add_suite ("Compton Processes", NULL, NULL);

  if (suite == NULL)
  {
    fprintf (stderr, "Failed to create `Compton Processes` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }

  /* Add CPU tests to suite */
  if ((CU_add_test (suite, "Klein-Nisina Formula", test_klein_nishina) == NULL) ||
      (CU_add_test (suite, "Compton Alpha - heating cross section ", test_compton_alpha) == NULL) ||
      (CU_add_test (suite, "Compton Beta - cooling cross section", test_compton_beta) == NULL) ||
      (CU_add_test (suite, "Compton Formula", test_compton_func) == NULL))
  {
    fprintf (stderr, "Failed to add tests to `Compton Processes` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
