/** ********************************************************************************************************************
 *
 *  @file test_translation.c
 *  @author James Matthews (james.matthews@physics.ox.ac.uk)
 *  @date October 2023
 *
 *  @brief Unit tests for translating photons and ds routines
 *
 * ****************************************************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "CUnit/CUnit.h"

#include "../assert.h"

/* we want to use the same photon and plane structures as defined in sirocco.h */
#include "../../atomic.h"
#include "../../sirocco.h"

double
test_ds_one_photon (struct photon *ptest, int force_positive_z)
{
  double ds1, ds2, ds3, rsphere;
  struct plane windplane;
  double r = 10.0;

  windplane.x[0] = windplane.x[1] = 0;
  windplane.x[2] = r;
  windplane.lmn[0] = windplane.lmn[1] = 0;
  windplane.lmn[2] = 1;

  ds1 = ds_to_plane (&windplane, ptest, force_positive_z);
  ds2 = ds_to_cylinder (r, ptest);
  rsphere = sqrt(2) * r; 
  ds3 = ds_to_sphere(rsphere, ptest);
  // printf("ds1 ds2 ds3 %8.4e %8.4e %8.4e\n", ds1, ds2, ds3);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (ds1, ds2, EPSILON * r);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (ds1, ds3, EPSILON * r); 
  return ds1;
}


/** *******************************************************************************************************************
 *
 * @brief Test case for ds routines
 *
 * @param [in] test_name the test name
 *
 * @return an error status
 *
 * @details
 * ****************************************************************************************************************** */

void
test_ds (void)
{
  struct photon ptest;
  double ds1, ds2;
  double r = 10.0;
  int idir, jdir;

  /* we loop over a variable that swaps x and y around, since the 
     ds routines should work in cylindrical symmetry */
  for (idir = 0; idir < 2; idir++)
  {
    if (idir == 0)
    {
      jdir = 1;
    }
    else
    {
      jdir = 0;
    }
      
    /* fire off a photon at 45 degrees, starting from a small positive distance from the origin */
    ptest.lmn[idir] = ptest.lmn[2] = 0.70710678118;
    ptest.lmn[jdir] = 0.0;
    ptest.x[0] = ptest.x[1] = ptest.x[2] = EPSILON; 
    test_ds_one_photon (&ptest, FALSE);

    /* move the photon to the -ve hemisphere and fire downwards at 45 deg */
    ptest.x[2] = -EPSILON;
    ptest.lmn[2] = -0.70710678118; 
    test_ds_one_photon (&ptest, TRUE);

    /* move the photon just into the +ve hemisphere and fire downwards at 45 deg */
    ptest.x[2] = EPSILON;
    ptest.lmn[2] = -0.70710678118; 
    test_ds_one_photon (&ptest, TRUE);

    /* check that if you move the photon halfway to its destination the
        distance remaining is half */
    ptest.lmn[2] = 0.70710678118; 
    ds1 = test_ds_one_photon (&ptest, FALSE);
    move_phot(&ptest, ds1 * 0.5);
    ds2 = test_ds_one_photon (&ptest, FALSE);
    CU_ASSERT_DOUBLE_EQUAL_FATAL (ds1 * 0.5, ds2, EPSILON * r);

    /* check that moving the photon 50% away from its destination 
        then reversing the direction makes the distance increase by 1.5 */
    ptest.lmn[idir] = ptest.lmn[2] = -0.70710678118; 
    ptest.x[0] = ptest.x[1] = ptest.x[2] = EPSILON; 
    move_phot(&ptest, ds1 * 0.5);
    ptest.lmn[idir] = ptest.lmn[2] = 0.70710678118; 
    ds2 = test_ds_one_photon (&ptest, FALSE);
    CU_ASSERT_DOUBLE_EQUAL_FATAL (ds1 * 1.5, ds2, EPSILON * r);
  }
}
/** *******************************************************************************************************************
 *
 * @brief Initialise the translate test suite
 *
 * @return An error code status
 *
 * @details
 * ****************************************************************************************************************** */


int
translate_suite_init (void)
{
  int error = EXIT_SUCCESS;
  return error;
}

/** *******************************************************************************************************************
 *
 * @brief Clean up after the translate test suite
 *
 * @return An error code status
 *
 * @details
 * ****************************************************************************************************************** */

int
translate_suite_teardown (void)
{
  int error = EXIT_SUCCESS;
  return error;
}

void
create_translate_test_suite (void)
{
  char *suite_name = "Translate Functions";

  CU_pSuite suite = CU_add_suite (suite_name, translate_suite_init, translate_suite_teardown);

  if (suite == NULL)
  {
    fprintf (stderr, "Failed to create `Translate Functions' suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }

  /* Add CPU tests to suite */
  if ((CU_add_test (suite, "Photon ds routines", test_ds) == NULL))
  {
    fprintf (stderr, "Failed to add tests to  `Translate Functions'\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
