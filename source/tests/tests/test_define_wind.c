/** ********************************************************************************************************************
 *
 *  @file test_define_wind.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date November 2023
 *
 *  @brief Unit tests for `define_wind`
 *
 *  TODO: checking actual against test data should be modular, rather than copy and paste
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
#include "../assert.h"
#include "../unit_test.h"

char *SIROCCO_ENV;
char TEST_CWD[LINELENGTH];
char ATOMIC_DATA_TARGET[LINELENGTH];
char ATOMIC_DATA_DEST[LINELENGTH];
char ATOMIC_DATA_TARGET_DEVELOPER[LINELENGTH];
char ATOMIC_DATA_DEST_DEVELOPER[LINELENGTH];

#define TEST_DATA_LENGTH 2056

/** *******************************************************************************************************************
 *
 * @brief Test an SV wind model for an AGN system with macro atoms.
 *
 * @details
 *
 * This uses the data from $SIROCCO/source/tests/test_data/define_wind/agn_macro.pf and
 * $SIROCCO/source/tests/test_data/define_wind/agn_macro.grid.txt. The latter was created using the Python script in the
 * $SIROCCO/source/tests/test_data/define_wind directory.
 *
 * TODO: we don't check anything to do with macro atoms (we don't output this from windsave2table)
 *
 * ****************************************************************************************************************** */

static void
test_sv_agn_macro_wind (void)
{
  int n;
  FILE *fp;
  char test_data_line[TEST_DATA_LENGTH];
  char test_data_filename[LINELENGTH];

  WindPtr wind_cell;
  PlasmaPtr plasma_cell;

  const int init_error = setup_model_grid ("agn_macro", ATOMIC_DATA_DEST);
  if (init_error)
  {
    CU_FAIL_FATAL ("Unable to initialise AGN Macro model");
  }

  define_wind ();

  snprintf (test_data_filename, LINELENGTH, "%s/source/tests/test_data/define_wind/agn_macro.grid.txt", SIROCCO_ENV);
  fp = fopen (test_data_filename, "r");
  if (fp == NULL)
  {
    CU_FAIL_FATAL ("Unable to open test data for AGN Macro");
  }

  /* For the CV model, we'll be looking at these parameters */
  int i, j, inwind;
  double x, z, xcen, zcen;
  double v_x, v_y, v_z;
  double vol, rho, ne, h1, c4;
  double t_e, t_r;
  double dv_x_dx, dv_x_dy, dv_x_dz;
  double dv_y_dx, dv_y_dy, dv_y_dz;
  double dv_z_dx, dv_z_dy, dv_z_dz;
  double div_v, dvds_max;
  double gamma;

  /* Skip the first line */
  if (fgets (test_data_line, TEST_DATA_LENGTH, fp) == NULL)
  {
    CU_FAIL_FATAL ("Unable to read first line of test data");
  }

  while (fgets (test_data_line, TEST_DATA_LENGTH, fp) != NULL)
  {
    /* Here's what the header of the file should look like:
     * # i j x z xcen zcen inwind v_x v_y v_z vol rho ne t_e t_r h1 c4 dv_x_dx dv_y_dx dv_z_dx dv_x_dy dv_y_dy dv_z_dy
     * dv_x_dz dv_y_dz dv_z_dz div_v dvds_max gamma
     */
    const short n_read = sscanf (test_data_line,
                                 "%d %d %le %le %le %le %d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le",
                                 &i, &j, &x, &z, &xcen, &zcen, &inwind, &v_x, &v_y, &v_z, &vol, &rho, &ne, &t_e, &t_r, &h1, &c4, &dv_x_dx,
                                 &dv_y_dx, &dv_z_dx, &dv_x_dy, &dv_y_dy, &dv_z_dy, &dv_x_dz, &dv_y_dz, &dv_z_dz, &div_v, &dvds_max, &gamma);
    if (n_read != 29)
    {
      CU_FAIL_FATAL ("Test data is in an invalid format");
    }

    /* Convert wind indices into an n in 1d n */
    wind_ij_to_n (0, i, j, &n);
    wind_cell = &wmain[n];
    plasma_cell = &plasmamain[wind_cell->nplasma];

    /* cell positions */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->x[0], x, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->x[2], z, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->xcen[0], xcen, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->xcen[2], zcen, TEST_TOLERANCE);
    CU_ASSERT_EQUAL_FATAL (wind_cell->inwind, inwind);

    /* The default behaviour of Python's output tools (e.g. windsave2table) is
     * to ignore file which are not fully in the wind. So we shall also ignore
     * them here */
    if (wind_cell->inwind != W_ALL_INWIND)
    {
      continue;
    }

    /* velocities */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[0], v_x, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[1], v_y, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[2], v_z, TEST_TOLERANCE);
    /* velocity gradients */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][0], dv_x_dx, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][1], dv_x_dy, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][2], dv_x_dz, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][0], dv_y_dx, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][1], dv_y_dy, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][2], dv_y_dz, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][0], dv_z_dx, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][1], dv_z_dy, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][2], dv_z_dz, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->div_v, div_v, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->dvds_max, dvds_max, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->xgamma, gamma, TEST_TOLERANCE);

    /* Some things (plasma properties) are stored in plasma cells */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->rho, rho, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->ne, ne, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->t_e, t_e, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->t_r, t_r, TEST_TOLERANCE);

    /* Ion abundances are tested in their number density relative to Hydrogen.
     * This is the default output option in windsave2table */
    const double n_h = rho2nh * plasma_cell->rho;
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->density[0] / (n_h * ele[0].abun), h1, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->density[8] / (n_h * ele[2].abun), c4, TEST_TOLERANCE);
  }

  fclose (fp);
  cleanup_model ("agn_macro");
}

/** *******************************************************************************************************************
 *
 * @brief Test a SV wind model for a CV system.
 *
 * @details
 *
 * This uses the data from $SIROCCO/source/tests/test_data/define_wind/cv.pf and
 * $SIROCCO/source/tests/test_data/define_wind/cv.grid.txt. The latter was created using the Python script in the
 * $SIROCCO/source/tests/test_data/define_wind directory.
 *
 * ****************************************************************************************************************** */

static void
test_sv_cv_wind (void)
{
  int n;
  FILE *fp;
  char test_data_line[TEST_DATA_LENGTH];
  char test_data_filename[LINELENGTH];
  char windsave_filename[LINELENGTH];

  WindPtr wind_cell;
  PlasmaPtr plasma_cell;

  const int init_error = setup_model_grid ("cv", ATOMIC_DATA_DEST);
  if (init_error)
  {
    cleanup_model ("cv");
    CU_FAIL_FATAL ("Unable to initialise CV model");
  }

  /* With the defined, we can try and create the wind */
  define_wind ();

  /* And now we can compare our created grid to the "ground truth" grid */
  snprintf (test_data_filename, LINELENGTH, "%s/source/tests/test_data/define_wind/cv.grid.txt", SIROCCO_ENV);
  fp = fopen (test_data_filename, "r");
  if (fp == NULL)
  {
    cleanup_model ("cv");
    CU_FAIL_FATAL ("Unable to open test data for CV model");
  }

  /* For the CV model, we'll be looking at these parameters */
  int i, j, inwind;
  double x, z, xcen, zcen;
  double v_x, v_y, v_z;
  double vol, rho, ne, h1, c4;
  double t_e, t_r;
  double dv_x_dx, dv_x_dy, dv_x_dz;
  double dv_y_dx, dv_y_dy, dv_y_dz;
  double dv_z_dx, dv_z_dy, dv_z_dz;
  double div_v, dvds_max;
  double gamma;

  /* Skip the first line */
  if (fgets (test_data_line, TEST_DATA_LENGTH, fp) == NULL)
  {
    cleanup_model ("cv");
    CU_FAIL_FATAL ("Unable to read first line of test data");
  }

  while (fgets (test_data_line, TEST_DATA_LENGTH, fp) != NULL)
  {
    /* Here's what the header of the file should look like:
     * # i j x z xcen zcen inwind v_x v_y v_z vol rho ne t_e t_r h1 c4 dv_x_dx dv_y_dx dv_z_dx dv_x_dy dv_y_dy dv_z_dy
     * dv_x_dz dv_y_dz dv_z_dz div_v dvds_max gamma
     */
    const short n_read = sscanf (test_data_line,
                                 "%d %d %le %le %le %le %d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le",
                                 &i, &j, &x, &z, &xcen, &zcen, &inwind, &v_x, &v_y, &v_z, &vol, &rho, &ne, &t_e, &t_r, &h1, &c4, &dv_x_dx,
                                 &dv_y_dx, &dv_z_dx, &dv_x_dy, &dv_y_dy, &dv_z_dy, &dv_x_dz, &dv_y_dz, &dv_z_dz, &div_v, &dvds_max, &gamma);
    if (n_read != 29)
    {
      CU_FAIL_FATAL ("Test data is in an invalid format");
    }

    /* Convert wind indices into an n in 1d wmain */
    wind_ij_to_n (0, i, j, &n);
    wind_cell = &wmain[n];
    plasma_cell = &plasmamain[wind_cell->nplasma];

    /* cell positions */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->x[0], x, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->x[2], z, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->xcen[0], xcen, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->xcen[2], zcen, TEST_TOLERANCE);
    CU_ASSERT_EQUAL_FATAL (wind_cell->inwind, inwind);

    /* The default behaviour of Python's output tools (e.g. windsave2table) is
     * to ignore file which are not fully in the wind. So we shall also ignore
     * them here */
    if (wind_cell->inwind != W_ALL_INWIND)
    {
      continue;
    }

    /* velocities */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[0], v_x, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[1], v_y, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[2], v_z, TEST_TOLERANCE);
    /* velocity gradients */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][0], dv_x_dx, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][1], dv_x_dy, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][2], dv_x_dz, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][0], dv_y_dx, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][1], dv_y_dy, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][2], dv_y_dz, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][0], dv_z_dx, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][1], dv_z_dy, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][2], dv_z_dz, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->div_v, div_v, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->dvds_max, dvds_max, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->xgamma, gamma, TEST_TOLERANCE);

    /* Some things (plasma properties) are stored in plasma cells */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->rho, rho, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->ne, ne, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->t_e, t_e, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->t_r, t_r, TEST_TOLERANCE);

    /* Ion abundances are tested in their number density relative to Hydrogen.
     * This is the default output option in windsave2table */
    const double n_h = rho2nh * plasma_cell->rho;
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->density[0] / (n_h * ele[0].abun), h1, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->density[8] / (n_h * ele[2].abun), c4, TEST_TOLERANCE);
  }

  /* For the CV model, we want to save the wind_save to use in another test */
  snprintf (windsave_filename, LINELENGTH, "%s/source/tests/test_data/define_wind/restart_cv.wind_save", SIROCCO_ENV);
  const int err = wind_save (windsave_filename);
  if (err == 0)
  {
    CU_FAIL ("Failed to produce wind_save for CV test case");
  }

  fclose (fp);
  cleanup_model ("cv");
}

/** *******************************************************************************************************************
 *
 * @brief Test a shell wind
 *
 * @details
 *
 * ****************************************************************************************************************** */

static void
test_shell_wind (void)
{
  FILE *fp;
  char test_data_line[TEST_DATA_LENGTH];
  char test_data_filename[LINELENGTH];

  WindPtr wind_cell;
  PlasmaPtr plasma_cell;

  const int init_error = setup_model_grid ("shell", ATOMIC_DATA_DEST_DEVELOPER);
  if (init_error)
  {
    cleanup_model ("shell");
    CU_FAIL_FATAL ("Unable to initialise shell model");
  }

  /* With the defined, we can try and create the wind */
  define_wind ();

  /* And now we can compare our created grid to the "ground truth" grid */
  snprintf (test_data_filename, LINELENGTH, "%s/source/tests/test_data/define_wind/shell.grid.txt", SIROCCO_ENV);
  fp = fopen (test_data_filename, "r");
  if (fp == NULL)
  {
    cleanup_model ("shell");
    CU_FAIL_FATAL ("Unable to open test data for shell model");
  }

  int i, inwind;
  double r, rcen;
  double v_x, v_y, v_z;
  double vol, rho, ne, h1, c4;
  double t_e, t_r;
  double dv_x_dx, dv_x_dy, dv_x_dz;
  double dv_y_dx, dv_y_dy, dv_y_dz;
  double dv_z_dx, dv_z_dy, dv_z_dz;
  double div_v, dvds_max;
  double gamma;

  /* Skip the first line */
  if (fgets (test_data_line, TEST_DATA_LENGTH, fp) == NULL)
  {
    CU_FAIL_FATAL ("Unable to read first line of test data");
  }

  while (fgets (test_data_line, TEST_DATA_LENGTH, fp) != NULL)
  {
    /*
     * We'll read in the following properties:
     * i r rcen inwind v_x v_y v_z vol rho ne t_e t_r h1 c4 dv_x_dx dv_y_dx dv_z_dx dv_x_dy dv_y_dy dv_z_dy
     * dv_x_dz dv_y_dz dv_z_dz div_v dvds_max gamma
     */

    const short n_read = sscanf (test_data_line,
                                 "%d %le %le %d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le",
                                 &i, &r, &rcen, &inwind, &v_x, &v_y, &v_z, &vol, &rho, &ne, &t_e, &t_r, &h1, &c4,
                                 &dv_x_dx,
                                 &dv_y_dx, &dv_z_dx, &dv_x_dy, &dv_y_dy, &dv_z_dy, &dv_x_dz, &dv_y_dz, &dv_z_dz, &div_v,
                                 &dvds_max, &gamma);
    if (n_read != 26)
    {
      cleanup_model ("shell");
      CU_FAIL_FATAL ("Test data is in an invalid format");
    }

    /* Convert wind indices into an n in 1d wmain */
    wind_cell = &wmain[i];
    plasma_cell = &plasmamain[wind_cell->nplasma];

    /* cell positions */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->r, r, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->rcen, rcen, TEST_TOLERANCE);
    CU_ASSERT_EQUAL_FATAL (wind_cell->inwind, inwind);

    /* The default behaviour of Python's output tools (e.g. windsave2table) is
     * to ignore file which are not fully in the wind. So we shall also ignore
     * them here */
    if (wind_cell->inwind != W_ALL_INWIND)
    {
      continue;
    }

    /* velocities */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[0], v_x, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[1], v_y, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[2], v_z, TEST_TOLERANCE);
    /* velocity gradients */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][0], dv_x_dx, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][1], dv_x_dy, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][2], dv_x_dz, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][0], dv_y_dx, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][1], dv_y_dy, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][2], dv_y_dz, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][0], dv_z_dx, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][1], dv_z_dy, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][2], dv_z_dz, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->div_v, div_v, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->dvds_max, dvds_max, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->xgamma, gamma, TEST_TOLERANCE);

    /* Some things (plasma properties) are stored in plasma cells */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->rho, rho, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->ne, ne, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->t_e, t_e, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->t_r, t_r, TEST_TOLERANCE);

    /* Ion abundances are tested in their number density relative to Hydrogen.
     * This is the default output option in windsave2table */
    const double n_h = rho2nh * plasma_cell->rho;
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->density[0] / (n_h * ele[0].abun), h1, TEST_TOLERANCE);
  }

  fclose (fp);
  cleanup_model ("shell");
}

/** *******************************************************************************************************************
 *
 * @brief Test a spherical stellar wind
 *
 * @details
 *
 * ****************************************************************************************************************** */

static void
test_spherical_star_wind (void)
{
  FILE *fp;
  char test_data_line[TEST_DATA_LENGTH];
  char test_data_filename[LINELENGTH];

  WindPtr wind_cell;
  PlasmaPtr plasma_cell;

  const int init_error = setup_model_grid ("star", ATOMIC_DATA_DEST);
  if (init_error)
  {
    cleanup_model ("star");
    CU_FAIL_FATAL ("Unable to initialise star model");
  }

  /* With the defined, we can try and create the wind */
  define_wind ();

  /* And now we can compare our created grid to the "ground truth" grid */
  snprintf (test_data_filename, LINELENGTH, "%s/source/tests/test_data/define_wind/star.grid.txt", SIROCCO_ENV);
  fp = fopen (test_data_filename, "r");
  if (fp == NULL)
  {
    cleanup_model ("star");
    CU_FAIL_FATAL ("Unable to open test data for star model");
  }

  int i, inwind;
  double r, rcen;
  double v_x, v_y, v_z;
  double vol, rho, ne, h1, c4;
  double t_e, t_r;
  double dv_x_dx, dv_x_dy, dv_x_dz;
  double dv_y_dx, dv_y_dy, dv_y_dz;
  double dv_z_dx, dv_z_dy, dv_z_dz;
  double div_v, dvds_max;
  double gamma;

  /* Skip the first line */
  if (fgets (test_data_line, TEST_DATA_LENGTH, fp) == NULL)
  {
    CU_FAIL_FATAL ("Unable to read first line of test data");
  }

  while (fgets (test_data_line, TEST_DATA_LENGTH, fp) != NULL)
  {
    /*
     * We'll read in the following properties:
     * i r rcen inwind v_x v_y v_z vol rho ne t_e t_r h1 c4 dv_x_dx dv_y_dx dv_z_dx dv_x_dy dv_y_dy dv_z_dy
     * dv_x_dz dv_y_dz dv_z_dz div_v dvds_max gamma
     */

    const short n_read = sscanf (test_data_line,
                                 "%d %le %le %d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le",
                                 &i, &r, &rcen, &inwind, &v_x, &v_y, &v_z, &vol, &rho, &ne, &t_e, &t_r, &h1, &c4, &dv_x_dx,
                                 &dv_y_dx, &dv_z_dx, &dv_x_dy, &dv_y_dy, &dv_z_dy, &dv_x_dz, &dv_y_dz, &dv_z_dz, &div_v, &dvds_max,
                                 &gamma);
    if (n_read != 26)
    {
      cleanup_model ("star");
      CU_FAIL_FATAL ("Test data is in an invalid format");
    }

    /* Convert wind indices into an n in 1d wmain */
    wind_cell = &wmain[i];
    plasma_cell = &plasmamain[wind_cell->nplasma];

    /* cell positions */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->r, r, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->rcen, rcen, TEST_TOLERANCE);
    CU_ASSERT_EQUAL_FATAL (wind_cell->inwind, inwind);

    /* The default behaviour of Python's output tools (e.g. windsave2table) is
     * to ignore file which are not fully in the wind. So we shall also ignore
     * them here */
    if (wind_cell->inwind != W_ALL_INWIND)
    {
      continue;
    }

    /* velocities */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[0], v_x, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[1], v_y, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[2], v_z, TEST_TOLERANCE);
    /* velocity gradients */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][0], dv_x_dx, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][1], dv_x_dy, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][2], dv_x_dz, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][0], dv_y_dx, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][1], dv_y_dy, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][2], dv_y_dz, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][0], dv_z_dx, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][1], dv_z_dy, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][2], dv_z_dz, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->div_v, div_v, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->dvds_max, dvds_max, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->xgamma, gamma, TEST_TOLERANCE);

    /* Some things (plasma properties) are stored in plasma cells */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->rho, rho, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->ne, ne, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->t_e, t_e, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->t_r, t_r, TEST_TOLERANCE);

    /* Ion abundances are tested in their number density relative to Hydrogen.
     * This is the default output option in windsave2table */
    const double n_h = rho2nh * plasma_cell->rho;
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->density[0] / (n_h * ele[0].abun), h1, TEST_TOLERANCE);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->density[8] / (n_h * ele[2].abun), c4, TEST_TOLERANCE);
  }

  fclose (fp);
  cleanup_model ("star");
}

/** *******************************************************************************************************************
 *
 * @brief Clean up after the `define_wind` tests have run.
 *
 * @details
 *
 * This will deallocate the dynamic memory allocated for the domain structure and will also remove the symbolic link
 * to the atomic data from the test directory.
 *
 * ****************************************************************************************************************** */

static int
suite_teardown (void)
{
  if (unlink (ATOMIC_DATA_DEST) != EXIT_SUCCESS)
  {
    perror ("Unable to unlink test data symbolic link");        /* We won't worry too hard about this */
  }

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Initialise the testing suite for testing `define_wind`
 *
 * @details
 *
 * To be able to create a plasma grid, we need atomic data. The main purpose of this function is to figure out where
 * the tests are being run and to create a symbolic link to the atomic data and to   create file paths to the atomic
 * data. This function will also initialise some global properties in Python, allocating space for the domain structure
 * and, initialise the geometry structure and initialising MPI if that is required. Note that the unit tests should only
 * ever be run using a single MPI rank.
 *
 * Also sets some global variables (to this file) associated with environment variables.
 *
 * ****************************************************************************************************************** */

static int
suite_init (void)
{
  struct stat sb;

  /* Find the SIROCCO env var and the directory where tests are called. We need
   * these to create a symbolic link to the atomic data required for the
   * tests */
  SIROCCO_ENV = getenv ("SIROCCO");
  if (SIROCCO_ENV == NULL)
  {
    fprintf (stderr, "Failed to find SIROCCO environment variable\n");
    return EXIT_FAILURE;
  }
  if (getcwd (TEST_CWD, LINELENGTH) == NULL)
  {
    perror ("Failed to find current working directory for tests");
    return EXIT_FAILURE;
  }

  /* Set global variables for atomic data */
  snprintf (ATOMIC_DATA_TARGET, LINELENGTH, "%s/xdata", SIROCCO_ENV);
  if (!(stat (ATOMIC_DATA_TARGET, &sb) == EXIT_SUCCESS && S_ISDIR (sb.st_mode)))
  {
    perror ("Unable to find atomic data directory");
    return EXIT_FAILURE;
  }

  snprintf (ATOMIC_DATA_DEST, LINELENGTH, "%s/data", TEST_CWD);
  if (symlink (ATOMIC_DATA_TARGET, ATOMIC_DATA_DEST) != EXIT_SUCCESS)
  {
    /* If the symlink exists, we'll try not worry about it as if something is
     * wrong with the atomic data it'll be caught later */
    if (errno != EEXIST)
    {
      perror ("Unable to created symbolic link for atomic data for test case");
      return EXIT_FAILURE;
    }
  }

  /* Set global variables for atomic data for developers */
  snprintf (ATOMIC_DATA_TARGET_DEVELOPER, LINELENGTH, "%s/zdata", SIROCCO_ENV);
  if (!(stat (ATOMIC_DATA_TARGET_DEVELOPER, &sb) == EXIT_SUCCESS && S_ISDIR (sb.st_mode)))
  {
    perror ("Unable to find atomic data directory");
    return EXIT_FAILURE;
  }

  snprintf (ATOMIC_DATA_DEST_DEVELOPER, LINELENGTH, "%s/zdata", TEST_CWD);
  if (symlink (ATOMIC_DATA_TARGET_DEVELOPER, ATOMIC_DATA_DEST_DEVELOPER) != EXIT_SUCCESS)
  {
    /* If the symlink exists, we'll try not worry about it as if something is
     * wrong with the atomic data it'll be caught later */
    if (errno != EEXIST)
    {
      perror ("Unable to created symbolic link for atomic data for test case");
      return EXIT_FAILURE;
    }
  }

  /* Now initialise the things Python will need to work, and things which aren't
   * specific to the model being tested such as the domain allocation */
  rel_mode = REL_MODE_FULL;

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Create a CUnit test suite for define wind
 *
 * @details
 *
 * This function will create a test suite for functions related to creating the wind and plasma grid
 *
 * ****************************************************************************************************************** */

void
create_define_wind_test_suite (void)
{
  CU_pSuite suite = CU_add_suite ("Define Wind", suite_init, suite_teardown);

  if (suite == NULL)
  {
    fprintf (stderr, "Failed to create `Define Wind` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }

  if ((CU_add_test (suite, "Shell wind", test_shell_wind) == NULL) ||
      (CU_add_test (suite, "Spherical: Supernova", test_spherical_star_wind) == NULL) ||
      (CU_add_test (suite, "SV: Cataclysmic Variable", test_sv_cv_wind) == NULL) ||
      (CU_add_test (suite, "SV: AGN Macro", test_sv_agn_macro_wind) == NULL))
  {
    fprintf (stderr, "Failed to add tests to `Define Wind` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
