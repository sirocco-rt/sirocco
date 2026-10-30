/** ********************************************************************************************************************
 *
 *  @file tests.h
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief Function prototypes for test suites
 *
 * ****************************************************************************************************************** */

#ifndef TESTS_H
#define TESTS_H

/* test_compton.c */
void create_compton_test_suite (void);

/* test_sv_cv_wind */
void create_define_wind_test_suite (void);

/* test_matrix.c */
void create_matrix_test_suite (void);

/* test_run_mode.c */
void create_run_mode_test_suite (void);

/* test_matrix.c */
void create_translate_test_suite (void);

/* test_macro_estimators.c */
void create_macro_estimators_test_suite (void);

/* test_coord_fraction.c */
void create_coord_fraction_test_suite (void);

/* test_atomic_data.c */
void create_atomic_data_test_suite (void);

/* test_alpha_sp.c */
void create_alpha_sp_test_suite (void);

#endif
