/** ********************************************************************************************************************
 *
 *  @file unit_test.h
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date Jan 2024
 *
 *  @brief
 *
 * ****************************************************************************************************************** */

#define TEST_TOLERANCE 0.01

#ifndef SIROCCO_UNIT_TEST_H
#define SIROCCO_UNIT_TEST_H

int cleanup_model (const char *root_name);
int setup_model_grid (const char *root_name, const char *atomic_data_location);
const char *get_sirocco_env_variable (void);

#endif
