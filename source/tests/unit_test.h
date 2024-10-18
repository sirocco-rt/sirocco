/** ********************************************************************************************************************
 *
 *  @file load_model.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date Jan 2024
 *
 *  @brief
 *
 * ****************************************************************************************************************** */

#ifndef PYTHON_UNIT_TEST_H
#define PYTHON_UNIT_TEST_H

int cleanup_model (const char *root_name);
int setup_model_grid (const char *root_name, const char *atomic_data_location);
const char *get_sirocco_env_variable (void);

#endif
