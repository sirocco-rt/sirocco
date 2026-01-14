/***********************************************************/
/** @file  windtable_vars.c
 * @author Long, Matthews
 * @date   January 2026
 *
 * @brief  Variable registry for windtable - THE ONLY PLACE to add new variables
 *
 * To add a new variable that get_one() already supports:
 *   1. Add one entry to the windtable_vars[] array below
 *   2. That's it. No other code changes needed.
 *
 * To add a new variable that get_one() doesn't support:
 *   1. Add the variable to get_one() in windsave2table_sub.c
 *   2. Add one entry to the windtable_vars[] array below
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atomic.h"
#include "sirocco.h"
#include "windtable.h"

/**
 * @brief Static variable registry
 *
 * Each entry defines one output column:
 *   {column_name, unit, get_one_name, is_int}
 *
 * - column_name: Name that appears in FITS table
 * - unit: Physical unit string (empty string if dimensionless)
 * - get_one_name: Variable name to pass to get_one()
 * - is_int: 1 if the value is an integer count, 0 for double
 */
static windtable_var_t windtable_vars[] = {
  /* === MASTER TABLE VARIABLES === */
  {"vol", "cm^3", "vol", 0},
  {"rho", "g/cm^3", "rho", 0},
  {"ne", "cm^-3", "ne", 0},
  {"t_e", "K", "t_e", 0},
  {"t_r", "K", "t_r", 0},
  {"dmo_dt_x", "g.cm/s^2", "dmo_dt_x", 0},
  {"dmo_dt_y", "g.cm/s^2", "dmo_dt_y", 0},
  {"dmo_dt_z", "g.cm/s^2", "dmo_dt_z", 0},
  {"ip", "", "ip", 0},
  {"xi", "erg.cm/s", "xi", 0},
  {"ntot", "", "ntot", 1},
  {"nrad", "", "nrad", 1},
  {"nioniz", "", "nioniz", 1},
  {"nscat_es", "", "nscat_es", 1},
  {"nscat_res", "", "nscat_res", 1},
  {"nscat_ff", "", "nscat_ff", 1},
  {"nscat_bf", "", "nscat_bf", 1},
  {"converge", "", "converge", 0},

  /* === HEAT TABLE VARIABLES === */
  {"w", "", "w", 0},
  {"ave_freq", "Hz", "ave_freq", 0},
  {"J", "erg/cm^2/s/Hz", "J", 0},
  {"J_direct", "erg/cm^2/s/Hz", "J_direct", 0},
  {"J_scatt", "erg/cm^2/s/Hz", "J_scatt", 0},
  {"lum_tot", "erg/s", "lum_tot", 0},
  {"heat_tot", "erg/s", "heat_tot", 0},
  {"heat_comp", "erg/s", "heat_comp", 0},
  {"heat_lines", "erg/s", "heat_lines", 0},
  {"heat_ff", "erg/s", "heat_ff", 0},
  {"heat_photo", "erg/s", "heat_photo", 0},
  {"heat_auger", "erg/s", "heat_auger", 0},
  {"cool_tot", "erg/s", "cool_tot", 0},
  {"cool_comp", "erg/s", "cool_comp", 0},
  {"lum_lines", "erg/s", "lum_lines", 0},
  {"cool_dr", "erg/s", "cool_dr", 0},
  {"lum_ff", "erg/s", "lum_ff", 0},
  {"lum_rr", "erg/s", "lum_rr", 0},
  {"cool_rr", "erg/s", "cool_rr", 0},
  {"cool_adiab", "erg/s", "cool_adiab", 0},
  {"heat_shock", "erg/s", "heat_shock", 0},
  {"heat_lines_macro", "erg/s", "heat_lines_macro", 0},
  {"heat_photo_macro", "erg/s", "heat_photo_macro", 0},
  {"cool_lines_macro", "erg/s", "cool_lines_macro", 0},
  {"cool_bf_macro", "erg/s", "cool_bf_macro", 0},

  /* === CONVERGENCE TABLE VARIABLES === */
  {"t_e_old", "K", "t_e_old", 0},
  {"dt_e", "K", "dt_e", 0},
  {"dt_e_old", "K", "dt_e_old", 0},
  {"t_r_old", "K", "t_r_old", 0},
  {"heat_tot_old", "erg/s", "heat_tot_old", 0},
  {"gain", "", "gain", 0},
  {"macro_bf_in", "", "macro_bf_in", 0},
  {"macro_bf_out", "", "macro_bf_out", 0},

  /* === VELOCITY GRADIENT TABLE VARIABLES === */
  {"dv_x_dx", "1/s", "dv_x_dx", 0},
  {"dv_y_dx", "1/s", "dv_y_dx", 0},
  {"dv_z_dx", "1/s", "dv_z_dx", 0},
  {"dv_x_dy", "1/s", "dv_x_dy", 0},
  {"dv_y_dy", "1/s", "dv_y_dy", 0},
  {"dv_z_dy", "1/s", "dv_z_dy", 0},
  {"dv_x_dz", "1/s", "dv_x_dz", 0},
  {"dv_y_dz", "1/s", "dv_y_dz", 0},
  {"dv_z_dz", "1/s", "dv_z_dz", 0},
  {"div_v", "1/s", "div_v", 0},
  {"dvds_max", "1/s", "dvds_max", 0},
  {"gamma", "", "gamma", 0},
  {"dfudge", "cm", "dfudge", 0},
};

static int n_windtable_vars = sizeof (windtable_vars) / sizeof (windtable_vars[0]);


/**********************************************************/
/**
 * @brief  Return pointer to the variable registry array
 *
 * @param [out] count   Number of variables in the registry
 * @return              Pointer to the windtable_vars array
 *
 ***********************************************************/

windtable_var_t *
windtable_get_vars (int *count)
{
  *count = n_windtable_vars;
  return windtable_vars;
}


/**********************************************************/
/**
 * @brief  Generate ion column descriptors dynamically
 *
 * @param [out] ions    Pointer to allocated array of ion descriptors
 * @param [out] count   Number of ion columns generated
 * @return              0 on success, -1 on error
 *
 * @details
 *
 * This function generates ion fraction column descriptors based on
 * the atomic data that was loaded. Column names follow the pattern:
 *   {element}_i{state}
 *
 * For example:
 *   H_i0  = neutral hydrogen fraction
 *   H_i1  = H+ fraction
 *   He_i0 = neutral helium
 *   He_i1 = He+
 *   He_i2 = He++
 *   C_i0  = neutral carbon
 *   ...
 *   C_i6  = C6+ (fully ionized)
 *
 * The elements included depend on the atomic data file, typically:
 * H, He, C, N, O, Na, Si, Ca, Fe (Z = 1, 2, 6, 7, 8, 11, 14, 20, 26)
 *
 ***********************************************************/

int
windtable_generate_ion_columns (windtable_ion_t ** ions, int *count)
{
  int i, j;
  int first_ion, number_ions;
  int total_ions;
  windtable_ion_t *ion_array;

  /* First pass: count total number of ion columns needed */
  total_ions = 0;
  for (i = 0; i < nelements; i++)
  {
    total_ions += ele[i].nions;
  }

  if (total_ions == 0)
  {
    Error ("windtable: No ions found in atomic data\n");
    *ions = NULL;
    *count = 0;
    return -1;
  }

  /* Allocate array */
  ion_array = (windtable_ion_t *) calloc (total_ions, sizeof (windtable_ion_t));
  if (ion_array == NULL)
  {
    Error ("windtable: Could not allocate memory for ion columns\n");
    *ions = NULL;
    *count = 0;
    return -1;
  }

  /* Second pass: fill in ion descriptors */
  total_ions = 0;
  for (i = 0; i < nelements; i++)
  {
    first_ion = ele[i].firstion;
    number_ions = ele[i].nions;

    for (j = 0; j < number_ions; j++)
    {
      windtable_ion_t *v = &ion_array[total_ions];

      /* Column name: element_iN where N is ionization state - 1 */
      snprintf (v->name, WINDTABLE_VAR_NAME_LEN, "%s_i%d", ele[i].name, ion[first_ion + j].istate - 1);

      v->z = ele[i].z;
      v->istate = ion[first_ion + j].istate;

      total_ions++;
    }
  }

  Log ("windtable: Generated %d ion fraction columns\n", total_ions);

  *ions = ion_array;
  *count = total_ions;
  return 0;
}


/**********************************************************/
/**
 * @brief  Free ion column array allocated by windtable_generate_ion_columns
 *
 * @param [in] ions   Pointer to ion array (can be NULL)
 *
 ***********************************************************/

void
windtable_free_ion_columns (windtable_ion_t * ions)
{
  if (ions != NULL)
  {
    free (ions);
  }
}
