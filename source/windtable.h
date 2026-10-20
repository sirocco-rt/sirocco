/***********************************************************/
/** @file  windtable.h
 * @author Long, Matthews
 * @date   January 2026
 *
 * @brief  Header for windtable program - exports wind data to FITS binary table
 *
 * This file defines the data structures used for the table-driven
 * variable registry that makes adding new output variables trivial.
 *
 ***********************************************************/

#ifndef WINDTABLE_H
#define WINDTABLE_H

/* Note: Include sirocco.h before this header - it provides atomic.h definitions */

#define WINDTABLE_VAR_NAME_LEN  32
#define WINDTABLE_VAR_UNIT_LEN  16

/**
 * @brief Variable descriptor for windtable columns
 *
 * Each entry in the variable registry describes one output column.
 * To add a new variable, add one entry to the windtable_vars[] array
 * in windtable_vars.c - no other code changes required (unless
 * the variable isn't already handled by get_one()).
 */
typedef struct windtable_var
{
  char name[WINDTABLE_VAR_NAME_LEN];        /**< Column name (FITS TTYPE) */
  char unit[WINDTABLE_VAR_UNIT_LEN];        /**< Physical unit (FITS TUNIT), "" if dimensionless */
  char get_one_name[WINDTABLE_VAR_NAME_LEN]; /**< Variable name passed to get_one() */
  int is_int;                                /**< 1 if integer type (TINT), 0 for double (TDOUBLE) */
} windtable_var_t;

/**
 * @brief Ion column descriptor (generated dynamically)
 */
typedef struct windtable_ion
{
  char name[WINDTABLE_VAR_NAME_LEN];  /**< Column name, e.g. "H_i0", "C_i4" */
  int z;                               /**< Atomic number */
  int istate;                          /**< Ionization state */
} windtable_ion_t;

/* Function prototypes - windtable_vars.c */
windtable_var_t *windtable_get_vars (int *count);
int windtable_generate_ion_columns (windtable_ion_t ** ions, int *count);
void windtable_free_ion_columns (windtable_ion_t * ions);

/* Function prototypes - windtable.c */
int windtable_write_fits (int ndom, char *rootname, windtable_var_t * vars,
                          int nvars, windtable_ion_t * ions, int nions);

#endif /* WINDTABLE_H */
