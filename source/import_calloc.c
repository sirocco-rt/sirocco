/* ************************************************************************** */
/**
 * @file  import_calloc.c
 * @author EJP
 * @date   Feb 2019
 *
 * @brief    Functions for memory management for the import structures.
 *
 * ************************************************************************** */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/* ************************************************************************** */
/**
 * @brief   Allocate memory for the import_model structures.
 *
 * @param[in]   int coord_type     The coordinate system of the imported grid
 *
 * @details
 *
 * This function is simply a wrapper around a lot of calloc () calls to allocate
 * memory for the different arrays in the structure. It was opted to do this way
 * rather than static memory to avoid storing the raw imported grid details
 * when it isn't required.
 *
 * NDIM_MAX (or NDIM_MAX2D) elements are allocated for each array, as
 * at this point we do not know how many elements comprise the grid. Note that
 *
 * ************************************************************************** */

void
calloc_import (int coord_type, int ndom)
{
  if (imported_model == NULL)
  {
    imported_model = calloc (MAX_DOM, sizeof *imported_model);
    if (imported_model == NULL)
    {
      Error ("calloc_import: unable to allocate %i domains for imported_model\n", MAX_DOM);
      Exit (1);
    }
  }

  imported_model[ndom].init_temperature = FALSE;

  if (coord_type == SPHERICAL)
  {
    imported_model[ndom].i = calloc (sizeof *imported_model[ndom].i, NDIM_MAX);
    imported_model[ndom].inwind = calloc (sizeof *imported_model[ndom].inwind, NDIM_MAX);
    imported_model[ndom].r = calloc (sizeof *imported_model[ndom].r, NDIM_MAX);
    imported_model[ndom].v_r = calloc (sizeof *imported_model[ndom].v_r, NDIM_MAX);
    imported_model[ndom].mass_rho = calloc (sizeof *imported_model[ndom].mass_rho, NDIM_MAX);
    imported_model[ndom].t_r = calloc (sizeof *imported_model[ndom].t_r, NDIM_MAX);
    imported_model[ndom].t_e = calloc (sizeof *imported_model[ndom].t_e, NDIM_MAX);
  }
  else if (coord_type == CYLIND || coord_type == RTHETA)
  {
    imported_model[ndom].i = calloc (sizeof *imported_model[ndom].i, NDIM_MAX2D);
    imported_model[ndom].j = calloc (sizeof *imported_model[ndom].j, NDIM_MAX2D);
    imported_model[ndom].inwind = calloc (sizeof *imported_model[ndom].inwind, NDIM_MAX2D);
    imported_model[ndom].v_x = calloc (sizeof *imported_model[ndom].v_x, NDIM_MAX2D);
    imported_model[ndom].v_y = calloc (sizeof *imported_model[ndom].v_y, NDIM_MAX2D);
    imported_model[ndom].v_z = calloc (sizeof *imported_model[ndom].v_z, NDIM_MAX2D);
    imported_model[ndom].mass_rho = calloc (sizeof *imported_model[ndom].mass_rho, NDIM_MAX2D);
    imported_model[ndom].t_r = calloc (sizeof *imported_model[ndom].t_r, NDIM_MAX2D);
    imported_model[ndom].t_e = calloc (sizeof *imported_model[ndom].t_e, NDIM_MAX2D);
    imported_model[ndom].wind_x = calloc (sizeof *imported_model[ndom].wind_x, NDIM_MAX2D);
    imported_model[ndom].wind_z = calloc (sizeof *imported_model[ndom].wind_z, NDIM_MAX2D);
    imported_model[ndom].wind_midx = calloc (sizeof *imported_model[ndom].wind_midx, NDIM_MAX2D);
    imported_model[ndom].wind_midz = calloc (sizeof *imported_model[ndom].wind_midz, NDIM_MAX2D);

    if (coord_type == CYLIND)
    {
      imported_model[ndom].x = calloc (sizeof *imported_model[ndom].x, NDIM_MAX2D);
      imported_model[ndom].z = calloc (sizeof *imported_model[ndom].z, NDIM_MAX2D);
    }
    else
    {
      imported_model[ndom].r = calloc (sizeof *imported_model[ndom].r, NDIM_MAX2D);
      imported_model[ndom].theta = calloc (sizeof *imported_model[ndom].theta, NDIM_MAX2D);
    }
  }
  else
  {
    Error ("calloc_import: Unknown coord_type %i\n", coord_type);
    Exit (1);
  }
}



/* ************************************************************************** */
/**
 * @brief    Free the memory allocated for the import structures.
 *
 * @param[in]   int coord_type     The coordinate system of the imported grid
 *
 * @details
 *
 * When the raw imported grid is no longer needed as the relevant data has been
 * assigned to the wind or domain structure, then we should free the memory
 * so as not to have extra data hanging around.
 *
 * ************************************************************************** */

void
free_import (int coord_type, int ndom)
{
  if (coord_type == SPHERICAL)
  {
    free (imported_model[ndom].i);
    free (imported_model[ndom].inwind);
    free (imported_model[ndom].r);
    free (imported_model[ndom].v_r);
    free (imported_model[ndom].mass_rho);
    free (imported_model[ndom].t_r);
    free (imported_model[ndom].t_e);
  }
  else if (coord_type == CYLIND || coord_type == RTHETA)
  {
    free (imported_model[ndom].i);
    free (imported_model[ndom].j);
    free (imported_model[ndom].inwind);
    free (imported_model[ndom].v_x);
    free (imported_model[ndom].v_y);
    free (imported_model[ndom].v_z);
    free (imported_model[ndom].mass_rho);
    free (imported_model[ndom].t_r);
    free (imported_model[ndom].t_e);
    free (imported_model[ndom].wind_x);
    free (imported_model[ndom].wind_z);
    free (imported_model[ndom].wind_midx);
    free (imported_model[ndom].wind_midz);
    free (imported_model[ndom].x);
    free (imported_model[ndom].z);
    free (imported_model[ndom].r);
    free (imported_model[ndom].theta);
  }
  else
  {
    Error ("free_import: Unknown coord_type %i\n", coord_type);
    Exit (1);
  }

  free (imported_model);
}
