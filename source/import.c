
/***********************************************************/
/** @file  import.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief   general purpose routines reading in model
 * grids
 *
 * The routines contained here are basically steering
 * routines. The real works is done in import_spherical,
 * etc
 *
 * For importing models, we first read in the data from a file.
 * We assume all of the data, positions, velocities and importantly
 * densities are given at the grid points of the imported model.
 *
 * We then map these models into the structures that Sirocco uses.
 * Most of the mapping is one-to-one, but Sirocco wants the densities
 * to be a the cell centers and not at the corners.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/** 
 * @brief      Import a gridded model
 *
 * @param [in] int  ndom   The domain for the model
 * @return   Always returns 0  
 *
 * @details
 *
 * This is a steering routine.  It reads the name of the file
 * to import and depending on the pre-established coordinate
 * system calls one of several coordinate system specific
 * routines to actually read in the model
 *
 * ### Notes ###
 *
 **********************************************************/

int
import_wind (int ndom)
{
  char filename[LINELENGTH];

  sprintf (filename, "%s", "e.g. foo.txt");

  rdstr ("Wind.model2import", filename);

  import_wind2 (ndom, filename);

  return (0);
}

int
import_wind2 (int ndom, char *filename)
{

  calloc_import (zdom[ndom].coord_type, ndom);

  if (zdom[ndom].coord_type == SPHERICAL)
  {
    import_1d (ndom, filename);
  }
  else if (zdom[ndom].coord_type == CYLIND)
  {
    import_cylindrical (ndom, filename);
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {
    import_rtheta (ndom, filename);
  }
  else if (zdom[ndom].coord_type == CYLIND3D)
  {
    import_cylindrical3d (ndom, filename);
  }
  else if (zdom[ndom].coord_type == SPH3D)
  {
    import_sph3d (ndom, filename);
  }
  else
  {
    Error ("%s : %i : Do not know how to import a model of coord_type %d\n", __FILE__, __LINE__, zdom[ndom].coord_type);
    Exit (0);
  }

  if (zdom[ndom].coord_type == CYLIND3D || zdom[ndom].coord_type == SPH3D)
    Log ("The imported model for domain %i has dimensions %d x %d x %d\n",
         ndom, imported_model[ndom].ndim, imported_model[ndom].mdim, imported_model[ndom].pdim);
  else
    Log ("The imported model for domain %i has dimensions %d x %d\n", ndom, imported_model[ndom].ndim, imported_model[ndom].mdim);

  return (0);
}


/**********************************************************/
/**
 * @brief      get parameters associated with an imported wind model
 *
 * @param[in] int  ndom               The domain associated with an imported
 *                                    model
 * @return    int  init_temperature   Will return TRUE if the temperature is to
 *                                    be set to the init temperature for this
 *                                    domain
 *
 * @details
 * At present this routine is simply a place holder and simply returns
 * to the call in setup_domains. At the moment, there is a flag returned which
 * indicates if a temperature has been provided with the imported model
 * data.
 *
 * ### Notes ###
 *
 **********************************************************/

int
import_set_wind_boundaries (int ndom)
{
  if (zdom[ndom].coord_type == SPHERICAL)
  {
    import_spherical_setup_boundaries (ndom);
  }
  else if (zdom[ndom].coord_type == CYLIND)
  {
    import_cylindrical_setup_boundaries (ndom);
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {
    import_rtheta_setup_boundaries (ndom);
  }
  else if (zdom[ndom].coord_type == CYLIND3D)
  {
    import_cylindrical3d_setup_boundaries (ndom);
  }
  else if (zdom[ndom].coord_type == SPH3D)
  {
    import_sph3d_setup_boundaries (ndom);
  }
  else
  {
    Error ("get_import_wind_params: unknown coord_type %d\n", zdom[ndom].coord_type);
    Exit (1);
  }

  return imported_model[ndom].init_temperature;
}




/**********************************************************/
/**
 * @brief  Check that boundary slices of an imported grid carry guard cells.
 *
 * @param [in] ndom   Domain index
 * @param [in] w      Wind array (cells already populated by make_grid_import)
 * @return  void — calls Exit(1) if any boundary violation is found
 *
 * @details
 * Guard cells (inwind == W_IGNORE) must line the outer edges of the grid so
 * that photons exiting the wind pass through a transparent cell rather than
 * leaving directly from an active wind cell.  Missing guard cells cause
 * photon trapping or incorrect path lengths near the boundary.
 *
 * After the import make_grid functions run, user-supplied W_NOT_INWIND cells
 * have been converted to W_IGNORE.  This routine checks that the outermost
 * boundary slice in each required direction contains no W_ALL_INWIND cells.
 * All violations are reported before the program exits.
 *
 * Boundaries checked by coordinate system:
 *   SPHERICAL  -- skipped (import_1d forces guard cells programmatically)
 *   CYLIND     -- outer rho (i=ndim-1), outer z (j=mdim-1)
 *   RTHETA     -- outer r  (i=ndim-1), pole (j=0)
 *   CYLIND3D   -- outer rho (i=ndim-1), outer z (j=mdim-1)
 *   SPH3D      -- outer r  (i=ndim-1), pole (j=0)
 * The phi direction never requires guard cells.
 *
 **********************************************************/

static void
check_import_guard_cells (int ndom, WindPtr w)
{
  int i, j, k;
  int nstart, ndim, mdim, pdim;
  int nbad, nbad_total;
  DomainPtr one_dom;

  one_dom = &zdom[ndom];
  if (one_dom->coord_type == SPHERICAL)
    return;

  nstart = one_dom->nstart;
  ndim = one_dom->ndim;
  mdim = one_dom->mdim;
  pdim = one_dom->pdim;
  nbad_total = 0;

  /* Outer radial boundary: i = ndim-1.
   * 2D index: nstart + i*mdim + j
   * 3D index: nstart + i*mdim*pdim + j*pdim + k          */
  nbad = 0;
  i = ndim - 1;
  if (one_dom->coord_type == CYLIND || one_dom->coord_type == RTHETA)
  {
    for (j = 0; j < mdim; j++)
      if (w[nstart + i * mdim + j].inwind == W_ALL_INWIND)
        nbad++;
    if (nbad)
      Error ("check_import_guard_cells: domain %d outer radial boundary (i=%d) has %d/%d active cells; guard cells (inwind=-1) required\n",
             ndom, i, nbad, mdim);
  }
  else
  {
    for (j = 0; j < mdim; j++)
      for (k = 0; k < pdim; k++)
        if (w[nstart + i * mdim * pdim + j * pdim + k].inwind == W_ALL_INWIND)
          nbad++;
    if (nbad)
      Error ("check_import_guard_cells: domain %d outer radial boundary (i=%d) has %d/%d active cells; guard cells (inwind=-1) required\n",
             ndom, i, nbad, mdim * pdim);
  }
  nbad_total += nbad;

  /* Outer z boundary: j = mdim-1  (CYLIND and CYLIND3D only) */
  if (one_dom->coord_type == CYLIND || one_dom->coord_type == CYLIND3D)
  {
    nbad = 0;
    j = mdim - 1;
    if (one_dom->coord_type == CYLIND)
    {
      for (i = 0; i < ndim; i++)
        if (w[nstart + i * mdim + j].inwind == W_ALL_INWIND)
          nbad++;
      if (nbad)
        Error ("check_import_guard_cells: domain %d outer z boundary (j=%d) has %d/%d active cells; guard cells (inwind=-1) required\n",
               ndom, j, nbad, ndim);
    }
    else
    {
      for (i = 0; i < ndim; i++)
        for (k = 0; k < pdim; k++)
          if (w[nstart + i * mdim * pdim + j * pdim + k].inwind == W_ALL_INWIND)
            nbad++;
      if (nbad)
        Error ("check_import_guard_cells: domain %d outer z boundary (j=%d) has %d/%d active cells; guard cells (inwind=-1) required\n",
               ndom, j, nbad, ndim * pdim);
    }
    nbad_total += nbad;
  }

  /* Polar boundary: j = 0  (RTHETA and SPH3D only) */
  if (one_dom->coord_type == RTHETA || one_dom->coord_type == SPH3D)
  {
    nbad = 0;
    j = 0;
    if (one_dom->coord_type == RTHETA)
    {
      for (i = 0; i < ndim; i++)
        if (w[nstart + i * mdim + j].inwind == W_ALL_INWIND)
          nbad++;
      if (nbad)
        Error ("check_import_guard_cells: domain %d polar boundary (j=0, theta~0) has %d/%d active cells; guard cells (inwind=-1) required\n",
               ndom, nbad, ndim);
    }
    else
    {
      for (i = 0; i < ndim; i++)
        for (k = 0; k < pdim; k++)
          if (w[nstart + i * mdim * pdim + j * pdim + k].inwind == W_ALL_INWIND)
            nbad++;
      if (nbad)
        Error ("check_import_guard_cells: domain %d polar boundary (j=0, theta~0) has %d/%d active cells; guard cells (inwind=-1) required\n",
               ndom, nbad, ndim * pdim);
    }
    nbad_total += nbad;
  }

  if (nbad_total > 0)
  {
    Error ("check_import_guard_cells: domain %d has %d boundary cells missing guard status -- import file must supply inwind=-1 at all outer edges\n",
           ndom, nbad_total);
    Exit (1);
  }
}


/**********************************************************/
/**
 * @brief      Make the Sirocco grid
 *
 * @param [in] WindPtr  w  The entire wind structure
 * @param [in] int  ndom   The domain for the imported model
 * @return     Always returns 0
 *
 * @details
 * This is merely a steering routine for calling one
 * of the coordinate-system specific routines for creating
 * a grid from one of the imported models
 *
 * ### Notes ###
 * The fact that w is provided to this routine is for consistency.
 *
 **********************************************************/

int
import_make_grid (int ndom, WindPtr w)
{
  if (zdom[ndom].coord_type == SPHERICAL)
  {
    spherical_make_grid_import (w, ndom);
  }
  else if (zdom[ndom].coord_type == CYLIND)
  {
    cylindrical_make_grid_import (w, ndom);
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {
    rtheta_make_grid_import (w, ndom);
  }
  else if (zdom[ndom].coord_type == CYLIND3D)
  {
    cylindrical3d_make_grid_import (w, ndom);
  }
  else if (zdom[ndom].coord_type == SPH3D)
  {
    sph3d_make_grid_import (w, ndom);
  }
  else
  {
    Error ("import_wind: Do not know how to import a model of coord_type %d\n", zdom[ndom].coord_type);
    Exit (0);
  }

  check_import_guard_cells (ndom, w);

  return (0);
}




/**********************************************************/
/** 
 * @brief      Get the velocity at a position for an imported model
 *
 * @param [in] int  ndom   The domain in which the velocity is to be determined
 * @param [in] double *  x   The position where the velocity is to be determined
 * @param [out] double *  v   The velocity in cartesian coordiantes
 * @return     The speed       
 *
 * @details
 * The routine simply calls one of several coordinate-system specific routines
 * for obtaining the velocity of the wind when an imported model is involved.
 *
 * ### Notes ###
 * The routine is used to set up velocities in wmain
 *
 **********************************************************/

double
import_velocity (int ndom, double *x, double *v)
{
  double speed = 0.0;

  if (zdom[ndom].coord_type == SPHERICAL)
  {
    speed = velocity_1d (ndom, x, v);
  }
  else if (zdom[ndom].coord_type == CYLIND)
  {
    speed = velocity_cylindrical (ndom, x, v);
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {
    speed = velocity_rtheta (ndom, x, v);
  }
  else if (zdom[ndom].coord_type == CYLIND3D)
  {
    speed = velocity_cylindrical3d (ndom, x, v);
  }
  else if (zdom[ndom].coord_type == SPH3D)
  {
    speed = velocity_sph3d (ndom, x, v);
  }
  else
  {
    Error ("import_velocity: unknown coord_type %d\n", zdom[ndom].coord_type);
    Exit (1);
  }

  return (speed);
}




/**********************************************************/
/** 
 * @brief      Get the density at an arbitrary position in an imported
 * model
 *
 * @param [in] int  ndom   The domain associated with an imported model
 * @param [in] double *  x   The position where we desire rho
 * @return     rho              
 *
 * @details
 * The routine simply calls coordinate system specific routines to get the
 * density for imported models
 *
 * ### Notes ###
 * The routine is used to map densities from an imported model 
 * where we assume that the density is given at the grid points.
 * In Sirocco, we want map the grid points to the edges of wind cells,
 * but we expect the densities to be given at the centers of the cells.
 *
 **********************************************************/

double
import_rho (int ndom, double *x)
{
  double rho = 0.0;

  if (zdom[ndom].coord_type == SPHERICAL)
  {
    rho = rho_1d (ndom, x);
  }
  else if (zdom[ndom].coord_type == CYLIND)
  {
    rho = rho_cylindrical (ndom, x);
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {
    rho = rho_rtheta (ndom, x);
  }
  else if (zdom[ndom].coord_type == CYLIND3D)
  {
    rho = rho_cylindrical3d (ndom, x);
  }
  else if (zdom[ndom].coord_type == SPH3D)
  {
    rho = rho_sph3d (ndom, x);
  }
  else
  {
    Error ("import_rho: unknown coord_type %d\n", zdom[ndom].coord_type);
    Exit (1);
  }

  return (rho);
}




/* ************************************************************************** */
/**
 * @brief  Get the temperature of an imported model at the given position x.
 *
 * @param[in]    int ndom       The domain of interest
 *
 * @param[in]    double x[3]    The position of interest
 *
 * @return       t_r            The radiation temperature at the position x
 *
 * @details
 *
 * The purpose of this function is to simply look up the temperature at a given
 * grid cell for an imported wind model. In some cases this will be the
 * temperature which is given by the model, however, we also allow one to not
 * provide a cell temperature. In these cases, a default temperature value is
 * used which is set to the domain's initial wind temperature.
 *
 * ************************************************************************** */

double
import_temperature (int ndom, double *x, int return_t_e)
{
  double temperature = 0;

  if (zdom[ndom].coord_type == SPHERICAL)
  {
    temperature = temperature_1d (ndom, x, return_t_e);
  }
  else if (zdom[ndom].coord_type == CYLIND)
  {
    temperature = temperature_cylindrical (ndom, x, return_t_e);
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {
    temperature = temperature_rtheta (ndom, x, return_t_e);
  }
  else if (zdom[ndom].coord_type == CYLIND3D)
  {
    temperature = temperature_cylindrical3d (ndom, x, return_t_e);
  }
  else if (zdom[ndom].coord_type == SPH3D)
  {
    temperature = temperature_sph3d (ndom, x, return_t_e);
  }
  else
  {
    Error ("import_temperature: unknown coord_type %d\n", zdom[ndom].coord_type);
    Exit (1);
  }

  return temperature;
}
