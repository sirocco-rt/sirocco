
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
 * We then map these models into the structures that Python uses.
 * Most of the mapping is one-to-one, but Python wants the densities
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
import_wind (ndom)
     int ndom;
{
  char filename[LINELENGTH];

  sprintf (filename, "%s", "e.g. foo.txt");

  rdstr ("Wind.model2import", filename);

  import_wind2 (ndom, filename);

  return (0);
}

int
import_wind2 (ndom, filename)
     int ndom;
     char *filename;
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
  else
  {
    Error ("%s : %i : Do not know how to import a model of coord_type %d\n", __FILE__, __LINE__, zdom[ndom].coord_type);
    Exit (0);
  }

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
import_set_wind_boundaries (ndom)
     int ndom;
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
  else
  {
    Error ("get_import_wind_params: unknown coord_type %d\n", zdom[ndom].coord_type);
    Exit (1);
  }

  return imported_model[ndom].init_temperature;
}




/**********************************************************/
/** 
 * @brief      Make the Python grid 
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
  else
  {
    Error ("import_wind: Do not know how to import a model of coord_type %d\n", zdom[ndom].coord_type);
    Exit (0);
  }

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
import_velocity (ndom, x, v)
     int ndom;
     double *x, *v;
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
 * In Python, we want map the grid points to the edges of wind cells,
 * but we expect the densities to be given at the centers of the cells.
 *
 **********************************************************/

double
import_rho (ndom, x)
     int ndom;
     double *x;
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
  else
  {
    Error ("import_temperature: unknown coord_type %d\n", zdom[ndom].coord_type);
    Exit (1);
  }

  return temperature;
}
