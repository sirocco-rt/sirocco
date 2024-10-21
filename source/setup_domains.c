
/***********************************************************/
/** @file   setup_domains.c
 * @author ksl
 * @date   January, 2018
 * @brief  Get the parameters needed to describe winds
 *
 * File containing several routines that collectively
 * define the components to a wind in sirocco. Each component
 * of the wind is said to be a domain, and the information
 * for each domain is stored in the elements of zdom
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/** 
 * @brief       Get inputs that describe a particular component of the wind
 *
 * @param [in] ndom  The number (begining with 0) of this particular domain
 * @return  0 
 *
 * Sets up the one domain, which includes defining wind type, e.g whether
 * it is a shell, or a biconical flow, or an imported model, as well
 * as the type of coordinate system and its simensions.
 *
 * If the wind is to be imported from a file, it is imported in this
 * routine.
 *
 * ###Notes###
 * 1801 -   Refactored into this file in 1801.  Updates in the
 *          fall of 17 were made to allow for importing models.
 *          Note that cyl_var coordinates are not currently 
 *          working
***********************************************************/


int
get_domain_params (ndom)
     int ndom;
{
  char answer[LINELENGTH];

  if (ndom >= geo.ndomain)
  {
    Error ("Trying to get grid params for a non-existent domain!\n");
    Exit (0);
  }

  strcpy (answer, "SV");

  if (geo.system_type == SYSTEM_TYPE_STAR)
  {
    strcpy (answer, "star");
  }

  zdom[ndom].wind_type = rdchoice ("Wind.type(SV,star,hydro,corona,kwd,homologous,shell,imported)", "0,1,3,4,5,6,9,10", answer);

  /* For a shell wind, we define the coordinate system in shell_wind */
  if (zdom[ndom].wind_type == SHELL)
  {
    return (0);
  }

  strcat (zdom[ndom].name, "Wind");

  /* Define the coordinate system for the grid and allocate memory for the wind structure */

  if (geo.system_type == SYSTEM_TYPE_STAR)
  {
    strcpy (answer, "spherical");
  }
  else
  {
    strcpy (answer, "cylindrical");
  }

  zdom[ndom].coord_type = rdchoice ("Wind.coord_system(spherical,cylindrical,polar,cyl_var)", "0,1,2,3", answer);

  if (zdom[ndom].wind_type == IMPORT)   // Do not define dimensions for imported model
  {
    import_wind (ndom);
  }
  else
  {
    rdint ("Wind.dim.in.x_or_r.direction", &zdom[ndom].ndim);
    if (zdom[ndom].ndim > NDIM_MAX)
    {
      Error
        ("get_domain_params: %d x-cells for domain %d is greater than the maximum allowed cells of %d. Try increasing NDIM_MAX in sirocco.h\n",
         zdom[ndom].ndim, ndom, NDIM_MAX);
      Exit (EXIT_FAILURE);
    }
    if (zdom[ndom].coord_type > SPHERICAL)
    {
      rdint ("Wind.dim.in.z_or_theta.direction", &zdom[ndom].mdim);
      if (zdom[ndom].mdim > NDIM_MAX)
      {
        Error
          ("get_domain_params: %d z-cells for domain %d is greater than the maximum allowed cells of %d. Try increasing NDIM_MAX in sirocco.h\n",
           zdom[ndom].mdim, ndom, NDIM_MAX);
        Exit (EXIT_FAILURE);
      }
      else if (zdom[ndom].mdim < 4)
      {
        Error ("sirocco: domain mdim must be at least 4 to allow for boundaries\n");
        Exit (EXIT_FAILURE);
      }
    }
    else
    {
      zdom[ndom].mdim = 1;
    }
  }

  allocate_domain_wind_coords (ndom);

  /* If we are in advanced then allow the user to modify scale lengths */
  if (modes.iadvanced)
  {
    strcpy (answer, "no");
    modes.adjust_grid = rdchoice ("@Diag.adjust_grid(yes,no)", "1,0", answer);

    if (modes.adjust_grid)
    {
      Log ("You have opted to adjust the grid scale lengths\n");
      rddoub ("@geo.xlog_scale", &zdom[ndom].xlog_scale);
      if (zdom[ndom].coord_type > SPHERICAL)
        rddoub ("@geo.zlog_scale", &zdom[ndom].zlog_scale);
    }
  }

  zdom[ndom].ndim2 = zdom[ndom].ndim * zdom[ndom].mdim;

  return (0);
}

/**********************************************************/
/**
 * @brief Allocate memory for the wind coordinate arrays and for special
 *        cases for certain coordinate systems.
 *
 * @param [in] ndom The domain number to allocate for
 *
 * @details
 *
 * Memory is allocated for the `wind_x`, `wind_z`, `wind_midx` and `wind_midz`
 * arrays. For a polar coordinate system an additional wind cone is allocated.
 * For cylindrical variable coordinate system `wind_z_var` and `wind_midz_var`
 * are allocated.
 *
 * If allocating fails for anything, the program will exit.
 *
***********************************************************/

void
allocate_domain_wind_coords (int ndom)
{
  /* Allocate x dimensions */
  zdom[ndom].wind_x = calloc (zdom[ndom].ndim, sizeof (double));
  zdom[ndom].wind_midx = calloc (zdom[ndom].ndim, sizeof (double));
  if (zdom[ndom].wind_x == NULL || zdom[ndom].wind_midx == NULL)
  {
    Error ("allocate_domain_wind_coords: Problem allocating memory for x-dim for domain %d\n", ndom);
    Exit (EXIT_FAILURE);
  }

  /* Allocate z dimensions */
  zdom[ndom].wind_z = calloc (zdom[ndom].mdim, sizeof (double));
  zdom[ndom].wind_midz = calloc (zdom[ndom].mdim, sizeof (double));
  if (zdom[ndom].wind_z == NULL || zdom[ndom].wind_midz == NULL)
  {
    Error ("allocate_domain_wind_coords: Problem allocating memory for z-dim for domain %d\n", ndom);
    Exit (EXIT_FAILURE);
  }

  /* allocate any special cases */
  if (zdom[ndom].coord_type == RTHETA)
  {
    zdom[ndom].cones_rtheta = calloc (zdom[ndom].mdim, sizeof (cone_dummy));
    if (zdom[ndom].cones_rtheta == NULL)
    {
      Error ("allocate_domain_wind_coords: Problem allocating memory for the r-theta cone in domain %d\n", ndom);
      Exit (EXIT_FAILURE);
    }
  }
  else if (zdom[ndom].coord_type == CYLVAR)
  {
    cylvar_allocate_domain (ndom);
  }
}

/**********************************************************/
/**   
 * @brief       Get detailed particular component of the wind
 *
 * @param [in] ndom  The number (begining with 0) of this particular domain
 * @return  0 
 *
 * Continues the setup of a single domain begun in get_domain_params.    
 *
 * Much of this routine is a steering routine that calls other
 * subroutines depending on the type of wind, e.g sv, for this
 * particular wind domain.
 *
 *
 * ###Notes###
 * 1801 -   Refactored into this file in 1801.  

***********************************************************/

int
get_wind_params (ndom)
     int ndom;
{
  int import_t_init = FALSE;

  /* Now get parameters that are specific to a given wind model

     Note: When one adds a new model, the only things that should be read in and modified
     are parameters in geo.  This is in order to preserve the ability to continue a calculation
     with the same basic wind geometry, without reading in all of the input parameters.  
   */

  zdom[ndom].rmax = 0;

  if (zdom[ndom].wind_type == STAR)
  {
    get_stellar_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == SV)
  {
    get_sv_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == HYDRO)
  {
    get_hydro_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == CORONA)
  {
    get_corona_params (ndom);
  }
  else if (zdom[ndom].wind_type == KNIGGE)
  {
    get_knigge_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == HOMOLOGOUS)
  {
    get_homologous_params (ndom);
  }
  else if (zdom[ndom].wind_type == SHELL)       //NSH 18/2/11 This is a new wind type to produce a thin shell.
  {
    get_shell_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == IMPORT)
  {
    /*
     * At this stage, we should set the boundaries for the wind imported wind
     * as this will avoid the situation where one has to provide Wind.radmax
     * or Wind.t.init when these have already been set by the data which
     * has been read in.
     *
     * Note there is a horrid spaghetti variable which is a flag for if the
     * temperature has been provided or not.
     */

    import_t_init = import_set_wind_boundaries (ndom);
  }
  else
  {
    Error ("get_wind_params(%i): unknown wind type %d\n", __LINE__, zdom[ndom].wind_type);
    Exit (0);
  }

  /* For many models, zdom[ndom].rmax is defined within the get_parameter file, e.g
   * for stellar wind, but for others zdom[ndom].rmax can be set independently of
   * the details of the model for thesse, we clean up at the end. Ultimately
   * we may want to push all of this into ghe various get_..._params routines
   * but for now we just check if zdom[ndom].rmax has been set already and if not
   * we ask a generic quesiton here.  Because we ultimately may want to change this
   * we throw an error at this point
   *
   * These are models known not to have zdom[ndom].rmax defined as part of the model
   *
   * sv, knigge, 
   *
   * Others currently have it defined within get_whatever_params
   *
   * stellar, homologous, shell, corona
   *
   * Note that we define the maximum wind radius elsewhere in Python, so we should
   * not check for the error condition here, as there is no reason to supply Wind.radmax
   * for an imported model.
   */

  Log ("Checking wind boundaries for domain %i\n", ndom);

  if (zdom[ndom].rmax == 0)
  {
    if (zdom[ndom].wind_type != SV || zdom[ndom].wind_type != KNIGGE)   // Not an error for these models, see above
      Error ("get_wind_params: zdom[ndom].rmax = 0 for wind type %d\n", zdom[ndom].wind_type);

    if (geo.system_type == SYSTEM_TYPE_AGN || geo.system_type == SYSTEM_TYPE_BH)
    {
      zdom[ndom].rmax = 50. * geo.rstar;
    }
    else
    {
      zdom[ndom].rmax = 1e12;
    }

    rddoub ("Wind.radmax(cm)", &zdom[ndom].rmax);
    if (zdom[ndom].rmax == 0)
    {
      Error ("get_wind_params: wind.radmax == 0. Should be non-zero number!\n");
      Exit (1);
    }
    if (zdom[ndom].rmax < 0)
    {
      Error ("get_wind_params: Wind.rmax %e < 0 -- multiplying by -1\n", zdom[ndom].rmax);
      zdom[ndom].rmax *= -1;
    }
  }

  /* Next lines are to assure that we have the largest possible value of the
   * sphere surrounding the system
   * JM 1710 -- if this is the first domain, then initialise geo.rmax see #305
   */
  if ((ndom == 0) || (zdom[ndom].rmax > geo.rmax))
  {
    geo.rmax = zdom[ndom].rmax;
  }
  geo.rmax_sq = geo.rmax * geo.rmax;

  if (zdom[ndom].rmax <= zdom[ndom].rmin)
  {
    Error ("get_wind_params: rmax (%10.4e) less than or equal to rmin %10.4e in domain %d\n", zdom[ndom].rmax, zdom[ndom].rmin, ndom);
    Exit (1);
  }

  /*
   * Only query for Wind.t.init if the model is not imported, or if it is an
   * imported model and no temperature has been provided with the model
   */

  zdom[ndom].twind = 40000;
  if (zdom[ndom].wind_type != IMPORT || (zdom[ndom].wind_type == IMPORT && import_t_init))
  {
    rddoub ("Wind.t.init", &zdom[ndom].twind);
    if (zdom[ndom].twind == 0)
    {
      Error ("get_wind_params: Wind.t.init == 0. Increase the temperature!\n", zdom[ndom].twind);
      Exit (1);
    }
    if (zdom[ndom].twind < 0)
    {
      Error ("get_wind_params: Wind.t.init %e < 0 -- multiplying by -1\n", zdom[ndom].twind);
      zdom[ndom].twind *= -1;
    }
  }

  /* Get the filling factor of the wind */
  /* JM 1606 -- the filling factor is now specified on a domain by domain basis. See #212
     XXX allows any domain to be allowed a filling factor but this should be modified when
     we know what we are doing with inputs for multiple domains. Could create confusion */

  zdom[ndom].fill = 1.;
  rddoub ("Wind.filling_factor(1=smooth,<1=clumped)", &zdom[ndom].fill);
  if (zdom[ndom].fill > 1)
  {
    Error ("get_wind_params: filling factor f = %e > 1\n", zdom[ndom].fill);
    Exit (1);
  }
  if (zdom[ndom].fill < 0)
  {
    Error ("get_wind_params: negative filling factor f = %e -- multiplying by -1\n", zdom[ndom].fill);
    zdom[ndom].fill *= -1;
  }

  return (0);
}




/**********************************************************/
/** 
 * @brief      sets up the windcones for each domain
 *
 * @return     Always returns 0
 *
 * @details
 * 
 * This routine takes input variables, the minimim and maximum
 * radius of the wind at the disk, and the flow angles bounding
 * the wind, and calculates the variables used to define the cones
 * that bound each wind domain involving biconical flows.
 *
 * ### Notes ###
 *
 * The routine cycles through all of the existing domains, and
 * uses variables which have been read in or entered previously.
 *
 * The input variables that are used are typicall, wind_rhomin_at_disk, wind_rhomax_at_disk
 * and wind_thetamin and max.  They are defined in routines like,
 * get_sv_parameters.
 *
 * The angles thetamin and
 * thetamax are all defined from the z axis, so that an angle of 0
 * is a flow that is perpindicular to to the disk and one that is
 * close to 90 degrees will be parallel to the plane of the disk
 *
 * The routine files a structure of type cone, that contain, an element
 * z,  the place where the windcone intercepts the z axis, and an element
 * dzdr is the slope of the cone
 *
 * Wind cones are defined azimuthally around the z axis.
 *
 * For models that are not biconical flows, the windcones are set 
 * to include the entire domain.
 * 
 *
 **********************************************************/

int
setup_windcone ()
{
  int ndom;
  double dzdr, z;

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {

    dzdr = 1. / tan (zdom[ndom].wind_thetamin);
    z = zdisk (zdom[ndom].wind_rhomin_at_disk);
    init_windcone (zdom[ndom].wind_rhomin_at_disk, z, dzdr, FALSE, &zdom[ndom].windcone[0]);

    dzdr = 1. / tan (zdom[ndom].wind_thetamax);
    z = zdisk (zdom[ndom].wind_rhomax_at_disk);
    init_windcone (zdom[ndom].wind_rhomax_at_disk, z, dzdr, FALSE, &zdom[ndom].windcone[1]);

  }
  return (0);
}




/***********************************************************/
/** 
 * @brief      initialize a windcone structure 
 * 
 * @param [in] double r  rho position where the windcone structure is being defined
 * @param [in] double z  the z position (not necessarily zero where the cone 
 *                     being defined.
 * @param [in] double dzdr the slope of the cone measured by dzdr
 * @param [in] int allow_negative_dzdr  A switch which allows windcones to be defined
 *                     with negative dz_dr, that is to say which causes the cone
 *                     to intersect the z axis above z 
 * @param [out] ConePtr one_windcone  A pointer to the wind cone structe being 
 *                       initilized.
 *
 * @return     Returns 0 for the normal case, 1 in special cases where the
 *             cone the z axis at a posiition
 *             greater than z.  
 *
 * @details
 *
 * windcones are surfaces of revolution, cones that are defined by where they enconter
 * the z axis and a slope.  They are used to set boundaries of certain types of
 * domains and in setting up certain types of coordinate systems.
 * 
 * ### Notes ###
 *
 * The allow_nagative_dz_dr is used because their are circumstances 
 * e.g for a KWD model, where one really does not want to define 
 * a cone that converges to a point along the positive z axix. (Very
 * likely this will be an error if a windcone of this type is created..)
 * 
 *
 **********************************************************/


int
init_windcone (r, z, dzdr, allow_negative_dzdr, one_windcone)
     double r, z, dzdr;
     int allow_negative_dzdr;
     ConePtr one_windcone;

{


  if (dzdr == 0 || (allow_negative_dzdr == FALSE && dzdr < 0))
  {
    Log ("Why am I here\n");
    one_windcone->dzdr = VERY_BIG;
    one_windcone->z = -VERY_BIG;
    if (dzdr == 0)
      return (0);
    else
      return (1);

  }

  one_windcone->z = z - dzdr * r;
  one_windcone->dzdr = dzdr;

  return (0);


}
