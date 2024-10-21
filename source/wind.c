
/***********************************************************/
/** @file  wind.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  Except for where_in_wind, these are convenience routines
 * which make it easier to incorporate wind models into Python.  
 *
 * where_in_wind determines whether a position is x is in the wind
 * and if so in what domain.  Since grid can overlap, the Python
 * convention is that the domains are overlaid on one another, and
 * the last one is the one that matters at a particular postion.  
 * where_in_wind is therefore one of the basic utilites to determine
 * where one is at any time.
 *
 * In Python, the actual radiative transfer is carried on various
 * grids (domains).  During radiative transfer quatnities like
 * velocity and density are obtained by interpolating on values
 * in the grid, and not from the formulae that defince for, e. g.
 * a SV model.  The routines, here model_velocity, model_vgrad,
 * model_rho, etc are used to populate the grid.
 *
 *
 * ### Notes ###
 *
 * The file also contains a routine wind_check that probably
 * could be removed altogether, as it has not revealed a problem
 * in a very long time.  
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/** 
 * @brief      determines whether a position x is in a wind region.
 *
 * @param [in] double  x[]    a position
 * @param [out] int *  ndomain   the wind domain that the position lies in, if it does
 * @return   W_ALL_INWIND if the position is in a wind, W_IN_DISK if inside a vertically
 * extended disk, W_IN_STAR if inside the star, W_NOT_INWIND otherwise 
 * 
 *
 * @details
 *
 * where_in_wind searches the wind domains in reverse order to decide whether a position
 * in a particular wind domain or not.  For most types of wind, where the edges of
 * the wind are defined by wind_cones and by min and max radii etc.  For imported models,
 * where the wind_cones are insufficient, we check whether the position is in a cell
 * that is in the wind.
 *
 *
 * ### Notes ###
 * Since domains can overlap in principle, e.g one can have say a bipolar wind which starts
 * at the disk surface and goes to a very large radius, and then a coronal just above the disk,
 * one needs a convention to know which domain to use when there are two possibilities.  The
 * Python convention is to lay the grids on top of one another so that the last one matters when 
 * there is a conflict. For this reason the grids are searched in reverse order.
 *
 * Where_in_wind does not tell you whether a position is in the grid or not, just whether 
 * he photon is between the two wind cones and/or inside minimum or maximum radius of the wind.     
 * 
 *
 * The routine first checks to see if the photon is inside the disk or the star. This can
 * happen because transphot pushes into a boundary by a small difference so that the
 * routine walls operates properly
 *
 **********************************************************/

int
where_in_wind (x, ndomain)
     double x[];
     int *ndomain;
{
  double rho, rad, z;
  int ireturn;
  int ndom, n;
  DomainPtr one_dom;


  /* First check for items, like the disk that are not domain
   * related, like the disk.  We check for the star first
   * because the way we have defined the vertically extended
   * disk the height does not go to zero at rstar */

  if (length (x) < geo.rstar)
  {
    *ndomain = -1;
    return (W_IN_STAR);
  }


  z = fabs (x[2]);              /* Necessary to get correct answer above
                                   and below plane */
  rho = sqrt (x[0] * x[0] + x[1] * x[1]);       /* This is distance from z axis */


  /* Check if position is inside the disk for a vertically extended disk */
  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    if (rho < geo.disk_rad_max && z < zdisk (rho))
    {
      *ndomain = -1;
      return (W_IN_DISK);
    }
  }

  /* Now check whether position is a wind region of any of the domains.
   * This is done in reverse order on the assumption that our domains
   * are layered on top of one another.  */

  ireturn = W_NOT_INWIND;
  *ndomain = -1;

  rad = length (x);

  for (ndom = geo.ndomain - 1; ndom > -1; ndom--)
  {

    one_dom = &zdom[ndom];

    /* First check to see if photon is inside or outside wind */

    if (rad < one_dom->rmin)
    {
      continue;                 /*x is inside the wind  radially */
    }
    if (rad > one_dom->rmax)
    {
      continue;                 /*the position is beyond the wind radially */
    }

    if (z > one_dom->zmax)
    {
      continue;                 /*the position is beyond the wind radially */
    }

    /* Check if one is inside the inner windcone */
    if (rho < (one_dom->wind_rhomin_at_disk + z * tan (one_dom->wind_thetamin)))
    {
      continue;
    }

    /* Finally check if positon is outside the outer windcone */
    /* NSH 130401 - The check below was taking a long time if geo.wind_thetamax was very close to pi/2.
       check inserted to simply return INWIND if geo.wind_thetamax is within machine precision of pi/2. */

    if (fabs (one_dom->wind_thetamax - PI / 2.0) > 1e-6)        /* Only perform the next check if thetamax is not equal to pi/2 */
    {
      if (rho > (one_dom->wind_rhomax_at_disk + z * tan (one_dom->wind_thetamax)))
      {
        continue;
      }

    }

    /* At this point global constraints (however poorly defined) have 
     * been applied to an arbitrary imported model, but we must still check 
     * whether this particular point is in the grid and whether that point is
     * in the wind or not.  We follow the usual practice of allowing the grid to
     * define whether it is in the grid or not.  We also do this for the
     * case where we want to exclude cells that are partially in the wind.
     */

    if (one_dom->wind_type == IMPORT)
    {
      n = where_in_grid (ndom, x);
      if (n >= 0)
      {
        *ndomain = ndom;
        ireturn = wmain[n].inwind;
        break;
      }
    }


    /* At this point we have passed all of the tests for being in the wind */

    *ndomain = ndom;
    ireturn = W_ALL_INWIND;
    break;

  }

  return (ireturn);
}




/**********************************************************/
/** 
 * @brief      Calculate the wind velocity at a specific point in space from the original 
 * usually analytic expressions
 *
 * @param [out] int  ndom   The domain of interest
 * @param [out] double  x[]   A position (nominally in the domain)
 * @param [out] double  v[]   The resulting 3 velicity 
 * @return  The speed at that position  
 *
 * @details
 * 
 * This is a convenience function that allows one to obtain the velocity at
 * any point in space for a given domain. The routine uses the domain number
 * to decide what velocity law is appropriate.  For models which are defined from 
 * a set of parametric equations the velocity will sill be calculated. For models
 * that are read in as a grid, the models will be interpolated from the imported
 * grid
 *
 * ### Notes ###
 * 
 * This routine is used to set up the grid of velocities at the corners of
 * grid cells.  It is not used later on.
 *
 * The routine works for imported models as well, even in the case
 * the model velocity is actually one of the inputs.
 *
 * This routine calculates velocities in the observer frame.  
 * 
 * The routine has a check to see that the velocities do not exceed
 * the speed of light; if they do the velocity is rescaled to 0.99
 * C
 **********************************************************/

double
model_velocity (ndom, x, v)
     double x[], v[];
     int ndom;
{
  double speed = 0;

  if (zdom[ndom].wind_type == SV)
  {
    speed = sv_velocity (x, v, ndom);
  }
  else if (zdom[ndom].wind_type == STAR)
  {
    speed = stellar_velocity (ndom, x, v);
  }
  else if (zdom[ndom].wind_type == HYDRO)
  {
    speed = hydro_velocity (ndom, x, v);
  }
  else if (zdom[ndom].wind_type == CORONA)
  {
    speed = corona_velocity (ndom, x, v);
  }
  else if (zdom[ndom].wind_type == KNIGGE)
  {
    speed = kn_velocity (ndom, x, v);
  }
  else if (zdom[ndom].wind_type == HOMOLOGOUS)
  {
    speed = homologous_velocity (ndom, x, v);
  }
  else if (zdom[ndom].wind_type == SHELL)
  {
    speed = stellar_velocity (ndom, x, v);
  }
  else if (zdom[ndom].wind_type == IMPORT)
  {
    speed = import_velocity (ndom, x, v);
  }
  else
  {
    Error ("wind: Unknown windtype %d for doman %d\n", zdom[ndom].wind_type, ndom);
    Exit (0);
  }

  if (speed > 0.99 * VLIGHT)
  {
    rescale (v, 0.99 * VLIGHT / speed, v);
  }

  return (speed);
}





/**********************************************************/
/** 
 * @brief      calculate  the co-moving frame velocity gradient at 
 * positions in the flow based on the analytic wind models and imported models.  
 *
 * @param [in]  int  ndom   The domain of interest
 * @param [in]  double  x[]   A position in the domain
 * @param [out] double  v_grad[][3]   The velocity gradient tensor at the position
 * @return   Always returns 0  
 *
 * @details
 * 
 * The routine calls model_velocity multiple times to calculate the velocity
 * gradient tensor in the co-moving frame at a particular position (in the
 * observer frame..
 *
 * ### Notes ###
 *
 * This routine is normally used only during the initialization
 * of the wind
 *
 * This uses a symmetric calculation of the derivative but does 
 * not adaptively adjust the length of the steps to check the
 * accuracy of the calculation.  See issue #782
 *
 * Since the routine calls model_velocity directly, and not vwind_xyz
 * it can be called before wmain.v has been populated.
 *
 **********************************************************/

int
model_vgrad (ndom, x, v_grad)
     double x[], v_grad[][3];
     int ndom;
{

  double v[3], v_forward[3], v_reverse[3];
  double dx_forward[3], dx_reverse[3];
  double dv[3];
  double ds_cmf;
  int i, j;
  double zero_vector[3], dx[3], dx_obs[3];

  zero_vector[0] = zero_vector[1] = zero_vector[2] = 0.0;




  ds_cmf = 0.00001 * length (x);
  if (ds_cmf == 0)
  {
    stuff_v (zero_vector, v_grad[0]);
    stuff_v (zero_vector, v_grad[1]);
    stuff_v (zero_vector, v_grad[2]);
    return (0);
  }

  if (ds_cmf < 1.e6)
    ds_cmf = 1.e6;

  model_velocity (ndom, x, v);

  for (i = 0; i < 3; i++)
  {
    /* first create vectors  which are offset by +-ds.  Note that
       we want the observer frame velocities at a point which is ds away
       in the cmf frame.  To do this, we have to find out how far away
       that point would be in the observer frame. */

    stuff_v (zero_vector, dx);
    dx[i] = ds_cmf;
    local_to_observer_frame_ruler_transform (v, dx, dx_obs);
    vadd (x, dx_obs, dx_forward);
    vsub (x, dx_obs, dx_reverse);

    /* calculate the velocity at these positions */
    model_velocity (ndom, dx_reverse, v_reverse);
    model_velocity (ndom, dx_forward, v_forward);



    observer_to_local_frame_velocity (v_forward, v_reverse, dv);


    for (j = 0; j < 3; j++)
      dv[j] /= 2 * ds_cmf;

    if (sane_check (dv[0]) || sane_check (dv[1]) || sane_check (dv[2]))
    {
      Error ("model_vgrad: x %12.4e %12.4e %12.4e dv %12.4e %12.4e %12.4e %12.4e\n", x[0], x[1], x[2], dv[0], dv[1], dv[2]);
    }


    stuff_v (dv, v_grad[i]);
  }

  return (0);


}



/**********************************************************/
/** 
 * @brief      calculate  the diverernce of  the velocity at positions 
 * in the flow based on  the analytic wind models and imported models.  
 *
 * @param [in] int  ndom   The domain of interest
 * @param [in] double  x[]   A position in the domain
 * @return   Always returns the divergence of the velocity in the co-moving
 * frame
 *
 * @details
 * 
 * The routine calls model_velocity multiple times to calculate the velocity
 * divergence at a particular postion
 *
 * ### Notes ###
 *
 * This routine is normally used only during the initialization
 * of the wind
 *
 * Since the routine calls model_velocity directly, and not vwind_xyz
 * it can be called before wmain.v has been populated.
 *
 **********************************************************/

double
get_div_v_in_cmf_frame (ndom, x)
     int ndom;
     double *x;
{
  int i;
  double v[3][3];
  double div_v = 0;

  model_vgrad (ndom, x, v);

  /* the trace of the velocity gradient tensor is the divergence */
  for (i = 0; i < 3; i++)
  {
    div_v += v[i][i];
  }

  return (div_v);
}



/**********************************************************/
/** 
 * @brief      calculate the density of the wind in from the flow based on the
 * analytic or imported wind models.
 *
 * @param [in] int  ndom   The domain of interest
 * @param [in] double  x[]   A position in the domain
 * @return     rho at the position
 *
 * @details
 * This routine simple checks the wind type, e. g. SV, for the domain, and then calls the
 * appropriate routine to that model, e.g sv_rho, to get the density
 *
 * ### Notes ###
 * This routine is used during the initialization
 * of the wind.
 *
 **********************************************************/

double
model_rho (ndom, x)
     int ndom;
     double x[];
{
  double rho = 0;
  int n = 0;


  if (modes.partial_cells == PC_ZERO_DEN)
  {
    n = where_in_grid (ndom, x);

    if (wmain[n].inwind != W_ALL_INWIND)
    {
      return (0);
    }
  }

  if (zdom[ndom].wind_type == SV)
  {
    rho = sv_rho (ndom, x);
  }
  else if (zdom[ndom].wind_type == STAR)
  {
    rho = stellar_rho (ndom, x);
  }
  else if (zdom[ndom].wind_type == HYDRO)
  {
    rho = hydro_rho (x);
  }
  else if (zdom[ndom].wind_type == CORONA)
  {
    rho = corona_rho (ndom, x);
  }
  else if (zdom[ndom].wind_type == KNIGGE)
  {
    rho = kn_rho (ndom, x);
  }
  else if (zdom[ndom].wind_type == HOMOLOGOUS)
  {
    rho = homologous_rho (ndom, x);
  }
  else if (zdom[ndom].wind_type == SHELL)
  {
    rho = stellar_rho (ndom, x);
  }
  else if (zdom[ndom].wind_type == IMPORT)
  {
    rho = import_rho (ndom, x);
  }
  else
  {
    Error ("wind2d: Unknown windtype %d for domain %d\n", zdom[ndom].wind_type, ndom);
    Exit (0);
  }

  return (rho);

}


/**********************************************************/
/** 
 * @brief      Simple checks of the wind structure for reasonability
 *
 * @param [in] WindPtr  www   The entire wind
 * @param [in] int  n   n >= 0  then an element of the array will be checked
 * @return     Always returns 0
 *
 * @details
 * The routine should be called in a number of ways
 * 
 * * wind_check(w,-1);  to check the entire structure
 * * wind_check(w,50);   to check element 50 of the structure
 * * wind_check(w[50],0) to check element 50 of the structure
 *
 * ### Notes ###
 * 
 * These checks are basic, just NaN (sane_checks) and a check that
 * the wind does not contain velocities that exceed the speed of light. 
 * 
 * The program will stop if any of the checks failed.
 *
 * The checks are made on the wind, without reference to domains
 *
 **********************************************************/

int
wind_check (www, n)
     WindPtr www;
     int n;
{
  int i, j, k, istart, istop;
  int ierr = 0;
  int ndom, ndim, mdim;
  double dxmin, dzmin;
  double drmin, dtmin;
  int outer_n, outer_m;
  double delta;
  double frac = 0.01;

  if (n < 0)
  {
    istart = 0;
    istop = NDIM2;
  }
  else
  {
    istart = n;
    istop = istart + 1;
  }

  for (i = istart; i < istop; i++)
  {
    if (length (www[i].v) > VLIGHT)
    {
      Error ("wind_check: greater than light speed velocity %e in wind element %d\n", length (www[i].v), i);
      ierr++;
    }
    for (j = 0; j < 3; j++)
    {
      if (sane_check (www[i].x[j]))
      {
        Error ("wind_check:sane_check www[%d].x[%d] %e\n", i, j, www[i].x[j]);
        ierr++;
      }
      if (sane_check (www[i].xcen[j]))
      {
        Error ("wind_check:sane_check www[%d].xcen[%d] %e\n", i, j, www[i].xcen[j]);
        ierr++;
      }
      if (sane_check (www[i].v[j]))
      {
        Error ("wind_check:sane_check www[%d].v[%d] %e\n", i, j, www[i].v[j]);
        ierr++;
      }
    }
    for (j = 0; j < 3; j++)
    {
      for (k = 0; k < 3; k++)
      {
        if (sane_check (www[i].v_grad[j][k]))
        {
          Error ("wind_check:sane_check www[%d].v_grad[%d][%d] %e\n", i, j, k, www[i].v_grad[j][k]);
          ierr++;
        }
      }

    }
  }

  if (ierr)
  {
    Error ("wind_check: Something is very seriously wrong with the wind.  %d problems Exiting\n", ierr);
    Exit (0);
  }


/* Now perform some checks to ensure DFUDGE is unlikely to punch through any cells  */
/* This versions does not change DFUDGE but simply logs places where problems might arise */


  delta = frac * DFUDGE;

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    ndim = zdom[ndom].ndim;
    mdim = zdom[ndom].mdim;
    if (zdom[ndom].coord_type == RTHETA)
    {
      drmin = 1e99;
      dtmin = 1e99;
      for (i = 0; i < ndim; i++)
      {
        for (j = 0; j < mdim; j++)
        {
          wind_ij_to_n (ndom, i, j, &n);
          if (wmain[n].vol > 0.0)
          {
            wind_ij_to_n (ndom, i + 1, j, &outer_n);
            wind_ij_to_n (ndom, i, j + 1, &outer_m);
            drmin = fabs (wmain[outer_n].r - wmain[n].r);
            dtmin = fabs (wmain[n].r * (wmain[outer_m].theta - wmain[n].theta) / RADIAN);
            if (drmin < delta || dtmin < delta)
            {
              Error ("wind_check: DFUDGE may be large in cell %d %d (%.1e %.1e)\n", i, j, drmin, dtmin);
            }
          }
        }
      }
    }
    else if (zdom[ndom].coord_type == CYLIND || zdom[ndom].coord_type == CYLVAR)
    {
      dxmin = 1e99;
      dzmin = 1e99;
      for (i = 0; i < ndim; i++)
      {
        for (j = 0; j < mdim; j++)
        {
          wind_ij_to_n (ndom, i, j, &n);
          if (wmain[n].vol > 0.0)
          {
            wind_ij_to_n (ndom, i + 1, j, &outer_n);
            wind_ij_to_n (ndom, i, j + 1, &outer_m);
            dxmin = fabs (wmain[outer_n].x[0] - wmain[n].x[0]);
            dzmin = fabs (wmain[outer_m].x[2] - wmain[n].x[2]);

            if (dxmin < delta || dzmin < delta)
            {
              Error ("wind_check: DFUDGE may be large in cell %d %d (%.1e %.1e)\n", i, j, dxmin, dzmin);
            }
          }
        }
      }
    }
    else if (zdom[ndom].coord_type == SPHERICAL)
    {
      drmin = 1e99;
      for (i = 0; i < ndim; i++)
      {
        if (wmain[i].vol > 0.0)
        {
          drmin = fabs (wmain[i + 1].r - wmain[i].r);
          if (drmin < delta)
          {
            Error ("wind_check: DFUDGE may be large in cell %d (%.1e)\n", i, drmin);
          }
        }
      }
    }
    else
    {
      Error ("wind_check: Disaster - unknown wind type\n");
      Exit (0);
    }
  }

  Log ("Wind_check: Punchthrough distance DFUDGE %e www[1].x[2] %e\n", DFUDGE, www[1].x[2]);
  return (0);
}
