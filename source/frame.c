/***********************************************************/
/** @file  frame.c
 * @author ksl
 * @date   May, 2020
 *
 * @brief
 * These routines have to do with switching from the observer 
 * to the local frame and back
 *
 * ###Notes###
 *
 * These routines are written so that if p and p_obs point 
 * to the same variable then everything occurs in place
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief  Check that a photon is in the desired frame.
 *
 * @param [in] PhotPtr  p      The photon whose frame is to be checked
 * @param [in] frame    frame  The desired frame of the photon
 * @param [in] char    *msg    A comment to print as an error message if the
 *                             photon is in the incorrect frame
 * @return 0 (False) if the photons is in the desired frame, 1 (True) if the
 *  photon is NOT in the desired frame.
 *
 * @details
 *
 * The purpose of this function is to avoid situations where photons are being
 * transformed incorrectly between the local and observer frame.
 *
 *
 **********************************************************/

int
check_frame (p, desired_frame, msg)
     PhotPtr p;
     enum frame desired_frame;
     char *msg;
{
  if (p->frame == desired_frame)
  {
    return EXIT_SUCCESS;
  }
  else
  {
    Error ("check_frame: %s :Photon (%5d) of istat (%2d) and origin (%2d) not in desired frame %d (0=Loc,1=Obs)\n", msg, p->np, p->istat,
           p->origin, desired_frame);
    return EXIT_FAILURE;
  }
}

/**********************************************************/
/**
 * @brief  Calculate the gamma factor for a given velocity
 *
 * @param [in]  double vel[3]  the velocity vector
 *
 * @return double the value of gamma
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

double
calculate_gamma_factor (double vel[3])
{
  return 1.0 / sqrt (1.0 - dot (vel, vel) / (VLIGHT * VLIGHT));
}

/**********************************************************/
/**
 * @brief      carries out the transformation of a all the quantities
 *      in a photon structure
 *      from the observer (or global) frame to the local (or co-moving)
 *      frame
 *
 * @param [in] PhotPtr  p_in   The photon in the observer frame
 * @param [out] PhotPtr  p_out   The photon in the local frame
 *
 * @return    The routine returns 0, unless there was an error
 *
 *
 * @details
 * The routine copies all of the quantities contain in the p_in to
 * p_out, obtains the velocity at that point in the wind, and then
 * performs an in place transformation from the * global to the local
 * frame of the photon allowing for special
 * relativity.
 *
 *
 * ### Notes ###
 *
 * It p_in and p_out are the same then the transformation will
 * be performed in place.
 *
 **********************************************************/

int
observer_to_local_frame (p_in, p_out)
     PhotPtr p_in, p_out;
{
  WindPtr one;
  int ndom;
  double v[3];
  int ierr;


  ierr = check_frame (p_in, F_OBSERVER, "Observer_to_local_frame");

  /* Initialize the output photon */
  stuff_phot (p_in, p_out);
  p_out->frame = F_LOCAL;

  /* Calculate the local velocity of the wind at this position */
  one = &wmain[p_in->grid];
  ndom = one->ndom;
  vwind_xyz (ndom, p_in, v);

  ierr = lorentz_transform (p_in, p_out, v);



  return (ierr);
}



/**********************************************************/
/**
 * @brief      carries out the transformation of all the quantities
 *      in a photon structure 
 *      from the local (or co-moving) frame to the observer (or global)
 *      frame
 *
 * @param [in] PhotPtr  p_in   The photon in the local frame                   
 * @param [out] PhotPtr  p_out   The photon in the global frame                  
 *
 * @return    The routine returns 0 unless there was an error       
 *
 *
 * @details
 * The routine copies all of the quantities contain in the p_in to 
 * p_out, obtains the velocity at that point in the wind, and then 
 * performs an in place transformation from the local to the observer
 * frame of the photon allowing for special 
 * relativity.
 *
 *
 * ### Notes ###
 *
 * It p_in and p_out are the same then the transformation will
 * be performed in place.
 *
 **********************************************************/


int
local_to_observer_frame (p_in, p_out)
     PhotPtr p_in, p_out;
{
  WindPtr one;
  int ndom;
  double v[3];
  int ierr;



  ierr = check_frame (p_in, F_LOCAL, "local_to_observer_frame");

  /* Initialize the output photon */
  stuff_phot (p_in, p_out);
  p_out->frame = F_OBSERVER;

  /* Calculate the local velocity of the wind at this position */
  one = &wmain[p_in->grid];
  ndom = one->ndom;
  vwind_xyz (ndom, p_in, v);
  rescale (v, -1, v);

  ierr = lorentz_transform (p_in, p_out, v);


  return (ierr);
}



/**********************************************************/
/**
 * @brief      carries out the transformation of a all the quantities
 *      in a photon structure
 *      from the observer (or global) frame to the local (or co-moving)
 *      frame
 *
 * @param [in] PhotPtr  p_in   The photon in the observer frame
 * @param [out] PhotPtr  p_out   The photon in the local frame
 *
 * @return    The routine returns 0, unless there was an error
 *
 *
 * @details
 * The routine copies all of the quantities contain in the p_in to
 * p_out, obtains the velocity at that point in the wind, and then
 * performs an in place transformation from the * global to the local
 * frame of the photon allowing for special
 * relativity.
 *
 *
 * ### Notes ###
 *
 * It p_in and p_out are the same then the transformation will
 * be performed in place.
 *
 **********************************************************/




int
observer_to_local_frame_disk (p_in, p_out)
     PhotPtr p_in, p_out;
{
//  WindPtr one;
//int ndom;
  double v[3];
  int ierr;


  ierr = check_frame (p_in, F_OBSERVER, "Observer_to_local_frame_disk");

  /* Initialize the output photon */
  stuff_phot (p_in, p_out);
  p_out->frame = F_LOCAL;


  /* Calculate the local velocity of the disk at this position */
  vdisk (p_in->x, v);

  ierr = lorentz_transform (p_in, p_out, v);



  return (ierr);
}



/**********************************************************/
/**
 * @brief      carries out the transformation of a all the quantities
 *      in a photon emitted from the surface of the disk
 *      from the local (or co-moving) frame to the observer (or global)
 *      frame
 *
 * @param [in] PhotPtr  p_in   The photon in the local frame                   
 * @param [out] PhotPtr  p_out   The photon in the global frame                  
 *
 * @return    The routine returns 0 unless the was an error         
 *
 *
 * @details
 * The routine copies all of the quantities contain in the p_in to 
 * p_out, obtains the velocity at that point in the wind, and then 
 * performs an in place transformation from the local to the observer
 * frame of the photon allowing for special 
 * relativity.
 *
 *
 * ### Notes ###
 *
 * It p_in and p_out are the same then the transformation will
 * be performed in place.  This is normally the way the 
 * routine is used.   One could have omitted p_out in this case
 * but to make the interface coisistent, both photon structures
 * are included.  
 *
 **********************************************************/



int
local_to_observer_frame_disk (p_in, p_out)
     PhotPtr p_in, p_out;
{
  double v[3];
  int ierr;



  ierr = check_frame (p_in, F_LOCAL, "local_to_observer_frame");

  /* Initialize the output photon */
  stuff_phot (p_in, p_out);
  p_out->frame = F_OBSERVER;

  /* Calculate the local velocity of the wind at this position */
  vdisk (p_in->x, v);
  rescale (v, -1, v);

  ierr = lorentz_transform (p_in, p_out, v);


  return (ierr);
}


/**********************************************************/
/**
 * @brief      calculate the distance a photon will travel
 *      in the local frame given a distance in the observer
 *      frame
 *
 * @param [in] PhotPtr  p_obs    The photon in the observer frame                
 * @param [in] double   ds_obs   The distance from a starting point
                                 for the photon to travel
 *
 * @return    The distance in the co-moving or local frame         
 *
 *
 * @details
 *
 *
 * ### Notes ###
 *
 *
 **********************************************************/

double
observer_to_local_frame_ds (p_obs, ds_obs)
     PhotPtr p_obs;
     double ds_obs;
{
  WindPtr one;
  int ndom;
  double v[3];
  double gamma;
  double ds_cmf;

  if (rel_mode == REL_MODE_LINEAR)
  {
    return (ds_obs);
  }

  /* Calculate the local velocity of the wind at this position */
  one = &wmain[p_obs->grid];
  ndom = one->ndom;
  vwind_xyz (ndom, p_obs, v);

  gamma = 1. / sqrt (1 - (dot (v, v) / (VLIGHT * VLIGHT)));

  ds_cmf = ds_obs;


  ds_cmf *= gamma * (1 - dot (p_obs->lmn, v) / VLIGHT);


  return (ds_cmf);



}

/**********************************************************/
/**
 * @brief      calculate the distance a photon will travel
 *      in the observer frame given a distance in the local
 *      frame
 *
 * @param [in] PhotPtr  p_obs    The photon in the observer frame                
 * @param [in] double   ds_cmf   The distance from a starting point
                                 for the photon to travel
 *
 * @return    The distance in the observer frame         
 *
 *
 * @details
 *
 * Note that the photon MUST BE in the observer frame here. 
 *
 *
 * ### Notes ###
 *
 *
 **********************************************************/

double
local_to_observer_frame_ds (p_obs, ds_cmf)
     PhotPtr p_obs;
     double ds_cmf;
{
  WindPtr one;
  int ndom;
  double v[3];
  double gamma;
  double ds_obs;

  if (rel_mode == REL_MODE_LINEAR)
  {
    return (ds_cmf);
  }

  /* Calculate the local velocity of the wind at this position */
  one = &wmain[p_obs->grid];
  ndom = one->ndom;
  vwind_xyz (ndom, p_obs, v);

  gamma = 1. / sqrt (1 - (dot (v, v) / (VLIGHT * VLIGHT)));

  ds_obs = ds_cmf;
  ds_obs /= gamma * (1 - dot (p_obs->lmn, v) / VLIGHT);


  return (ds_obs);



}

/**********************************************************/
/**
 * @brief      calculate a velocity in the local frame given
 *      a velocity in the observer frame
 *      
 *
 * @param [in] double   *v_obs       A velocity in the observer frame                
 * @param [in] double   *v           The velocity of the local frame as 
 *                                   measured in the observers frame
 * @param [out] double  *v_cmf       The peculiar velocity in the local frame
 *
 * @return    The speed in the local frame               
 *
 *
 * @details
 * This uses the standard special relativitistic velocity addition
 * law to calculate a peculiar velocity in the local frame.
 *
 *
 * ### Notes ###
 *
 *
 **********************************************************/

double
observer_to_local_frame_velocity (v_obs, v, v_cmf)
     double *v_obs;
     double *v;
     double *v_cmf;
{
  double gamma, c1, c2;
  double a[3], b[3];
  double vdotv;

  if (rel_mode == REL_MODE_LINEAR)
  {
    vsub (v_obs, v, v_cmf);
    return (length (v_cmf));
  }

  gamma = 1. / sqrt (1 - (dot (v, v) / (VLIGHT * VLIGHT)));

  vdotv = dot (v_obs, v) / (VLIGHT * VLIGHT);

  c1 = 1. / (1 - vdotv);
  c2 = (gamma / (1 + gamma) * vdotv);

  rescale (v_obs, 1. / gamma, a);
  vsub (a, v, a);
  rescale (v, c2, b);
  vadd (a, b, a);
  rescale (a, c1, v_cmf);

/*  The transformation formula is
   v_cmf=c1*(v_obs/xgamma-v+c2*v);
*/
  return length (v_cmf);


}

/**********************************************************/
/**
 * @brief      calculate a velocity in the observer frame
 *      a velocity in the local frame
 *      
 *
 * @param [in] double   *v_cmf       A peculiar velocity in the local frame                
 * @param [in] double   *v           The velocity of the local frame in
 *                                   the observers frame
 * @param [out] double  *v_obs       The velocity in the global frame
 *
 * @return    The speed in the observer frame               
 *
 *
 * @details
 * This uses the standard special relativitistic velocity addition
 * law to calculate a velocity in the observer frame from a "peculiar"
 * velocity in the local frame.
 *
 *
 *
 * ### Notes ###
 *
 *
 **********************************************************/

double
local_to_observer_frame_velocity (v_cmf, v, v_obs)
     double *v_cmf;
     double *v;
     double *v_obs;
{
  double gamma, c1, c2;
  double a[3], b[3];
  double vdotv;

  if (rel_mode == REL_MODE_LINEAR)
  {
    vadd (v_cmf, v, v_obs);
    return (length (v_obs));
  }

  gamma = 1. / sqrt (1 - (dot (v, v) / (VLIGHT * VLIGHT)));

  vdotv = dot (v_obs, v) / (VLIGHT * VLIGHT);

  c1 = 1. / (1 + vdotv);
  c2 = (gamma / (1 + gamma) * vdotv);

  rescale (v_cmf, 1. / gamma, a);
  vadd (a, v, a);
  rescale (v, c2, b);
  vadd (a, b, a);
  rescale (a, c1, v_obs);

/*  The transformation formula is
   v_obs=c1*(v+v_cmf/xgamma+c2*v)
*/
  return length (v_obs);


}



/**********************************************************/
/**
 * @brief      calculate how a 3-vector would look in the observer
 *      frame given a 3-vector in the local frame.
 *      
 *
 * @param [in]  double  V           the v  of the co-moving frame 
 * @param [in]  double  dx_cmf      the 3 vector in the local frame
 * @param [out] double  dx_obs     the resulting 3 vector in the obsevr frame
 *
 * @return    Always returns 0               
 *
 *
 * @details
 *
 * ### Notes ###
 *
 * Except for the fact that the distance along the velocity vector
 * gets contracted in the observer frame, this is identical to
 * the routine observer_to_local_frame_ruler_transform
 *
 **********************************************************/

int
local_to_observer_frame_ruler_transform (v, dx_cmf, dx_obs)
     double v[], dx_cmf[], dx_obs[];
{

  double beta, gamma, speed;
  double lmn[3];
  double ds_cmf_par, dx_cmf_par[3];
  double dx_cmf_perp[3];
  double dx_obs_par[3];

  if (rel_mode == REL_MODE_LINEAR)
  {
    stuff_v (dx_cmf, dx_obs);
    return (0);
  }



  speed = length (v);
  beta = speed / VLIGHT;
  gamma = 1. / (1. - beta * beta);


  rescale (v, 1. / speed, lmn);
  ds_cmf_par = dot (lmn, dx_cmf);
  rescale (lmn, ds_cmf_par, dx_cmf_par);

  vsub (dx_cmf, dx_cmf_par, dx_cmf_perp);

  rescale (dx_cmf_par, 1. / gamma, dx_obs_par);

  vadd (dx_obs_par, dx_cmf_perp, dx_obs);

  return (0);


}

/**********************************************************/
/**
 * @brief      calculate how a 3-vector would look in the observer
 *      frame given a 3-vector in the cmf frame.
 *      
 *
 * @param [in]  double  V           the v  of the co-moving frame 
 * @param [in] double   dx_obs      the 3 vector in the observer frame
 * @param [out] double  dx_cmf     the resulting 3 vector in the local frame
 *
 * @return    Always returns 0               
 *
 *
 * @details
 *
 * ### Notes ###
 *
 * Except for the fact that the distance along the velocity vector
 * gets larger in the local frame, this is identical to
 * the routine local_to_obseerver_frame_ruler_transform
 *
 **********************************************************/


int
observer_to_local_frame_ruler_transform (v, dx_obs, dx_cmf)
     double v[], dx_obs[], dx_cmf[];
{

  double beta, gamma, speed;
  double lmn[3];
  double ds_obs_par, dx_obs_par[3];
  double dx_obs_perp[3];
  double dx_cmf_par[3];

  if (rel_mode == REL_MODE_LINEAR)
  {
    stuff_v (dx_obs, dx_cmf);
    return (0);
  }



  speed = length (v);
  beta = speed / VLIGHT;
  gamma = 1. / (1. - beta * beta);


  rescale (v, 1. / speed, lmn);
  ds_obs_par = dot (lmn, dx_obs);
  rescale (lmn, ds_obs_par, dx_obs_par);

  vsub (dx_obs, dx_obs_par, dx_obs_perp);

  rescale (dx_obs_par, gamma, dx_cmf_par);

  vadd (dx_cmf_par, dx_obs_perp, dx_cmf);

  return (0);


}


/**********************************************************/
/**
 * @brief      carries out the transformation of all the quantities
 *      in a photon structure from one frame to another
 *
 * @param [in] PhotPtr  p_in   The photon in the observer frame
 * @param [out] PhotPtr  p_out   The photon in the local frame
 *
 * @return    The routine returns 0, unless there was an error
 *
 *
 * @details
 *
 *  This does a transfrom from the current frame into a frame moving with velocity v
 *  with respect to the current frame.
 *
 * ### Notes ###
 *
 * One needs to be careful about the direction of the transformation
 * As stated this tranforms into a frame that has a velocity v with
 * respect to the current frame.
 *
 * It p_in and p_out are the same then the transformation will
 * be performed in place.
 *
 **********************************************************/


int
lorentz_transform (p_in, p_out, v)
     PhotPtr p_in, p_out;
     double v[];
{
  double f_out, f_in;
  double x;
  double vel;
  double gamma;
  int i, ierr;

  ierr = FALSE;


  vel = dot (p_in->lmn, v);


  f_in = p_in->freq;

  if (rel_mode == REL_MODE_LINEAR)
  {
    if (p_in->frame == F_LOCAL)
    {
      f_out = p_out->freq = f_in * (1. - vel / VLIGHT);
    }
    else
    {
      f_out = p_out->freq = f_in / (1. + vel / VLIGHT);
    }
    return (ierr);
  }



  gamma = 1. / (sqrt (1 - (dot (v, v) / (VLIGHT * VLIGHT))));

  f_out = p_out->freq = f_in * gamma * (1. - vel / VLIGHT);

  x = gamma / VLIGHT * (1.0 - (gamma * vel / ((gamma + 1) * VLIGHT)));

  for (i = 0; i < 3; i++)
  {
    p_out->lmn[i] = f_in / f_out * (p_in->lmn[i] - x * v[i]);
  }

  p_out->w *= (f_out / f_in);


  return (ierr);
}
