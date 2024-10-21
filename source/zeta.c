
/***********************************************************/
/** @file  zeta.c
 * @author ksl,nsh
 * @date   April, 2018
 *
 * @brief  The routines in this file are to do with computing the 
 * correction to the Saha equation for a modified on the spot
 * approach to ionization.  
 *
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "sirocco.h"
//OLD #include "recipes.h"


/**********************************************************/
/** 
 * @brief      calculates the correction term used in in the L+M and 
 * Sim equations as a modfication to the Saha equation
 *
 * @param [in] double  temp   The temperature for the calculation
 * @param [in] int  nion  the ion we are recombining to
 * @param [in] int  mode  The mode  allow that controls how zeta is calcuated
 * @return     zeta, the correction factor to the recombination rates.
 *
 * @details
 *
 * zeta is part of correction factor suggested in LM93 in a modified version
 * of the on-the spot approximation.  Specifically it is the
 * fraction of recombinations that go directly to the ground state
 * Saha equation.  
 *
 * * mode 1 is the original version, it takes ilow, ihi and interpfrac and 
 *  interpolates in the ground state fraction table.
 * * mode 2 attempts to modify zeta by application of the dielectronic recombination
 *  rate - this calculates zeta from first principles using ground state recombination
 * 	rates and total recombination rates. It defaults to mode 1 if is it missing data
 * 	but it returns an error.
 *
 * 	The program will exit if the mode is unknown.
 *
 * ### Notes ###
 *
 * Values for zeta were tabulated by CK (circa 2000)
 *
 * When compute_zeta is called, nion is the lower ion in the pair whose abundances are being calculated.
 * Calls to the various recombination rate coefficient subroutines are made with nion+1, since these
 * rates are tabulated for the recombining ion.
 *
 * The data are read in by get_atomicdata.c (around line 2100)
 **********************************************************/

double
compute_zeta (temp, nion, mode)
     double temp;
     int mode, nion;
{
  double zeta, interpfrac, dummy;
  int ihi, ilow;

#define MIN_FUDGE  1.e-10
#define MAX_FUDGE  10.
#define TMIN_ZETA    2000.      /* This minimum temperature  TMIN_ZETA is specific to CK's ground fraction tables */

  zeta = 0.0;


  /* now get the right place in the ground_frac tables  CK */
  dummy = temp / TMIN_ZETA - 1.;
  ilow = dummy;                 /* have now truncated to integer below */
  ihi = ilow + 1;               /*these are the indeces bracketing the true value */
  interpfrac = (dummy - ilow);  /*this is the interpolation fraction */
  if (ilow < 0)
  {
    ilow = 0;
    ihi = 0;
    interpfrac = 1.;
  }
  if (ihi > 19)
  {
    ilow = 19;
    ihi = 19;
    interpfrac = 1.;
  }


  if (mode == 1)
  {
    /* This is the original method of getting zeta - just from the ground state tables computed by Christian */

    zeta = ground_frac[nion].frac[ilow] + interpfrac * (ground_frac[nion].frac[ihi] - ground_frac[nion].frac[ilow]);
  }


  /* Best try at full blown zeta including DR, if we have the data, else default to the old way of doing things */
  else if (mode == 2)
  {
    /* We call this with the lower ion in a pair, so if that ion has only one electron, 
       (ie. Carbon 6) then we cannot have DR into this ion, so there will be no DR rate 
       associated with it. NSH 140317 We also do this if we dont have DR data at all. THE DR 
       data is associated with the ion doing the recombining */
    if (ion[nion].istate == ion[nion].z || ion[nion + 1].drflag == 0)   //Either we dont have DR data, or DR cannot happen       
    {
      if (ion[nion + 1].total_rrflag == 1)      //We have total RR data
      {
        zeta = gs_rrate (nion + 1, temp) / total_rrate (nion + 1, temp);
      }

      else                      //We do not have total RR data, so we must default back to the old way of doing things
      {
        Error ("Compute zeta: total RR rate missing for element %i state %i\n", ion[nion + 1].z, ion[nion + 1].istate);
        zeta = ground_frac[nion].frac[ilow] + interpfrac * (ground_frac[nion].frac[ihi] - ground_frac[nion].frac[ilow]);
      }
    }
    else                        //We do have DR data, so we can include this in zeta
    {
      if (ion[nion + 1].drflag > 0 && ion[nion + 1].total_rrflag == 1)  //We have total RR data
      {
        compute_dr_coeffs (temp);
        zeta = gs_rrate (nion + 1, temp) / (total_rrate (nion + 1, temp) + dr_coeffs[nion + 1]);
      }
      else                      //We dont have total RR data
      {
        Error ("Compute zeta: total RR rate missing for element %i state %i\n", ion[nion + 1].z, ion[nion + 1].istate);
        zeta = ground_frac[nion].frac[ilow] + interpfrac * (ground_frac[nion].frac[ihi] - ground_frac[nion].frac[ilow]);
      }
    }
  }
  else
  {
    Error ("Compute zeta: Unkown mode %i \n", mode);
    Exit (0);
  }


  return (zeta);
}
