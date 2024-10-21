
/***********************************************************/
/** @file  wind_sum.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  Just prints a subset of quantities to the screen
 *
 * ### Notes ###
 *
 * At some level, a separate file for this does not seem
 * merited.  One should consider either adding other routines
 * to this file, or adding it to wind_updates2d.c 
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/** 
 * @brief      Prints out t_r and t_e and the total number of photons
 * in a subsample of the the grid that fits on a typical output terminal
 *
 * @param [in] WindPtr  w   The entire wind
 * @return     Always returns 0
 *
 * @details
 *
 * This is just a logging routine that prints out a few key parameters
 * at the end of each ionization cycle so that one can get a crude overview
 * of what is happening.  
 *
 * ### Notes ###
 * 
 * This routine is adapted from a very similar routine in swind 
 *
 * If the number of grid cells in the x direction is less more than 30 
 * the grid gets subsampled in the x direction
 * 
 * The routine won't fail but the output is less useful
 * when coodinates systems other than cylindrical are
 * used.
 *
 * 180526 - ksl - Fixed #384
 *
 *
 **********************************************************/

int
xtemp_rad (w)
     WindPtr w;
{
  int i, j, n;
  double x;
  int swind_min, swind_max, swind_delta;
  int nplasma;
  int ndom, ndim, mdim, nstart;
  int ntot;



  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {

    Log ("Results for Domain %d:\n", ndom);
    ndim = zdom[ndom].ndim;
    mdim = zdom[ndom].mdim;
    nstart = zdom[ndom].nstart;



    swind_min = 0;
    swind_max = ndim;

    /* swind_delta can be used to subsample the array */
    swind_delta = 1;
    if (mdim > 30)
      swind_delta = 1 + mdim / 30;


    if (zdom[ndom].coord_type == SPHERICAL)
    {
      Log ("\n R\n");
      j = 1;
      for (i = swind_min; i < swind_max; i += 1)
      {
        Log ("%8.2e ", w[nstart + i * mdim].r);
        if (j % 10 == 0)
        {
          Log ("\n");
        }
        j++;
      }
      Log ("\n T rad\n");

      j = 1;
      for (i = swind_min; i < swind_max; i += 1)
      {
        n = nstart + i;
        nplasma = w[n].nplasma;
        Log ("%8.2e ", plasmamain[nplasma].t_r);
        if (j % 10 == 0)
        {
          Log ("\n");
        }
        j++;
      }

      Log ("\n T e\n");

      j = 1;
      for (i = swind_min; i < swind_max; i += 1)
      {
        n = nstart + i;
        nplasma = w[n].nplasma;
        Log ("%8.2e ", plasmamain[nplasma].t_e);
        if (j % 10 == 0)
        {
          Log ("\n");
        }
        j++;
      }

      Log ("\n nphot\n");

      j = 1;
      for (i = swind_min; i < swind_max; i += 1)
      {
        n = nstart + i;
        nplasma = w[n].nplasma;
        if (wmain[n].inwind >= 0)
        {
          ntot = plasmamain[nplasma].ntot;
        }
        else
          ntot = 0;
        Log ("%8d ", ntot);
        if (j % 10 == 0)
        {
          Log ("\n");
        }
        j++;
      }

      Log ("\n");

      return (0);
    }
    /* 2d coord systems */
    if (zdom[ndom].coord_type != CYLIND)
    {
      Log ("Warning: Since coord type is not cylindrical, next print out may look odd\n");
    }


    Log ("\n T rad\n");

    Log ("   z\\x   ");
    for (i = swind_min; i < swind_max; i += swind_delta)
      Log ("%8.2e ", w[nstart + i * mdim].x[0]);
    Log ("\n");

    for (j = 0; j < mdim; j++)
    {
      Log ("%8.2e ", w[j].x[2]);
      for (i = swind_min; i < swind_max; i += swind_delta)
      {
        n = nstart + i * mdim + j;
        if (w[n].inwind >= 0)
        {
          nplasma = w[n].nplasma;
          x = plasmamain[nplasma].t_r;
        }
        else
          x = 0.0;
        Log ("%8.2g ", x);
      }
      Log ("\n");
    }

    Log ("\n T e\n");

    Log ("   z\\x   ");
    for (i = swind_min; i < swind_max; i += swind_delta)
      Log ("%8.2e ", w[nstart + i * mdim].x[0]);
    Log ("\n");

    for (j = 0; j < mdim; j++)
    {
      Log ("%8.2e ", w[j].x[2]);
      for (i = swind_min; i < swind_max; i += swind_delta)
      {
        n = nstart + i * mdim + j;
        nplasma = w[n].nplasma;
        if (wmain[n].inwind >= 0)
        {
          x = plasmamain[nplasma].t_e;
        }
        else
          x = 0.0;
        Log ("%8.2g ", x);
      }
      Log ("\n");
    }

    Log ("\n ntot \n");

    Log ("   z\\x   ");
    for (i = swind_min; i < swind_max; i += swind_delta)
      Log ("%8.2e ", w[nstart + i * mdim].x[0]);
    Log ("\n");

    for (j = 0; j < mdim; j++)
    {
      Log ("%8.2e ", w[j].x[2]);
      for (i = swind_min; i < swind_max; i += swind_delta)
      {
        n = nstart + i * mdim + j;
        nplasma = w[n].nplasma;
        if (wmain[n].inwind >= 0)
        {
          ntot = plasmamain[nplasma].ntot;
        }
        else
          ntot = -99;
        Log ("%8d ", ntot);
      }
      Log ("\n");
    }
  }


  return (0);

}
