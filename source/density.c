
/***********************************************************/
/** @file  density.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  Contains a rroutine to find the density of an ion at a position
 * designamted by a Pothon pointner
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
 * @brief      finds the density of a specific ion designated by nion at a position 
 *
 * @param [in] int  ndom   The domain number where the density is to be obtained
 * @param [in] double  x[] A 3 vector contaiing a postion
 * @param [in] int  nion   The ion number of the ion for which a density is needed
 * @return     The density of the ion in cgs units
 *
 * @details
 * The routine interpolates a density from the densitiees ad the centers of cells.
 *
 * ### Notes ###
 *
 * The routine does not determine what domain a parituclar position refers to, so
 * this must be provided.
 *
 * Densities are defined at the centers of cells. The routine depends upon the
 * fact that densities are defined in the edge cells to work out what to do
 * at the edges of the wind.  The routine uses corrd_frac to determine what
 * cells have to be interpolated.
 *
 **********************************************************/

double
get_ion_density (ndom, x, nion)
     int ndom;
     double x[];
     int nion;
{
  double dd;
  int nn, nnn[4], nelem;
  double frac[4];
  int nplasma;


/* nnn is an array which contains the elements of the Wind (coordinate)
structure that must be summed. We need to convert that to a position
in the plasma structure*/

  if ((coord_fraction (ndom, 1, x, nnn, frac, &nelem)) > 0)
  {

    dd = 0;

    for (nn = 0; nn < nelem; nn++)
    {
      nplasma = wmain[nnn[nn]].nplasma;
      dd += plasmamain[nplasma].density[nion] * frac[nn];
    }
  }
  else
  {
    dd = 0;
  }


  return (dd);
}
