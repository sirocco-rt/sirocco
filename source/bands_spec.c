
/***********************************************************/
/** @file  bands_spec.c
 * @author ksl
 * @date   April 2025   
 *
 * @brief  Setup the frequency bands used for and
 * for characterisin  spectra in wind cells.
 *
 * The subroutines here deal with * how to record the spectrum 
 * that pass through individual cells.  This is distinct from
 * the banding used in creating the photons 
 *
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

/* Actual structures are in sirocco.h.  Here for reference only.

#define NBANDS 10
struct xbands
{
	double f1[NBANDS],f2[NBANDS];
	double min_fraction[NBANDS];
	double nat_fraction[NBANDS];          // The fraction of the accepted luminosity in this band
	double used_fraction[NBANDS];
	double flux[NBANDS];                  //The "luminosity" within a band
	double weight[NBANDS];
	int nphot[NBANDS];
	int nbands;           // Actual number of bands in use
}
xband;
*/


/**********************************************************/
/** 
 * @brief copy frequencies from the from the photon generation
 * banding structure to the spectrum charcaterization bands
 * in each plasma cell.
 *
 *
 * @details
 * 
 * This is just a simple routine to initialize the spectrl
 * charateriztion band frequences to the same bands that
 * are being used for photon generation.  
 **********************************************************/

void
band_copy()
{
    int n,nband;
    for (n=0;n<NPLASMA;n++){
        plasmamain[n].nbands=xband.nbands;
        for(nband=0;nband<xband.nbands;nband++){
            plasmamain[n].f1[nband]=xband.f1[nband];
            plasmamain[n].f2[nband]=xband.f2[nband];

        }
    }
}
