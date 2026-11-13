
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
 * @brief copy the ionization bands to the spectrum
 * characterization bands in each plasma cell.
 *
 *
 * @details
 * 
 * The bands in which the spectrum in each cell is
 * characterized (the band sums and the power-law or
 * exponential models fitted to them) are the ionization
 * bands set up by ion_bands_init (geo.xfreq, geo.nxfreq),
 * not the photon generation bands (xband).  The two
 * sets are the same only with band model photon_bands
 * and at least MIN_N_IONBANDS generation bands.
 *
 * Each cell has its own copy of the band edges, so that
 * cells could in future have different bands, but at
 * present every cell gets the same ones.
 **********************************************************/

void
band_copy ()
{
  int n, nband;

  if (geo.nxfreq < 1 || geo.nxfreq > NXBANDS)
  {
    Error ("band_copy: number of ionization bands %d is not between 1 and NXBANDS (%d)\n", geo.nxfreq, NXBANDS);
    Exit (1);
  }

  for (n = 0; n < NPLASMA; n++)
  {
    plasmamain[n].state.nbands = geo.nxfreq;
    for (nband = 0; nband < geo.nxfreq; nband++)
    {
      plasmamain[n].state.f1[nband] = geo.xfreq[nband];
      plasmamain[n].state.f2[nband] = geo.xfreq[nband + 1];
    }
    /* Initialize remaining unused band frequencies to safe values */
    for (nband = geo.nxfreq; nband < NXBANDS + 1; nband++)
    {
      plasmamain[n].state.f1[nband] = 0.0;
      plasmamain[n].state.f2[nband] = 0.0;
    }
  }
}
