
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
band_copy ()
{
  int n, nband;

  /* banding26: unless the old behaviour is requested, the cells get the ionization
     bands (geo.xfreq) rather than the photon generation bands */
  if (modes.ion_band_rule != ION_BANDS_GENERATION)
  {
    for (n = 0; n < NPLASMA; n++)
    {
      plasmamain[n].state.nbands = geo.nxfreq;
      for (nband = 0; nband < geo.nxfreq; nband++)
      {
        plasmamain[n].state.f1[nband] = geo.xfreq[nband];
        plasmamain[n].state.f2[nband] = geo.xfreq[nband + 1];
      }
      for (nband = geo.nxfreq; nband < NXBANDS + 1; nband++)
      {
        plasmamain[n].state.f1[nband] = 0.0;
        plasmamain[n].state.f2[nband] = 0.0;
      }
    }
    Log ("band_copy: cells use %d ionization bands (rule %d)\n", geo.nxfreq, modes.ion_band_rule);
    return;
  }

  for (n = 0; n < NPLASMA; n++)
  {
    plasmamain[n].state.nbands = xband.nbands;
    for (nband = 0; nband < xband.nbands; nband++)
    {
      plasmamain[n].state.f1[nband] = xband.f1[nband];
      plasmamain[n].state.f2[nband] = xband.f2[nband];
    }
    /* Initialize remaining unused band frequencies to safe values */
    for (nband = xband.nbands; nband < NXBANDS + 1; nband++)
    {
      plasmamain[n].state.f1[nband] = 0.0;
      plasmamain[n].state.f2[nband] = 0.0;
    }
  }
}
