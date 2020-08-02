/***********************************************************/
/** @file  estimators_simple.c
 * @author ksl
 * @date   July, 2020   
 *
 * @brief  Routines related to updating estimators for the
 * case of simple atoms.  The routine is primarily intended
 * for simple atoms, but also may be called when a macro-atom
 * calculation is carried out.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/** 
 * @brief      updates the estimators required for determining crude
 * spectra in each Plasma cell
 *
 * @param [in,out] PlasmaPtr  xplasma   PlasmaPtr for the cell of interest
 * @param [in] PhotPtr  p   Photon pointer
 * @param [in] double  ds   ds travelled
 * @param [in] double  w_ave   the weight of the photon in the cell. 
 *
 * @return  Always returns 0
 *
 *
 *
 * @details
 * 
 * Increments the estimators that allow one to construct a crude
 * spectrum in each cell of the wind.  The frequency intervals
 * in which the spectra are constructed are in geo.xfreq. This information
 * is used in different ways (or not at all) depending on the ionization mode.
 *
 * It also records the various parameters intended to describe the radiation field, 
 * including the IP.
 *
 * ### Notes ###
 *
 * The term direct refers to photons that have not been scattered by the wind.
 * 
 * In non macro atom mode, w_ave
 * this is an average weight (passed as w_ave), but 
 * in macro atom mode weights are never reduced (so p->w 
 * is used).
 *
 * This routine is called from bf_estimators in macro_atom modes
 * and from radiation (above).  Although the historical documentation
 * suggests it is only called for certain ionization modes, it appears
 * to be called in all cases, though clearly it is only provides diagnostic
 * information in some of them.
 *
 * 
 **********************************************************/


/* A couple of external variables to improve the counting of ionizing
   photons coming into a cell
*/
int nioniz_nplasma = -1;
int nioniz_np = -1;

/* A couple of external variables to improve the counting of photons
   in a cell
*/

int plog_nplasma = -1;
int plog_np = -1;

int
update_banded_estimators (xplasma, p, ds, w_ave, ndom)
     PlasmaPtr xplasma;
     PhotPtr p;
     double ds;
     double w_ave;
     int ndom;
{
  int i;
  double flux[3];
  double p_dir_cos[3];
  struct photon phot_mid;

  /*photon weight times distance in the shell is proportional to the mean intensity */

  xplasma->j += w_ave * ds;

  if (p->nscat == 0)
  {
    xplasma->j_direct += w_ave * ds;
  }
  else
  {
    xplasma->j_scatt += w_ave * ds;
  }



/* frequency weighted by the weights and distance in the shell .  See eqn 2 ML93 */
  xplasma->mean_ds += ds;
  xplasma->n_ds++;
  xplasma->ave_freq += p->freq * w_ave * ds;


/* The lines below compute the flux element of this photon */

  stuff_phot (p, &phot_mid);    // copy photon ptr
  move_phot (&phot_mid, ds / 2.);       // get the location of the photon mid-path 
  stuff_v (p->lmn, p_dir_cos);  //Get the direction of the photon packet

  renorm (p_dir_cos, w_ave * ds);       //Renormalise the direction into a flux element
  project_from_xyz_cyl (phot_mid.x, p_dir_cos, flux);   //Go from a direction cosine into a cartesian vector

  if (p->x[2] < 0)              //If the photon is in the lower hemisphere - we need to reverse the sense of the z flux
    flux[2] *= (-1);


  /*Deal with the special case of a spherical geometry */

  if (zdom[ndom].coord_type == SPHERICAL)
  {
    renorm (phot_mid.x, 1);     //Create a unit vector in the direction of the photon from the origin
    flux[0] = dot (p_dir_cos, phot_mid.x);      //In the spherical geometry, the first comonent is the radial flux
    flux[1] = flux[2] = 0.0;    //In the spherical geomerry, the theta and phi compnents are zero    
  }

/* We now update the fluxes in the three bands */


  if (p->freq < UV_low)
  {
    vadd (xplasma->F_vis, flux, xplasma->F_vis);
    xplasma->F_vis[3] += length (flux);
  }
  else if (p->freq > UV_hi)
  {
    vadd (xplasma->F_Xray, flux, xplasma->F_Xray);
    xplasma->F_Xray[3] += length (flux);
  }
  else
  {
    vadd (xplasma->F_UV, flux, xplasma->F_UV);
    xplasma->F_UV[3] += length (flux);
  }


  /* The next loop updates the banded versions of j and ave_freq, analogously to routine inradiation
     nxfreq refers to how many frequencies we have defining the bands. So, if we have 5 bands, we have 6 frequencies, 
     going from xfreq[0] to xfreq[5] Since we are breaking out of the loop when i>=nxfreq, this means the last loop 
     will be i=nxfreq-1 */

  /* note that here we can use the photon weight and don't need to calculate anm attenuated average weight
     as energy packets are indisivible in macro atom mode */


  for (i = 0; i < geo.nxfreq; i++)
  {
    if (geo.xfreq[i] < p->freq && p->freq <= geo.xfreq[i + 1])
    {
      xplasma->xave_freq[i] += p->freq * w_ave * ds;    /* frequency weighted by weight and distance */
      xplasma->xsd_freq[i] += p->freq * p->freq * w_ave * ds;   /* input to allow standard deviation to be calculated */
      xplasma->xj[i] += w_ave * ds;     /* photon weight times distance travelled */
      xplasma->nxtot[i]++;      /* increment the frequency banded photon counter */
      /* work out the range of frequencies within a band where photons have been seen */
      if (p->freq < xplasma->fmin[i])
      {
        xplasma->fmin[i] = p->freq;
      }
      if (p->freq > xplasma->fmax[i])
      {
        xplasma->fmax[i] = p->freq;
      }

    }
  }

  /* NSH 131213 slight change to the line computing IP, we now split out direct and scattered - this was 
     mainly for the progha_13 work, but is of general interest */
  /* 70h -- nsh -- 111004 added to try to calculate the IP for the cell. Note that 
   * this may well end up not being correct, since the same photon could be counted 
   * several times if it is rattling around.... */

  /* 1401 JM -- Similarly to the above routines, this is another bit of code added to radiation
     which originally did not get called in macro atom mode. */

  /* NSH had implemented a scattered and direct contribution to the IP. This doesn't really work 
     in the same way in macro atoms, so should instead be thought of as 
     'direct from source' and 'reprocessed' radiation */

  if (xplasma->nplasma != plog_nplasma || p->np != plog_np)
  {
    xplasma->ntot++;

    /* NSH 15/4/11 Lines added to try to keep track of where the photons are coming from, 
     * and hence get an idea of how 'agny' or 'disky' the cell is. */
    /* ksl - 1112 - Fixed this so it records the number of photon bundles and not the total
     * number of photons.  Also used the PTYPE designations as one should as a matter of 
     * course
     */

    if (p->origin == PTYPE_STAR)
      xplasma->ntot_star++;
    else if (p->origin == PTYPE_BL)
      xplasma->ntot_bl++;
    else if (p->origin == PTYPE_DISK)
      xplasma->ntot_disk++;
    else if (p->origin == PTYPE_WIND)
      xplasma->ntot_wind++;
    else if (p->origin == PTYPE_AGN)
      xplasma->ntot_agn++;
    plog_nplasma = xplasma->nplasma;
    plog_np = p->np;
  }







  if (HEV * p->freq > 13.6)     // only record if above H ionization edge
  {

    /*
     * Calculate the number of H ionizing photons, see #255
     * EP 11-19: moving the number of ionizing photons counter into this
     * function so it will be incremented for both macro and non-macro modes
     */
    if (xplasma->nplasma != nioniz_nplasma || p->np != nioniz_np)
    {
      xplasma->nioniz++;
      nioniz_nplasma = xplasma->nplasma;
      nioniz_np = p->np;
    }

    /* IP needs to be radiation density in the cell. We sum contributions from
       each photon, then it is normalised in wind_update. */
    xplasma->ip += ((w_ave * ds) / (PLANCK * p->freq));

    if (HEV * p->freq < 13600)  //Tartar et al integrate up to 1000Ryd to define the ionization parameter
    {
      xplasma->xi += (w_ave * ds);
    }

    if (p->nscat == 0)
    {
      xplasma->ip_direct += ((w_ave * ds) / (PLANCK * p->freq));
    }
    else
    {
      xplasma->ip_scatt += ((w_ave * ds) / (PLANCK * p->freq));
    }
  }




  return (0);
}

