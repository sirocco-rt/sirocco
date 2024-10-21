
/***********************************************************/
/** @file  pi_rates.c
 * @author nsh
 * @date   Aug, 2014
 *
 * @brief  Routines associated with calculating photionization rates
 *
 * ???
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

struct topbase_phot *xtop;      //Topbase description of a photoionization x-section - this is the only type we use - tabulated.

double qromb_temp;              //Temperature used in integrations - has to be an external variable so qromb can use it

//Model parameters for integrations, also passed externally so qromb can use them

double xpl_alpha, xpl_w, xpl_logw;
double xexp_temp, xexp_w;


/**********************************************************/
/**
 * @brief      calculate the photoionization rate coefficient for an ion in a given cell
 *
 * @param [in] int  nion   The ion being ionized (or the index to an inner shell cross section)
 * @param [in] PlasmaPtr  xplasma   The cell for which rates are calulated
 * @param [in] int  mode Use either a piecewise model (1) or a blackbody model (2) for J_nu
 * @param [in] type Compute either an outer shell (1) or inner shell (2) rate
 * @return     The photioinization rate coefficient for the ion.
 *
 * @details
 *  calc_pi_rate  calculates the photoionization rate coefficient for ion
 * 	nion, based upon the mean intensity stored in cell xplasma
 * 	The mode tells the subroutine wether we are modelling the
 * 	mean intensity as a dilute blackbody (mode2) or as a series
 * 	of power laws and exponentials (mode1). The type allows for the calculation
 *      of inner shell rates - in this case nion is the index into an inner
 * 	shell cross section record - this links to the relevant ion internally.
 *      the reason for the difference is that there is exactly one outer rate
 * 	per ion, but there can be many inner shell rates. Most of the information
 * 	needed for the calculations is stored in the xplasma structure (e.g. temperature
 * 	and spectral model)
 *
 * ### Notes ###
 * This was created in Summer 2014 in preparation for matrix ionization solver. Previously, this code
 * was contained in two subroutines bb_correct_2 and pl_correct_2.The functionality of these two
 * have been combined into one - hence the requirement for the mode parameter. It was further extended
 * to deal with inner shell rates - hence the type parameter
 *
 **********************************************************/

double
calc_pi_rate (nion, xplasma, mode, type)
     int nion;
     PlasmaPtr xplasma;
     int mode;
     int type;
{
  int j;
  double pi_rate;
  int ntmin, nvmin;
  double fthresh, fmax, fmaxtemp;
  double f1, f2;
  double exp_qromb, pl_qromb;


  exp_qromb = 1e-4;             /*These are the two tolerance values for the rhomberg integrals */
  pl_qromb = 1e-4;

  ntmin = nvmin = -1;           /* Initialize these to an unreasonable number. We dont use them all the time */

  if (type == 1)                //We are computing a normal outer shell rate
  {
    if (-1 < nion && nion < nions)      //Get cross section for this specific ion_number
    {
      ntmin = ion[nion].ntop_ground;    /*We only ever use the ground state cross sections. This is for topbase */
      nvmin = ion[nion].nxphot; /* this gets any verner cross section */
    }
    else
    {
      Error ("calc_pi_rate: %d is unacceptable value of nion\n", nion);
      Exit (0);
      return (1.0);
    }
    if (ion[nion].phot_info > 0)        //topbase or hybrid VFKY (GS)+TB excited
    {
      xtop = &phot_top[ntmin];
    }
    else if (ion[nion].phot_info == 0)  // verner
    {
      xtop = &phot_top[nvmin];  //If we don't have a topbase or hybrid cross section we fall back to VFKY
    }
    else
    {
      Error ("calc_pi_rate: No photoionization xsection for ion %d (element %d, ion state %d)\n", nion, ion[nion].z, ion[nion].istate);
      Exit (0);
    }
  }
  else if (type == 2)           //We are computing an inner shell rate - nion refers to an actual cross section.
  {
    if (-1 < nion && nion < n_inner_tot)        //We have a reasonable value for nion_in
    {
      xtop = &inner_cross[nion];        //Set the cross sections for the integral to the inner shell cross section
    }
    else
    {
      Error
        ("calc_pi_rate: No inner shell xsection for record %d (element %d, ion state %d)\n",
         nion, inner_cross[nion].z, inner_cross[nion].istate);
      Exit (0);
    }
  }
  else
  {
    Error ("calc_pi_rate: unknown mode %i\n", type);
  }



/* At this stage, xtop points to either an outer or inner shell photoionization cross section
	We now prepare to explicitly integrate J x csec over modelled bands */

  fthresh = xtop->freq[0];      //The first frequency for which we have a cross section
  fmax = xtop->freq[xtop->np - 1];      //The last frequency
  pi_rate = 0;                  //Initialise the pi rate - it will be computed through several integrals



  if (mode == 1)                //Modelled version of J
  {
    for (j = 0; j < geo.nxfreq; j++)    //We loop over all the bands
    {
      xpl_alpha = xplasma->pl_alpha[j]; //set the various model parameters to those for this model
      xpl_logw = xplasma->pl_log_w[j];
      xexp_temp = xplasma->exp_temp[j];
      xexp_w = xplasma->exp_w[j];
      if (xplasma->spec_mod_type[j] != SPEC_MOD_FAIL)   //Only bother doing the integrals if we have a model in this band
      {
        f1 = xplasma->fmin_mod[j];      //NSH 131114 - Set the low frequency limit to the lowest frequency that the model applies to
        f2 = xplasma->fmax_mod[j];      //NSH 131114 - Set the high frequency limit to the highest frequency that the model applies to
        if (f1 < fthresh && fthresh < f2 && f1 < fmax && fmax < f2)     //Case 1-
        {
          if (xplasma->spec_mod_type[j] == SPEC_MOD_PL)
          {
            pi_rate += num_int (tb_logpow, fthresh, fmax, pl_qromb);

          }
          else
          {
            pi_rate += num_int (tb_exp, fthresh, fmax, exp_qromb);
          }
        }
        else if (f1 < fthresh && fthresh < f2 && f2 < fmax)     //case 2
        {
          if (xplasma->spec_mod_type[j] == SPEC_MOD_PL)
          {
            pi_rate += num_int (tb_logpow, fthresh, f2, pl_qromb);
          }
          else
          {
            pi_rate += num_int (tb_exp, fthresh, f2, exp_qromb);
          }
        }
        else if (f1 > fthresh && f1 < fmax && fmax < f2)        //case 3
        {
          if (xplasma->spec_mod_type[j] == SPEC_MOD_PL)
          {
            pi_rate += num_int (tb_logpow, f1, fmax, pl_qromb);
          }
          else
          {
            pi_rate += num_int (tb_exp, f1, fmax, exp_qromb);
          }
        }
        else if (f1 > fthresh && f2 < fmax)     // case 4
        {
          if (xplasma->spec_mod_type[j] == SPEC_MOD_PL)
          {
            pi_rate += num_int (tb_logpow, f1, f2, pl_qromb);
          }
          else
          {
            pi_rate += num_int (tb_exp, f1, f2, exp_qromb);
          }
        }
        else                    //case 5 - should only be the case where the band is outside the range for the integral.
        {
          pi_rate += 0;         // Add nothing - bit of a null statement, but makes the code look nice.
        }
      }                         //End of loop to only integrate in this band if there is power
    }

  }
  else if (mode == 2)           //blackbody mode
  {
    fmaxtemp = xtop->freq[xtop->np - 1];        //Set the maximum frequency temporarily to the maximum cross section frequency
    fmax = check_freq_max (fmaxtemp, xplasma->t_r);     /*Check that the requested maximum frequency is sensible - if it is way
                                                           off the end of the wien tail then the integration can fail - reset if necessary. */
    if (fthresh > fmax)         //The threshold for PI is above the maximum frequency of the radiation
    {
      pi_rate = 0.0;
    }
    else                        //We are OK - do the integral
    {
      qromb_temp = xplasma->t_r;
      pi_rate = xplasma->w * num_int (tb_planck, fthresh, fmax, 1.e-4);
    }
  }

  pi_rate = (4 * PI * pi_rate) / PLANCK;        //We multiply by 4pi and divide by photon energy - the division by nu is done in the integrands




  return (pi_rate);
}


/**********************************************************/
/**
 * @brief      The integrand for working out the PI rate in a BB radiation field
 *
 * @param [in] double  freq   the frequency
 * @param [in] void  params   An extra (unused) variable to make it paletable for the gsl integrator
 * @return     value of sigma B_nu/nu at nu
 *
 * @details
 * tb_planck is the function a\nuB\nu and is called by Qromb in order to integrate over the frequency range where the
 *    ionisation cross section is significant. This is the function for ions with a topbase cross section NSH
 *
 * ### Notes ###
 *
 **********************************************************/

double
tb_planck (double freq, void *params)
{
  double answer, bbe;
  bbe = exp ((PLANCK * freq) / (BOLTZMANN * qromb_temp));
  answer = (2. * PLANCK * pow (freq, 3.)) / (pow (VLIGHT, 2));
  answer *= (1 / (bbe - 1));
//      answer*=weight;
  answer *= sigma_phot (xtop, freq);

  answer /= freq;

  return (answer);
}


/**********************************************************/
/**
 * @brief      The integrand for working out the PI rate in a radiation field modelled by a power law
 *
 * @param [in] double  freq   the frequency
 * @param [in] void  params   An extra (unused) variable to make it paletable for the gsl integrator
 * @return     value of sigma J_nu/nu at nu
 *
 * @details
 * This is the integrand for computing the PI rate coefficient when j_nu is described by a power
 * law. The returned value is J_nu x cross section / nu. The power law is computed in log space to
 * allow a larger range of possible values.
 *
 * ### Notes ###
 *
 **********************************************************/

double
tb_logpow (double freq, void *params)
{
  double answer;
  answer = pow (10, xpl_logw + (xpl_alpha - 1.0) * log10 (freq));       //NB - the alpha-1.0 appears because we divide by nu
  answer *= sigma_phot (xtop, freq);    // and finally multiply by the cross section.

  return (answer);
}


/**********************************************************/
/**
 * @brief      The integrand for working out the PI rate in a radiation field modelled by an exponential
 *
 * @param [in] double  freq   the frequency
 * @param [in] void  params   An extra (unused) variable to make it paletable for the gsl integrator
 * @return     value of sigma J_nu/nu at nu
 *
 * @details
 * This is the integrand for computing the PI rate coefficient when j_nu is described by an exponential
 * law. The returned value is J_nu x cross section / nu.
 *
 * ### Notes ###
 *
 **********************************************************/

double
tb_exp (double freq, void *params)
{
  double answer;

  answer = xexp_w * exp ((-1.0 * PLANCK * freq) / (BOLTZMANN * xexp_temp));
  answer *= sigma_phot (xtop, freq);    // and finally multiply by the cross section.

  answer /= freq;               //then finally finally divide by the frequency
  return (answer);
}
