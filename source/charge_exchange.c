
/***********************************************************/
/** @file  charge_exchanger.c
 * @author nsh
 * @date   May, 2020
 *
 * @brief  charge_exchange contains routines to compute charge exchange rates
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
 * @brief      computes the rates due to charge exchange
 *
 * @param [in] double  T  temperature in K of the cell in question
 * @return     0 on success - the rates are stored in a global variable
 *
 * @details
 * This routine computes the charge exchange recombination coeffcieicnts for
 * most ions, and the charge exchange ionization coefficients for a few
 * ions. 
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

int
compute_ch_ex_coeffs (double T)
{
  int n, m;
  double temp;
  double T_used;

  /*first compute the (metal) recombination rates. There is an approximation to use for ions without a rate (generally those
     above 4x ionized) so we compute a rate for all ions other thasn the neutral ions and hydrogen */

  for (n = 0; n < nions; n++)
  {
    if (ion[n].z != 1 && ion[n].istate > 1)     //No rates for hydrogen or neutral ions
    {
      if (ion[n].n_ch_ex < 0)   //We dont have a proper rate, so use the approximation
      {
        charge_exchange_recomb_rates[n] = 1.92e-9 * (ion[n].istate - 1);        //Approxmation - multiplied by charge 
      }
      else
      {
        m = ion[n].n_ch_ex;
        if (T < charge_exchange[m].tmin)        /*Follow cloudy here - outside the applicable temperature range, use the rate for the last temperature */
          T_used = charge_exchange[m].tmin;
        else if (T > charge_exchange[m].tmax)
          T_used = charge_exchange[m].tmax;
        else
          T_used = T;
        temp = charge_exchange[m].c * exp (charge_exchange[m].d * T_used / 1e4);
        temp = temp + 1.0;
        temp = temp * pow (T_used / 1e4, charge_exchange[m].b);
        charge_exchange_recomb_rates[n] = temp * charge_exchange[m].a * 1e-9;
      }
    }
    else
    {
      charge_exchange_recomb_rates[n] = 0.0;    //Set the rate to zero if hydrogen or neutral ion
    }
  }

  //Now we make ionization rates - do it in a slightly different way because there are only a small number of rates and we dont approximate
  for (n = 0; n < n_charge_exchange; n++)
  {
    if (ion[charge_exchange[n].nion2].z == 1)   //The recombining ion is hydrogen so we have an ionization rate to make
    {
      if (T < charge_exchange[n].tmin)
        T_used = charge_exchange[n].tmin;
      else if (T > charge_exchange[n].tmax)
        T_used = charge_exchange[n].tmax;
      else
        T_used = T;

      temp = charge_exchange[n].c * exp (charge_exchange[n].d * T_used / 1e4);
      temp = temp + 1.0;
      temp = temp * pow (T_used / 1e4, charge_exchange[n].b);
      charge_exchange_ioniz_rates[n] = temp * charge_exchange[n].a * 1e-9;
      charge_exchange_ioniz_rates[n] = charge_exchange_ioniz_rates[n] * exp (-1. * charge_exchange[n].delta_e_ovr_k * 1e4 / T); //Apply the boltzman factor for the actual temperature
    }
    else
    {
      charge_exchange_ioniz_rates[n] = 0.0;     //Set the rate to zero          
    }
  }
  return (0);
}

/**********************************************************/
/**
 * @brief      computes the heating effect due to charge exchange
 *
 * @param [in] WindPtr  one   pointer to cell
 * @param [in] double  t_e  temperature in K of the cell in question
 * @return     the heating rate in erg/s/cm^3
 *
 * @details
 * This routine computes heating effect of charge exchange ionization
 * and recombination. It first computes the coefficients, then the 
 * recombination heating then ionization. There will be some 'cooling'
 * terms but predominantly this is a heating effect. 
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
ch_ex_heat (one, t_e)
     WindPtr one;               // Pointer to the current wind cell - we need the cell volume, this is not in the plasma structure
     double t_e;                //Current electron temperature of the cell

{
  double x;                     //The returned variable
  int nplasma;                  //The cell number in the plasma array
  PlasmaPtr xplasma;            //pointer to the relevant cell in the plasma structure
  int n;                        //loop pointers
  double nh1, nh2;              //the neutral and ionizaed hydrogen number densities.


  nplasma = one->nplasma;       //Get the correct plasma cell related to this wind cell
  xplasma = &plasmamain[nplasma];       //copy the plasma structure for that cell to local variable
  x = 0;                        //zero the luminosity
  nh1 = nh2 = 0;


  compute_ch_ex_coeffs (t_e);   //Calculate the charge exchange rate coefficients for this cell

  for (n = 0; n < nions; n++)   //We need the hydrogen densities, get them here via a rather belt and braces approach
  {
    if (ion[n].z == 1 && ion[n].istate == 1)
    {
      nh1 = xplasma->density[n];
    }
    if (ion[n].z == 1 && ion[n].istate == 2)
    {
      nh2 = xplasma->density[n];
    }
  }

  /* first compute the effect for the recombination process */


  for (n = 0; n < nions; n++)
  {
    if (ion[n].z != 1 && ion[n].istate > 1)     //No rates for hydrogen or neutral ions
    {
      if (ion[n].n_ch_ex < 0)   //We dont have a proper rate, so use the approximation
      {
        x += xplasma->vol * charge_exchange_recomb_rates[n] * nh1 * xplasma->density[n] * 2.86 * (ion[n].istate - 1) * EV2ERGS;
      }
      else
      {
        x += xplasma->vol * charge_exchange_recomb_rates[n] * nh1 * xplasma->density[n] * charge_exchange[ion[n].n_ch_ex].energy_defect;
      }
    }
  }

  /* and now the ionization effect */

  for (n = 0; n < n_charge_exchange; n++)
    if (ion[charge_exchange[n].nion2].z == 1)   //A hydrogen recomb - metal ionization rate
    {
      x +=
        xplasma->vol * charge_exchange_ioniz_rates[n] * nh2 * xplasma->density[charge_exchange[n].nion1] * charge_exchange[n].energy_defect;
    }
  return (x);
}
