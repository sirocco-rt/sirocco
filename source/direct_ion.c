
/***********************************************************/
/** @file  direct_ion.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  direct_ion contains the routines relating to direct (collisional) ionization and
 * threebody recombination
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_expint.h>  //We need this gsl library to evaluate the first exponential integral

#include "atomic.h"
#include "sirocco.h"
//OLD #include "recipes.h"

/* Ratio of hnu / kT beyond which we don't bother calculating
   see #197 */
#define ALPHABIG_DIRECT_ION 100.


/**********************************************************/
/**
 * @brief      computes the rates due to direct (collisional)
 *   ionization from Dere data.
 *
 * @param [in] double  T  temperature in K of the cell in question
 * @return     0 on success
 *
 * @details
 * This routine controls the calculation of the direct ionization
 * coefficients for each ionic species at a temperature T. They are
 * stored in a global array di_coeffs - usually for use in calculating
 * the ionization balance. The actual calculations are done in a
 * seperate routine, q_ioniz_dere
 * Data is from  Ionization rate coefficients for elements H to Zn (Dere+, 2007)
 * Table  J_A_A_466_771_table29:
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

int
compute_di_coeffs (T)
     double T;
{
  int n;

  for (n = 0; n < nions; n++)   //Loop over all ions
  {
    if (ion[n].dere_di_flag == 0)       //If there isnt DI data for the ion
    {
      di_coeffs[n] = 0.0;       //Set the DI rate to zero
    }
    else
    {
      di_coeffs[n] = q_ioniz_dere (n, T);       //Otherwise compute the rate
    }
  }                             //End of loop over ions

  return (0);
}

/**********************************************************/
/**
 * @brief      Returns the collisional ionization rate coefficient for a given ion/temperature
 *
 * @param [in] int  nion   Index into ion structure for current ion under analysis
 * @param [in] double  t_e   Electron temperature
 * @return     the collisional ionization rate coefficient
 *
 * @details
 * This uses data from Dere 07 to calculate a collisional ionization rate coefficient.
 * Data is provided as spline fits over a range of temperatures and so an interpolation is
 * required to arrive at the correct data for a given temperature. There is also an
 * exponential integral which is done using a gsl routine.
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
q_ioniz_dere (nion, t_e)
     int nion;
     double t_e;
{
  double coeff, t, scaled_t;
  double exp_int, dt, drdt, rate;
  int nrec, imax, imin, i;

  imin = imax = -1;             //Set to avoid compiler warning

  nrec = ion[nion].nxderedi;    // find the correct coefficient via an index in the ion structure


  t = (BOLTZMANN * t_e) / (dere_di_rate[nrec].xi * EV2ERGS);
  scaled_t = 1.0 - ((log (2.0)) / (log (2.0 + t)));     //The fit is in 'reduced temperature' units

  /* 1/t is (hnu) / (k t_e). when this ratio is >100 we return 0, see #197 The integral will fail if we try to do it */
  if ((1.0 / t) > ALPHABIG_DIRECT_ION)
    return 0.0;

  if (scaled_t < dere_di_rate[nrec].temps[0])   //we are below the range of DI data data
  {

    Log_silent ("compute_di_coeffs: Requested temp %e is below limit of data for ion %i(Tmin= %e)\n",
                scaled_t, nion, dere_di_rate[nrec].temps[0]);
    imax = 1;
    imin = 0;                   //We will try to extrapolate
  }

  else if (scaled_t >= dere_di_rate[nrec].temps[dere_di_rate[nrec].nspline - 1])        //we are above the range of GS data
  {
    Log_silent ("compute_di_coeffs: Requested temp %e is above limit (%e) of data for ion %i\n",
                scaled_t, dere_di_rate[nrec].temps[dere_di_rate[nrec].nspline - 1], nion);
    imax = dere_di_rate[nrec].nspline - 1;
    imin = dere_di_rate[nrec].nspline - 2;      //We will try to extrapolate.
  }

  else                          //We are within the range of tabulated data
  {
    for (i = 0; i < dere_di_rate[nrec].nspline - 1; i++)
    {
      if (dere_di_rate[nrec].temps[i] <= scaled_t && scaled_t < dere_di_rate[nrec].temps[i + 1])        //We have bracketed the correct temperature
      {
        imin = i;
        imax = i + 1;           //Set the indices to straddle the required temperature
      }
    }
  }

  //Interpolate in log space

  drdt =
    ((dere_di_rate[nrec].rates[imax]) - (dere_di_rate[nrec].rates[imin])) / ((dere_di_rate[nrec].temps[imax]) -
                                                                             (dere_di_rate[nrec].temps[imin]));
  dt = ((scaled_t) - (dere_di_rate[nrec].temps[imin]));
  rate = dere_di_rate[nrec].rates[imin] + drdt * dt;

  coeff = pow (t, -0.5) * pow (dere_di_rate[nrec].xi, -1.5) * rate;     //perform the interpolation

  if (exp (-1.0 / t) < (1.0 / VERY_BIG))        //Check to make sure the exponential integral will give a sensible result
  {
    exp_int = 0.0;              //If its too small, it will just be zero
  }
  else
  {
    exp_int = gsl_sf_expint_E1 (1.0 / t);       //Evaluate the first exponential integral using gsl library.
  }

  coeff *= exp_int;             //Finalise the coefficient.

  return (coeff);
}

/**********************************************************/
/**
 * @brief      computes the total cooling due collisional ionization
 *
 * @param [in] WindPtr  one   pointer to cell
 * @param [in] double  t_e   electron temperature
 * @return     The total di cooling
 *
 * @details
 * This routine computes the total CMF cooling rate in a cell due to collisional
 * ionization. This is just the rate, times the ionization thereshold of the
 * state being recombined into.
 *
 * ### Notes ###
 *
 * This calculation is a CMF calculation.
 *
 **********************************************************/

double
total_di (one, t_e)
     WindPtr one;
     double t_e;

{
  double cooling_rate;
  int nplasma;
  PlasmaPtr xplasma;
  int n;


  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  cooling_rate = 0;


  compute_di_coeffs (t_e);      //Calculate the DI rate coefficients for this cell


  for (n = 0; n < nions; n++)
  {
    if (ion[n].dere_di_flag)
    {
      cooling_rate += xplasma->vol * xplasma->ne * xplasma->density[n] * di_coeffs[n] * dere_di_rate[ion[n].nxderedi].xi * EV2ERGS;
    }
  }
  return (cooling_rate);
}





/**********************************************************/
/**
 * @brief      computes the rate coefficient for three body recombination from Dere data.
 *
 * @param [in] double  T   Temperature in K
 * @return     0 on success
 *
 * @details
 * This routine calculates the three body recombination rate coefficient for all
 * ions in a cell at temperature T. The resulting rate coefficients are stored
 * in a global array for use by the calling routine - normally matrix_ion as
 * part of the search for an ionization balance.
 * It works by essentially applying the Milne relation to collisional ionization
 * data.
 * This is a driver routine, the hard work is done in a seperate routine (q_recomb_dere)
 *
 * ### Notes ###
 * 1508 JM coded
 *
 **********************************************************/

int
compute_qrecomb_coeffs (T)
     double T;
{
  int n, nvmin, ntmin;
  struct topbase_phot *xtop;

  xtop = NULL;                  //Avoid compiler warning, though it seem odd we need to do this at all.

  for (n = 0; n < nions; n++)   //We need to generate data for the ions doing the recombining.
  {
    if (ion[n].istate > 1)      //We cannot recombine to anything from the neutral iom
    {
      if (ion[n - 1].dere_di_flag == 0) //If there isn't a collisional ionization rate for the lower ion
      {
        qrecomb_coeffs[n] = 0.0;        //We can't do anything - so set the rate to zero
      }

      else
      {
        /* we need to know about the bound-free jump, so we need the details
           from the ground state cross-section for for the ion below this one */

        ntmin = ion[n - 1].ntop_ground; /* We only ever use the ground state cont_ptr.
                                           This is for topbase */
        nvmin = ion[n - 1].nxphot;      //VFKY

        if (ion[n - 1].phot_info > 0)   //topbase or hybrid VFKY (GS)+TB excited
        {
          xtop = &phot_top[ntmin];
        }
        else if (ion[n - 1].phot_info == 0)     // verner
        {                       //just the ground state ionization fraction.
          xtop = &phot_top[nvmin];
        }
        else                    //We dont have any information about the bound-free jump
        {
          Error ("compute_qrecomb_coeffs: no coll ionization data for ion %i\n", n - 1);
        }

        /* this will return 0 if there aren't any coeffs */
        qrecomb_coeffs[n] = q_recomb_dere (xtop, T);    //We have everything we need, so calculate the coeff
      }

    }                           //End of if statement for neutral ions
    else                        //We are a neutral ion - so there can be no recombination
    {
      qrecomb_coeffs[n] = 0.0;
    }
  }                             //End of loop over ions

  return (0);
}


/**********************************************************/
/**
 * @brief      Returns the collisional recombination co-efficient
 *
 * @param [in] struct topbase_phot a pointer to the topbdae style PI cross section
 * @param [in] double  electron_temperature
 * @return     ??? RETURNS ???
 *
 * @details
 * This equation comes from considering TE and getting the expression
 * q_recomb = 2.07e-16 * gl/gu * exp(E/kT) * q_ioniz
 * then substituting the above expression for q_ioniz. We get the
 * ionization threshold from the supplied PI cross section.
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
q_recomb_dere (cont_ptr, electron_temperature)
     struct topbase_phot *cont_ptr;
     double electron_temperature;
{
  int nion;
  double u0;
//OLD  double gaunt, coeff;
  double coeff;
  double root_etemp;

//OLD  gaunt = 0.1;                  //for now - from Mihalas for hydrogen
  u0 = cont_ptr->freq[0] * H_OVER_K / electron_temperature;
  nion = cont_ptr->nion;

  /* when hnu / kT ratio is >100 we return 0, see #197 */
  if (u0 > ALPHABIG_DIRECT_ION)
    return (0.0);

  /* if ion[n].dere_di_flag == 1 then we have direct ionization data for this ion
     only do this if it is the ground state */
  if (ion[nion].dere_di_flag == 1 && xconfig[cont_ptr->nlev].ilv == 1)
  {
    root_etemp = sqrt (electron_temperature);
    coeff = 2.07e-16 / (root_etemp * root_etemp * root_etemp);

    /* We now multiply by the ratio of the multiplicity of the ground states.
       Dere's data comes from a vast range of sources, some PI cross sections
       include excited states, others dont. We will make the approximation
       of using juast the ground state multiplicities. If there is an error
       introduced here it will be small */
    coeff *= ion[nion].g / ion[nion + 1].g;

    coeff *= exp (u0);
    coeff *= q_ioniz_dere (nion, electron_temperature);


    /* do a sane check here, as there's an exponential which could blow up */
    if (sane_check (coeff))
    {
      Error ("q_recomb is %8.4e for ion %i at temperature %8.4e, setting to zero\n", coeff, nion, electron_temperature);
      coeff = 0.0;
    }
  }
  else
    coeff = 0.0;

  return (coeff);
}



/**********************************************************/
/**
 * @brief      returns the collisional ionization rate co-efficient
 *
 * @param [in] struct topbase_phot  - pointer to a Topbase style cross section - u
 * @param [in] double  electron_temperature - electron temperature
 * @return     collisional ionization rate coefficieient
 *
 * @details
 * This returns the collisional ionization co-efficient
 * Calculated following equation 5-79 of Mihalas or from data from
 * Dere 2007.
 *
 * ### Notes ###
 * This is extrememly similar to the routine compute_di_coeffs
 * apart from the fact that it computes just one rate coefficient,
 * and can use an approximate formula in the absence of data. It
 * is an earlier routine used only in macroatom calculations.
 *
 **********************************************************/

double
q_ioniz (cont_ptr, electron_temperature)
     struct topbase_phot *cont_ptr;
     double electron_temperature;
{
  double coeff;
  double gaunt;
  double u0;
  int nion;

  /* these next two quantities only used in Hydrogen, no Dere data case */
  u0 = cont_ptr->freq[0] * H_OVER_K / electron_temperature;
  nion = cont_ptr->nion;
  gaunt = 0.1 * ion[nion].z;    //for now - from Mihalas for hydrogen and Helium

  /* when hnu / kT ratio is >100 we return 0, see #197 */
  if (u0 > ALPHABIG_DIRECT_ION)
    return (0.0);


  /* if ion[n].dere_di_flag == 1 then we have direct ionization data for this ion
     only do this if it is the ground state */
  if (ion[nion].dere_di_flag == 1 && xconfig[cont_ptr->nlev].ilv == 1 && ion[nion].macro_info == 0)
  {
    coeff = q_ioniz_dere (nion, electron_temperature);
  }

  /* let's only apply the approximation from Mihalas for Hydogen and Helium */
  else if (ion[nion].z < 3 && ion[nion].macro_info == 1)
  {
    coeff = 1.55e13 / sqrt (electron_temperature) * gaunt * cont_ptr->x[0] * exp (-1. * u0) / u0;
  }

  /* otherwise return 0 */
  else
    coeff = 0.0;


  return (coeff);
}


/**********************************************************/
/**
 * @brief      Returns a collisional recombination rate co-efficient
 *
 * @param [in out] struct topbase_phot - pointer to a topbase style PI cross section
 * @param [in out] double  electron_temperature - electron temperature
 * @return     Collisional recombination rate coefficient
 *
 * @details
 * This equation comes from considering TE and getting the expression
 * q_recomb = 2.07e-16 * gl/gu * exp(E/kT) * (T_e**-3/2) * q_ioniz
 * then substituting the above expression for q_ioniz.
 *
 * ### Notes ###
 * This routine is very similar to compute_qrecomb_coeffs apart from
 * the fact that this routine only computes one recomb coefficient
 * and applies an approximate calculation where dere data is
 * missing. It is used in macroatom calculations.
 *
 **********************************************************/

double
q_recomb (cont_ptr, electron_temperature)
     struct topbase_phot *cont_ptr;
     double electron_temperature;
{
  double coeff;
  double gaunt, u0;
  int nion;

  nion = cont_ptr->nion;
  u0 = cont_ptr->freq[0] * H_OVER_K / electron_temperature;
  gaunt = 0.1 * ion[nion].z;    //for now - from Mihalas for hydrogen and Helium

  /* when hnu / kT ratio is >100 we return 0, see #197 */
  if (u0 > ALPHABIG_DIRECT_ION)
    return (0.0);

  /* if ion[n].dere_di_flag == 1 then we have direct ionization data for this ion
     only do this if it is the ground state */
  /* Still use the Mihalas approximation for macro-atoms */
  if (ion[nion].dere_di_flag == 1 && xconfig[cont_ptr->nlev].ilv == 1 && ion[nion].macro_info == 0)
  {
    coeff = q_recomb_dere (cont_ptr, electron_temperature);
  }

  /* let's only apply the approximation from Mihalas for Hydogen and Helium macro-atoms */
  else if (ion[nion].z < 3 && ion[nion].macro_info == 1)
  {
    coeff = 3.2085e-3 / electron_temperature * gaunt * cont_ptr->x[0];  // normal constants * 1/T times gaunt * cross section

    coeff /= cont_ptr->freq[0] * H_OVER_K;      // divide by h nu / k
    coeff *= xconfig[cont_ptr->nlev].g / xconfig[cont_ptr->uplev].g;
  }

  /* otherwise return 0 */
  else
    coeff = 0.0;


  return (coeff);
}
