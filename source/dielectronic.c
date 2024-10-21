
/***********************************************************/
/** @file  dielectronic.c
 * @author Nick
 * @date   January, 2018
 *
 * @brief  The routines in this file all have to do with dielectronic recombination
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
 * @brief      returns the volumetric dielectronic rate
 *  	coefficients for a given temperature.
 *
 * @param [in] double  temp   the temperature at which to compute the coefficients
 * @return     nothing, but populates the array dr_coeffs
 *
 * @details
 * Takes as its input a temperature, and it populates the dr_coeffs array.
 * called every time a cell needs a set of DR coefficients.
 *
 * ### Notes ###
 * the rates are associated with the ion being recombined into.
 *
 **********************************************************/

int
compute_dr_coeffs (temp)
     double temp;
{
  int n, n1, n2;
  double Adi, Bdi, T0, T1;
  for (n = 1; n < nions + 1; n++)
  {
    if (ion[n].drflag == 0)     //There are no dielectronic coefficients relating to this ion
    {
      dr_coeffs[n] = 0.0;       //So set the coefficients to zero
    }
    else
    {
      n1 = ion[n].nxdrecomb;    //Obtain an index into the drecomb structure for this ion's data
      dr_coeffs[n] = 0.0;       //Zero the coefficient
      if (drecomb[n1].type == DRTYPE_BADNELL)   //It is a badnell type coeefficient
      {
        for (n2 = 0; n2 < drecomb[n1].nparam; n2++)     //Loop over the parameters
        {
          dr_coeffs[n] += (drecomb[n1].c[n2] * exp (-1 * (drecomb[n1].e[n2] / temp)));
        }
        dr_coeffs[n] *= pow (temp, -1.5);
      }
      else if (drecomb[n1].type == DRTYPE_SHULL)        //A schiull type parameter
      {
        Adi = drecomb[n1].shull[0];
        Bdi = drecomb[n1].shull[1];
        T0 = drecomb[n1].shull[2];
        T1 = drecomb[n1].shull[3];
        dr_coeffs[n] = Adi * pow (temp, -1.5) * exp ((-1.0 * T0) / temp);
        dr_coeffs[n] *= (1 + Bdi * exp ((-1.0 * T1) / temp));
      }
      else
      {
        Error ("Compute_dr_coeffs: Unknown DR data rtype for ion %i\n", n);
        dr_coeffs[n] = 0.0;
      }
    }
  }
  return (0);
}



/**********************************************************/
/**
 * @brief      calculates the total luminosity from DR.
 *
 * @param [in] WindPtr  one   Pointer to the wind cell currently under analysis
 * @param [in] double  t_e   The temperature of the cell - this is from the plasma structure
 * @return     the total luminosity of this cell due to dielectronic recombinaions - currently set to zero
 *
 * @details
 * This routine is a dummy routine at the moment - it mirrors the other
 * routines which compute the luminosity due to various cooling routines, but
 * currently always returns zero, because
 * as of the current time we do not know how to work out the luminosity
 * from dielectronic recombination. We know the rate, which is used in ionization
 * calculations but the energy released by each recombination is complex and
 * depends on what happens to the excited ion.
 *
 * ### Notes ###
 * This is issue #186
 *
 **********************************************************/

double
total_dr (one, t_e)
     WindPtr one;               // Pointer to the current wind cell - we need the cell volume, this is not in the plasma structure
     double t_e;                //Current electron temperature of the cell
{
  double x;                     //The returned variable
//OLD  double meanv, meanke;         //The mean velocity and kinetic energy of electrons in the cell
//OLD  int nplasma;                  //The cell number in the plasma array
//OLD  PlasmaPtr xplasma;            //pointer to the relevant cell in the plasma structure
  int n;                        //loop pointers


//OLD  nplasma = one->nplasma;       //Get the correct plasma cell related to this wind cell
//OLD  xplasma = &plasmamain[nplasma];       //copy the plasma structure for that cell to local variable
  x = 0;                        //zero the luminosity


  compute_dr_coeffs (t_e);      //Calculate the DR coefficients for this cell

/* The next two lines calculate the mean K.E. of electrons - it was an early attempt to
  estimate the energy released in a DR. It hugely over estimates the energy release */
//OLD  meanv = pow ((2 * BOLTZMANN * t_e / MELEC), 0.5);
//OLD  meanke = 0.5 * MELEC * meanv * meanv;

/* In ordet to compute tht total luminosity for a cell, we loop over all ions */

  for (n = 1; n < nions + 1; n++)
  {
    if (ion[n].drflag == 0)     //We have no DR for this ion.
    {
      x += 0.0;                 //Add nothing to the sum of coefficients
    }
    else                        //This is the place where one would put code to compute the DR luminosity for a given ion
    {
/* One possibility is to use the mean k.e. of electrons this is based on the idea that an eleectron
		is absorbed, so a guess would be to just say its energy is re-radiated. This is clearly an over estimate
		since some of the enegy is lost to binding energy 	*/
//              x += xplasma->vol * xplasma->ne * xplasma->density[n + 1] * dr_coeffs[n] * meanke;
/* A different estimate is to use the ionization potential of the ion - this is order of magnitude at best */
//              x += xplasma->vol * xplasma->ne * xplasma->density[n + 1] * dr_coeffs[n] * ion[n].ip;
      x += 0.0;                 //At present, we just set it to zero - obviously an underestimate!
    }
  }
  return (x);
}
