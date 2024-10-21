
/***********************************************************/
/** @file  lines.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  Subroutines associated with with resonance line radiation
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"
//OLD #include "recipes.h"


/**********************************************************/
/**
 * @brief
 * Calculate the total line emission from a wind cell
 *
 *
 * @param [in] PlasmaPtr  xplasma   A specific plasma cell in the wind
 * @param [in] double  f1   A minimum frequncy
 * @param [in] double  f2   A maximum frequency
 * @return     The line lumionisity of the cell (in erg/s)
 *
 * @details
 * Using t_e and the densities of ions in a cell the routine
 * simply calculates the total amount of line emission in
 * the cell
 *
 * ### Notes ###
 *
 * The routine simply use f1 and f2 to define the
 * range of lines to calculate the luminosity for and then
 * calls lum_lines
 *
 * @bug It is not exactly clear why two routines total_line_emission
 * and lum_lines is needed, as lum_lines appears to be called
 * only from total_line_emission.  ksl - 180509
 * 
 * Total line emission in this function is calculated CMF. If called in
 * the cooling routines then this is left in the CMF. If called in the photon
 * generation routines, this is converted to the observer frame.
 *
 **********************************************************/

double
total_line_emission (xplasma, f1, f2)
     PlasmaPtr xplasma;         /* WindPtr to a specific cell in the wind */
     double f1, f2;             /* Minimum and maximum frequency */
{

  double lum;
  double t_e;

  t_e = xplasma->t_e;

  if (t_e <= 0 || f2 < f1)
    return (0);

  /* Update nline_min and nline_max in atomic.h which define which
   * lines lie in the frequency range f1-f2 in the frepeuncy ordered
   * version of the lines
   */

  limit_lines (f1, f2);

  lum = lum_lines (xplasma, nline_min, nline_max);



  return (lum);

}


/**********************************************************/
/**
 * @brief      Calculate the line lumiosity of lines between
 * nmin and nmax of the frequency ordered list of lines
 *
 * @param [in] PlasmaPtr  xplasma   A wind cell
 * @param [in] int  nmin   The minimum number of a line in the frequency ordered list
 * @param [in] int  nmax   The maximum number of a line in the frequency ordered list
 * @return     The total line luminosity between element nmin and nmax of the frequncy
 * ordered list of lines
 *
 * @details
 * Using densities in the plasma cell associated with this wind cell, calculate the
 * total line luminosity.
 *
 * ### Notes ###
 * The individual line luminosities are stored in lin_ptr[n]->pow
 *
 * This is a co-moving frame calculation.                            
 *
 **********************************************************/

double
lum_lines (xplasma, nmin, nmax)
     PlasmaPtr xplasma;
     int nmin, nmax;            /* The min and max index in lptr array for which the power is to be calculated */
{
  int n;
  double lum, x, z;
  double dd, d1, d2;
  double q;
  double t_e;
  double foo1, foo2, foo3, foo4;
  t_e = xplasma->t_e;
  lum = 0;
  for (n = nmin; n < nmax; n++)
  {
    dd = xplasma->density[lin_ptr[n]->nion];

    if (dd > LDEN_MIN)
    {                           /* potentially dangerous step to avoid lines with no power */
      two_level_atom (lin_ptr[n], xplasma, &d1, &d2);
      x = foo1 = lin_ptr[n]->gu / lin_ptr[n]->gl * d1 - d2;

      z = exp (-H_OVER_K * lin_ptr[n]->freq / t_e);


      //Next lines required if want to use escape probabilities

      q = 1. - scattering_fraction (lin_ptr[n], xplasma);

      x *= foo2 = q * a21 (lin_ptr[n]) * z / (1. - z);

      x *= foo3 = PLANCK * lin_ptr[n]->freq * xplasma->vol;
      if (geo.line_mode == LINE_MODE_ESC_PROB)
        x *= foo4 = p_escape (lin_ptr[n], xplasma);     // Include effects of line trapping
      else
      {
        foo4 = 0.0;             // Added to prevent compilation warning
      }

      lum += lin_ptr[n]->pow = x;
      if (x < 0)
      {
        Log
          ("lum_lines: foo %10.3g (%10.3g %10.3g %10.3g) %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g\n",
           foo1, d1, d2, dd, foo2, foo3, foo4, lin_ptr[n]->el, xplasma->t_r, t_e, xplasma->w);
      }
      if (sane_check (x) != 0)
      {
        Error ("total_line_emission:sane_check %e %e\n", x, z);
      }
    }
    else
    {
      lin_ptr[n]->pow = 0;
    }
  }


  return (lum);
}



struct lines *old_line_ptr;
double old_ne, old_te, old_w, old_tr, old_dd;
double old_d1, old_d2, old_n2_over_n1;

/**********************************************************/
/**
 * @brief      calculates the ratio n2/n1 and gives the individual
 * densities for the states of a two level atom.
 *
 * @param [in] struct lines *  line_ptr   The line of interest
 * @param [in] PlasmaPtr  xplasma   The plasma cell of interest
 * @param [out] double *  d1   The calculated density of the lower level for the line of interest
 * @param [out] double *  d2   The calculated density of the upper levl
 * @return     The density ratio d2/d1
 *
 * @details
 * In the two level approximation,
 *
 * n2 / n1 = ( c12 + g2/g1 c*c / (2 h nu * nu * nu )  * A21 J12) /
 * (c21 + A21 + c*c / (2 h nu * nu * nu )  * A21 J21)
 *
 * In the on-the-spot approx. we assume,.
 *
 * J21= W  (2 h nu * nu * nu )  / (exp (h nu / k T(rad)) -1)
 *
 * and this is what is calculated here
 *
 * ### Notes ###
 * This routine is not (should not be) called for macro atoms.
 * The program will exit if this happens
 *
 **********************************************************/


double
two_level_atom (line_ptr, xplasma, d1, d2)
     struct lines *line_ptr;
     PlasmaPtr xplasma;
     double *d1, *d2;
{
  double a, a21 ();
  double q, q21 (), c12, c21;
  double freq;
  double g2_over_g1;
  double n2_over_n1;
  double n1_over_ng;
  double n2_over_ng;
  double z;
  int gg;
  double xw;
  double ne, te, w, tr, dd;
  int nion;
  double J;                     //Model of the specific intensity



  if (line_ptr->macro_info == TRUE && geo.rt_mode == RT_MODE_MACRO && geo.macro_simple == FALSE)
  {
    Error ("Calling two_level_atom for macro atom line. Abort.\n");
    Exit (0);
  }

/* Move variables used in the calculation from the xplasma structure into subroutine variables */
  ne = xplasma->ne;
  te = xplasma->t_e;
  tr = xplasma->t_r;
  w = xplasma->w;
  nion = line_ptr->nion;
  dd = xplasma->density[nion];

  /* Calculate the number density of the lower level for the transition using the partition function */
  ;
  if (ion[nion].nlevels > 0)
  {
    dd *= xconfig[ion[nion].firstlevel].g / xplasma->partition[nion];
  }

  if (old_line_ptr == line_ptr && old_ne == ne && old_te == te && old_w == w && old_tr == tr && old_dd == dd)
  {                             // Then there is no need to recalculate eveything
    *d1 = old_d1;
    *d2 = old_d2;
  }
  else
  {
    if (line_ptr->el == 0.0)
    {                           // Then the lower level is the ground state

/* For a ground state connected transition we correct for the partition
function in calculating the density of the lower level, and then we
make an on-the-spot approximation for the upper level.  The only reason
this is a "improvement" over the numbers available from the partition
function directly is the allowance for collisions, and also for the
possibility that not all lines have upper levels that are included
in the configuration structure. 01dec ksl */

      a = a21 (line_ptr);
      q = q21 (line_ptr, te);
      freq = line_ptr->freq;
      g2_over_g1 = line_ptr->gu / line_ptr->gl;


      c21 = ne * q;
      c12 = c21 * g2_over_g1 * exp (-H_OVER_K * freq / te);


      z = (VLIGHT * VLIGHT) / (2. * PLANCK * freq * freq * freq);       //This is the factor which relates the A coefficient to the b coefficient


      /* we call mean intensity with mode 1 - this means we are happy to use the
         dilute blackbody approximation even if we havent run enough spectral cycles
         to have a model for J */
      J = mean_intensity (xplasma, freq, MEAN_INTENSITY_BB_MODEL);

      /* this equation is equivalent to equation 4.29 in NSH's thesis with the
         einstein b coefficients replaced by a multiplied by suitable conversion
         factors from the einstein relations. */
      n2_over_n1 = (c12 + g2_over_g1 * a * z * J) / (c21 + a * (1. + (J * z)));




      *d1 = dd;
      *d2 = *d1 * n2_over_n1;

    }
    else
    {                           // The transition has both levels above the ground state

/*
 * In the event that both levels are above the ground state, we assume
 * that the upper level population is given by an on-the-spot approximation.
 * We make the same assumption for the lower level, unless the lower level
 * is matastable in which case we set the weight to 1 and force equlibrium
*/

      gg = ion[line_ptr->nion].g;
      z = w / (exp (line_ptr->eu / (BOLTZMANN * tr)) + w - 1.);
      n2_over_ng = line_ptr->gu / gg * z;

/* For lower level, use an on the spot approximation if the lower level has a short radiative lifetive;
Othewise, assert that the lower level is metastable and set the radiative weight to 1
ERROR -- At present I don't believe what one should do with metastable lower levels is
ERROR -- worked out, either in the program... where are we getting radiative rates
ERROR -- or conceptually
07mar - ksl - We still need to determine whether this makes sense at all !!
*/

      xw = w;                   // Assume all lower levels are allowed at present

      z = xw / (exp (line_ptr->el / (BOLTZMANN * tr)) + xw - 1.);
      n1_over_ng = line_ptr->gl / gg * z;

      *d1 = dd * n1_over_ng;
      *d2 = dd * n2_over_ng;
      n2_over_n1 = n2_over_ng / n1_over_ng;
    }


    old_line_ptr = line_ptr;
    old_ne = ne;
    old_te = te;
    old_w = w;
    old_tr = tr;
    old_dd = dd;
    old_d1 = (*d1);
    old_d2 = (*d2);
    old_n2_over_n1 = n2_over_n1;
  }

  return (old_n2_over_n1);

}



/**********************************************************/
/**
 * @brief
 * Calculate the total line absorption x-section for a specific transition
 * allowing for stimulated emission
 *
 * @param [in] struct lines *  line_ptr   The line of interest
 * @param [int] PlasmaPtr  xplasma   The plasma cell of interest
 * @return     The total x-section for the line
 *
 * @details
 * The total x-section for a line in the two level approximation
 *
 * ### Notes ###
 * The effects of stimulated emission are included
 *
 **********************************************************/

double
line_nsigma (line_ptr, xplasma)
     struct lines *line_ptr;
     PlasmaPtr xplasma;
{
  double d1, d2, x;

  two_level_atom (line_ptr, xplasma, &d1, &d2);

  x = (d1 - line_ptr->gl / line_ptr->gu * d2);
  x *= PI_E2_OVER_MC * line_ptr->f;
  return (x);
}







/**********************************************************/
/**
 * @brief      calculate the fraction of excited
state atoms which correspond to scattered photons, i.e. the portion which are
excited by radiation and return to the ground state via spontaneous emission.
 *
 * @param [in] struct lines *  line_ptr   The line of interest
 * @param [in] PlasmaPtr  xplasma   The cell of interest
 * @return     The fraction of excitations which result in a scattering
 * event
 *
 * @details
 * In the radiative transfer equation for lines, we separate the fraction of photons
 * which are "absorbed" and the fraction which are "scattered".   The results
 * depend on the line mode
 *
 * * If line_mode==0, the atomosphere is a completely absorbing and no photons
 * 		will be scattered.  In this mode, assuming the wind is a source
 * 		of emission, the emissivity will be the Einstein A coefficient
 * * If line_mode==1, the atmosphere is completely scattering, there will be no
 * 		interchange of energy between the photons and the electrons
 * 		as a result of radiation transfer
 * * If line_mode==2, then a simple single scattering approximation is applied in which
 * 		case the scattered flux is just  A21/(C21+A21*(1-exp(-h nu/k T_e).
 * * If line_mode==3, then radiation trapping is included as well.  The basic idea
 * 		is to calculate the average number of scatters in a single
 * 		interaction and there is heat lost in each of these scatters.
 *
 * ### Notes ###
 * It may be more efficient to combine several of these routines in view
 * of the fact that the exp is calculated several times
 *
 * @bug The modes used by scattering fraction are hardwired, which is not
 * our standard approach.  Note however there is a lot about the line_mode
 * which is quite complicated, not the least of which being that there is
 * a related variable scatter_mode
 *
 **********************************************************/

double
scattering_fraction (line_ptr, xplasma)
     struct lines *line_ptr;
     PlasmaPtr xplasma;
{
  double q, escape;
  double a, c, z;
  double sf;
  double ne, te;
  double w;                     /* the radiative weight, and radiation tempeature */

  if (geo.line_mode == LINE_MODE_ABSORB)
    return (0.0);               //purely absorbing atmosphere

  else if (geo.line_mode == LINE_MODE_SCAT)
    return (1.);                //purely scattering atmosphere

  //Populate variable from previous calling structure
  ne = xplasma->ne;
  te = xplasma->t_e;
  w = xplasma->w;

  c = (-H_OVER_K * line_ptr->freq / te);
  a = exp (c);
  z = 1.0 - a;
  a = a21 (line_ptr);
  c = ne * q21 (line_ptr, te) * z;
  q = c / (a + c);              //q == epsilon in Rybicki and elsewhere

  if (geo.line_mode == LINE_MODE_SINGLE_SCAT)
    return (1 - q);             //single scattering atmosphere

  else if (geo.line_mode == LINE_MODE_ESC_PROB)
  {                             // atmosphere with  line trapping

    escape = p_escape (line_ptr, xplasma);
    //The following is exact
    sf = (1. - q) * escape / (q + escape * (1. - q));
    if (sane_check (sf))
    {
      Error ("scattering fraction:sane_check sf %8.2e q %8.2e escape %8.2e w %8.2e\n", sf, q, escape, w);
    }
    return (sf);
  }

  else
  {                             // Unknown treatment of line radiation

    Error ("scattering_fraction: Cannot handle %d line_mode\n", geo.line_mode);
    Exit (0);
    return (0);
  }
}



struct lines *pe_line_ptr;
double pe_ne, pe_te, pe_dd, pe_dvds, pe_w, pe_tr;
double pe_escape;

/**********************************************************/
/**
 * @brief      Estimate the esccapte probability for a line
 * in a plasma cell
 *
 * @param [in] struct lines *  line_ptr   The element in the lines structure of a line
 * @param [in] PlasmaPtr  xplasma   An element of the plasma structure
 * @return     The escape probability of the line in that plasma cell
 *
 * @details
 * Estimate the escape probability using the Sobolev approximation
 *
 * ### Notes ###
 *
 **********************************************************/

double
p_escape (line_ptr, xplasma)
     struct lines *line_ptr;
     PlasmaPtr xplasma;
{
  double tau, two_level_atom ();
  double escape;
  double ne, te;
  double dd;                    /* density of the relevent ion */
  double dvds;
  double w, tr;                 /* the radiative weight, and radiation tempeature */
  WindPtr one;

  ne = xplasma->ne;
  te = xplasma->t_e;
  tr = xplasma->t_r;
  w = xplasma->w;

  dd = xplasma->density[line_ptr->nion];

  one = &wmain[xplasma->nwind];
  dvds = one->dvds_ave;

  if (dvds <= 0.0)
  {
    Error ("Warning: p_escape: dvds <=0 \n");
    return (0.0);
  }


  if (pe_line_ptr != line_ptr || pe_ne != ne || pe_te != te || pe_dd != dd || pe_dvds != dvds || pe_w != w || pe_tr != tr)
  {

    tau = sobolev (one, one->x, dd, line_ptr, dvds);

    escape = p_escape_from_tau (tau);


    pe_line_ptr = line_ptr;
    pe_ne = ne;
    pe_te = te;
    pe_dd = dd;
    pe_dvds = dvds;
    pe_w = w;
    pe_tr = tr;

    pe_escape = escape;
  }


  return (pe_escape);
}




/**********************************************************/
/**
 * @brief      Given an optical depth estimate the escape
 * probability
 *
 * @param [in] double  tau   An optical depth
 * @return     The escape probability
 *
 * @details
 * p_escape_from_tau calculates the probability of escape
 *   given an actual tau. It simply returns the equation
 *
 *   (1. - exp (-tau)) / tau;
 *
 *   Except for high and low tau where it returns 1/tau and
 *   1.0 respectively. This is used by p_escape above, which
 *   calculates the sobolev escape probability, and
 *   also by the anisotropic scattering routines.
 *
 *
 * ### Notes ###
 *
 **********************************************************/

double
p_escape_from_tau (tau)
     double tau;
{
  double escape;

  /* TAU_MIN is defined in sirocco.h */

  if (tau < TAU_MIN)
    escape = 1.;
  else if (tau < 10.0)
    escape = (1. - exp (-tau)) / tau;

  else
    escape = 1. / tau;

  return (escape);
}




/**********************************************************/
/**
 * @brief
 * calculates the amount of line heating that during a resonance.
 *
 * @param [in out] PlasmaPtr  xplasma   The plasma cell where the resonance
 * occurs
 * @param [in,out] PhotPtr  pp   The photon bundle associated with the event
 * @param [in] int  nres   The number of the resonance
 * @return   Alway returns 0  f
 *
 * xplasma->heat_lines and heat_total are updated.  The weight of photon
 * is decreased by the amount of its energy that goes into heating
 *
 * @details
 * The routine calls scttering_fraction to determine the fraction of the
 * energy that is scattered, and adds to the heating and decrements
 * the energy of the photon bundle.
 *
 * ### Notes ###
 * It is called by trans_phot in sirocco
 *
 **********************************************************/

int
line_heat (xplasma, pp, nres)
     PlasmaPtr xplasma;
     PhotPtr pp;
     int nres;
{
  double x, sf;


  if (check_plasma (xplasma, "line_heat"))
  {
    Error ("line_heat: Attempting heat dummy plasma cell\n");
    return (0);
  }

  sf = scattering_fraction (lin_ptr[nres], xplasma);

  if (sane_check (sf))
  {
    Error ("line_heat:sane_check scattering fraction %g\n", sf);
  }
  x = pp->w * (1. - sf);
  xplasma->heat_lines += x;
  xplasma->heat_tot += x;

  // Reduce the weight of the photon bundle


  pp->w *= sf;

  return (0);

}
