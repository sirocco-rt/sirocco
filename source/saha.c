
/***********************************************************/
/** @file  saha.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  Routines to (re)calculate the ion densities
 *
 * The driving routine for calculating the ion desnsities is here, 
 * as well as the routines for carrying out many of the older ionization
 * schemes, but more modern matrix based approaches re elsewhere.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/** 
 * @brief   steering routine for    
 * modifying the densities of ions, levels, and
 * partition functions of ions within a cell of the wind based upon the mode, 
 * and other data contained within the WindPtr itself, such as t_r, t_e, w, 
 * based on the "mode".  
 * 
 *
 * @param [in,out] PlasmaPtr  xplasma   A single plasma cell      
 * @param [in] int  mode   A mode which describes which of the various ionization schemes to use
 * @return   A status message (returned from the individual routines that calculate the ionization. 
 * 0 generally signifies success.
 *
 * @details
 *
 * ###Notes####
 **********************************************************/

int
nebular_concentrations (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;
{
  int m;

  // LTE all the way -- uses tr or te for partition functions depending on nebular mode
  if (mode == NEBULARMODE_TR || mode == NEBULARMODE_TE)
  {
    partition_functions (xplasma, mode);
    m = concentrations (xplasma, mode);
  }
  else if (mode == NEBULARMODE_ML93)    // This is the standard LM method
  {
    partition_functions (xplasma, mode);        // calculate partition functions, using t_r with weights (dilute blackbody W)

    /* JM 1308 -- concentrations then populates xplasma with saha abundances. 
       in macro atom mode it also call macro_pops, which is done incorrectly, 
       the escape probabilities are calculated with saha ion densities each time,
       because saha() repopulates xplasma in each cycle. */

    concentrations (xplasma, NEBULARMODE_TR);   // Saha equation using t_r

    /* JM 1308 -- lucy then applies the lucy mazzali correction factors to the saha abundances. 
       in macro atom mode it also call macro_pops which is done correctly in this case, as lucy_mazzali, 
       which is called in lucy, does repopulate xplasma like saha does. This means that in this mode,
       it doesn't actually matter that concentrations does macro level populations wrong, as that is 
       corrected here. It should be sorted very soon, however. */

    m = lucy (xplasma);         // Main routine for running Lucy Mazzali
  }
  else if (mode == NEBULARMODE_MATRIX_BB || mode == NEBULARMODE_MATRIX_SPECTRALMODEL || mode == NEBULARMODE_MATRIX_ESTIMATORS)
  {
    /*
     * Use rate matrices based on pairwise bb approximations (matrix_bb) or
     * power law approximations for the SED in a cell (matrix_est, matrix_pow)
     */

    m = matrix_ion_populations (xplasma, mode);
  }
  else
  {
    Error ("nebular_concentrations: Unknown mode %d\n", mode);
    Exit (EXIT_FAILURE);
    exit (EXIT_FAILURE);        // avoids compiler warnings about return being uninitialized
  }


  return (m);
}





//#define SAHA 4.82907e15               /* 2* (2.*PI*MELEC*k)**1.5 / h**3  (Calculated in constants) */
//#define MAXITERATIONS 200
//#define FRACTIONAL_ERROR 0.03
//#define THETAMAX       1e4
//#define MIN_TEMP         100. NSH 0712 - moved into sirocco.h because this is used in severalpaces



/**********************************************************/
/** 
 * @brief      Calculate concentrations in a plasma cell using the Saha 
 * equation with various assumptions about what to use for the temperature
 *
 * @param [in,out] PlasmaPtr  xplasma   A singll plasma cell
 * @param [in] int  mode   A swich which defines what to use for the temperature
 * @return     concentrations returns 0 if it converged; -1 if it did not converge.  If it does
 * 	not converge, both density and *ne are likely to be garbage.
 *
 * 	The densities are updated in xplasma
 *
 * 	The modes are as follows:
 * 	* NEBULAR_TR -- use the radiation temperature
 * 	* NEBULaR_TE -- use the electron temperature
 * 	* NEBULAR_ML93 - From Mazzli & Lucy 1993, use the geometric mean of the electron and radiaion temperature
 * 
 *
 * @details
 *
 * This routine iterates to find the densities of individual ions using the Saha 
 * equation.
 *
 *
 * The Saha eqn. is  (e.g. Motz, p 110, eqn 4.14
 * 
 * *  N_r+1 N_e / N_r =  (2*pi*m*kT)**1.5/h**3 2 Z_r+1(T)/Z_r * exp (-X/kT)
 *    
 * 	where N_r+1 is the number density of the r+1 ionization state
 * 	and   Z_r+1 is the partition function of that state.  
 * 
 *The program works by perfoming a crude iteration on ne.  It could undoubtedly 
 *be improved, if one wanted better accuracy or higher speed.
 *
 * ### Notes ###
 * So that one doesn't get zero-divides in other programs, e.g. nebular_concentrations, 
 * a floor is set on the density of any ion or the electron density.
 *
 **********************************************************/

int
concentrations (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;
{
  int nion, niterate;
  double xne, xxne, xnew, xsaha;
  double theta, x;
  double get_ne ();
  double t, nh;
  int saha ();


  // This needs to be moved up into nebular_concentrations given that we
  // do not want these routines to do mode shifting.
  //
  if (mode == NEBULARMODE_TR)
  {
    t = xplasma->t_r;
  }
  else if (mode == NEBULARMODE_TE)
  {
    t = xplasma->t_e;
  }
  else if (mode == NEBULARMODE_ML93)
  {
    t = sqrt (xplasma->t_e * xplasma->t_r);
  }
  else
  {
    Error ("Concentrations: Unknown mode %d\n", mode);
    Exit (0);
    return (0);
  }

  nh = xplasma->rho * rho2nh;

  /* make an initial estimate of ne based on H alone,  Our guess
     assumes ion[0] is H1.  Note that x below is the fractional
     ionization of H and should vary from 0 to 1

     Note that in a pure H plasma, the left hand side of the Saha 
     equation can be written as (x*x*nh*nh)/((1-x)*nh)  and can
     then be converted into a quadratic.  That is what is done
     below.  

     Since this is an initial estimate g factors have been ignored
     in what follows.
   */

  if (t < MIN_TEMP)
    t = MIN_TEMP;

  xsaha = SAHA * pow (t, 1.5);

  theta = xsaha * exp (-ion[0].ip / (BOLTZMANN * t)) / nh;

  if (theta < THETAMAX)
  {
    x = (-theta + sqrt (theta * theta + 4 * theta)) / 2.;
    xne = xxne = x * nh;
  }
  else
    xne = xxne = nh;

  if (xne < 1.e-6)
    xne = 1.e-6;

  /* At this point we have an initial estimate of ne. */


  niterate = 0;
  while (niterate < MAXITERATIONS)
  {

    /* Assuming a value of ne calculate the relative densities of each ion for each element */
    /* JM 1308 -- Here we  actually populate the plasma structure with saha abundances. In Lucy Mazzali mode, 
       we apply a correction to these abundances. Note that in macro atom mode this call should 
       really only happen before the first ionization cycle */

    saha (xplasma, xne, t);


    /* New (SS Apr 04) call the routine that will replace the macro atoms ionization fractions and
       level populations with those computed with the Monte Carlo estimators. geo.macro_iniz_mode
       controls the use of macro_pops. */
    /* JM 1308 -- At the moment, the escape probabilities in macro_pops are calculated using saha ion 
       densities which is wrong. Fortunately, macro_pops is called again in lucy() and converges on the 
       correct value, but this should be fixed. I am not sure if it even needs to be called at all. */

    if (geo.macro_ioniz_mode == MACRO_IONIZ_MODE_ESTIMATORS)
    {
      macro_pops (xplasma, xne);
    }

    /*Set some floor so future divisions are sensible */
    for (nion = 0; nion < nions; nion++)
    {
      if (xplasma->density[nion] < DENSITY_MIN)
        xplasma->density[nion] = DENSITY_MIN;
    }

    /* Now determine the new value of ne from the ion abundances */
    xnew = get_ne (xplasma->density);

    if (xnew < DENSITY_MIN)
      xnew = DENSITY_MIN;

    if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)
      break;

    xne = (xnew + xne) / 2.;    /* Make a new estimate of xne */
    niterate++;
  }

  if (niterate == MAXITERATIONS)
  {
    Error ("concentrations: failed to converge t %.2g nh %.2g xnew %.2g\n", t, nh, xnew);
    Error ("concentrations: xxne %e theta %e\n", xxne, theta);
    return (-1);
  }

  xplasma->ne = xnew;
  return (0);
}




/**********************************************************/
/** 
 * @brief      Calculate the  densities for all of the ions in a single
 *    	cell given ne and t.
 *
 * @param [in,out] PlasmaPtr  xplasma   The cell
 * @param [in] double  ne   A trial value for ne
 * @param [in] double  t   A value of the temperature
 * @return     Always returns 0
 *
 * The calculate abundances are stored in xplasma
 *
 * @details
 * This routine simply uses the Saha equantion to calculate the densities
 * of all of the ions.  It assumes that ne is known, and so in general
 * if one calcualted ne from the concentrations, one would not get the
 * input value.  Iteration to calculate ne takes place in concentrations.
 *
 * ### Notes ###
 * All of the partition
 * functions are assumed to have been calculated before entering
 * the routine.
 * 	
 *
 **********************************************************/

int
saha (xplasma, ne, t)
     PlasmaPtr xplasma;
     double ne, t;

{
  double nh;
  int nelem;
  double *density, *partition;

  int first, last, nion;
  double xsaha;
  double sum, a, b;
  double big;

  density = xplasma->density;
  partition = xplasma->partition;

  nh = xplasma->rho * rho2nh;   //LTE
  xsaha = SAHA * pow (t, 1.5);

  for (nelem = 0; nelem < nelements; nelem++)
  {
    first = ele[nelem].firstion;        /*first and last identify the position in the array */

    last = first + ele[nelem].nions;    /*  So for H which has 2 ions, H1 and H2, first will generally
                                           be 0 and last will be 2 so the for loop below will just be done once for nion = 1 */

    if (first < 0 || first >= nions || last < 0 || last > nions)
    {
      Error ("saha: Confusion for element %d with first ion %d and last ion %d\n", nelem, first, last);
      Exit (0);
    }

/*    These lines were put in to make sim work properly, ideally there should be a switch so 
      if we are doing things the old way, we keep the old numbers. But, saha doesnt know the mode....
      sum = density[first] = 1e-250;
      big = pow (10., 250. / (last - first));
      big=big*1e6;   */

    sum = density[first] = 1.0;
    big = pow (10., 250. / (last - first));

    for (nion = first + 1; nion < last; nion++)
    {

      /* JM 1309 -- this next if statement is to ensure that saha densities are only calculated if the 
         ion in question is not a macro ion. Otherwise, this will affect the escape 
         probabilities that are calculated in macro_pops. The exception to this is prior
         to the first ionization cycle when we need to populate saha densities as a first guess */

      if ((ion[nion].macro_info == FALSE) || (geo.macro_ioniz_mode == MACRO_IONIZ_MODE_NO_ESTIMATORS))
      {
        b = xsaha * partition[nion] * exp (-ion[nion - 1].ip / (BOLTZMANN * t)) / (ne * partition[nion - 1]);
        if (b > big)
          b = big;              //limit step so there is no chance of overflow

        a = density[nion - 1] * b;

        sum += density[nion] = a;

        if (density[nion] < 0.0)
          Error ("saha: ion %i has negative density %8.4e", nion, density[nion]);

        if (sane_check (sum))
          Error ("saha:sane_check failed for density summation\n");

      }
    }

    a = nh * ele[nelem].abun / sum;

    for (nion = first; nion < last; nion++)
    {

      /* JM 1309 -- this next if statement is to ensure that saha densities are only calculated if the 
         ion in question is not a macro ion. Otherwise, this will affect the escape 
         probabilities that are calculated in macro_pops. The exception to this is prior
         to the first ionization cycle when we need to populate saha densities as a first guess */

      if ((ion[nion].macro_info == FALSE) || (geo.macro_ioniz_mode == MACRO_IONIZ_MODE_NO_ESTIMATORS))
      {

        density[nion] *= a;
        if (sane_check (density[nion]))
          Error ("saha:sane_check failed for density of ion %i\n", nion);
      }
    }
  }


  return (0);
}




#define MIN_FUDGE  1.e-10
#define MAX_FUDGE  10.


/**********************************************************/
/** 
 * @brief      Calculates ion abundances allowing for the Lucy Mazzali correction factors to the saha abundances
 * 	contained in xplasma. It thus assumes that saha() has been called before it, 
 * 	in the concentrations routine.
 * 
 *
 * @param [in,out] PlasmaPtr  xplasma   A single plasma cell
 * @return   Always retuns 0
 *
 *
 * @details
 * This routine solves for ne (and the concentrations) using the Lucy Mazzli correction factors.  
 *
 * ### Notes ###
 * The routine needs the Saha abundances before the correction factors can be found.  As a result,
 * the sequence is to call parttion functions, concentrations and then this routine.
 * 
 * Procedurally, the routine is analogous to concentrations()
 *
 * Note that  also applied the macro_pops routine, and required that the 
 * ne convergence criteria is rigorous enough that the macro atom level populations
 * converge on the correct values.
 **********************************************************/

int
lucy (xplasma)
     PlasmaPtr xplasma;
{
  int nelem, nion, niterate;
  double xne, xnew;
  double newden[NIONS];
  double t_r, nh;
  double t_e, www;

  t_r = xplasma->t_r;
  t_e = xplasma->t_e;
  www = xplasma->w;


  /* Initally assume electron density from the LTE densities */

  xne = xplasma->ne;
  if (xne < DENSITY_MIN)
  {
    Error ("nebular_concentrations: Very low ionization: ne initially %8.2e\n", xne);
    xne = DENSITY_MIN;
  }

  nh = xplasma->rho * rho2nh;   //LTE -- Not clear needed at this level

  /* Begin iteration loop to find ne */
  niterate = 0;
  while (niterate < MAXITERATIONS)
  {
    for (nelem = 0; nelem < nelements; nelem++)
    {

      /* JM1308 -- Here we apply the lucy/mazzali correction factors to the saha abundances
         which are contained in xplasma. These corrected abundances are copied to newden, which
         is transferred over to xplasma when we converge on ne */

      lucy_mazzali1 (nh, t_r, t_e, www, nelem, xplasma->ne, xplasma->density, xne, newden);
    }

    /* Re solve for the macro atom populations with the current guess for ne */
    /* JM1308 -- note that unlike lucy mazzali above, here we actually modify the xplasma
       structure for those ions which are being treated as macro ions. This means that the
       the newden array will contain wrong values for these particular macro ions, but due
       to the if loop at the end of this subroutine they are never passed to xplasma */
    if (geo.macro_ioniz_mode == MACRO_IONIZ_MODE_ESTIMATORS)
    {
      macro_pops (xplasma, xne);
    }


    for (nion = 0; nion < nions; nion++)
    {

      /* if the ion is being treated by macro_pops then use the populations just computed */
      if ((ion[nion].macro_info == TRUE) && (geo.macro_simple == FALSE) && (geo.macro_ioniz_mode == MACRO_IONIZ_MODE_ESTIMATORS))
      {
        newden[nion] = xplasma->density[nion];
      }

      /*Set some floor so future divisions are sensible */
      if (newden[nion] < DENSITY_MIN)
        newden[nion] = DENSITY_MIN;
    }
    xnew = get_ne (newden);
    if (xnew < DENSITY_MIN)
      xnew = DENSITY_MIN;

    /* Check to see whether the search for xne has converged and if so exit loop */
    if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)
      break;

    /* else start another iteration of the main loop */
    xne = (xnew + xne) / 2.;    /* Make a new estimate of xne */
    niterate++;
  }
  /* End of main iteration loop */

  if (niterate == MAXITERATIONS)
  {
    Error ("nebular_concentrations: failed to converge:nh %8.2e www %8.2e t_e %8.2e  t_r %8.2e \n", nh, www, t_e, t_r);
    return (-1);
  }

/* Finally transfer the calculated densities to the real density array */

  xplasma->ne = xnew;
  for (nion = 0; nion < nions; nion++)
  {
    /* If statement added here to suppress interference with macro populations (SS Apr 04) */
    if (ion[nion].macro_info == FALSE || geo.macro_ioniz_mode == MACRO_IONIZ_MODE_NO_ESTIMATORS || geo.macro_simple == TRUE)
    {
      xplasma->density[nion] = newden[nion];
    }
  }
  return (0);
}






/**********************************************************/
/** 
 * @brief      calculates ion
 * densities of a single element (defined by nelem) according to a modified on
 * the spot approximation, in this case equation 11 of Mazzali and Lucy which
 * attempts to consider ionizations from excieted and ground state levels.
 *
 * @param [in] double  nh   Number density of H
 * @param [in] double  t_r   The radiation temeperture
 * @param [in] double  t_e   The electorn temperature                                        
 * @param [in] double  www   dilution factor for the radiation field
 * @param [in] int  nelem   the element number for which the calculation is carried out
 * @param [in] double  ne   electron density
 * @param [in] double  density[]   densities of individual ions under pure Saha conditions
 * @param [in out] double  xne   ???
 * @param [out] double  newden[]   Array containing the modified densities
 * @return  0 if it converged to an answer, -1 otherwise.  On an
 *  	abnormal return the density array and ne are left unchanged.
 *
 * @details
 * This is a first attempt to fudge the ionization equation in the wind to account 
 * for the fact that one is not sitting in a BB radiation field.  The routine assumes one
 * has already put the concentrations of of the ions under LTE conditions for a temperature
 * t_r into density.  It uses the same dumb procedure as for "concentrations."
 *
 * ### Notes ###
 * It's probably OK but it would be worthwhile worrying whether I have introduced 
 * any problems by setting the floor to any species and ne of DENSITY_MIN
 *
 **********************************************************/

int
lucy_mazzali1 (nh, t_r, t_e, www, nelem, ne, density, xne, newden)
     double nh, t_r, t_e, www;
     int nelem;
     double ne, density[], xne, newden[];
{
  double fudge;
  double fudge2, q;
  double sum, a;
  int first, last, nion;
  double numerator, denominator;

  if (ele[nelem].z == 26)
  {
    Debug ("Working on Fe\n");
  }

  if (t_r > MIN_TEMP)
  {
    fudge = www * sqrt (t_e / t_r);
  }

  else
  {
    Error_silent ("lucy_mazzali1: t_r too low www %8.2e t_e %8.2e  t_r %8.2e \n", www, t_e, t_r);
    fudge = 0.0;
  }

  if (fudge < MIN_FUDGE || MAX_FUDGE < fudge)
  {
    Error_silent ("lucy_mazzali1: fudge>10 www %8.2e t_e %8.2e  t_r %8.2e \n", www, t_e, t_r);
  }

  /* Initialization of fudges complete */

  first = ele[nelem].firstion;  /*identify the position of the first and last ion in the array */
  last = first + ele[nelem].nions;      /*  So for H which has 2 ions, H1 and H2, first will generally
                                           be 0 and last will be 2 so the for loop below will just be done once for nion = 1 */


  /* JM 1411 -- the next two while loops check if there are ions low down
     or high up in the sequence of ions that have negligible density-
     in which case we don't bother applying the mazzali lucy 
     scheme to those ions */

  /* JM 1411 -- we only want to increment first if 
     we don't go off the end of the array, see #93 */
  while (density[first] < 1.1 * DENSITY_MIN && first < last)
  {
    newden[first] = DENSITY_MIN;

    first++;

    /* if first has been incremented all the way up to last 
       then this means all ions have negligible densities - 
       so throw an error */
    if (first == last)
      Error ("lucy_mazzali1: elem %i has all ions < %8.4e density\n", nelem, DENSITY_MIN);

  }

  /* JM 1411 -- we only want to decrement last if 
     we don't go off the end of the array, see #93 */
  while (density[last - 1] < 1.1 * DENSITY_MIN && last > first)
  {
    newden[last - 1] = DENSITY_MIN;

    if (last > first)
      last--;
  }



  sum = newden[first] = 1.;

  for (nion = first + 1; nion < last; nion++)
  {
    numerator = newden[nion - 1] * fudge * (ne) * density[nion];
    denominator = density[nion - 1] * xne;
    q = numerator / denominator;

    /* now apply the Mazzali and  Lucy fudge, i.e. zeta + w (1-zeta)
       first find fraction of recombinations going directly to
       the ground state for this temperature and ion */

    /* NSH 1207 - call external function, mode 2 uses chianti and badnell 
       data to try to take account of DR in zeta - if atomic data is not read in, 
       then the old interpolated zeta will be returned */

    fudge2 = compute_zeta (t_e, nion - 1, 2);

    /* Since nion-1 points at state i-1 (saha: n_i/n_i-1) we want ground_frac[nion-1].
       Note that we NEVER access the last ion for any element in that way
       which is just as well since you can't have recombinations INTO
       a fully ionized state -- The corresponding lines in recombination.out 
       are just dummies to make reading the data easier */

    fudge2 = fudge2 + www * (1.0 - fudge2);
    newden[nion] = fudge2 * q;
    sum += newden[nion];

    // This is the equation being calculated
    // sum+=newden[nion]=newden[nion-1]*fudge*(*ne)*density[nion]/density[nion-1]/xne;
  }

  a = nh * ele[nelem].abun / sum;

  for (nion = first; nion < last; nion++)
    newden[nion] *= a;



  return (0);
}


int fix_con_start = 0;
struct force_con
{
  int z, istate;
  float frac;
}
con_force[1000];
int nforce;


/**********************************************************/
/** 
 * @brief      sets the concentrations of individual ions
 *  assuming hard coded ionization fractions that are read from a file
 *
 * @param [in,out] PlasmaPtr  xplasma   A plasma cell
 * @param [in] int  mode   UNUSED
 * @return     Always returns 0
 *
 * @details
 * The first time the routine is called it reads the file specified by
 * geo.con_file which populates the force_con array.
 * 
 * On subsequent calls the values of one of an xplasma element are 
 * modified.
 * 
 * The format of the file is fairly obvious
 * 
 * Element.z   Ion   ion_fraction
 * 
 * 6            4       1
 * 
 * would set the CIV fraction to 1.  
 * 	
 * Note there is nothing forcing the totals to be anything sensible.
 *
 * ### Notes ###
 * o change the ions that you want to fix, simply update con_force,
 * AND MAKE SURE NFORCE IS SET TO THE NUMBER OF IONS
 *
 * ### Programming Comment ### 
 *
 * It is not obvious that this is really what we want for fixed ionization.  This routine was
 * written when we only were concerned about scattering and does not really
 * take into account emission (for which one needs n_e.  One possiblity is to
 * calculated ion abundances for H and He properly, in which case we would need
 * the mode variable.  This is issue #291.
 *
 **********************************************************/

int
fix_concentrations (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;                  // 0=saha using tr, 1=saha using te, 2= Lucy & Mazzali

{
  int nelem, nion;
  int n;
  double get_ne ();
  FILE *fopen (), *cptr;
  char line[LINELENGTH];

  double nh;

  nh = xplasma->rho * rho2nh;


  /* Define the element and ion which will be present in the wind */


  if (fix_con_start == 0)
  {
    if ((cptr = fopen (geo.fixed_con_file, "r")) == NULL)
    {
      Error ("fix_concentrations: Could not open %s to read concentrations\n", geo.fixed_con_file);
      Exit (0);
    }

    nforce = 0;
    while (fgets (line, LINELENGTH, cptr) != NULL && nforce < 1000)
    {
      if ((n = sscanf (line, "%d %d %f", &con_force[nforce].z, &con_force[nforce].istate, &con_force[nforce].frac)) != 3)
      {
        Error ("fix_concentrations: Read %d items from %s\n", n, line);
      }
      Log ("Fixing element %d istate %d to fraction %f\n", con_force[nforce].z, con_force[nforce].istate, con_force[nforce].frac);
      nforce++;
    }

    fix_con_start = 1;
    fclose (cptr);
  }

  /* Set all the ion abundances to 0 and *ne to 0 */
  for (nion = 0; nion < nions; nion++)
    xplasma->density[nion] = 0;

  /* Search for matches and if one is found set the density of that ion to be the density of that
     atom */
  for (n = 0; n < nforce; n++)
  {
    nion = 0;
    while (nion < nions && !(ion[nion].z == con_force[n].z && ion[nion].istate == con_force[n].istate))
      nion++;
    if (nion < nions)
    {                           /* There was a match */
      /* Find the element associated with the ion */
      nelem = 0;
      while (nelem < nelements && ele[nelem].z != con_force[n].z)
        nelem++;
      /* Increment the ion density and the electron density */
      xplasma->density[nion] = nh * ele[nelem].abun * con_force[n].frac;
    }
  }

  xplasma->ne = get_ne (xplasma->density);

  partition_functions (xplasma, NEBULARMODE_TR);

  return (0);
}



/**********************************************************/
/** 
 * @brief      calculates the electron density given the 
 * 	densities of each ion.
 *
 * @param [in] double  density[]   A trial set of densities for ions
 * @return   An electron density
 *
 * @details
 * The routine is used in situations where we have a trial set 
 * of densities for a cell, but do not want to populate the 
 * plasma structure with them (yet).
 *
 * ### Notes ###
 * This makes use of the fact that the densities are stored
 * 	in ion order.
 *
 **********************************************************/

double
get_ne (density)
     double density[];
{
  int n;
  double ne;
  ne = 0;
  for (n = 0; n < nions; n++)
  {
    ne += density[n] * (ion[n].istate - 1);
  }
  return (ne);
}
