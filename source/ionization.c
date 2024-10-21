
/***********************************************************/
/** @file  ionization.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  Routines used to calculate and update ion densities
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

/**********************************************************/
/**
 * @brief
 *
 * @param [in,out] PlasmaPtr  xplasma   The plasma cell to update
 *
 * @return
 *
 * @details
 *
 **********************************************************/

void
update_old_plasma_variables (PlasmaPtr xplasma)
{
  xplasma->dt_e_old = xplasma->dt_e;
  xplasma->dt_e = xplasma->t_e - xplasma->t_e_old;
  xplasma->t_e_old = xplasma->t_e;
  xplasma->lum_tot_old = xplasma->lum_tot;
  xplasma->heat_tot_old = xplasma->heat_tot;
}


/**********************************************************/
/**
 * @brief      ionization routines for the wind one cell at a time
 *
 * @param [in,out] PlasmaPtr  xplasma   The cell in which the ioniztion is to be calculated
 * @param [in] int  mode   A parameter describing how to calculate the ionization
 * @return     The routine returns status messages dereived from the individual routines
 * used to calculate the abundances.
 *
 * @details
 * The intent is that the routine ion_abundances is the steering routine for
 * all calculations of the abundances
 *
 * ### Notes ###
 *
 * EP 15 Jun 2020: It appears to me that this function can only ever return 0
 * as the functions called in here only return 0.
 *
 **********************************************************/

int
ion_abundances (PlasmaPtr xplasma, int mode)
{
  int ireturn;

  if (mode == IONMODE_ML93_FIXTE)
  {
    /* on-the-spot approximation using existing t_e.   This routine does not attempt
       to match heating and cooling in the wind element! */
    if ((ireturn = nebular_concentrations (xplasma, NEBULARMODE_ML93)))
    {
      Error ("ionization_abundances: nebular_concentrations failed to converge\n");
      Error ("ionization_abundances: j %8.2e t_e %8.2e t_r %8.2e w %8.2e\n", xplasma->j, xplasma->t_e, xplasma->w);
    }
  }
  else if (mode == IONMODE_LTE_TR)
  {
    /* LTE using t_r */
    ireturn = nebular_concentrations (xplasma, NEBULARMODE_TR);
  }
  else if (mode == IONMODE_LTE_TE)
  {
    /* LTE using t_e */
    ireturn = nebular_concentrations (xplasma, NEBULARMODE_TE);
  }
  else if (mode == IONMODE_FIXED)
  {                             //  Hardwired concentrations

    ireturn = fix_concentrations (xplasma, 0);
  }
  else if (mode == IONMODE_ML93)
  {
    /* On the spot, setting te to 0.9 t_r before calculating densities */

    xplasma->dt_e_old = xplasma->dt_e;
    xplasma->dt_e = xplasma->t_e - xplasma->t_e_old;
    xplasma->t_e_old = xplasma->t_e;
    xplasma->lum_tot_old = xplasma->lum_tot;
    xplasma->heat_tot_old = xplasma->heat_tot;
    ireturn = 0;
    xplasma->t_e = 0.9 * xplasma->t_r;
    if ((ireturn = nebular_concentrations (xplasma, NEBULARMODE_ML93)))
    {
      Error ("ionization_abundances: nebular_concentrations failed to converge\n");
      Error ("ionization_abundances: j %8.2e t_e %8.2e t_r %8.2e w %8.2e\n", xplasma->j, xplasma->t_e, xplasma->w);
    }
    convergence (xplasma);

  }
  else if (mode == IONMODE_MATRIX_BB)
  {

    update_old_plasma_variables (xplasma);
    ireturn = one_shot (xplasma, mode);

    convergence (xplasma);
  }
  else if (mode == IONMODE_MATRIX_SPECTRALMODEL || mode == IONMODE_MATRIX_ESTIMATORS)
  {
/*  spectral_estimators does the work of getting banded W and alpha. Then oneshot gets called. */

    spectral_estimators (xplasma);
    update_old_plasma_variables (xplasma);
    ireturn = one_shot (xplasma, mode);

    convergence (xplasma);
  }
  else
  {
    Error ("ion_abundances: Could not calculate abundances for mode %d\n", mode);
    Exit (EXIT_FAILURE);
    exit (EXIT_FAILURE);        // avoids compiler warnings about return being uninitialized
  }

  return (ireturn);

}



/**********************************************************/
/**
 * @brief      checks to see whether a single cell
 * 	is or is not converging
 *
 * @param [in out] PlasmaPtr  xplasma   The cell of interest
 * @return    A number which describes the degree to which
 * the cell is converging
 *
 * @details
 * The routine attempts to determine whether a cell is
 * on the track towards a final solution by checking whether
 * the electron and radiation temperatures are getting
 * smaller with each cycle and whether the difference between
 * heating and cooling is decreasing.
 *
 * The routine returns a number between 0 and 3,
 * depending on the number of convergence checks that
 * are passed.  If all convergence tests are passed
 * then the number returned will be 0.
 *
 * The routine also adjusts the gain which controls how
 * much the electron temperature can change in a cycle.
 *
 * ### Notes ###
 *
 **********************************************************/

int
convergence (PlasmaPtr xplasma)
{
  int trcheck, techeck, hccheck, whole_check;
  double min_gain, gain_damp, max_gain, gain_amp, cyc_frac;
  double epsilon;

  // TODO: are these values optimal?
  min_gain = 0.1;
  gain_damp = 0.7;
  epsilon = 0.05;

  trcheck = techeck = hccheck = CONVERGENCE_CHECK_PASS;
  xplasma->trcheck = xplasma->techeck = xplasma->hccheck = CONVERGENCE_CHECK_PASS;      // NSH 70g - zero the global variables

  /*
   * Check the convergence of the radiation temperature
   */

  xplasma->converge_t_r =       // Radiation temperature check
    fabs (xplasma->t_r_old - xplasma->t_r) / (xplasma->t_r_old + xplasma->t_r);
  if (xplasma->converge_t_r > epsilon)
    xplasma->trcheck = trcheck = CONVERGENCE_CHECK_FAIL;

  /*
   * Check the convergence for electron temperature and heat + cooling rates
   * CHANGES:
   * --------
   * - 110919 nsh modified line below to include the adiabatic cooling in the check that heating equals cooling
   * - 111004 nsh further modification to include DR and compton cooling, now moved out of lum_tot
   * - 130722 added a fabs to the bottom, since it is now conceivable that this could be negative if
   *   cool_adiabatic is large and negative - and hence heating
   * - NSH 130711 - also changed to have fabs on top and bottom, since heating can now be negative!)
   * - NSH 130725 - moved the hc check to be within the if statement about overtemp - we cannot expect hc to
   *   converge if we are hitting the maximum temperature
   */

  if (xplasma->t_e < TMAX)
  {
    xplasma->converge_t_e =     // Electron temperature check
      fabs (xplasma->t_e_old - xplasma->t_e) / (xplasma->t_e_old + xplasma->t_e);
    if (xplasma->converge_t_e > epsilon)
      xplasma->techeck = techeck = CONVERGENCE_CHECK_FAIL;

    xplasma->converge_hc =      // Heating and cooling rates check
      fabs (xplasma->heat_tot + xplasma->heat_shock - xplasma->cool_tot) / fabs (xplasma->heat_tot + xplasma->heat_shock +
                                                                                 xplasma->cool_tot);
    if (xplasma->converge_hc > epsilon)
      xplasma->hccheck = hccheck = CONVERGENCE_CHECK_FAIL;
  }
  else                          // If the cell has reached the maximum temperature we mark it as over-limit
  {
    xplasma->techeck = techeck = xplasma->hccheck = hccheck = CONVERGENCE_CHECK_OVER_TEMP;
  }

  /*
   * whole_check is the sum of the temperature checks and the heating check - the higher this is, the more convergence checks have failed.
   */

  xplasma->converge_whole = whole_check = trcheck + techeck + hccheck;

  /*
   * Now we check to see if a cell is converging:
   * Converging is a situation where the change in electron temperature is dropping with time and the cell is
   * oscillating around a temperature. If that is the case, we drop the amount by which the temperature can change in
   * this cycle. Else if the cell is not converging, we increase the amount by which the temperature can change in this
   * cycle.
   */

  /*
   * The cell is converging as the electron temperature is oscillating and the change in temperature is decreasing
   */
  if (xplasma->dt_e_old * xplasma->dt_e < 0 && fabs (xplasma->dt_e) < fabs (xplasma->dt_e_old))
  {
    xplasma->converging = CELL_CONVERGING;

    // TODO: is this optimal for converging cells? See discussion on Bug #631
    xplasma->gain *= gain_damp;
    if (xplasma->gain < min_gain)
      xplasma->gain = min_gain;
  }
  /*
   * The cell is not converging, which means either that the temperature is consistently moving in one direction or
   * that the oscillations of the temperature have increased in the past two cycles.
   */
  else
  {
    /*
     * EP: allow the gain to increase more for the first cyc_frac * cycles to
     * allow the plasma to change temperature more rapidly -- right now this
     * is controlled by some magic numbers and should probably be fine tuned
     * to find the best numbers
     */

    xplasma->converging = CELL_NOT_CONVERGING;

    cyc_frac = 0.5;

    if (geo.wcycle <= floor (cyc_frac * geo.wcycles))
    {
      gain_amp = 1.5;
      max_gain = 0.999;
    }
    else
    {
      gain_amp = 1.1;
      max_gain = 0.8;
    }

    xplasma->gain *= gain_amp;
    if (xplasma->gain > max_gain)
      xplasma->gain = max_gain;
  }

  return (whole_check);
}



/**********************************************************/
/**
 * @brief      The routine summarizes the how well the wind is converging
 * to a solution as a whole
 *
 * @return     Always returns 0
 *
 * @details
 * Basically the routine just looks at the numbers of cells which haave passed/failed
 * the various convergence tests and writes this to the log file
 *
 * ### Notes ###
 *
 * All of the checks are done in the routine convergence, which also uses
 * the checks to set the gain for the the next cycle.  When a value for a
 * check is 0, then a particular convergence check has been passed.
 * Non-zero values represent conditions which resulted in the check being
 * failed, and may be different for different checks.
 *
 **********************************************************/

int
check_convergence (void)
{
  int n;
  int nconverge, nconverging, ntot;
  int nte, ntr, nhc;
  int nmax;
  double xconverge, xconverging;

  nconverge = nconverging = ntot = 0;
  ntr = nte = nhc = nmax = 0;

  for (n = 0; n < NPLASMA; n++)
  {
    if (wmain[plasmamain[n].nwind].inwind == W_ALL_INWIND || modes.partial_cells == PC_INCLUDE)
    {
      ntot++;
      if (plasmamain[n].converge_whole == CONVERGENCE_CHECK_PASS)
        nconverge++;
      if (plasmamain[n].trcheck == CONVERGENCE_CHECK_PASS)
        ntr++;
      if (plasmamain[n].techeck == CONVERGENCE_CHECK_PASS)
        nte++;
      if (plasmamain[n].hccheck == CONVERGENCE_CHECK_PASS)
        nhc++;
      if (plasmamain[n].techeck == CONVERGENCE_CHECK_OVER_TEMP)
        nmax++;
      if (plasmamain[n].converging == CELL_CONVERGING)
        nconverging++;
    }
  }

  xconverge = ((double) nconverge) / ntot;
  xconverging = ((double) nconverging) / ntot;
  geo.fraction_converged = xconverge;

  Log ("!!Check_convergence: %4d (%.3f) converged and %4d (%.3f) converging of %d cells actually in the wind\n",
       nconverge, xconverge, nconverging, xconverging, ntot);
  Log ("!!Check_convergence: t_r %4d t_e(real) %4d t_e(maxed) %4d hc(real) %4d\n", ntr, nte, nmax, nhc);
  Log_flush ();

  return (0);
}



/* An externall pointer reference used by zero_emit.  */
PlasmaPtr xxxplasma;


/**********************************************************/
/**
 * @brief      calculates new densities of ions in a single element of the wind
 * 	after (usually) having found the
 * 	temperature which matches heating and cooling for the previous
 * 	densities
 *
 * @param [in,out] PlasmaPtr  xplasma   The plasma cell of interest
 * @param [in] int  mode   A switch describing what approximation to use in determinging the
 * densities
 * @return     Always returns 0
 *
 * @details
 * This routine attempts to match heating and cooling in the wind element!
 * To do this it calls calc_te.  Based on the returned value of te, the
 * routine then calculates densities for various ions in the cell.  The densities
 * in xplasma are updated.
 *
 * ### Notes ###
 *
 *
 * Special exceptions are made for Zeus; it is not clear why this is necessary
 *
 **********************************************************/

int
one_shot (PlasmaPtr xplasma, int mode)
{
  double te_old, te_new;
  double gain;



  gain = xplasma->gain;

  te_old = xplasma->t_e;

  if (modes.zeus_connect == 1 || modes.fixed_temp == 1)
  {
    te_new = te_old;            //We don't want to change the temperature
    xxxplasma = xplasma;
    zero_emit (te_old);         //But we do still want to compute all heating and cooling rates
  }
  else                          //Do things to normal way - look for a new temperature
  {
    te_new = calc_te (xplasma, 0.7 * te_old, 1.3 * te_old);     //compute the new t_e - no limits on where it can go
    xplasma->t_e = (1 - gain) * te_old + gain * te_new; /*Allow the temperature to move by a fraction gain towards
                                                           the equilibrium temperature */

    /* NSH 130722 - NOTE - at this stage, the cooling terms are still those computed from
     * the 'ideal' t_e, not the new t_e - this may be worth investigatiing. */
    if (xplasma->t_e > TMAX)    //check to see if we have maxed out the temperature.
    {
      xplasma->t_e = TMAX;
    }
    zero_emit (xplasma->t_e);   //Get the heating and cooling rates correctly for the new temperature
  }



/* Modes in the driving routines are not identical to those in nebular concentrations.
The next lines are an attempt to mediate this problem.  It might be better internally
at least to define a flag for using one shot, and have the modes take on the
meaning in nebular concentrations.
*/

  if (mode == IONMODE_ML93)
    mode = NEBULARMODE_ML93;    // This is weird, why not continue
  else if (mode <= 1 || mode == 5 || mode > 10)
  {
    /* There is no mode 5 at present  - SIM + two new modes in Feb 2012  + mode 5 now removed */

    Error ("one_shot: Sorry, Charlie, don't know how to process mode %d\n", mode);
    Exit (0);
  }

  if (xplasma->t_r > 10.)
  {                             /* Then modify to an on the spot approx */
    if (nebular_concentrations (xplasma, mode))
    {
      Error ("ionization_on_the_spot: nebular_concentrations failed to converge\n");
      Error ("ionization_on_the_spot: j %8.2e t_e %8.2e t_r %8.2e w %8.2e nphot %i\n", xplasma->j, xplasma->t_e, xplasma->t_r, xplasma->w,
             xplasma->ntot);
    }
    if (xplasma->ne < 0 || VERY_BIG < xplasma->ne)
    {
      Error ("ionization_on_the_spot: ne = %8.2e out of range\n", xplasma->ne);
    }
  }
  else
  {
    Error ("ionization_on_the_spot: t_r exceptionally small %g\n", xplasma->t_r);
    Exit (0);
  }


  return (0);
}





/**********************************************************/
/**
 * @brief  find the electron termperature for a cell where the heating and cooling 
 * would match, assuming that is between tmin and tmax.
 *
 * @param [in] PlasmaPtr  xplasma   A plasma cell in the wind
 * @param [in] double  tmin   A bracketing minimum temperature
 * @param [in] double  tmax   A bracketing mxximum temperature
 * @return     The temperature where heating and cooling match
 *
 * @details
 * The routine iterates to find the temperature in a cell, where heating and cooling are matched.
 *
 * calc_te does not modify any abundances.  It simply takes the current value of the heating in the
 * cell and adjusts the electron temperature to a value that matches the pre-calculated heating.
 *
 * ### Notes ###
 * Ion densities are NOT updated in this process.
 *
 * xxxplasma is just a way to tranmit information to zero_emit
 *
 **********************************************************/

double
calc_te (PlasmaPtr xplasma, double tmin, double tmax)
{
  double z1, z2;
  int ierr = FALSE;


  /* we assign a plasma pointer here to a fixed structure because
   * we need to call zbrent and we cannot pass the xplasma ptr directly
   */

  xxxplasma = xplasma;

  xxxplasma->heat_tot += xxxplasma->heat_ch_ex;

  xplasma->t_e = tmin;
  z1 = zero_emit (tmin);
  xplasma->t_e = tmax;
  z2 = zero_emit (tmax);

  /* The way this works is that if we have a situation where the cooling
   * at tmax and tmin brackets the heating, then we use zbrent to improve
   * the estimated temperature, but if not we chose the best direction
   */

  if ((z1 * z2 < 0.0))
  {                             // Then the interval is bracketed
    xplasma->t_e = zero_find (zero_emit2, tmin, tmax, 50., &ierr);
    if (ierr)
    {
      Error ("calc_te: zero_find failed to find a temperature\n");
    }

  }
  else if (fabs (z1) < fabs (z2))
  {
    xplasma->t_e = tmin;
  }
  else
  {
    xplasma->t_e = tmax;
  }
  /* With the new temperature in place for the cell, get the correct value of heat_tot.
     SS June  04 */

  /* At this point we know the temperature that balances heating and cooling
   * within the constraints set by tmin and tmax.
   */


  /* Update heat_tot and heat_lines for macro_bb_heating at the new temperature. 
   * We subtract the current value and then compute at the new temperature and
   * add this back */

  xplasma->heat_tot -= xplasma->heat_lines_macro;
  xplasma->heat_lines -= xplasma->heat_lines_macro;

  xplasma->heat_lines_macro = macro_bb_heating (xplasma, xplasma->t_e);

  xplasma->heat_tot += xplasma->heat_lines_macro;
  xplasma->heat_lines += xplasma->heat_lines_macro;

  /* Similaryly for macro_atom_bf_heating */

  xplasma->heat_tot -= xplasma->heat_photo_macro;
  xplasma->heat_photo -= xplasma->heat_photo_macro;

  xplasma->heat_photo_macro = macro_bf_heating (xplasma, xplasma->t_e);

  xplasma->heat_tot += xplasma->heat_photo_macro;
  xplasma->heat_photo += xplasma->heat_photo_macro;


  return (xplasma->t_e);

}




/**********************************************************/
/**
 * @brief      Compute the cooling for a cell given a temperature t, and compare it
 * to the heating seen in the cell in the previous ionization cycle
 *
 * @param [in] double  t   A trial temperature
 * @return     The difference between the recorded heating and the cooling calculated
 * at a specifc temperature.  
 *
 * @details
 * This routine is used in the process of estimating a new temperature for a cell
 * given the
 * heating of the cell in the previous ionization cycle.  When 0 the heating
 * and cooling are matched.
 *
 * For the non-macro atom case heating is calculated as photons go through
 * the wind and is not altered by a change in temperature, however for macro
 * atoms a change in temperature affects the amount of heating as well as the 
 * cooling
 *
 *
 * ### Notes ###
 *
 * This routine is poorly named; what it does is simply calculate the
 * sum the heating that occurred during the current cycle, and calculate
 * the cooling for the input temperature given the current abundances
 * it returns the difference between these two numbers.  It does update
 * t_e in the plasma cell.
 *
 * The abundances of ions in the cell are not modified.  Results are stored
 * in the cell of interest.  This routine is used in connection with a zero
 * finding routine
 *
 * T1he equation now includes a term for non-radiative heating (heat_shock)
 * that was used for FU Ori
 *
 **********************************************************/

double
zero_emit (double t)
{
  double difference;

  /*Original method */
  xxxplasma->t_e = t;


  /* Correct heat_tot for the change in temperature. SS June 04. */
  xxxplasma->heat_tot -= xxxplasma->heat_lines_macro;
  xxxplasma->heat_lines -= xxxplasma->heat_lines_macro;

  xxxplasma->heat_lines_macro = macro_bb_heating (xxxplasma, t);

  xxxplasma->heat_tot += xxxplasma->heat_lines_macro;
  xxxplasma->heat_lines += xxxplasma->heat_lines_macro;

  xxxplasma->heat_tot -= xxxplasma->heat_photo_macro;
  xxxplasma->heat_photo -= xxxplasma->heat_photo_macro;

  xxxplasma->heat_photo_macro = macro_bf_heating (xxxplasma, t);

  xxxplasma->heat_tot += xxxplasma->heat_photo_macro;
  xxxplasma->heat_photo += xxxplasma->heat_photo_macro;

  /* Finished macro atom corrections */


  cooling (xxxplasma, t);

  difference = xxxplasma->heat_tot + xxxplasma->heat_shock - xxxplasma->cool_tot;



  return (difference);
}

/**********************************************************/
/**
 * @brief     A wrapper function used by zero_find as part the 
 * calculation of a new temperature for a plasma cell
 *
 * @param [in] double  t   A trial temperature
 * @param [in] void *params
 * @return     The difference between the recorded heating and the cooling calculated
 * at a specifc temperature.  
 *
 * @details
 *
 * See zero_emit
 *
 * ### Notes ###
 *
 **********************************************************/


double
zero_emit2 (double t, void *params)
{
  return (zero_emit (t));
}
