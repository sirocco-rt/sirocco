
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
 * This routine, ion_abundances, is the steering routine for
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
    /* LTE using t_r - force Saha equation for all ions including macro atoms */
    int save_macro_ioniz_mode = geo.macro_ioniz_mode;
    geo.macro_ioniz_mode = MACRO_IONIZ_MODE_NO_ESTIMATORS;
    ireturn = nebular_concentrations (xplasma, NEBULARMODE_TR);
    geo.macro_ioniz_mode = save_macro_ioniz_mode;
  }
  else if (mode == IONMODE_LTE_TE)
  {
    /* LTE using t_e - force Saha equation for all ions including macro atoms */
    int save_macro_ioniz_mode = geo.macro_ioniz_mode;
    geo.macro_ioniz_mode = MACRO_IONIZ_MODE_NO_ESTIMATORS;
    ireturn = nebular_concentrations (xplasma, NEBULARMODE_TE);
    geo.macro_ioniz_mode = save_macro_ioniz_mode;
  }
  else if (mode == IONMODE_LTE_ITERATE)
  {
    /* LTE with heating/cooling balance.  The approach is:
       1) Find t_e where heating = cooling, with Saha ionization
       recalculated at each trial temperature AND MC heating scaled
       to reflect the new densities.  This makes heating and cooling
       self-consistent at every trial temperature.
       2) Apply gain damping to the temperature.
       3) Blend the LTE densities with the old densities using gain^2
       to damp the ionization-opacity feedback loop.
       Force Saha equation for all ions including macro atoms. */
    int save_macro_ioniz_mode = geo.macro_ioniz_mode;
    geo.macro_ioniz_mode = MACRO_IONIZ_MODE_NO_ESTIMATORS;

    update_old_plasma_variables (xplasma);

    double te_old = xplasma->t_e;
    double gain = xplasma->gain;

    /* Save the current ion densities for blending later */
    double density_old[nions];
    int nion;
    for (nion = 0; nion < nions; nion++)
    {
      density_old[nion] = xplasma->density[nion];
    }

    /* Find t_e where heating = cooling with coupled LTE ionization.
       calc_te_lte uses zero_emit_lte which recalculates Saha ionization
       and scales MC heating at each trial temperature. */
    double te_new = calc_te_lte (xplasma, 0.7 * te_old, 1.3 * te_old);

    /* Apply gain damping */
    xplasma->t_e = (1 - gain) * te_old + gain * te_new;

    if (xplasma->t_e > TMAX)
    {
      xplasma->t_e = TMAX;
    }
    if (xplasma->t_e < MIN_TEMP)
    {
      xplasma->t_e = MIN_TEMP;
    }

    /* Recompute ionization and heating/cooling at the gain-damped temperature */
    nebular_concentrations (xplasma, NEBULARMODE_TE);

    /* Blend the new LTE densities with the old densities.  Because the
       Saha equation is exponentially sensitive to temperature, even modest t_e
       changes produce huge density changes (e.g. He2 shifting by factors of
       1000).  We therefore use a much smaller blending fraction for the
       densities than for the temperature, so the opacities evolve slowly
       and the heating changes smoothly between cycles.  The blending preserves
       total element abundances since both old and new densities individually
       sum to the correct element density. */
    double density_gain = gain * gain;
    if (density_gain < 0.01)
      density_gain = 0.01;
    for (nion = 0; nion < nions; nion++)
    {
      xplasma->density[nion] = (1 - density_gain) * density_old[nion] + density_gain * xplasma->density[nion];
    }
    xplasma->ne = get_ne (xplasma->density);

    convergence (xplasma);

    geo.macro_ioniz_mode = save_macro_ioniz_mode;
    ireturn = 0;
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
  else if (mode == IONMODE_MATRIX_MULTISHOT)
  {
/*  This is a new development mode
    spectral_estimators does the work of getting banded W and alpha. Then oneshot gets called. */

    spectral_estimators (xplasma);
    update_old_plasma_variables (xplasma);
    int kkk, jjj;
    double xte[MAX_MULTISHOT + 1];
    double delta[MAX_MULTISHOT + 1];

    for (kkk = 0; kkk < MAX_MULTISHOT; kkk++)
    {
      ireturn = one_shot (xplasma, NEBULARMODE_MATRIX_SPECTRALMODEL);
      xte[kkk] = xplasma->t_e;
      if (kkk > 1)
      {
        delta[kkk] = (xte[kkk] - xte[kkk - 1]) / (0.5 * (xte[kkk] + xte[kkk - 1]));
        if (fabs (delta[kkk]) < DELTA_MULTISHOT)
          break;
      }
      else
      {
        delta[kkk] = 1.0;
      }
    }

    for (jjj = 0; jjj < kkk; jjj++)
    {
      Log ("XXXXX %5d %10.3e %8.3f\n", jjj, xte[jjj], delta[jjj]);
    }

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
   * The cell is converging as the electron temperature is oscillating and the change in temperature is decreasing.
   * For LTE_ITERATE mode, any oscillation triggers damping because the ionization-opacity feedback
   * can cause self-reinforcing oscillations where increasing the gain makes things worse.
   */
  if (xplasma->dt_e_old * xplasma->dt_e < 0 && (fabs (xplasma->dt_e) < fabs (xplasma->dt_e_old) || geo.ioniz_mode == IONMODE_LTE_ITERATE))
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

/* Storage for original MC-phase values used by zero_emit_lte to scale
   heating when densities change during the temperature search */
static double *lte_density_orig = NULL;
static double lte_ne_orig;
static double lte_heat_photo_orig, lte_heat_ff_orig;
static double lte_heat_comp_orig, lte_heat_ind_comp_orig;
static double lte_heat_lines_orig, lte_heat_auger_orig;
static double lte_heat_lines_macro_orig, lte_heat_photo_macro_orig;
static double lte_heat_qrecomb_macro_orig;
static double lte_heat_ch_ex_orig;


/**********************************************************/
/**
 * @brief      calculates new densities of ions in a single element of the wind
 * 	after (usually) having found the
 * 	temperature which matches heating and cooling for the previous
 * 	densities
 *
 * @param [in,out] PlasmaPtr  xplasma   The plasma cell of interest
 * @param [in] int  mode   A switch describing what ionization mode to use in determinging the
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
 * Some of the complication in this routine reflects the fact that we have
 * two sets of definitions, one for IONMODES and one for NEBULARMODES. 
 * The later governs how the routine nebular_concentrations works.
 *
 **********************************************************/

int
one_shot (PlasmaPtr xplasma, int mode)
{
  double te_old, te_new;
  double gain;


  gain = xplasma->gain;

  te_old = xplasma->t_e;

  if (modes.zeus_connect == TRUE || modes.fixed_temp == TRUE)
  {
    // Special handling for rad_hydro wher we don not want temperature uptdated.
    te_new = te_old;
    xxxplasma = xplasma;
    zero_emit (te_old);
  }
  else                          //Find a new teperature where heating and cooling match
  {
    te_new = calc_te (xplasma, 0.7 * te_old, 1.3 * te_old);     //compute the new t_e - no limits on where it can go
    xplasma->t_e = (1 - gain) * te_old + gain * te_new; /*Allow the temperature to move by a fraction gain towards
                                                           the equilibrium temperature */

    if (xplasma->t_e > TMAX)    //check to see if we have maxed out the temperature.
    {
      xplasma->t_e = TMAX;
    }
    zero_emit (xplasma->t_e);   //Get the heating and cooling rates correctly for the new temperature
  }


  if (nebular_concentrations (xplasma, mode))
  {
    Error ("one_shot: nebular_concentrations failed to converge\n");
    Error ("one_shot: j %8.2e t_e %8.2e t_r %8.2e w %8.2e nphot %i\n", xplasma->j, xplasma->t_e, xplasma->t_r, xplasma->w, xplasma->ntot);
  }
  if (xplasma->ne < 0 || VERY_BIG < xplasma->ne)
  {
    Error ("one_shot: ne = %8.2e out of range\n", xplasma->ne);
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

  /* Similarly for macro_atom_bf_heating, for both photoionization 
     and three body recombination components */

  xplasma->heat_tot -= xplasma->heat_photo_macro;
  xplasma->heat_photo -= xplasma->heat_photo_macro;
  xplasma->heat_tot -= xplasma->heat_qrecomb_macro;
  xplasma->heat_photo -= xplasma->heat_qrecomb_macro;

  xplasma->heat_photo_macro = macro_photo_heating (xplasma, xplasma->t_e);
  xplasma->heat_qrecomb_macro = macro_qrecomb_heating (xplasma, xplasma->t_e);

  xplasma->heat_tot += xplasma->heat_photo_macro;
  xplasma->heat_photo += xplasma->heat_photo_macro;
  xplasma->heat_tot += xplasma->heat_qrecomb_macro;
  xplasma->heat_photo += xplasma->heat_qrecomb_macro;


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

  /* Original method */
  xxxplasma->t_e = t;


  /* Correct heat_tot for the change in temperature. SS June 04. */
  xxxplasma->heat_tot -= xxxplasma->heat_lines_macro;
  xxxplasma->heat_lines -= xxxplasma->heat_lines_macro;

  xxxplasma->heat_lines_macro = macro_bb_heating (xxxplasma, t);

  xxxplasma->heat_tot += xxxplasma->heat_lines_macro;
  xxxplasma->heat_lines += xxxplasma->heat_lines_macro;

  /* Similarly for macro_atom_bf_heating, recompute for both photoionization 
     and three body recombination components */

  xxxplasma->heat_tot -= xxxplasma->heat_photo_macro;
  xxxplasma->heat_photo -= xxxplasma->heat_photo_macro;
  xxxplasma->heat_tot -= xxxplasma->heat_qrecomb_macro;
  xxxplasma->heat_photo -= xxxplasma->heat_qrecomb_macro;

  xxxplasma->heat_photo_macro = macro_photo_heating (xxxplasma, t);
  xxxplasma->heat_qrecomb_macro = macro_qrecomb_heating (xxxplasma, t);

  xxxplasma->heat_tot += xxxplasma->heat_photo_macro;
  xxxplasma->heat_photo += xxxplasma->heat_photo_macro;
  xxxplasma->heat_tot += xxxplasma->heat_qrecomb_macro;
  xxxplasma->heat_photo += xxxplasma->heat_qrecomb_macro;


  /* Finished macro atom corrections. Now compute the total cooling and heating - cooling */

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


/**********************************************************/
/**
 * @brief      Calculate heating - cooling for a trial temperature,
 * updating LTE ionization at each trial temperature.
 *
 * @param [in] double  t   A trial temperature
 * @return     The difference between heating and cooling at temperature t
 * with LTE ionization calculated at t.
 *
 * @details
 * This is similar to zero_emit, but recalculates LTE ionization at each
 * trial temperature using the Saha equation. This properly couples the
 * temperature and ionization for LTE calculations.
 *
 * The routine assumes geo.macro_ioniz_mode has been set to
 * MACRO_IONIZ_MODE_NO_ESTIMATORS by the caller so that Saha is applied
 * to all ions including macro atoms.
 *
 **********************************************************/

double
zero_emit_lte (double t)
{
  double difference;
  double scaled_heat_photo, scaled_heat_auger;
  double ne_ratio;
  int nion;
  double density_min = 1.e-30;

  xxxplasma->t_e = t;

  /* Update ionization to LTE at this trial temperature */
  nebular_concentrations (xxxplasma, NEBULARMODE_TE);

  /* Scale MC-phase heating to reflect the new densities from Saha.
     The MC heating was accumulated with the original densities, but
     the Saha equation has now changed them.  We scale per-ion heating
     by the density ratio and ne-dependent heating by the ne ratio. */

  ne_ratio = (lte_ne_orig > 0) ? xxxplasma->ne / lte_ne_orig : 1.0;

  /* Scale per-ion photoionization heating */
  scaled_heat_photo = 0.0;
  for (nion = 0; nion < nions; nion++)
  {
    if (lte_density_orig[nion] > density_min)
    {
      scaled_heat_photo += xxxplasma->heat_ion[nion] * xxxplasma->density[nion] / lte_density_orig[nion];
    }
    else
    {
      scaled_heat_photo += xxxplasma->heat_ion[nion];
    }
  }
  xxxplasma->heat_photo = scaled_heat_photo;

  /* Scale per-ion Auger (inner shell) heating */
  scaled_heat_auger = 0.0;
  for (nion = 0; nion < nions; nion++)
  {
    if (lte_density_orig[nion] > density_min)
    {
      scaled_heat_auger += xxxplasma->heat_inner_ion[nion] * xxxplasma->density[nion] / lte_density_orig[nion];
    }
    else
    {
      scaled_heat_auger += xxxplasma->heat_inner_ion[nion];
    }
  }
  xxxplasma->heat_auger = scaled_heat_auger;

  /* Scale ne-dependent heating */
  xxxplasma->heat_ff = lte_heat_ff_orig * ne_ratio;
  xxxplasma->heat_comp = lte_heat_comp_orig * ne_ratio;
  xxxplasma->heat_ind_comp = lte_heat_ind_comp_orig * ne_ratio;

  /* Scale non-macro line heating by ne ratio */
  double non_macro_lines = lte_heat_lines_orig - lte_heat_lines_macro_orig;
  xxxplasma->heat_lines = non_macro_lines * ne_ratio;

  /* Reconstruct heat_tot from scaled components (without macro yet) */
  xxxplasma->heat_tot = xxxplasma->heat_photo + xxxplasma->heat_auger
    + xxxplasma->heat_ff + xxxplasma->heat_comp + xxxplasma->heat_ind_comp + xxxplasma->heat_lines + lte_heat_ch_ex_orig;

  /* Now subtract the original macro contributions (which are included
     in the scaled heat_photo and heat_lines above) and replace with
     freshly computed macro heating at the new temperature and densities */
  /* note that we have to include collisional recombination (qrecomb) */
  xxxplasma->heat_photo -= lte_heat_photo_macro_orig;
  xxxplasma->heat_photo -= lte_heat_qrecomb_macro_orig;

  xxxplasma->heat_photo_macro = macro_photo_heating (xxxplasma, t);
  xxxplasma->heat_photo += xxxplasma->heat_photo_macro;

  xxxplasma->heat_qrecomb_macro = macro_qrecomb_heating (xxxplasma, t);
  xxxplasma->heat_photo += xxxplasma->heat_qrecomb_macro;

  xxxplasma->heat_lines_macro = macro_bb_heating (xxxplasma, t);

  xxxplasma->heat_lines += xxxplasma->heat_lines_macro;

  /* Reconstruct heat_tot with macro corrections */
  xxxplasma->heat_tot = xxxplasma->heat_photo + xxxplasma->heat_auger
    + xxxplasma->heat_ff + xxxplasma->heat_comp + xxxplasma->heat_ind_comp + xxxplasma->heat_lines + lte_heat_ch_ex_orig;

  cooling (xxxplasma, t);

  difference = xxxplasma->heat_tot + xxxplasma->heat_shock - xxxplasma->cool_tot;

  return (difference);
}

/**********************************************************/
/**
 * @brief     A wrapper function for zero_emit_lte used by zero_find
 *
 * @param [in] double  t   A trial temperature
 * @param [in] void *params   Not used
 * @return     The difference between heating and cooling
 *
 **********************************************************/

double
zero_emit_lte2 (double t, void *params)
{
  return (zero_emit_lte (t));
}


/**********************************************************/
/**
 * @brief  Find the electron temperature for LTE where heating and cooling
 * match, properly coupling temperature and ionization.
 *
 * @param [in] PlasmaPtr  xplasma   A plasma cell in the wind
 * @param [in] double  tmin   A bracketing minimum temperature
 * @param [in] double  tmax   A bracketing maximum temperature
 * @return     The temperature where heating and cooling match with LTE ionization
 *
 * @details
 * This is similar to calc_te, but uses zero_emit_lte which recalculates
 * LTE ionization at each trial temperature. This ensures that the temperature
 * found is self-consistent with LTE ionization.
 *
 * Unlike calc_te, this routine DOES modify ion densities as it searches
 * for the equilibrium temperature.
 *
 **********************************************************/

double
calc_te_lte (PlasmaPtr xplasma, double tmin, double tmax)
{
  double z1, z2;
  double heat1, cool1, heat2, cool2;
  int ierr = FALSE;
  int bracketed = TRUE;

  xxxplasma = xplasma;

  xxxplasma->heat_tot += xxxplasma->heat_ch_ex;

  /* Save original MC-phase heating values and ion densities.
     These are used by zero_emit_lte to scale heating when Saha
     ionization changes the densities at each trial temperature. */
  int nion;
  if (lte_density_orig == NULL)
  {
    lte_density_orig = calloc (nions, sizeof (double));
  }
  for (nion = 0; nion < nions; nion++)
  {
    lte_density_orig[nion] = xplasma->density[nion];
  }
  lte_ne_orig = xplasma->ne;
  lte_heat_photo_orig = xplasma->heat_photo;
  lte_heat_ff_orig = xplasma->heat_ff;
  lte_heat_comp_orig = xplasma->heat_comp;
  lte_heat_ind_comp_orig = xplasma->heat_ind_comp;
  lte_heat_lines_orig = xplasma->heat_lines;
  lte_heat_auger_orig = xplasma->heat_auger;
  lte_heat_lines_macro_orig = xplasma->heat_lines_macro;
  lte_heat_photo_macro_orig = xplasma->heat_photo_macro;
  lte_heat_qrecomb_macro_orig = xplasma->heat_qrecomb_macro;
  lte_heat_ch_ex_orig = xplasma->heat_ch_ex;

  /* Evaluate heating-cooling difference at bracket endpoints */
  xplasma->t_e = tmin;
  z1 = zero_emit_lte (tmin);
  heat1 = xxxplasma->heat_tot;
  cool1 = xxxplasma->cool_tot;

  xplasma->t_e = tmax;
  z2 = zero_emit_lte (tmax);
  heat2 = xxxplasma->heat_tot;
  cool2 = xxxplasma->cool_tot;

  if ((z1 * z2 < 0.0))
  {                             // Then the interval is bracketed
    xplasma->t_e = zero_find (zero_emit_lte2, tmin, tmax, 50., &ierr);
    if (ierr)
    {
      Error ("calc_te_lte: zero_find failed to find a temperature\n");
    }
  }
  else
  {
    /* Root not bracketed. This indicates a fundamental mismatch between the
       heating (from MC phase) and what LTE cooling can provide. */
    bracketed = FALSE;

    Error ("calc_te_lte: Root not bracketed for cell %d\n", xplasma->nplasma);
    Error ("calc_te_lte:   At tmin=%.0f K: heat=%.2e cool=%.2e diff=%.2e\n", tmin, heat1, cool1, z1);
    Error ("calc_te_lte:   At tmax=%.0f K: heat=%.2e cool=%.2e diff=%.2e\n", tmax, heat2, cool2, z2);

    /* Choose the temperature that minimizes |heat - cool| */
    if (fabs (z1) < fabs (z2))
    {
      xplasma->t_e = tmin;
      Error ("calc_te_lte:   Using tmin=%.0f K (smaller imbalance)\n", tmin);
    }
    else
    {
      xplasma->t_e = tmax;
      Error ("calc_te_lte:   Using tmax=%.0f K (smaller imbalance)\n", tmax);
    }
  }

  /* Ensure ionization and heating/cooling are set correctly at the final temperature.
     zero_emit_lte will call nebular_concentrations and cooling. */
  zero_emit_lte (xplasma->t_e);

  /* Store whether we found a proper solution */
  xplasma->trcheck = bracketed ? 0 : 1;

  return (xplasma->t_e);
}
