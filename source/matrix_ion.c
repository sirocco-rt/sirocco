
/***********************************************************/
/** @file  matrix_ion.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  Contains routines which compute the ionization state using a matrix inversion technique
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <float.h>
#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief      A matrix solver for the ionization state in a cell
 *
 * @param [in,out] PlasmaPtr  xplasma   The plasma cell we are working on
 * @param [in] int  mode   What type of model to use for J_nu - 1=fitted model 2=blackbody
 * @return     0 if successful
 *
 * The abundances contained in xplasma are updated witht he results of
 * the calculation.
 *
 * @details
 * modes:
 *   1 - use a power law and or exponential model for J
 *   2 - use a dilute BB model for J defined by t_r and w
 * The general structure of matrix_ion_populations is as follows:
 * First we compute all of the rates that we have data for linking all ionization stages
 * We make an initial guess at the electron density
 * We attempt to solve the ionization rates matrix, and produce a new guess at the electron density
 * We then calculate new rates, re-solve and proceed until the electron density converges.
 * We solve the matrix equation A x = b for the vector x, where A is a square rate matrix linking
 * each state with each other state (often containing zeros if the ions arent linked) b is a vector
 * containing the total density of each element arranged in a certain way (see later) and
 * x is a vector containing the relative ion populations. We invert A to get x=bA^-1
 *
 * ### Notes ###
 * Uses a relative abundance scheme, in order to reduce large number issues
 *
 * Various parameters for the calculation, and in particular the t_e are passed
 * via the PlasmaPtr
 *
 **********************************************************/

int
matrix_ion_populations (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;

{
  double elem_dens[NELEMENTS];  //The density of each element
  int nn, mm, nrows;
  double rate_matrix[nions][nions];     //The rate matrix that we are going to try and solve
  double newden[NIONS];         //A temporary array to hold our intermediate solutions
  double nh, nh1, nh2, t_e;
  double xne, xxne, xxxne;      //Various stores for intermediate guesses at electron density
  double b_temp[nions];         //The b matrix
  double *b_data, *a_data;      //These arrays are allocated later and sent to the matrix solver
  double *populations;          //This array is allocated later and is retrieved from the matrix solver
  int matrix_err, niterate;     //counters for errors and the number of iterations we have tried to get a converged electron density
  double xnew;
  int xion[nions];              // This array keeps track of what ion is in each line
//OLD  int xelem[nions];             // This array keeps track of the element for each ion
  double pi_rates[nions];       //photoionization rate coefficients
  double rr_rates[nions];       //radiative recombination rate coefficients
  double inner_rates[n_inner_tot];      //This array contains the rates for each of the inner shells. Where they go to requires the electron yield array

  nh1 = nh2 = 0;

  /* Copy some quantities from the cell into local variables */

  nh = xplasma->rho * rho2nh;   // The number density of hydrogen ions - computed from density
  t_e = xplasma->t_e;           // The electron temperature in the cell - used for collisional processes


  /* We now calculate the total abundances for each element to allow us to use fractional abundances */

  /* Belt and braces, we zero our array */
  for (mm = 0; mm < ion[nions - 1].z + 1; mm++)
  {
    elem_dens[mm] = 0.0;
  }

  /* Now we populate the elemental abundance array */
  for (mm = 0; mm < nions; mm++)
  {
    elem_dens[ion[mm].z] = elem_dens[ion[mm].z] + xplasma->density[mm];
  }

  /* Dielectronic recombination, collisional ionization coefficients, three body recombination and
     charge_exchange rate coefficients depend only on electron temperature, calculate them now -
     they will not change they are all stored in global arrays */

  compute_dr_coeffs (t_e);
  compute_di_coeffs (t_e);
  compute_qrecomb_coeffs (t_e);
  compute_ch_ex_coeffs (t_e);

  /* In the following loop, over all ions in the simulation, we compute the radiative recombination rates, and photionization
     rates OUT OF each ionization stage. The PI rates are calculated either using the modelled mean intensity in a cell, or
     using the dilute blackbody approximation, depending on which mode we are in. At the same time, we copy the ion densities
     from the plasma structure into a local array. We will only overwrite the numbers in the structure if we believe the
     results are an improvement on what is there. */

  for (mm = 0; mm < nions; mm++)
  {
    newden[mm] = xplasma->density[mm] / elem_dens[ion[mm].z];   // newden is our local fractional density array
    xion[mm] = mm;              // xion is an array we use to track which ion is in which row of the matrix
    if (mm != ele[ion[mm].nelem].firstion)      // We can recombine since we are not in the first ionization stage
    {
      rr_rates[mm] = total_rrate (mm, xplasma->t_e);    // radiative recombination rates
    }
    if (ion[mm].istate != ele[ion[mm].nelem].istate_max)        // we can photoionize, since we are not in the highest ionization state
    {
      if (mode == NEBULARMODE_MATRIX_BB)
      {
        pi_rates[mm] = calc_pi_rate (mm, xplasma, 2, 1);        // PI rate for the BB model
      }
      else if (mode == NEBULARMODE_MATRIX_SPECTRALMODEL)
      {
        pi_rates[mm] = calc_pi_rate (mm, xplasma, 1, 1);        // PI rate for an explicit spectral model
      }
      else if (mode == NEBULARMODE_MATRIX_ESTIMATORS)
      {
        pi_rates[mm] = xplasma->ioniz[mm] / xplasma->density[mm];       // PI rate logged during the photon passage
      }
      else
      {
        Error ("matrix_ion_populations: Unknown mode %d\n", mode);
        Exit (0);
      }
    }

/* ksl: This whole lope is not actually used, and produces warnings in gcc11 */
//OLD    for (nn = 0; nn < nelements; nn++)
//OLD    {
//OLD      if (ion[mm].z == ele[nn].z)
//OLD      {
//OLD        xelem[mm] = nn;         /* xelem logs which element each row in the arrays refers to. This is important because we need
//OLD                                   to know what the total density will be for a group of rows all representing the same
//OLD                                   element. */
//OLD      }
//OLD    }
  }


  /* The next loop generates the inner shell ionization rates, if they are present in the atomic data and
     we wish to compute auger ionizaion rates. This only computes the rates out of each ion, we also need to
     consult the electron yield array if we are to compute the change in electron number */

  if (geo.auger_ionization == 1)
  {
    for (mm = 0; mm < n_inner_tot; mm++)
    {
      if (mode == NEBULARMODE_MATRIX_BB)
      {
        inner_rates[mm] = calc_pi_rate (mm, xplasma, 2, 2);
      }
      else if (mode == NEBULARMODE_MATRIX_SPECTRALMODEL)
      {
        inner_rates[mm] = calc_pi_rate (mm, xplasma, 1, 2);
      }
      else if (mode == NEBULARMODE_MATRIX_ESTIMATORS)   //We are using estimators - so we will need to reorder the rates from freq order to cross section order to match the electron yields. This takes time, so we only do it if we need to.
      {
        for (nn = 0; nn < n_inner_tot; nn++)
        {
          if (inner_cross[mm].nion == inner_cross_ptr[nn]->nion && inner_cross[mm].freq[0] == inner_cross_ptr[nn]->freq[0])     //Check for a match
          {
            inner_rates[mm] = xplasma->inner_ioniz[nn] / xplasma->density[inner_cross_ptr[nn]->nion];
          }
        }
      }
    }
  }

  /* This next line sets the partition function for each ion. This has always been the place here sirocco calculates the
     partition functions and sets the level densities for each ion. It needs to be done, or other parts of the code which rely
     on sensible level populations don't work properly. In the case of the dilute blackbody, the code works well, however we do
     not currently (v78) have a procedure to calucate the levels for a spectral model case. We therefore call partition
     functions with mode 4 - this is a special mode which forces the ions to be in the ground state. This is reasonable for a
     radiation dominated plasma, since any excitations into higher states will quickly be followed by radative de-excitation.
     We should really do better though... */

  if (mode == NEBULARMODE_MATRIX_BB)
  {
    partition_functions (xplasma, NEBULARMODE_ML93);    // We use t_r and the radiative weight
  }
  else if (mode == NEBULARMODE_MATRIX_SPECTRALMODEL || mode == NEBULARMODE_MATRIX_ESTIMATORS)
  {
    partition_functions (xplasma, NEBULARMODE_LTE_GROUND);      // Set to ground state
  }

  /* Next we need to obtain an initial guess for the electron density. In the past this has been done by calculating the
     hydrogen density directly from the Saha equation at the current electron temperature. In initial testing of this mode -
     this seemed to be a little unstable. At the moment, we use the last ionization state to compute n_e. For the first
     iteraion, this will just be calculated in LTE from the initial temperature structure of the wind, which will give almost
     the same result as the original procedure, or for successive calculations, it should be a better guess. I've leftin the
     original code, commented out...  */

  xne = xxne = xxxne = get_ne (xplasma->density);       //Even though the abundances are fractional, we need the real electron density


  /* xne is the current working number xxne */


  /* We are now going to iterate on the electron density - MAXITERATIONS is set in sirocco.h and is currently (78) set to 200.
     We would normally expect to converge much fater than this */

  niterate = 0;
  while (niterate < MAXITERATIONS)
  {

    /* In order to compute charge exchange reactions, we need the number density of neutral and ionized hydrogen
       the current ethod of matrix solving does not allow is to link ions multiplicatively so we need to compute
       these outside the loop. If hydrogen were not dominant - this could cause iddues since the ionization(recombination)
       of a metal or helium via this process would cause an identical recombination(ionization) of hydrogen. The loop below
       is a bit belt and braces, since one would almost always expect hydrogen to be ion=0 and 1  */

    for (nn = 0; nn < nions; nn++)
    {
      if (ion[nn].z == 1 && ion[nn].istate == 1)
      {
        nh1 = newden[nn] * elem_dens[ion[nn].z];
      }
      if (ion[nn].z == 1 && ion[nn].istate == 2)
      {
        nh2 = newden[nn] * elem_dens[ion[nn].z];
      }
    }


    populate_ion_rate_matrix (rate_matrix, pi_rates, inner_rates, rr_rates, b_temp, xne, nh1, nh2);


    /* The array is now fully populated, and we can begin the process of solving it */

    nrows = nions;              /* This is a placeholder, we may end up removing rows and columns that have no rates (or very
                                   low rates) connecting them to other rows. This may improve stability but will need to be
                                   done carefully */

    /* Here we solve the matrix equation M x = b, where x is our vector containing level populations as a fraction w.r.t the
       whole element. The actual LU decomposition- the process of obtaining a solution - is done by the routine solve_matrix() */

    a_data = (double *) calloc (sizeof (double), nrows * nrows);
    populations = (double *) calloc (nrows, sizeof (double));

    /* This b_data column matrix is the total number density for each element, placed into the row which relates to the neutral
       ion. This matches the row in the rate matrix which is just 1 1 1 1 for all stages. NB, we could have chosen any line for
       this. */

    b_data = (double *) calloc (nrows, sizeof (double));

    /* We now copy our rate matrix into the prepared matrix */
    for (mm = 0; mm < nrows; mm++)
    {
      for (nn = 0; nn < nrows; nn++)
      {
        a_data[mm * nrows + nn] = rate_matrix[mm][nn];  /* row-major */
      }
    }

    for (nn = 0; nn < nrows; nn++)
    {
      b_data[nn] = b_temp[nn];
    }

    matrix_err = solve_matrix (a_data, b_data, nrows, populations, xplasma->nplasma);

    if (matrix_err)
    {
      Error ("matrix_ion_populations: %s\n", get_matrix_error_string (matrix_err));
    }

    /* free memory */
    free (a_data);
    free (b_data);

    if (matrix_err == 4)
    {
      free (populations);
      return (-1);
    }

    /* Calculate level populations for macro-atoms */
    if (geo.macro_ioniz_mode == MACRO_IONIZ_MODE_ESTIMATORS)
    {
      int mp_err = macro_pops (xplasma, xne);
      if (mp_err != EXIT_SUCCESS)
      {
        return -1;
      }
    }

    /* We now have the populations of all the ions stored in the matrix populations. We copy this data into the newden array
       which will temperarily store all the populations. We wont copy this to the plasma structure until we are sure thatwe
       have made things better. We just loop over all ions, and copy. The complexity in the loop is to future proof us against
       the possibility that there are some ions that are not included in the matrix scheme becuse there is no route into or
       out of it. */

    for (nn = 0; nn < nions; nn++)
    {
      newden[nn] = 0.0;

      /* if the ion is being treated by macro_pops then use the populations just computed */
      if ((ion[nn].macro_info == TRUE) && (geo.macro_simple == FALSE) && (geo.macro_ioniz_mode == MACRO_IONIZ_MODE_ESTIMATORS)
          && (modes.no_macro_pops_for_ions == FALSE))
      {
        newden[nn] = xplasma->density[nn] / elem_dens[ion[nn].z];
      }

      /* if the ion is "simple" then find it's calculated ionization state in populations array */
      else
      {

        for (mm = 0; mm < nrows; mm++)  // inner loop over the elements of the population array
        {
          if (xion[mm] == nn)   // if this element contains the population of the ion is question
          {
            newden[nn] = populations[mm];       // get the population
          }
        }
      }

      if (newden[nn] < DENSITY_MIN)     // this wil also capture the case where population doesnt have a value for this ion
        newden[nn] = DENSITY_MIN;
    }
    free (populations);


/* We need to get the 'true' new electron density so we need to do a little loop here to compute it */


    xnew = 0.0;
    for (nn = 0; nn < nions; nn++)
    {
      xnew += newden[nn] * (ion[nn].istate - 1) * elem_dens[ion[nn].z];
    }



    if (xnew < DENSITY_MIN)
      xnew = DENSITY_MIN;       /* fudge to keep a floor on ne */


    if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)        /* We have converged, or have the situation where we
                                                                                   have a neutral plasma */
    {
      break;                    /* Break out of the while loop - we have finished our iterations */
    }
    xne = xxxne = (xnew + xne) / 2.;    /* New value of ne */

    niterate++;


    if (niterate == MAXITERATIONS)
    {
      Error ("matrix_ion_populations: failed to converge for cell %i t %e nh %e xnew %e\n", xplasma->nplasma, t_e, nh, xnew);


      for (nn = 0; nn < geo.nxfreq; nn++)
      {
        Log
          ("numin= %e (%e) numax= %e (%e) Model= %2d PL_log_w= %e PL_alpha= %e Exp_w= %e EXP_temp= %e\n",
           xplasma->fmin_mod[nn], geo.xfreq[nn], xplasma->fmax_mod[nn],
           geo.xfreq[nn + 1], xplasma->spec_mod_type[nn],
           xplasma->pl_log_w[nn], xplasma->pl_alpha[nn], xplasma->exp_w[nn], xplasma->exp_temp[nn]);
      }

      Error ("matrix_ion_populations: xxne %e theta %e\n", xxne);

      return (-1);              /* If we get to MAXITERATIONS, we return without copying the new populations into plasma */
    }
  }                             /* This is the end of the iteration loop */


  xplasma->ne = xnew;
  for (nn = 0; nn < nions; nn++)
  {
    /* If statement added here to suppress interference with macro populations */
    if (ion[nn].macro_info == FALSE || geo.macro_ioniz_mode == MACRO_IONIZ_MODE_NO_ESTIMATORS || geo.macro_simple == TRUE)
    {
      xplasma->density[nn] = newden[nn] * elem_dens[ion[nn].z]; //We return to absolute densities here
    }
    if ((sane_check (xplasma->density[nn])) || (xplasma->density[nn] < 0.0))
      Error ("matrix_ion_populations: ion %i has population %8.4e in cell %i\n", nn, xplasma->density[nn], xplasma->nplasma);
  }

  xplasma->ne = get_ne (xplasma->density);

  if (n_charge_exchange > 0)
  {
    xplasma->heat_ch_ex = ch_ex_heat (&wmain[xplasma->nwind], xplasma->t_e);    //Compute the charge exchange heating
    xplasma->heat_tot += xplasma->heat_ch_ex;
  }

  /*We now need to populate level densities in order to later calculate line emission (for example).
     We call partition functions to do this. At present, we do not have a method for accurately computing
     level populations for modelled mean intensities. We get around this by setting all level populations
     other than ground state to zero - this is a pretty good situation for very dilute radiation fields,
     which is what we are normally looking at, however one should really do better in the long term.
   */

  partition_functions (xplasma, NEBULARMODE_LTE_GROUND);

  return (0);
}


/**********************************************************/
/**
 * @brief      populates a rate_matrix
 *
 * @param [out] double  rate_matrix[nions][nions] - the rate matrix that will be filled by this array
 * @param [in] double  pi_rates[nions] - vector of photionization rates
 * @param [in] double  inner_rates[n_inner_tot] - vector of inner shell photoionization rates
 * @param [in] double  rr_rates[nions] - vector of radiative recobination rates
 * @param [out] double  b_temp[nions] - the vector of total elemental densities that we also fill here
 * @param [in] double  xne - current electron density
 * @param [in] double  nh1 - current neutral hydrogen density
 * @param [in] double  nh2 - current ionized hydrogen density
 * @return - zero if successful
 *
 * @details
 * populate_ion_rate_matrix populates a rate_matrix of shape nions x nions
 * with the pi_rates and rr_rates supplied at the density xne in question.
 * It also populates the b matrix - that is the total elemental abundance -
 * we use relative abundances so this is just a set of 1s and 0s.
 *
 * ### Notes ###
 * This routine includes the process of replacing the first row of the matrix with
 *     1s in order to make the problem soluble.
 *
 **********************************************************/

int
populate_ion_rate_matrix (rate_matrix, pi_rates, inner_rates, rr_rates, b_temp, xne, nh1, nh2)
     double rate_matrix[nions][nions];
     double pi_rates[nions];
     double inner_rates[n_inner_tot];
     double rr_rates[nions];
     double xne;
     double b_temp[nions];
     double nh1, nh2;

{
//  int nn, mm, zcount;
  int nn, mm;
  int n_elec, d_elec, ion_out;  //The number of electrons left in a current ion



  /* First we initialise the matrix */
  for (nn = 0; nn < nions; nn++)
  {
    for (mm = 0; mm < nions; mm++)
    {
      rate_matrix[nn][mm] = 0.0;
    }
  }


  /* The next block of loops populate the matrix. For simplicity of reading the code each process has its own loop. Some rates
     actually dont change during each iteration, but those that depend on n_e will. All are dealt with together at the moment,
     but this could be streamlined if it turns out that there is a bottleneck. */

  /* Now we populate the elements relating to PI depopulating a state */


  for (mm = 0; mm < nions; mm++)
  {
    if (ion[mm].istate != ele[ion[mm].nelem].istate_max)        // we have electrons
    {
      rate_matrix[mm][mm] -= pi_rates[mm];
    }
  }


  /* Now we populate the elements relating to PI populating a state */

  for (mm = 0; mm < nions; mm++)
  {
    for (nn = 0; nn < nions; nn++)
    {

      if (mm == nn + 1 && ion[nn].istate != ele[ion[nn].nelem].istate_max && ion[mm].z == ion[nn].z)
      {
        rate_matrix[mm][nn] += pi_rates[nn];
      }
    }
  }

  /* Now we populate the elements relating to direct ionization depopulating a state */

  for (mm = 0; mm < nions; mm++)
  {
    if (ion[mm].istate != ele[ion[mm].nelem].istate_max && ion[mm].dere_di_flag > 0)    // we have electrons and a DI rate
    {
      rate_matrix[mm][mm] -= (xne * di_coeffs[mm]);
    }
  }

  /* Now we populate the elements relating to direct ionization populating a state - this does depend on the electron density */

  for (mm = 0; mm < nions; mm++)
  {
    for (nn = 0; nn < nions; nn++)
    {
      if (mm == nn + 1 && ion[nn].istate != ele[ion[nn].nelem].istate_max && ion[mm].z == ion[nn].z && ion[nn].dere_di_flag > 0)
      {
        rate_matrix[mm][nn] += (xne * di_coeffs[nn]);
      }
    }
  }


  /* Now we populate the elements relating to radiative recomb depopulating a state */

  for (mm = 0; mm < nions; mm++)
  {
    if (mm != ele[ion[mm].nelem].firstion)      // we have space for electrons
    {
      rate_matrix[mm][mm] -= xne * (rr_rates[mm] + xne * qrecomb_coeffs[mm]);
    }
  }


  /* Now we populate the elements relating to radiative recomb populating a state */

  for (mm = 0; mm < nions; mm++)
  {
    for (nn = 0; nn < nions; nn++)
    {
      if (mm == nn - 1 && nn != ele[ion[nn].nelem].firstion && ion[mm].z == ion[nn].z)
      {
        rate_matrix[mm][nn] += xne * (rr_rates[nn] + xne * qrecomb_coeffs[nn]);
      }
    }
  }

  /* Now we populate the elements relating to dielectronic recombination depopulating a state */

  for (mm = 0; mm < nions; mm++)
  {
    if (mm != ele[ion[mm].nelem].firstion && ion[mm].drflag > 0)        // we have space for electrons
    {
      rate_matrix[mm][mm] -= (xne * dr_coeffs[mm]);
    }
  }


  /* Now we populate the elements relating to dielectronic recombination populating a state */



  for (mm = 0; mm < nions; mm++)
  {
    for (nn = 0; nn < nions; nn++)
    {
      if (mm == nn - 1 && nn != ele[ion[nn].nelem].firstion && ion[mm].z == ion[nn].z && ion[mm].drflag > 0)
      {
        rate_matrix[mm][nn] += (xne * dr_coeffs[nn]);
      }
    }
  }


  /* Now we populate the elements relating to charge exchange recombination -  */
  if (n_charge_exchange > 0)
  {
    for (mm = 0; mm < nions; mm++)      //This is a loop over ions - rates are computed for all ions.
    {
      if (ion[mm].z != 1 && ion[mm].istate != 0)        //Only compute for helium and up, and not for neutral species
      {
        rate_matrix[mm][mm] -= charge_exchange_recomb_rates[mm] * nh1;  //This is the depopulation
        rate_matrix[mm - 1][mm] += charge_exchange_recomb_rates[mm] * nh1;      //This is the population
      }
    }
  }

  /*And now charge exchange ionization - only a very small nummber of ions */

  for (mm = 0; mm < n_charge_exchange; mm++)    //We loop over the charge exchange rates - most will not be ionization
    if (ion[charge_exchange[mm].nion2].z == 1)  //A hydrogen recomb - metal ionization rate
    {
      ion_out = charge_exchange[mm].nion1;      //This is the ion that is being depopulated
      rate_matrix[ion_out][ion_out] -= charge_exchange_ioniz_rates[mm] * nh2;   //This is the depopulation
      rate_matrix[ion_out + 1][ion_out] += charge_exchange_ioniz_rates[mm] * nh2;       //This is the population
    }


  for (mm = 0; mm < n_inner_tot; mm++)  //There mare be several rates for each ion, so we loop over all the rates
  {
    if (inner_cross[mm].n_elec_yield != -1)     //we only want to treat ionization where we have info about the yield
    {
      ion_out = inner_cross[mm].nion;   //this is the ion which is being depopulated

      rate_matrix[ion_out][ion_out] -= inner_rates[mm]; //This is the depopulation
      n_elec = ion[ion_out].z - ion[ion_out].istate + 1;
      if (n_elec > 11)
        n_elec = 11;
      for (d_elec = 1; d_elec < n_elec; d_elec++)       //We do a loop over the number of remaining electrons
      {
        nn = ion_out + d_elec;  //We will be populating a state d_elec stages higher
        rate_matrix[nn][ion_out] += inner_rates[mm] * inner_elec_yield[inner_cross[mm].n_elec_yield].prob[d_elec - 1];
      }
    }
  }





  /* Now, we replace the first line for each element with 1's and 0's. This is done because we actually have more equations
     than unknowns. We replace the array elements relating to each ion stage in this element with a 1, and all the other array
     elements (which relate to other elements (in the chemicalsense) with zeros. This is equivalent to the equation
     1*n1+1*n2+1*n3 = n_total - i.e. the sum of all the partial number densities adds up to the total number density for that
     element. This loop also produces the 'b matrix'. This is the right hand side of the matrix equation, and represents the
     total number density for each element. This can be a series of 1's for each row in the matrix relating to the ground
     state of the repsective element, however one can just use the total number density for that elements and then the
     densities which the solver computes are just the actual number densities for each ion. This does mean that different rows
     are orders of magnitude different, so one can imagine numerical issues. However each element is essentially solved for
     seperately.. Something to watch */
//OLD  zcount = 0;
  for (nn = 0; nn < nions; nn++)
  {
    if (nn == ele[ion[nn].nelem].firstion)
    {
      b_temp[nn] = 1.0;         //In the relative abundance schene this equals one.

      for (mm = 0; mm < nions; mm++)
      {
        if (ion[mm].z == ion[nn].z)
        {
          rate_matrix[nn][mm] = 1.0;
        }
        else
        {
          rate_matrix[nn][mm] = 0.0;
        }
      }
    }
    else
    {
      b_temp[nn] = 0.0;
    }
  }


  return (0);
}
