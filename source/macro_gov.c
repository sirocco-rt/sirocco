
/***********************************************************/
/** @file  macro_gov.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief
 * Contains the main functions to govern macro atoms and calculate densities.
 *
 * @details
 * This is the file which contains the functions which govern macro-atoms and
 * obtain their level populations and ion densities. The actual functions which
 * do the jumps inside an activated macro-atom are in matom.c. This is done to
 * prevent overly long files.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief macro_gov is a routine that governs the excitation and de-excitation of macro-atoms, kpkts and simple-atoms.

 * @param [in,out]  PhotPtr  p   the packet at the point of activation
 * @param [in,out]  int *  nres   the process which activates the Macro Atom
 * @param [in]      int  matom_or_kpkt   initially excite a matom (MATOM) or create a kpkt (KPKT)
 * @param [in,out]     int *  which_out   set to 1 if return is via macro atom and 2 if via kpkt
 * @return 0        Will return an r-packet after (possibly) several calls to matom and kpkt
 *
 * @details macro_gov sits at a higher level in the code than either matom or kpkt and governs the passage
 * of packets between these routines. At the moment, since matom and kpkt
 * call each other it can all get rather confusing and loads of nested subroutine calls ensue.
 * macro_gov removes this by calling matom and kpkt which on return tell  to either return
 * an r-packet to resonate or make another call to either kpkt or matom as appropriate.
 *
 * ### Notes ###
 * Stuart wrote this to be as general as possible so that if we want to improve the treatment of simple
 * ions it should not need to be re-written.
 *
 * During the spectrum calculation the emission of r-packets within the spectral region of interest
 * is done using emissivities obtained during the ionization cycles. Therefore whenever an r-packet
 * is converted into a k-packet or an excited macro atom that ends the need to follow the packet any
 * further. To deal with this, this routine sets the weights of such packets to zero. Upon return to
 * trans_phot these packets will then be thrown away.
 *
 **********************************************************/

int
macro_gov (p, nres, matom_or_kpkt, which_out)
     PhotPtr p;
     int *nres;
     int matom_or_kpkt;
     int *which_out;
{
  int escape;                   //this tells us when the r-packet is escaping
  int n_jump = 0;
  int n_jump_tot = 0;
  int n_loop = 0;
  int new_uplvl, uplvl;
  PlasmaPtr xplasma;
  MacroPtr mplasma;
  WindPtr one;

  one = &wmain[p->grid];
  xplasma = &plasmamain[one->nplasma];
  mplasma = &macromain[one->nplasma];

  /* before we do anything else we look to see if we are exciting
     simple/fake two-level ions */
  if (matom_or_kpkt == MATOM)
  {
    escape = FALSE;

    /*  it's a bb transition  without the full macro atom treatment. */
    if (*nres > (-1) && *nres < NLINES && (geo.macro_simple == TRUE || lin_ptr[*nres]->macro_info == FALSE))
    {
      fake_matom_bb (p, nres, &escape);

      /* at this point it has either generated an r-packet (escape == TRUE)
         or a k-packet (escape == FALSE) */
      if (escape == TRUE)
      {
        return (0);
      }
      else
      {
        matom_or_kpkt = KPKT;
      }
    }

    /* if it's bf continuum without the full macro atom treatment.

       In the pre-2018 approach, we process the photon in a way that makes it return a bf photon of the same type
       as caused the excitation.  In the old approach, escape will be set to 1, and we will escape.
       In the new "simple emissivity" approach, we should never satisfy the do loop, and so an error
       is thrown and we exit. */

    else if (*nres > NLINES && (phot_top[*nres - NLINES - 1].macro_info == FALSE || geo.macro_simple == TRUE))
    {
#if BF_SIMPLE_EMISSIVITY_APPROACH
      Error ("Macro_gov: Error - trying to access fake_matom_bf in alternate bf treatment.\n");
      Exit (0);
#endif
      fake_matom_bf (p, nres, &escape);

      /* at this point it has either generated an r-packet (escape == TRUE)
         or a k-packet (escape == FALSE) */
      if (escape == TRUE)
      {
        return (0);
      }
      else
      {
        matom_or_kpkt = KPKT;
      }
    }
  }

  /* if we got here then we know we are exciting a macro-atom or creating a kpkt */
  if (geo.matom_radiation == 1)
  {
    /* During the spectrum cycles we want to throw these photons away. */
    p->w = 0.0;
    escape = TRUE;
    return (0);
  }
  else
  {
    /* start with it not being ready to escape as an r-packet */
    escape = FALSE;
  }

  if (mplasma->matom_transition_mode == MATOM_MATRIX)
  {
    if (matom_or_kpkt == MATOM)
    {
      uplvl = 0;
      if (*nres < NLINES)
      {
        uplvl = lin_ptr[*nres]->nconfigu;
      }
      else if (*nres > NLINES)
      {
        uplvl = phot_top[*nres - NLINES - 1].uplev;
      }
      else
      {
        Error ("matom: upper level not identified. nres = %d in photon %d of cycle %d/%d in thread %d\n",
               *nres, p->np, geo.wcycle, geo.pcycle, rank_global);
        escape = TRUE;
        p->istat = P_ERROR_MATOM;
        return (-1);
      }
    }
    else
    {
      uplvl = nlevels_macro;
    }

    new_uplvl = matom_deactivation_from_matrix (xplasma, uplvl);

    if (xconfig[new_uplvl].nauger > 0)
    {
      Log ("AUGER: %d %d Jumped to Auger level %d from %d %d old level %d\n",
           xconfig[new_uplvl].z, xconfig[new_uplvl].istate, new_uplvl, xconfig[uplvl].z, xconfig[uplvl].istate, uplvl);
    }

    if (new_uplvl == nlevels_macro)
    {
      /* XMACRO improve this so that kpkt only deals with k->r in certain modes */
      kpkt (p, nres, &escape, KPKT_MODE_CONT_PLUS_ADIABATIC);

      *which_out = KPKT;
    }
    /* XMACRO -- what do we do about frequency boundaries here? */
    /* XMACRO -- change wmain to be one? */
    else
    {
      emit_matom (wmain, p, nres, new_uplvl, 0, VERY_BIG);
      *which_out = MATOM;
    }

    if (p->origin < 10)
      p->origin += 10;

    escape = TRUE;
    return (0);

  }

  /* using the old MATOM_MC_JUMPS scheme */
  else if (mplasma->matom_transition_mode == MATOM_MC_JUMPS)
  {
    /* Beginning of the main loop for processing a macro-atom */
    while (escape == FALSE)
    {
      if (matom_or_kpkt == MATOM)       //excite a macro atom
      {

        /* if it's a bb transition of a full macro atom  */
        if (*nres > (-1) && *nres < NLINES && geo.macro_simple == FALSE && lin_ptr[*nres]->macro_info == TRUE)
        {
          n_jump = matom (p, nres, &escape);

          if (escape == TRUE)
          {
            /* It escapes as a r-packet that was created by de-activation of a macro atom.
             */
            *which_out = MATOM;

            /* Update the the photon origin to indicate the packet has been processed
               by a macro atom */
            if (p->origin < 10)
              p->origin += 10;
            return (0);
          }
        }

        /* if it's bf transition of a full macro atom. */
        else if (*nres > NLINES && phot_top[*nres - NLINES - 1].macro_info == TRUE && geo.macro_simple == FALSE)
        {
          n_jump = matom (p, nres, &escape);

          if (escape == TRUE)
          {
            /* It  escapes as a r-packet that was created by de-activation of a macro atom.
             */
            *which_out = MATOM;
            /* Update the the photon origin to indicate the packet has been processed
               by a macro atom */
            if (p->origin < 10)
              p->origin += 10;

            //If reverb is on, and this is the last ionisation cycle, then track the photon path
            if (geo.reverb == REV_MATOM && geo.ioniz_or_extract == CYCLE_IONIZ && geo.fraction_converged > geo.reverb_fraction_converged)
            {
              line_paths_add_phot (&(wmain[p->grid]), p, nres);
            }

            return (0);
          }
        }

        /* If it did not escape then it must have had a
           de-activation by collision processes, and so we label it a kpkt.
         */

        matom_or_kpkt = KPKT;
      }


      /* This the end of the section of the loop that deals with matom excitations. next domes the
         section of the loop that deals with kpts */
      else if (matom_or_kpkt == KPKT)
      {
        kpkt (p, nres, &escape, KPKT_MODE_ALL); // 1 implies include the possibility of deactivation due to non-thermal processes

        /* if it did not escape then the k-packet must have been
           destroyed by collisionally exciting a macro atom so...
         */
        matom_or_kpkt = MATOM;
      }
      else
      {
        Error ("macro_gov: Unknown choice for next action. Abort.\n");
        Exit (0);
      }


      /*XXXX test */
      if (n_jump > -1)
      {
        //XXXX This is a test that fails many times for the agn_macro model.  Just set to n_mump it is the same as in matom
        n_jump_tot += n_jump;
        if (n_jump > MAXJUMPS)
        {
          Error ("macro_gov: Exceed MAXJUMPS (last %d tot %d) in n_loops %d for phot %d in cell %d\n", n_jump, n_jump_tot, n_loop, p->np,
                 p->grid);
          escape = TRUE;
          p->istat = P_ERROR_MATOM;
        }
      }
      n_loop++;
    }

    *which_out = KPKT;
  }

  /* End of main matom processing loop that began with while (escape==FALSE)

     If it gets here, the escape is as a KPKT
   */

  /* Update the the photon origin to indicate the packet has been processed
     by a macro atom */

  if (p->origin < 10)
    p->origin += 10;

  return (0);
}

/**********************************************************/
/**
 * @brief      uses the Monte Carlo estimators to compute a set
 *             of level populations for levels of macro atoms.
 *
 * @param [in out] PlasmaPtr  xplasma   Plasma pointer of cell in question
 * @param [in out] double  xne   -> current value for electron density in this shell
 * @return   Should compute the fractional level populations for
 *           macro atoms and store them in "levden" array. The ion fractions
 *           are also computed and stored in w[n].density[nion]
 *
 * @details This routine uses a matrix inversion method to get the level populations.
 * For now the matrix solver used is the LU decomposition method provided
 * by the Gnu Scientific Library (GSL, which is free). This requires the
 * include files that I've added to the top of this file and access to the
 * GSL "library" file (I've added the library into the Makefile too). I found
 * GSL to be very easy to install but if there are problems in the future we
 * may need to switch to another matrix solver.
 *
 * We also clean for population inversion in this routine.
 *
 * The details are in Matthews' thesis.
 *
 **********************************************************/

int
macro_pops (xplasma, xne)
     PlasmaPtr xplasma;
     double xne;
{
  int i, j, index_element, index_lvl;
  int matrix_err, numerical_error, populations_ok;
  int n_macro_lvl;
  int n_iterations, n_inversions;
  double *a_data, *b_data;
  double *populations;
  double rate_matrix[NLEVELS_MACRO][NLEVELS_MACRO];
  int radiative_flag[NLEVELS_MACRO][NLEVELS_MACRO];     // array to flag if two levels are radiatively linked
  int conf_to_matrix[NLEVELS_MACRO];    // links config number to elements in arrays
  MacroPtr mplasma = &macromain[xplasma->nplasma];

  /*
   * In cells where there are no photons, we don't have any estimators yet
   * so we should first use dilute estimators to avoid problems further
   * down the road -- this should be a sensible choice, but could potentially
   * hide errors in very extreme cases.
   */

  if (xplasma->ntot == 0)
  {
    get_dilute_estimators (xplasma);
  }

  /* Start with an outer loop over elements: there are no rates that couple
     levels of different elements so we can always separate them out. */

  for (index_element = 0; index_element < nelements; index_element++)
  {
    /* Zero all elements of the matrix before doing anything else. */

    for (i = 0; i < NLEVELS_MACRO; i++)
    {
      for (j = 0; j < NLEVELS_MACRO; j++)
      {
        rate_matrix[j][i] = 0.0;
        radiative_flag[j][i] = 0;
      }
    }

    /* See if this element uses a macro atom treatment or is a simple element.
       For now I'm assuming that either all ions of a given element are
       treated using the macro atom method, or else none are (mixing and
       matching is probably a bad idea because of the way in which bf
       processes couple different ionisation stages). */

    if (ion[ele[index_element].firstion].macro_info == TRUE && geo.macro_simple == FALSE)       /* The check is against the first ion of the element. */
    {
      n_iterations = 0;
      populations_ok = FALSE;
      while (populations_ok == FALSE)
      {
        n_iterations++;

        /* Having established that the ion requires a macro atom treatment we
           are going to construct a matrix of rates between the levels and
           invert that matrix to get the level populations. The first thing we need
           to do is work out how many levels we are dealing with in total. This is
           easily done by summing up the number of levels of each ion. */

        n_macro_lvl = macro_pops_fill_rate_matrix (mplasma, xplasma, xne, index_element, rate_matrix, radiative_flag, conf_to_matrix);

        /* The rate matrix is now filled up. Since the problem is not closed as it stands, the next
           thing is to replace one of the rows of the matrix (say the first row) with the constraint
           that the sum of all the populations is 1.0 (this will let us get population fractions). */

        for (index_lvl = 0; index_lvl < n_macro_lvl; index_lvl++)
        {
          rate_matrix[0][index_lvl] = 1.0;
        }

        /* Now we can just invert the matrix to get the fractional level populations.
           The block that follows (down to next line of ***s) is to do the
           matrix inversion. It uses LU decomposition - the code for doing this is
           taken from the GSL manual with very few modifications.
           Here we solve the matrix equation M x = b, where x is our vector containing
           level populations as a fraction w.r.t the whole element */

        /* 211101 - ksl - Check added to avoid gcc11 warning */
        if (n_macro_lvl > SIZE_MAX / sizeof (double) || n_macro_lvl * n_macro_lvl > SIZE_MAX / sizeof (double))
        {
          Error ("macro_pops: n_macro_lvl %d too large for memory allocation\n", n_macro_lvl);
          Exit (EXIT_FAILURE);
        }
        else
        {

          a_data = (double *) calloc (n_macro_lvl * n_macro_lvl, sizeof (double));

          for (i = 0; i < n_macro_lvl; i++)
          {
            for (j = 0; j < n_macro_lvl; j++)
            {
              a_data[i * n_macro_lvl + j] = rate_matrix[i][j];  /* row-major ordering */
            }
          }

          b_data = (double *) calloc (n_macro_lvl, sizeof (double));
          populations = (double *) calloc (n_macro_lvl, sizeof (double));

          /* replace the first entry with 1.0 - this is part of the normalisation constraint */
          b_data[0] = 1.0;

          /* this next routine is a general routine which solves the matrix equation
             via LU decomposition */
          matrix_err = solve_matrix (a_data, b_data, n_macro_lvl, populations, xplasma->nplasma);

          free (a_data);
          free (b_data);

          if (matrix_err)
          {
            Error ("macro_pops: %s\n", get_matrix_error_string (matrix_err));
          }

          /* Now we take the population array and check to see if anything is very
           * small and set it to zero. This is basically some pre-emptive cleaning
           * since we could clean this up later, I suppose. */

          for (i = 0; i < n_macro_lvl; i++)
          {
            if (populations[i] < DENSITY_MIN)
            {
              populations[i] = 0.0;
            }
          }

          n_inversions = macro_pops_check_for_population_inversion (index_element, populations, radiative_flag, conf_to_matrix);

          if (n_inversions > 0)
            Debug ("macro_pops: iteration %d: there were %d levels which were cleaned due to population inversions in plasma cell %d\n",
                   n_iterations, n_inversions, xplasma->nplasma);

          /* 1 - IF the variable numerical_error has been set to TRUE then that means we had either a negative or
             non-finite level population somewhere. If that is the case, then set all the estimators
             to dilute blackbodies instead and go through the solution again.
             2 - IF we didn't set numerical_error to TRUE then we have a realistic set of populations, so set
             populations_ok to 1 to break the while loop, and copy the populations into the arrays
           */

          numerical_error =
            macro_pops_check_densities_for_numerical_errors (xplasma, index_element, populations, conf_to_matrix, n_iterations);

          if (numerical_error)
          {
            Error
              ("macro_pops: iteration %d: unreasonable population(s) in plasma cell %i. Using dilute BBody excitation with w %8.4e t_r %8.4e\n",
               n_iterations, xplasma->nplasma, xplasma->w, xplasma->t_r);
            get_dilute_estimators (xplasma);
          }
          else
          {
            populations_ok = TRUE;
            macro_pops_copy_to_xplasma (xplasma, index_element, populations, conf_to_matrix);
          }

          free (populations);

          if (n_iterations == MAXITERATIONS)
          {
            Error ("macro_pops: failed to converge for plasma cell %d\n", xplasma->nplasma);
            return EXIT_FAILURE;
          }
        }
      }                         // end of populations_ok == FALSE sane loop
    }                           // end of if statement for macro-atoms
  }                             // end of elements loop


  return (0);
}

/**********************************************************/
/**
 * @brief  Populate the rate matrix for the given element's ions and levels.
 *
 * @param[in] MacroPtr mplasma        The macro atom quantities for the current cell
 * @param[in] PlasmaPtr xplasma       The current plasma cell
 * @param[in] double xne              The current value of the electron density
 * @param[in] int index_element       The index of the current element to populate
 * @param[out] double rate_matrix     The populated rate matrix for the current element
 * @param[out] double radiative_flag  Flags for if two levels are radiatively linked
 * @param[out] int conf_to_matrix     A map to link congfiruation number to elements in
 *                                    the populations matrix
 *
 * @return int n_macro_lvl   The number of macro atom levels for this element
 *
 * @details
 * Having established that the ion requires a macro atom treatment we are going
 * to construct a matrix of rates between the levels and invert that matrix to
 * get the level populations. The first thing we need to do is work out how many
 * levels we are dealing with in total. This is easily done by summing up the
 * number of levels of each ion.
 *
 **********************************************************/

int
macro_pops_fill_rate_matrix (MacroPtr mplasma, PlasmaPtr xplasma, double xne, int index_element,
                             double rate_matrix[NLEVELS_MACRO][NLEVELS_MACRO], int radiative_flag[NLEVELS_MACRO][NLEVELS_MACRO],
                             int conf_to_matrix[NLEVELS_MACRO])
{
  int index_bbu, index_bbd;
  int index_bfu, index_bfd;
  int upper, lower;
  int index_fast_col;
  int index_ion, index_lvl;
  double rate;
  int n_macro_lvl = 0;
  struct lines *line_ptr;
  struct topbase_phot *cont_ptr;

  /* I want to be able to easily go from knowing the index of a level in the
     configurations structure to its position in the rates matrix. So I'm making
     two arrays here that allow the mapping between these indices to be done easily.
   */

  for (index_ion = ele[index_element].firstion; index_ion < (ele[index_element].firstion + ele[index_element].nions); index_ion++)
  {
    for (index_lvl = ion[index_ion].first_nlte_level; index_lvl < ion[index_ion].first_nlte_level + ion[index_ion].nlte; index_lvl++)
    {
      conf_to_matrix[index_lvl] = n_macro_lvl;
      n_macro_lvl++;
    }
  }

  /* We now know how many levels there are and therefore how big the matrix we
     need to invert will be. */

  /* Now we want to populate the matrix with all the rates between the levels. */

  for (index_ion = ele[index_element].firstion; index_ion < (ele[index_element].firstion + ele[index_element].nions); index_ion++)
  {
    index_lvl = ion[index_ion].first_nlte_level;

    /* The next loop judges whether or not a level is to be fixed in population relative to ground
       star. The input radiative lifetime is used to judge this at the moment. If the lifetime was set
       to be long (essentially infite) then a very fast collisional transition is put in to dominate
       all other rates into and out of this level.

       Whether this is really the best thing to do I don't know, but it's an improvement over ignoring
       this issue altogether! SS Aug 05 */

    for (index_fast_col = index_lvl; index_fast_col < ion[index_ion].first_nlte_level + ion[index_ion].nlte - 1; index_fast_col++)
    {
      {
        if (xconfig[index_fast_col + 1].rad_rate > 1.e15)
        {
          fast_line.gl = xconfig[index_lvl].g;
          fast_line.gu = xconfig[index_fast_col + 1].g;
          fast_line.freq = (xconfig[index_fast_col + 1].ex - xconfig[index_lvl].ex) / PLANCK;
          fast_line.f = 1e4;
          rate = q12 (&fast_line, xplasma->t_e) * xne;
          lower = conf_to_matrix[index_lvl];
          upper = conf_to_matrix[index_fast_col + 1];
          rate_matrix[lower][lower] += -1. * rate;
          rate_matrix[upper][lower] += rate;
          rate = q21 (&fast_line, xplasma->t_e) * xne;
          rate_matrix[upper][upper] += -1. * rate;
          rate_matrix[lower][upper] += rate;
        }
      }
    }

    for (index_lvl = ion[index_ion].first_nlte_level; index_lvl < ion[index_ion].first_nlte_level + ion[index_ion].nlte; index_lvl++)
    {
      /* Now add contribution for all the bb and bf processes between the levels. This is
         done by looping over the numbers "bbu, bbd, bfu, bfd" which tell us how many
         processes there are. */

      for (index_bbu = 0; index_bbu < xconfig[index_lvl].n_bbu_jump; index_bbu++)
      {
        /* These are bb upwards jumps. The rate in such a jump is given by
           Jbar which has been computed as a Monte Carlo estimator. I'm also
           including a collisional term (which depends on ne). */

        line_ptr = &line[xconfig[index_lvl].bbu_jump[index_bbu]];
        rate = b12 (line_ptr) * mplasma->jbar_old[xconfig[index_lvl].bbu_indx_first + index_bbu];
        rate += q12 (line_ptr, xplasma->t_e) * xne;

        /* This is the rate out of the level in question. We need to add it
           to the matrix in two places: firstly as a -ve contribution to the
           diagonal and secondly as a +ve contribution for the off-diagonal
           corresponding to the level populated by this process. */

        /* Get the matix indices for the upper and lower states of the jump. */

        lower = conf_to_matrix[index_lvl];
        upper = conf_to_matrix[line_ptr->nconfigu];

        rate_matrix[lower][lower] += -1. * rate;
        rate_matrix[upper][lower] += rate;

        if (rate < 0.0 || sane_check (rate))
        {
          Error ("macro_pops: bbu rate is %8.4e in plasma cell/matom %i\n", rate, xplasma->nplasma);
        }

        /* There's a radiative jump between these levels, so we want to clean
           for popualtion inversions. Flag this jump */
        radiative_flag[index_lvl][line_ptr->nconfigu] = 1;
      }

      for (index_bbd = 0; index_bbd < xconfig[index_lvl].n_bbd_jump; index_bbd++)
      {
        /* These are bb downwards jumps. The rate in such a jump is given by
           the A-value. I'm also
           including a collisional term (which depends on ne). */

        line_ptr = &line[xconfig[index_lvl].bbd_jump[index_bbd]];
        rate = (a21 (line_ptr) * p_escape (line_ptr, xplasma));
        rate += q21 (line_ptr, xplasma->t_e) * xne;

        /* This is the rate out of the level in question. We need to add it
           to the matrix in two places: firstly as a -ve contribution to the
           diagonal and secondly as a +ve contribution for the off-diagonal
           corresponding to the level populated by this process. */

        /* Get the matix indices for the upper and lower states of the jump. */

        upper = conf_to_matrix[index_lvl];
        lower = conf_to_matrix[line_ptr->nconfigl];

        rate_matrix[upper][upper] += -1. * rate;
        rate_matrix[lower][upper] += rate;

        /* There's a radiative jump between these levels, so we want to clean
           for population inversions. Flag this jump */
        radiative_flag[line_ptr->nconfigl][index_lvl] = 1;

        if (rate < 0.0 || sane_check (rate))
        {
          Error ("macro_pops: bbd rate is %8.4e in plasma cell/matom %i\n", rate, xplasma->nplasma);
        }
      }

      for (index_bfu = 0; index_bfu < xconfig[index_lvl].n_bfu_jump; index_bfu++)
      {
        /* These are bf upwards jumps. The rate in such a jump is given by
           gamma which has been computed as a Monte Carlo estimator. */

        cont_ptr = &phot_top[xconfig[index_lvl].bfu_jump[index_bfu]];
        rate = mplasma->gamma_old[xconfig[index_lvl].bfu_indx_first + index_bfu];
        rate += q_ioniz (cont_ptr, xplasma->t_e) * xne;

        /* This is the rate out of the level in question. We need to add it
           to the matrix in two places: firstly as a -ve contribution to the
           diagonal and secondly as a +ve contribution for the off-diagonal
           corresponding to the level populated by this process. */

        /* Get the matrix indices for the upper and lower states of the jump. */

        lower = conf_to_matrix[index_lvl];
        upper = conf_to_matrix[cont_ptr->uplev];

        rate_matrix[lower][lower] += -1. * rate;
        rate_matrix[upper][lower] += rate;

        if (rate < 0.0 || sane_check (rate))
        {
          Error ("macro_pops: bfu rate is %8.4e in plasma cell/matom %i\n", rate, xplasma->nplasma);
        }

        /* Now deal with the stimulated emission. */
        /* Lower and upper are the same, but now it contributes in the
           other direction. */

        rate = mplasma->alpha_st_old[xconfig[index_lvl].bfu_indx_first + index_bfu] * xne;

        rate_matrix[upper][upper] += -1. * rate;
        rate_matrix[lower][upper] += rate;

        if (rate < 0.0 || sane_check (rate))
        {
          Error ("macro_pops: st. recomb rate is %8.4e in plasma cell/matom %i\n", rate, xplasma->nplasma);
        }
      }

      for (index_bfd = 0; index_bfd < xconfig[index_lvl].n_bfd_jump; index_bfd++)
      {
        /* These are bf downwards jumps. The rate in such a jump is given by
           the alpha value. */
        cont_ptr = &phot_top[xconfig[index_lvl].bfd_jump[index_bfd]];
        /* Get new values of the recombination rates and store them. */
        mplasma->recomb_sp[xconfig[index_lvl].bfd_indx_first + index_bfd] = alpha_sp (cont_ptr, xplasma, 0);
        mplasma->recomb_sp_e[xconfig[index_lvl].bfd_indx_first + index_bfd] = alpha_sp (cont_ptr, xplasma, 2);
        rate = mplasma->recomb_sp[xconfig[index_lvl].bfd_indx_first + index_bfd] * xne;
        rate += q_recomb (cont_ptr, xplasma->t_e) * xne * xne;

        /* This is the rate out of the level in question. We need to add it
           to the matrix in two places: firstly as a -ve contribution to the
           diagonal and secondly as a +ve contribution for the off-diagonal
           corresponding to the level populated by this process. */

        /* Get the matrix indices for the upper and lower states of the jump. */

        upper = conf_to_matrix[index_lvl];
        lower = conf_to_matrix[cont_ptr->nlev];

        rate_matrix[upper][upper] += -1. * rate;
        rate_matrix[lower][upper] += rate;

        if (rate < 0.0 || sane_check (rate))
        {
          Error ("macro_pops: bfd rate is %8.4e in plasma cell/matom %i\n", rate, xplasma->nplasma);
        }
      }
    }
  }

  return n_macro_lvl;
}

/**********************************************************/
/**
 * @brief  Check for population inversions and deal with them appropriately.
 *
 * @param[in]  int index_element     The index for the element
 * @param[in]  double *populations   The calculated population densities
 * @param[in]  int **radiative_flag  Flags for if two levels are radiatively linked
 * @param[in]  int *conf_to_matrix   A map to link congfiruation number to elements in
 *                                    the populations matrix
 *
 * @return  void
 *
 * @details
 * MC noise can cause population inversions (particularly amongst highly excited
 * states) which are never a good thing and are most likely unphysical. We
 * therefore follow Leon Lucy's procedure (Lucy 2003) and remove inversions.
 *
 **********************************************************/

int
macro_pops_check_for_population_inversion (int index_element, double *populations, int radiative_flag[NLEVELS_MACRO][NLEVELS_MACRO],
                                           int conf_to_matrix[NLEVELS_MACRO])
{
  int i, index_ion, index_lvl;
  double inversion_test;
  int n_total_inversions = 0;

  for (index_ion = ele[index_element].firstion; index_ion < (ele[index_element].firstion + ele[index_element].nions); index_ion++)
  {
    /* Start loop with lowest level of the ion. For each level in turn check to see if there's a population
       inversion i.e. is  upper_pop > lower_pop * g_upper / g_lower. If it is then replace upper_pop with
       lower_pop * g_upper / g_lower. We loop over all levels higher than the currently chosen lower level. */
    for (index_lvl = ion[index_ion].first_nlte_level; index_lvl < ion[index_ion].first_nlte_level + ion[index_ion].nlte; index_lvl++)
    {
      for (i = index_lvl + 1; i < (ion[index_ion].first_nlte_level + ion[index_ion].nlte); i++)
      {
        /* this if statement means we only clean if there's a radiative jump between the levels */
        if (radiative_flag[index_lvl][i])
        {
          // todo: learn why we have a correction factor at the end
          inversion_test = populations[conf_to_matrix[index_lvl]] * xconfig[i].g / xconfig[index_lvl].g * 0.999999;

          if (populations[conf_to_matrix[i]] > inversion_test)
          {
            n_total_inversions += 1;
            populations[conf_to_matrix[i]] = inversion_test;
          }
        }
      }
    }
  }

  return n_total_inversions;
}

/**********************************************************/
/**
 * @brief  Check for negative or non-finite values in the calculated ion
 *         densities and populations.
 *
 * @param[in] PlasmaPtr xplasma        The plasma cell in question
 * @param[in] int  index_element       The index of the element in question
 * @param[in,out] double *populations  The calculated population densities
 * @param[in] int *conf_to_matrix      ?
 *
 * @return  TRUE if there has been a numerical issue, FALSE otherwise
 *
 * @details
 *
 **********************************************************/

int
macro_pops_check_densities_for_numerical_errors (PlasmaPtr xplasma, int index_element, double *populations,
                                                 int conf_to_matrix[NLEVELS_MACRO], int n_iterations)
{
  int index_ion, index_lvl;
  double this_ion_density, ion_density_temp;
  double level_population;

  int numerical_error = FALSE;

  for (index_ion = ele[index_element].firstion; index_ion < (ele[index_element].firstion + ele[index_element].nions); index_ion++)
  {
    this_ion_density = 0.0;
    for (index_lvl = ion[index_ion].first_nlte_level; index_lvl < ion[index_ion].first_nlte_level + ion[index_ion].nlte; index_lvl++)
    {
      level_population = populations[conf_to_matrix[index_lvl]];
      this_ion_density += level_population;
    }

    /* Check that the ion density is positive and finite */

    ion_density_temp = this_ion_density * ele[index_element].abun * xplasma->rho * rho2nh;
    if (sane_check (ion_density_temp) || ion_density_temp < 0.0)
    {
      Error ("macro_pops: iteration %d: ion %i has calculated a frac. pop. of %8.4e in plasma cell %i\n", n_iterations, index_ion,
             ion_density_temp, xplasma->nplasma);
      numerical_error = TRUE;
    }

    /* Check that the level populations for this ion are positive and finite */

    for (index_lvl = ion[index_ion].first_nlte_level; index_lvl < ion[index_ion].first_nlte_level + ion[index_ion].nlte; index_lvl++)
    {
      if (populations[conf_to_matrix[index_lvl]] < 0.0 || sane_check (populations[conf_to_matrix[index_lvl]]))
      {
        Error ("macro_pops: iteration %d: level %i has a calculated pop. of %8.4e in plasma cell %i\n",
               n_iterations, index_lvl, populations[conf_to_matrix[index_lvl]], xplasma->nplasma);
        numerical_error = TRUE;
      }
    }
  }

  return numerical_error;
}

/**********************************************************/
/**
 * @brief  Copy the calculated populations into the xplasma arrays.
 *
 * @param[in, out] PlasmaPtr xplasma  The plasma cell to update
 * @param[in] int index_element       The index of the element
 * @param[in] double *populations     The calculated densities to copy
 * @param[in] int *conf_to_matrix     A map to link congfiruation number to elements in
 *                                    the populations matrix
 *
 * @return void
 *
 * @details
 * The populations are now known. The populations need to be stored firstly as
 * ion populations and secondly as fractional level populations within an ion.
 * Get the ion populations and write them to one->density[nion]. The level
 * populations are to be put in "levden".
 *
 **********************************************************/

void
macro_pops_copy_to_xplasma (PlasmaPtr xplasma, int index_element, double *populations, int conf_to_matrix[NLEVELS_MACRO])
{
  int index_ion, index_lvl;
  double this_ion_density, fractional_population;

  for (index_ion = ele[index_element].firstion; index_ion < (ele[index_element].firstion + ele[index_element].nions); index_ion++)
  {
    this_ion_density = 0.0;
    for (index_lvl = ion[index_ion].first_nlte_level; index_lvl < ion[index_ion].first_nlte_level + ion[index_ion].nlte; index_lvl++)
    {
      this_ion_density += populations[conf_to_matrix[index_lvl]];
    }

    xplasma->density[index_ion] = this_ion_density * ele[index_element].abun * xplasma->rho * rho2nh;

    /* to maintain consistency with the higher level routines, only allow density to drop to DENSITY_MIN */

    if (xplasma->density[index_ion] < DENSITY_MIN)
    {
      xplasma->density[index_ion] = DENSITY_MIN;
    }

    for (index_lvl = ion[index_ion].first_nlte_level; index_lvl < ion[index_ion].first_nlte_level + ion[index_ion].nlte; index_lvl++)
    {
      fractional_population = populations[conf_to_matrix[index_lvl]] / this_ion_density;
      if (this_ion_density < DENSITY_MIN || fractional_population < DENSITY_MIN)
      {
        xplasma->levden[xconfig[index_lvl].nden] = DENSITY_MIN;
      }
      else
      {
        xplasma->levden[xconfig[index_lvl].nden] = fractional_population;
      }
    }
  }
}
