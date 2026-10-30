/***********************************************************/
/** @file  matom.c
 * @author SS,JM,SWM
 * @date   February, 2004
 * @brief  Macro-atom functions
 *
 * File containing Macro-atom functions.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_integration.h>

#include "atomic.h"
#include "sirocco.h"

/* External variables for use in matom */

double jprbs_known[NLEVELS_MACRO][2 * (NBBJUMPS + NBFJUMPS)], eprbs_known[NLEVELS_MACRO][2 * (NBBJUMPS + NBFJUMPS)];
double pjnorm_known[NLEVELS_MACRO], penorm_known[NLEVELS_MACRO];
int prbs_known[NLEVELS_MACRO];
int matom_cell = -1;
int matom_z = -1;
int matom_cycle = -1;

/**********************************************************/
/**
 * @brief The core of the implementation of Macro Atoms in sirocco
 *
 * @param [in,out]  PhotPtr p   the packet at the point of activation and deactivation
 * @param [in,out]  int nres    the process which activates and deactivates the Macro Atom
 * @param [out]  int escape  flag to tell us whether the matom de-activation
 *                             is via an r-packet (TRUE) or a k-packet (FALSE)
 * @return returns the number of level jumps before exiting the routine, or -1
 * if an error was encountered.
 *
 * @details
 * Matom is the core of the implementation of Leon Lucy's Macro Atom
 * approach to dealing with radiation-matter interactions. It is called whenever
 * a photon packet activates a macro atom. As input it takes the process by
 * which activation occurred (identified by the label "nres") which allow it
 * to deduce the level that has been excited. It then calculates all the
 * Macro Atom jumping/deactivation probabilities following Lucy and so deduces
 * the process by which deactivation occurs.
 *
 * The two possibilities are that the macro-atom deactivates by an r-packet
 * (radiative deactivations) or a k-packet (putting the energy  into the
 * thermal pool) as indicated by the variable escape.
 *
 * At the point where the routine exits, nres idenifies the process of
 * deactivation and the packet information is also updated
 *
 *
 * ###Notes###
 *
 * in macro atom mode, nres = -1 indicates electron scattering,
 * nres = -2 indicates ff, and nres > NLINES indicates bound-free.
 * [nres == NLINES is never used. Note also that NLINES is the *max* number of lines, whereas nlines
 * is the *actual* number of lines. So, actually, it's not just nres = NLINES that's never used, but
 * the entire range of nlines <= nres <= NLINES]
 *
 *
 *
***********************************************************/

int
matom (PhotPtr p, int *nres, int *escape)
{
  struct lines *line_ptr;
  struct topbase_phot *cont_ptr;
  struct auger *auger_ptr;
  int uplvl, uplvl_old;
//OLD  int icheck;
  double jprbs[2 * (NBBJUMPS + NBFJUMPS)];
  double eprbs[NBBJUMPS + NBFJUMPS];
  double pjnorm, penorm;
  double threshold, run_tot;
  double sp_rec_rate;
  int n, m;
  int njumps;
  int nbbd, nbbu, nbfd, nbfu;
  double t_e, ne;
  double bb_cont, choice, bf_cont;
  WindPtr one;
  double rad_rate, coll_rate;
  PlasmaPtr xplasma;
  MacroPtr mplasma;
  int z;
  int nauger, iauger, target_level;
  double auger_rate;


  one = &wmain[p->grid];
  xplasma = &plasmamain[one->nplasma];
  if (check_plasma (xplasma, "matom"))
  {
    Error ("matom:Working in dummy plasma cell\n");
    return (-1);
  }

  mplasma = &macromain[xplasma->nplasma];

  t_e = xplasma->state.t_e;
  ne = xplasma->state.ne;



  /* The first step is to identify the configuration that has been excited.
   * If *nres < NLINES the level will have been excited by a bb transioion
   * but if it is greater than this, it will bave been produced by photoionization.
   */

  uplvl = 0;
  if (*nres < NLINES)
  {
    uplvl = lin_ptr[*nres]->nconfigu;
    z = lin_ptr[*nres]->z;

  }
  else if (*nres > NLINES)
  {
    uplvl = phot_top[*nres - NLINES - 1].uplev;
    z = phot_top[*nres - NLINES - 1].z;
  }
  else
  {
    Error ("matom: upper level not identified. nres = %d in photon %d of cycle %d/%d in thread %d\n",
           *nres, p->np, geo.wcycle, geo.pcycle, rank_global);
    *escape = TRUE;
    p->istat = P_ERROR_MATOM;
    return (-1);
  }

  if (z != matom_z || p->grid != matom_cell || geo.wcycle != matom_cycle)
  {
    for (n = 0; n < NLEVELS_MACRO; n++)
    {
      prbs_known[n] = FALSE;
    }
    matom_z = z;
    matom_cell = p->grid;
    matom_cycle = geo.wcycle;
  }

  /* Now follows the main loop to govern the macro atom jumps. Keeps jumping until
     an emission occurs (at which point it breaks). */


  for (njumps = 0; njumps < MAXJUMPS; njumps++)
  {
    /*  The excited configuration is now known. Now compute all the probabilities of deactivation
       /jumping from this configuration. Then choose one. */

    nbbd = xconfig[uplvl].n_bbd_jump;   // number of bb downward jumps
    nbbu = xconfig[uplvl].n_bbu_jump;   // number of bb upward jump from this configuration
    nbfd = xconfig[uplvl].n_bfd_jump;   // number of bf downward jumps from this transition
    nbfu = xconfig[uplvl].n_bfu_jump;   // number of bf upward jumps from this transiion
    nauger = xconfig[uplvl].nauger;     // number of auger jumps
    iauger = xconfig[uplvl].iauger;

    if (prbs_known[uplvl] == FALSE)
    {

      m = 0;                    //m counts the total number of possible ways to leave the level
      pjnorm = 0.0;             //stores the total jump probability
      penorm = 0.0;             //stores the total emission probability

      for (n = 0; n < nbbd + nbfd; n++)
      {
        eprbs_known[uplvl][n] = eprbs[n] = 0;   //zero the individual emission probabilities
        jprbs_known[uplvl][n] = jprbs[n] = 0;   //zero  the individual jump probabilities
      }
      for (n = nbbd + nbfd; n < nbbd + nbfd + nbbu + nbfu; n++)
      {
        jprbs_known[uplvl][n] = jprbs[n] = 0;   /*slots for upward jumps */
      }


      /* bb */

      /* First downward jumps. (I.e. those that have emission probabilities. */

      /* For bound-bound decays the jump probability is A-coeff * escape-probability * energy */
      /* At present the escape probability is only approximated (p_escape). This should be improved. */
      /* The collisional contribution to both the jumping and deactivation probabilities are now added (SS, Apr04) */

      for (n = 0; n < nbbd; n++)
      {
        line_ptr = &line[xconfig[uplvl].bbd_jump[n]];

        rad_rate = (a21 (line_ptr) * p_escape (line_ptr, xplasma));
        coll_rate = ne * q21 (line_ptr, t_e);


        bb_cont = rad_rate + coll_rate;
        jprbs_known[uplvl][m] = jprbs[m] = bb_cont * xconfig[line_ptr->nconfigl].ex;    //energy of lower state

        eprbs_known[uplvl][m] = eprbs[m] = bb_cont * (xconfig[uplvl].ex - xconfig[line[xconfig[uplvl].bbd_jump[n]].nconfigl].ex);       //energy difference


        pjnorm += jprbs[m];
        penorm += eprbs[m];
        m++;
      }

      /* bf */
      for (n = 0; n < nbfd; n++)
      {

        cont_ptr = &phot_top[xconfig[uplvl].bfd_jump[n]];       //pointer to continuum

        sp_rec_rate = mplasma->est.recomb_sp[xconfig[uplvl].bfd_indx_first + n];
        bf_cont = (sp_rec_rate + q_recomb (cont_ptr, t_e) * ne) * ne;

        jprbs_known[uplvl][m] = jprbs[m] = bf_cont * xconfig[phot_top[xconfig[uplvl].bfd_jump[n]].nlev].ex;     //energy of lower state
        eprbs_known[uplvl][m] = eprbs[m] = bf_cont * (xconfig[uplvl].ex - xconfig[phot_top[xconfig[uplvl].bfd_jump[n]].nlev].ex);       //energy difference
        pjnorm += jprbs[m];
        penorm += eprbs[m];
        m++;
      }

      /* Auger ionization */
      if (iauger >= 0)
      {
        auger_ptr = &auger_macro[iauger];

        for (n = 0; n < nauger; n++)
        {
          target_level = auger_ptr->nconfig_target[n];
          auger_rate = auger_ptr->Avalue_auger * auger_ptr->branching_ratio[n];

          jprbs_known[uplvl][m] = jprbs[m] = auger_rate * xconfig[target_level].ex;     //energy of lower state
          eprbs_known[uplvl][m] = eprbs[m] = auger_rate * (xconfig[uplvl].ex - xconfig[target_level].ex);       //energy difference

          pjnorm += jprbs[m];
          penorm += eprbs[m];
          m++;
        }
      }

      /* Now upwards jumps. */

      /* bb */
      /* For bound-bound excitation the jump probability is B-coeff times Jbar with a correction
         for stimulated emission. To avoid the need for recalculation all the time, the code will
         be designed to include the stimulated correction in Jbar - i.e. the stimulated correction
         factor will NOT be included here. (SS) */
      /* There is no emission probability for upwards transitions. */
      /* Collisional contribution to jumping probability added. (SS,Apr04) */

      for (n = 0; n < nbbu; n++)
      {
        line_ptr = &line[xconfig[uplvl].bbu_jump[n]];
        rad_rate = (b12 (line_ptr) * mplasma->state.jbar_old[xconfig[uplvl].bbu_indx_first + n]);

        coll_rate = ne * q12 (line_ptr, t_e);   // this is multiplied by ne below


        jprbs_known[uplvl][m] = jprbs[m] = (rad_rate + coll_rate) * xconfig[uplvl].ex;  //energy of lower state


        pjnorm += jprbs[m];
        m++;
      }

      /* bf */
      for (n = 0; n < nbfu; n++)
      {
        /* For bf ionization the jump probability is just gamma * energy
           gamma is the photoionisation rate. Stimulated recombination also included. */
        cont_ptr = &phot_top[xconfig[uplvl].bfu_jump[n]];       //pointer to continuum

        // jumping probability, so uses energy of lower state 
        jprbs_known[uplvl][m] = jprbs[m] =
          (mplasma->state.gamma_old[xconfig[uplvl].bfu_indx_first + n] -
           (mplasma->state.alpha_st_old[xconfig[uplvl].bfu_indx_first + n] * stim_recomb_factor (xplasma, cont_ptr)) +
           (q_ioniz (cont_ptr, t_e) * ne)) * xconfig[uplvl].ex;

        /* this error condition can happen in unconverged hot cells where T_R >> T_E.
           for the moment we set to 0 and hope spontaneous recombiantion takes care of things */
        /* note that we check and report this in check_stimulated_recomb() in estimators.c once a cycle */
        if (jprbs[m] < 0.)      //test (can be deleted eventually SS)
        {
          //Error ("Negative probability (matom, 6). Abort?\n");
          jprbs_known[uplvl][m] = jprbs[m] = 0.0;

        }
        pjnorm += jprbs[m];
        m++;
      }

      pjnorm_known[uplvl] = pjnorm;
      penorm_known[uplvl] = penorm;
      prbs_known[uplvl] = TRUE;

    }


    if ((pjnorm_known[uplvl] + penorm_known[uplvl]) <= 0.0)
    {
      Error ("matom: macro atom level has no way out: uplvl %d pj %g pe %g t_e %.3g  ne %.3g\n", uplvl, pjnorm_known[uplvl],
             penorm_known[uplvl], t_e, ne);
      Error ("matom: macro atom level has no way out: z %d istate %d nion %d ilv %d nbfu %d nbfd %d nbbu %d nbbd %d\n", xconfig[uplvl].z,
             xconfig[uplvl].istate, xconfig[uplvl].nion, xconfig[uplvl].ilv, nbfu, nbfd, nbbu, nbbd);
      *escape = TRUE;
      p->istat = P_ERROR_MATOM;
      return (-1);
    }

    /* Probabilities of jumping (j) and emission (e) are now known. The normalisation factor pnorm has
       also been recorded. The integer m now gives the total number of possibilities too.
       Now select what happens next. Start by choosing the random threshold value at which the
       event will occur. */


    threshold = random_number (0.0, 1.0);
    if (((pjnorm_known[uplvl] / (pjnorm_known[uplvl] + penorm_known[uplvl])) < threshold) || (pjnorm_known[uplvl] == 0))
      break;                    // A deactivation of the macro-atom has occurred and so we leave the for loop.

    /* Othewise, a transition/jump to another state of the macro-atom has occurred, so we need
       to decide what the new state is. We use a running total to decide the new upper level */

    uplvl_old = uplvl;

    threshold = random_number (0.0, 1.0);
    threshold = threshold * pjnorm_known[uplvl_old];

    n = 0;
    run_tot = 0;
    while (run_tot < threshold)
    {
      run_tot += jprbs_known[uplvl_old][n];
      n++;
    }
    /* This added to prevent case where threshold is essentially 0.
     */

    if (n > 0)
    {
      n = n - 1;
    }

    /* n now identifies the jump that occurs - now set the new level. */
    if (n < nbbd)
    {                           /* bb downwards jump */
      uplvl = line[xconfig[uplvl].bbd_jump[n]].nconfigl;
    }
    else if (n < (nbbd + nbfd))
    {                           /* bf downwards jump */
      uplvl = phot_top[xconfig[uplvl].bfd_jump[n - nbbd]].nlev;
    }
    else if (n < (nbbd + nbfd + nbbu))
    {                           /* bb upwards jump */
      uplvl = line[xconfig[uplvl].bbu_jump[n - nbbd - nbfd]].nconfigu;
    }
    else if (n < (nbbd + nbfd + nbbu + nbfu))
    {                           /* bf upwards jump */
      uplvl = phot_top[xconfig[uplvl].bfu_jump[n - nbbd - nbfd - nbbu]].uplev;
    }
    else if (n < (nbbd + nbfd + nbbu + nbfu + nauger))
    {                           /* auger ionization jump */
      uplvl = auger_macro[iauger].nconfig_target[n - nbbd - nbfd - nbbu - nbfu];
//OLD      icheck = 4;
    }
    else
    {
      Error ("Trying to jump but nowhere to go! Matom. Abort");
      *escape = TRUE;
      p->istat = P_ERROR_MATOM;
      return (-1);
    }

  }


  if (njumps == MAXJUMPS)
  {
    Error ("Matom: jumped %d times with no emission for photon %d from upper level %d  pjnorm %e penorm %e Abort.\n", MAXJUMPS, p->np,
           uplvl, pjnorm_known[uplvl], penorm_known[uplvl]);
    *escape = TRUE;
    p->istat = P_ERROR_MATOM;
    return (-1);
  }


  /* At this point we know the macro actom has deactivated and we know the level
   * at which the deactivation occurs.  We still need to decide the exact procees
   * by which the macro actom deactivates.
   */

  run_tot = 0;
  n = 0;

  threshold = random_number (0.0, 1.0);

  threshold = threshold * penorm_known[uplvl];  //normalise to total emission prob.

  while (run_tot < threshold)
  {
    run_tot += eprbs_known[uplvl][n];
    n++;
  }
  n = n - 1;
  /* n now identifies the jump that occurs - now set nres for the return value. */
  if (n < nbbd)
  {                             /* bb downwards jump */
    /* With collisions included we need to decide whether the deactivation is radiative (r-packet)
       or collisional (k-packet). Get a random number and then use the ratio of the collisional
       emission probability to the (already known) collisional+radiative probability to decide whether
       collisional or radiative deactivation occurs. */
    choice = random_number (0.0, 1.0);


    line_ptr = &line[xconfig[uplvl].bbd_jump[n]];       //pointer for the bb transition

    rad_rate = a21 (line_ptr) * p_escape (line_ptr, xplasma);

    coll_rate = ne * q21 (line_ptr, t_e);

    if (choice > (coll_rate / (rad_rate + coll_rate)))
    {
      /* It's a r-packet (radiative) deactivation */
      *escape = TRUE;
      *nres = line[xconfig[uplvl].bbd_jump[n]].where_in_list;
      p->freq = line[xconfig[uplvl].bbd_jump[n]].freq;
    }
    else
    {
      /* It's a k-packet (collisional deactivation. */
      *escape = FALSE;
    }
  }
  else if (n < (nbbd + nbfd))
  {                             /* bf downwards jump */
    /* With collisional recombination included we need to decide whether to deactivate
       radiatively or make a k-packet. */

    cont_ptr = &phot_top[xconfig[uplvl].bfd_jump[n - nbbd]];

    rad_rate = mplasma->est.recomb_sp[xconfig[uplvl].bfd_indx_first + n - nbbd];        //again using recomb_sp rather than alpha_sp (SS July 04)
    coll_rate = ne * q_recomb (cont_ptr, t_e);

    choice = random_number (0.0, 1.0);

    if (choice > (coll_rate / (rad_rate + coll_rate)))
    {
      /*It's a r-packet (radiative) deactivation */
      *escape = TRUE;
      *nres = xconfig[uplvl].bfd_jump[n - nbbd] + NLINES + 1;
      p->freq = matom_select_bf_freq (one, xconfig[uplvl].bfd_jump[n - nbbd]);

    }
    else
    {
      /* It's a k-packet (collisional) deactivation */
      *escape = FALSE;
    }

  }
  else if (n < (nbbd + nbfd + nauger))
  {                             /* auger ionization */
    /* Auger processes only undergo k-packet (collisional) deactivation */
    *escape = FALSE;
  }
  else
  {
    Error ("Trying to emitt from Macro Atom but no available route (matom). Abort.");
    *escape = TRUE;
    p->istat = P_ERROR_MATOM;
    return (-1);
  }

  return (njumps);
}





/************************************/
/**
 * @brief the b12 Einstein coefficient.
 *
 * @param [in] struct lines line_ptr line pointer to calculate
 *
 * @return  retruns the B12 Einstein coefficient.
 * ###Notes###
 * ksl OK, Probably should be moved to lines.c for consistency, but that can be done
 * later.
***********************************/

#define B12_CONSTANT 5.01983e25

struct lines *b12_line_ptr;
double b12_a;

double
b12 (struct lines *line_ptr)
{
  double freq;

  if (b12_line_ptr != line_ptr)
  {
    freq = line_ptr->freq;
    b12_a = B12_CONSTANT * line_ptr->f / freq;
    b12_line_ptr = line_ptr;
  }

  return (b12_a);
}

/**********************************************************/
/**
 * @brief the stimulated-recombination factor n_e n_u / n_l of a bf continuum
 *
 * @param [in] PlasmaPtr xplasma the plasma cell
 * @param [in] struct topbase_phot cont_ptr the continuum
 * @return n_e n_u / n_l, with n_l the lower (ionizing) and n_u the upper (recombined-from) level density; 0 if n_l <= 0
 *
 * @details
 * The net photoionization rate per lower-level ion is gamma - alpha_st * stim_recomb_factor. It is the one
 * definition used by matom, macro_accelerate and scatter, which clip a negative net rate to zero, in the same
 * way as the bb stimulated-emission correction in normalise_macro_estimators.
 **********************************************************/

double
stim_recomb_factor (PlasmaPtr xplasma, struct topbase_phot *cont_ptr)
{
  double lower_density = den_config (xplasma, cont_ptr->nlev);

  if (lower_density <= 0.0)
    return (0.0);
  return (xplasma->state.ne * den_config (xplasma, cont_ptr->uplev) / lower_density);
}

/**********************************************************/
/**
 * @brief the number of Gauss-Legendre points that integrate exp (g (t)) over a piece to BF_RATE_EPSILON
 *
 * @param [in] double h the half-width of the piece in t
 * @param [in] double d1 an upper bound of |g'| on the piece
 * @param [in] double d2 an upper bound of |g''| on the piece
 * @return the number of points, at most BF_RATE_MAX_POINTS
 *
 * @details
 * On [-1, 1] the n-point rule has the error C_n f^(2n) (xi), C_n = (n!)^4 2^(2n+1) / ((2n+1) ((2n)!)^3). For
 * f = exp (G) with |G'| <= h d1 and |G''| <= h^2 d2, |f^(2n)| <= f L^(2n) with L = h (d1 + sqrt (2 n d2)) (exactly so
 * where G is quadratic; the higher derivatives of the exponential factors carry further powers of h), and the mean
 * of f is at least max f / (1 + 2 L), so the relative error is at most C_n / 2 (1 + 2 L) L^(2n). bf_rate_reach[n]
 * is an L at which this is at most BF_RATE_EPSILON: with A = (2 BF_RATE_EPSILON / C_n)^(1 / 2n),
 * A (1 + 2 A)^(-1 / 2n) is just below the root.
***********************************************************/

#define BF_RATE_EPSILON DBL_EPSILON     /**< relative error of the integrals of bf_rate_integrals */
#define BF_RATE_MAX_POINTS 20   /**< the largest rule; GSL tabulates every rule up to 20 points */

static gsl_integration_glfixed_table *bf_rate_rule[BF_RATE_MAX_POINTS + 1];
static double bf_rate_reach[BF_RATE_MAX_POINTS + 1], bf_rate_sqrt_2n[BF_RATE_MAX_POINTS + 1];

/** @brief allocate the Gauss-Legendre rules and fill bf_rate_reach, once */

static void
bf_rate_rules_init (void)
{
  int n;

  for (n = 1; n <= BF_RATE_MAX_POINTS; n++)
  {
    const double log_c = 4 * lgamma (n + 1.) + (2 * n + 1) * M_LN2 - log (2 * n + 1.) - 3 * lgamma (2 * n + 1.);
    const double a = exp ((log (2 * BF_RATE_EPSILON) - log_c) / (2 * n));

    if ((bf_rate_rule[n] = gsl_integration_glfixed_table_alloc (n)) == NULL)
    {
      Error ("bf_rate_rules_init: could not allocate the %d-point Gauss-Legendre rule\n", n);
      Exit (EXIT_FAILURE);
    }
    bf_rate_reach[n] = a * pow (1 + 2 * a, -0.5 / n);
    bf_rate_sqrt_2n[n] = sqrt (2. * n);
  }
}

static int
bf_rate_points (double h, double d1, double d2)
{
  const double sqrt_d2 = sqrt (d2);
  int n;

  for (n = 1; n < BF_RATE_MAX_POINTS; n++)
  {
    if (h * (d1 + bf_rate_sqrt_2n[n] * sqrt_d2) <= bf_rate_reach[n])
      break;
  }
  return (n);
}

/** @brief add one node of the Gauss-Legendre sum of bf_rate_segment, at t with nu - nu_threshold = above */

static void
bf_rate_node (double t, double above, double weight, double slope, double beta, double beta_r, double fthresh, double sum[3])
{
  const double freq = fthresh + above;
  double f = weight * exp ((slope + 3) * t - beta * above - beta_r * freq);

  if (beta_r > 0)
    f /= -expm1 (-beta_r * freq);
  sum[0] += f;
  sum[1] += f * freq / fthresh;
  sum[2] += f * above / fthresh;
}

/**********************************************************/
/**
 * @brief add the integrals of bf_rate_integrals from lo to hi, inside table interval i
 *
 * @details
 * In t = ln (nu / nu_i) the integrand times dnu / dt is x_i nu_i^3 exp ((s + 3) t - beta (nu - nu_threshold) -
 * beta_r nu) / (1 - exp (-beta_r nu)), with s the slope of the interval. Its log g has |g'| <= |s + 3.5 - (beta +
 * beta_r) nu| + 0.5 (the weight nu / nu_threshold adds 1 to s + 3) + 1 (from the Planck factor, if any) and
 * |g''| <= (beta + beta_r) nu. The interval is split into the fewest equal pieces that the largest rule
 * integrates, and each piece gets the fewest points that bf_rate_points allows. The nodes come in mirror pairs
 * t_mid +- d, at nu_mid exp (+-d), so one expm1 gives the frequencies of both; frequencies are carried as
 * nu - nu_threshold, so that beta (nu - nu_threshold) keeps its precision near the threshold.
***********************************************************/

static void
bf_rate_segment (struct topbase_phot *cont_ptr, int i, double lo, double hi, double beta, double beta_r, double sum[3])
{
  const double fthresh = cont_ptr->freq[0], freq_i = cont_ptr->freq[i], big_b = beta + beta_r;
  /* the width from the ratio: a difference of the stored log frequencies loses digits in a narrow interval */
  const double t_end = log (cont_ptr->freq[i + 1] / freq_i);
  const double slope = (cont_ptr->log_x[i + 1] - cont_ptr->log_x[i]) / t_end;
  const double scale = cont_ptr->x[i] * freq_i * freq_i * freq_i, c = slope + 3.5, planck = beta_r > 0;
  const double t_lo = log (lo / freq_i), t_hi = hi < cont_ptr->freq[i + 1] ? log (hi / freq_i) : t_end;
  const double d1 = fmax (fabs (c - big_b * lo), fabs (c - big_b * hi)) + 0.5 + planck;
  const double reach = 0.5 * (t_hi - t_lo) * (d1 + bf_rate_sqrt_2n[BF_RATE_MAX_POINTS] * sqrt (big_b * hi));
  const int pieces = reach > bf_rate_reach[BF_RATE_MAX_POINTS] ? (int) ceil (reach / bf_rate_reach[BF_RATE_MAX_POINTS]) : 1;
  const double h = 0.5 * (t_hi - t_lo) / pieces;
  double freq0 = lo;
  int k, l;

  for (k = 0; k < pieces; k++)
  {
    const double t_mid = t_lo + (2 * k + 1) * h, freq1 = k == pieces - 1 ? hi : freq_i * exp (t_mid + h);
    const double above_mid = (freq_i - fthresh) + freq_i * expm1 (t_mid), freq_mid = fthresh + above_mid;
    const int n = bf_rate_points (h, fmax (fabs (c - big_b * freq0), fabs (c - big_b * freq1)) + 0.5 + planck, big_b * freq1);
    const gsl_integration_glfixed_table *rule = bf_rate_rule[n];

    /* GSL stores the (n + 1) / 2 non-negative nodes on [-1, 1], with 0 first for odd n */
    for (l = 0; l < (n + 1) / 2; l++)
    {
      const double d = h * rule->x[l], weight = h * rule->w[l] * scale, e = expm1 (d);

      bf_rate_node (t_mid + d, above_mid + freq_mid * e, weight, slope, beta, beta_r, fthresh, sum);
      if (d > 0)
        bf_rate_node (t_mid - d, above_mid - freq_mid * e / (1 + e), weight, slope, beta, beta_r, fthresh, sum);
    }
    freq0 = freq1;
  }
}

/** @brief add the integrals of bf_rate_integrals from lo to hi, table interval by table interval */

static void
bf_rate_range (struct topbase_phot *cont_ptr, double lo, double hi, double beta, double beta_r, double sum[3])
{
  int i;

  for (i = 0; i < cont_ptr->np - 1 && cont_ptr->freq[i] < hi; i++)
  {
    const double a = fmax (lo, cont_ptr->freq[i]), b = fmin (hi, cont_ptr->freq[i + 1]);

    /* sigma_phot is zero inside an interval with a zero cross section at either end */
    if (b > a && cont_ptr->x[i] > 0 && cont_ptr->x[i + 1] > 0)
      bf_rate_segment (cont_ptr, i, a, b, beta, beta_r, sum);
  }
}

/**********************************************************/
/**
 * @brief the frequency integrals of a photoionization continuum that the macro-atom rates need
 *
 * @param [in] struct topbase_phot cont_ptr the continuum
 * @param [in] double freq_min, freq_max the frequency range, clipped to the continuum's table
 * @param [in] double beta h / k T_e of the Maxwellian factor exp (beta (nu_threshold - nu)); 0 for none
 * @param [in] double beta_r h / k T_r of the Planck factor 1 / (exp (beta_r nu) - 1); 0 for none
 * @param [out] double sum[3] the integrals of sigma (nu) nu^2 exp (beta (nu_threshold - nu)) / (exp (beta_r nu) - 1)
 * times 1, nu / nu_threshold and (nu - nu_threshold) / nu_threshold, each to about BF_RATE_EPSILON relative
 *
 * @details
 * sigma_phot interpolates the cross section linearly in log-log between the table points, so between two points the
 * integrand is a power law times exponentials, smooth whatever the resonance structure of the table; bf_rate_segment
 * integrates it there.
 *
 * The exponential factors fall by exp (-B dnu), B = beta + beta_r, so the integrals are first taken up to
 * nu_1 = freq_min - ln (BF_RATE_EPSILON) / B. Beyond nu_1 the rest is at most
 *
 *    R = x_max / nu_threshold exp (beta (nu_threshold - nu_1) - beta_r nu_1) / (1 - exp (-beta_r nu_1))
 *        (nu_1^3 / B + 3 nu_1^2 / B^2 + 6 nu_1 / B^3 + 6 / B^4)
 *
 * (x_max the largest cross section of the table, nu^2 and the weights at most nu^3 / nu_threshold), and since nu^3
 * grows by at most exp (3 dnu / nu_1), with B nu_1 > 3, the rest beyond nu_1 + ln (R / F) / (B - 3 / nu_1) is below F
 * = BF_RATE_EPSILON times the smallest sum. The integrals are completed to there. Infinite beta or beta_r (zero
 * temperature) give zero integrals.
***********************************************************/

void
bf_rate_integrals (struct topbase_phot *cont_ptr, double freq_min, double freq_max, double beta, double beta_r, double sum[3])
{
  const double fthresh = cont_ptr->freq[0], big_b = beta + beta_r;
  const double lo = fmax (freq_min, fthresh), end = fmin (freq_max, cont_ptr->freq[cont_ptr->np - 1]);
  double nu_1, rest, floor, x_max = 0;
  int i;

  sum[0] = sum[1] = sum[2] = 0;
  if (!(beta >= 0 && beta_r >= 0))
  {
    Error ("bf_rate_integrals: exponential factors %e %e are not >= 0\n", beta, beta_r);
    return;
  }
  if (isinf (big_b))
    return;
  if (bf_rate_rule[1] == NULL)
    bf_rate_rules_init ();

  nu_1 = big_b > 0 ? fmin (end, lo - log (BF_RATE_EPSILON) / big_b) : end;
  bf_rate_range (cont_ptr, lo, nu_1, beta, beta_r, sum);
  if (nu_1 == end)
    return;

  for (i = 0; i < cont_ptr->np; i++)
    x_max = fmax (x_max, cont_ptr->x[i]);
  rest = x_max / fthresh * exp (beta * (fthresh - nu_1) - beta_r * nu_1)
    * (((nu_1 / big_b + 3 / (big_b * big_b)) * nu_1 + 6 / (big_b * big_b * big_b)) * nu_1 + 6 / (big_b * big_b * big_b * big_b));
  if (beta_r > 0)
    rest /= -expm1 (-beta_r * nu_1);
  floor = BF_RATE_EPSILON * fmin (sum[0], sum[2]);
  if (rest > floor)
    bf_rate_range (cont_ptr, nu_1, fmin (end, nu_1 + log (rest / floor) / (big_b - 3 / nu_1)), beta, beta_r, sum);
}



/**********************************************************/
/**
 * @brief the integrals of alpha_sp_all at temperature temp, remembering the last ALPHA_SP_MEMO temperatures of
 * each continuum
 *
 * @details
 * A cell's temperature search asks for the same integrals again and again: the cooling is evaluated at both ends
 * of a bracket, and the level populations at the final T_e in every iteration of the ionization solve. So each
 * continuum keeps its last ALPHA_SP_MEMO temperatures and their integrals and gives them back, bit for bit, for an
 * equal temperature. Only positive temperatures are kept, so any other still reaches bf_rate_integrals (and its
 * error message). alpha_sp_forget clears the memory when the cross sections are freed.
***********************************************************/

#define ALPHA_SP_MEMO 2         /**< temperatures kept per continuum: the two ends of a bracket */

static struct alpha_sp_memo
{
  double temp[ALPHA_SP_MEMO], sum[ALPHA_SP_MEMO][3];
  int last;                     /**< the entry used last; the other one is replaced next */
} *alpha_sp_memo = NULL;

/** @brief forget the integrals remembered by alpha_sp_sums */

void
alpha_sp_forget (void)
{
  free (alpha_sp_memo);
  alpha_sp_memo = NULL;
}

static void
alpha_sp_sums (struct topbase_phot *cont_ptr, double temp, double sum[3])
{
  const long ncont = cont_ptr - phot_top;
  struct alpha_sp_memo *memo;
  int k;

  if (!(temp > 0 && ncont >= 0 && ncont < nphot_total))
  {
    bf_rate_integrals (cont_ptr, cont_ptr->freq[0], cont_ptr->freq[cont_ptr->np - 1], H_OVER_K / temp, 0, sum);
    return;
  }
  if (alpha_sp_memo == NULL && (alpha_sp_memo = calloc (nphot_total, sizeof *alpha_sp_memo)) == NULL)
  {
    Error ("alpha_sp_sums: could not allocate the memory of the integrals\n");
    Exit (EXIT_FAILURE);
  }
  memo = &alpha_sp_memo[ncont];
  for (k = 0; k < ALPHA_SP_MEMO; k++)
  {
    if (memo->temp[k] == temp)
    {
      memcpy (sum, memo->sum[k], sizeof memo->sum[k]);
      memo->last = k;
      return;
    }
  }
  bf_rate_integrals (cont_ptr, cont_ptr->freq[0], cont_ptr->freq[cont_ptr->np - 1], H_OVER_K / temp, 0, sum);
  memo->last = k = (memo->last + 1) % ALPHA_SP_MEMO;
  memo->temp[k] = temp;
  memcpy (memo->sum[k], sum, sizeof memo->sum[k]);
}

/* The methods used before #1219, kept as options (@Diag.bf_integration sum or gsl).  The GSL integrand
   gets the continuum, the temperature and the choice through these externals. */

static struct topbase_phot *cont_ext_ptr;       //continuum pointer passed externally
static double temp_ext;         //temperature passed externally
static int temp_choice;         //choice of type of calculation for alpha_sp

/**********************************************************/
/**
 *  @brief the integrand of alpha_sp at a chosen frequency, for num_int
 *
 * @param [in] double freq
 * @param [in] void *params
 *
 * Note:
 *  The peculiar variable *params is needed by the gsl routine
 *  that does the integration.
***********************************************************/

static double
alpha_sp_integrand (double freq, void *params)
{
  double fthresh;
  double x;
  double integrand;
  double tt;

  fthresh = cont_ext_ptr->freq[0];
  tt = temp_ext;

  if (freq < fthresh)
    return (0.0);               // No recombination at frequencies lower than the threshold freq occur

  x = sigma_phot (cont_ext_ptr, freq);  //this is the cross-section

  integrand = x * freq * freq * exp (H_OVER_K * (fthresh - freq) / tt);


  if (temp_choice == 1)
    return (integrand * freq / fthresh);        //energy weighted case
  if (temp_choice == 2)
    return (integrand * (freq - fthresh) / fthresh);    // difference case
  return (integrand);           //spontaneous case
}

/**********************************************************/
/**
 * @brief the integrals of alpha_sp by the methods used before #1219
 *
 * @param [in] struct topbase_phot cont_ptr the continuum
 * @param [in] double temp the electron temperature
 * @param [in] int want[3] which of the integrals (ichoice 0, 1, 2) to calculate
 * @param [out] double integral[3] the integrals asked for
 *
 * @details
 * The integrals stop at ALPHA_MATOM_NUMAX_LIMIT k T / h above the threshold.  With BF_INTEG_SUM they are a
 * 10,000-point sum in log nu, all three from one pass; with BF_INTEG_GSL each is a call to num_int.
***********************************************************/

static void
alpha_sp_integrals_old (struct topbase_phot *cont_ptr, double temp, const int want[3], double integral[3])
{
  double fthresh, flast;
  int ichoice;

  fthresh = cont_ptr->freq[0];  //first frequency in list
  flast = cont_ptr->freq[cont_ptr->np - 1];     //last frequency in list
  if ((H_OVER_K * (flast - fthresh) / temp) > ALPHA_MATOM_NUMAX_LIMIT)
  {
    //flast is currently very far into the exponential tail: so reduce flast to limit value of h nu / k T.
    flast = fthresh + temp * ALPHA_MATOM_NUMAX_LIMIT / H_OVER_K;
  }

  integral[0] = integral[1] = integral[2] = 0;

  if (modes.bf_integration == BF_INTEG_GSL)
  {
    temp_ext = temp;
    cont_ext_ptr = cont_ptr;
    for (ichoice = 0; ichoice < 3; ichoice++)
    {
      if (want[ichoice])
      {
        temp_choice = ichoice;
        integral[ichoice] = num_int (alpha_sp_integrand, fthresh, flast, 1e-4);
      }
    }
    return;
  }

  double freq, dfreq;
  double freq1, freq2;
  double x, integrand;
  int n, nmax;
  double log_fmin, dlog_freq;

  nmax = 10000;
  log_fmin = log10 (cont_ptr->freq[0]);
  dlog_freq = log10 (flast / cont_ptr->freq[0]) / nmax;

  freq1 = cont_ptr->freq[0];
  for (n = 1; n < nmax; n++)
  {
    freq2 = pow (10., log_fmin + n * dlog_freq);
    freq = 0.5 * (freq1 + freq2);
    dfreq = freq2 - freq1;
    x = sigma_phot (cont_ptr, freq);
    freq1 = freq2;

    if (freq > flast)
    {
      break;
    }

    integrand = x * freq * freq * exp (H_OVER_K * (fthresh - freq) / temp);

    if (want[0])
      integral[0] += integrand * dfreq;
    if (want[1])
      integral[1] += integrand * (freq / fthresh) * dfreq;      //energy weighted case
    if (want[2])
      integral[2] += integrand * ((freq - fthresh) / fthresh) * dfreq;  // difference case
  }
}

/**********************************************************/
/**
 * @brief the alpha_sp rates asked for, by the methods used before #1219
 *
 * @details
 * The integrals of alpha_sp_integrals_old times the constants, in the same order of operations as before,
 * so that the results are those of the earlier code.
***********************************************************/

static void
alpha_sp_old (struct topbase_phot *cont_ptr, PlasmaPtr xplasma, const int want[3], double alpha[3])
{
  int ichoice;

  alpha_sp_integrals_old (cont_ptr, xplasma->state.t_e, want, alpha);
  for (ichoice = 0; ichoice < 3; ichoice++)
  {
    if (!want[ichoice])
      continue;
    if (cont_ptr->macro_info == TRUE && geo.macro_simple == FALSE)
    {
      alpha[ichoice] = alpha[ichoice] * xconfig[cont_ptr->nlev].g / xconfig[cont_ptr->uplev].g * pow (xplasma->state.t_e, -1.5);
    }
    else                        //case for simple element
    {
      alpha[ichoice] = alpha[ichoice] * xconfig[cont_ptr->nlev].g / ion[cont_ptr->nion + 1].g * pow (xplasma->state.t_e, -1.5); //g for next ion up used
    }

    alpha[ichoice] = alpha[ichoice] * ALPHA_SP_CONSTANT;

    if (sane_check (alpha[ichoice]))
    {
      Error ("alpha_sp:sane_check value is %e\n", alpha[ichoice]);
    }
  }
}

/**********************************************************/
/**
 * @brief the matom estimator for the spontaneous recombination rate.
 *
 * @param [in] struct topbase_phot cont_ptr pointer to calculate
 * @param [in] struct PlasmaPtr xplasma the plasma cell of interesest
 * @param [in] int ichoice one of several types or rates to calculate
 *
 * @return the recombination rate is returned
 *
 * @details
 * The rate is given by
 *
 *    (4 pi /c2) (gu/gl) (h2/2 pi m k T)^(3/2)
 * times the integral of   a(nu) nu2 exp [(chi- h nu)/kT].
 *
 * The choices are
 * * ichoice = 0   --> spontanous recombination
 * * ichoice = 1   --> energy weighted recombination
 * * ichoice = 2   --> the difference between energy_weighted and spontaneous
 *
 * ###Notes###
 *
 * This basically coumputs integrals similar to eqn 13 of Lucy 2003
 *
 *  Energy weighted means that the integrand has an extra factor nu/nu_threshold
 *  The difference case is (nu-nu_threshold)/nu_threhold
***********************************************************/

double
alpha_sp (struct topbase_phot *cont_ptr, PlasmaPtr xplasma, int ichoice)
{
  double alpha[3];

  if (modes.bf_integration != BF_INTEG_EXACT)
  {
    int want[3] = { 0, 0, 0 };

    want[ichoice] = 1;
    alpha_sp_old (cont_ptr, xplasma, want, alpha);
    return (alpha[ichoice]);
  }
  alpha_sp_all (cont_ptr, xplasma, alpha);
  return (alpha[ichoice]);
}

/**********************************************************/
/**
 * @brief all three matom spontaneous recombination estimators of one continuum at once
 *
 * @param [in] struct topbase_phot cont_ptr pointer to calculate
 * @param [in] struct PlasmaPtr xplasma the plasma cell of interesest
 * @param [out] double alpha[3] the rates for ichoice 0, 1 and 2 of alpha_sp
 *
 * @details
 * The integrals, from bf_rate_integrals, run over the whole cross-section table.  With @Diag.bf_integration
 * sum or gsl they are done as before #1219 instead, see alpha_sp_old.
***********************************************************/

void
alpha_sp_all (struct topbase_phot *cont_ptr, PlasmaPtr xplasma, double alpha[3])
{
  const double temp = xplasma->state.t_e;
  double g_upper;
  int ichoice;

  if (modes.bf_integration != BF_INTEG_EXACT)
  {
    const int want[3] = { 1, 1, 1 };

    alpha_sp_old (cont_ptr, xplasma, want, alpha);
    return;
  }

  if (cont_ptr->macro_info == TRUE && geo.macro_simple == FALSE)
  {
    g_upper = xconfig[cont_ptr->uplev].g;
  }
  else                          //case for simple element
  {
    g_upper = ion[cont_ptr->nion + 1].g;        //g for next ion up used
  }

  alpha_sp_sums (cont_ptr, temp, alpha);
  for (ichoice = 0; ichoice < 3; ichoice++)
  {
    alpha[ichoice] *= xconfig[cont_ptr->nlev].g / g_upper * pow (temp, -1.5) * ALPHA_SP_CONSTANT;
    if (sane_check (alpha[ichoice]))
    {
      Error ("alpha_sp:sane_check value is %e\n", alpha[ichoice]);
    }
  }
}

/**********************************************************/
/**
 * @brief the band-limited integral of alpha_sp by the method used before #1219 (num_int)
 *
 * @details
 * Used for both @Diag.bf_integration sum and gsl, as before.  The integral stops at ALPHA_MATOM_NUMAX_LIMIT
 * k T / h above the threshold.
***********************************************************/

static double
scaled_alpha_sp_integral_band_limited_old (struct topbase_phot *cont_ptr, PlasmaPtr xplasma, int ichoice, double freq_min, double freq_max)
{
  double alpha_sp_value;
  double fthresh, flast;

  temp_choice = ichoice;
  temp_ext = xplasma->state.t_e;        //external for use in alph_sp_integrand
  cont_ext_ptr = cont_ptr;      //"
  fthresh = cont_ptr->freq[0];  //first frequency in list
  flast = cont_ptr->freq[cont_ptr->np - 1];     //last frequency in list
  if ((H_OVER_K * (flast - fthresh) / temp_ext) > ALPHA_MATOM_NUMAX_LIMIT)
  {
    //flast is currently very far into the exponential tail: so reduce flast to limit value of h nu / k T.
    flast = fthresh + temp_ext * ALPHA_MATOM_NUMAX_LIMIT / H_OVER_K;
  }
  if (flast < freq_max)
  {
    freq_max = flast;
  }
  alpha_sp_value = num_int (alpha_sp_integrand, freq_min, freq_max, 1e-4);

  return (alpha_sp_value);
}

/**********************************************************/
/**
 *  @brief the matom estimator for the band-limited spontaneous recombination rate.
 *
 * @param [in] struct topbase_phot cont_ptr pointer to calculate
 * @param [in] struct PlasmaPtr xplasma the plasma cell of interesest
 * @param [in] int ichoice one of several types or rates to calculate
 * @param [in] double freq_min  the mimimum frequency
 * @param [in] double freq_max  the maximum frequencey
 *
 * @returns  The routine returns the spontanous recombination rate, in various forms
 *
 * @details
 *
 * This routine differs from alpha_sp in that the reate is calculated over
 * a specific wavelength range, and it returns the integral alone, without the constant factors of alpha_sp.
 *
 * The choices are:
 *
 * * ichoice = 0   --> spontanous recombination
 * * ichoice = 1   --> energy weighted recombination
 * * ichoice = 2   --> the difference between energy_weighted spontaneous
***********************************************************/

double
scaled_alpha_sp_integral_band_limited (struct topbase_phot *cont_ptr, PlasmaPtr xplasma, int ichoice, double freq_min, double freq_max)
{
  double sum[3];

  if (modes.bf_integration != BF_INTEG_EXACT)
    return (scaled_alpha_sp_integral_band_limited_old (cont_ptr, xplasma, ichoice, freq_min, freq_max));

  /* the whole continuum is the same integrals as alpha_sp_all's */
  if (freq_min <= cont_ptr->freq[0] && freq_max >= cont_ptr->freq[cont_ptr->np - 1])
    alpha_sp_sums (cont_ptr, xplasma->state.t_e, sum);
  else
    bf_rate_integrals (cont_ptr, freq_min, freq_max, H_OVER_K / xplasma->state.t_e, 0, sum);
  return (sum[ichoice]);
}



/**********************************************************/
/**
 * @brief deals with the elimination of k-packets.
 *
 *
 * @param [in,out]  PhotPtr p   the packet at the point of activation and deactivation
 * @param [in,out]  int nres    the process which activates and deactivates the Macro Atom
 * @param [out]  int escape  to tell us whether the matom de-activation
 *                             is via an r-packet (1 TRUE) or a k-packet (0 FALSE)
 * @param [in] int mode         switch which allows photon to be deactivated by a non-radiative
 * term.  (non_zero is true)
 * @return 0
 * @details
 *
 * Whenever a k-packet is created it is
 * immediately destroyed insitu by this routine. At output "nres" identified the process
 * that destroys the k-packet and the packet information has been updated in the same
 * way as at the end of matom
 *
 * ###Notes###
************************************************************/

int
kpkt (PhotPtr p, int *nres, int *escape, int mode)
{

  int i;
  double cooling_adiabatic;
  double cooling_normalisation;
  double destruction_choice;
  double upweight_factor;
  WindPtr xwind;
  PlasmaPtr xplasma;
  MacroPtr mplasma;
  double freqmin, freqmax;
  double cooling_bbtot, cooling_bf_coltot, cooling_bb_use;


  /* Idea is to calculated the cooling
     terms for all the processes and then choose one at random to destroy the k-packet and
     turn it back into a photon bundle.
     The routine considers bound-free, collision excitation and ff
     emission. */

  /* The loop here runs over all the bf processes, regardless of whether they are macro or simple
     ions. The reason that we can do this is that in the macro atom method stimulated recombination
     has already been considered before creating the k-packet (by the use of gamma-twiddle) in "scatter"
     and so we only have spontaneous recombination to worry about here for ALL cases. */


  xwind = &wmain[p->grid];
  xplasma = &plasmamain[xwind->nplasma];
  if (check_plasma (xplasma, "kpkt"))
  {
    Error ("kpkt:Photon appears to be having interaction in wind cell associated with dummy plasma structure\n");
  }

  mplasma = &macromain[xplasma->nplasma];

//OLD  electron_temperature = xplasma->state.t_e;

  /* Set maximum and minimum frequency limits.  See #187. We need band limits for free free packet
     generation (see call to one_ff below). A bandaid is applied if there is not enough
     separation between freqmin and max */

  freqmin = xband.f1[0];
  freqmax = ALPHA_FF * xplasma->state.t_e / H_OVER_K;
  if (freqmax < 1.1 * freqmin)
  {
    freqmax = 1.1 * freqmin;
  }

  /* If the kpkt destruction rates for this cell are not known they are calculated here.  This happens
   * every time the wind is updated. Note that most of the cooling rates do not include a factor
   * of ne, since one power of ne applies to all of the cooling rates and we are only concerned
   * with relative cooling rates */

  if (mplasma->derived.kpkt_rates_known != TRUE)
  {
    fill_kpkt_rates (xplasma, escape, p);
  }

/* This is the end of the cooling rate calculation, which is done only once for each cell
   and once for each cycle.

   The next little section deals whith handling adiabatic cooling and shock heating.
   */

  cooling_normalisation =
    mplasma->est.cooling_normalisation - mplasma->est.cooling_adiabatic - mplasma->est.cooling_bbtot - mplasma->est.cooling_bf_coltot;
  cooling_adiabatic = 0.0;

  /* if kpkt mode is all processes, or continuum + adiabatic, then include adiabatic cooling */
  if (mode == KPKT_MODE_ALL || mode == KPKT_MODE_CONT_PLUS_ADIABATIC)
  {
    if (KPKT_NET_HEAT_MODE && geo.nonthermal)
    {
      if (xplasma->derived.cool_adiabatic > xplasma->derived.heat_shock)
      {
        cooling_adiabatic = (xplasma->derived.cool_adiabatic - xplasma->derived.heat_shock) / xplasma->state.vol / xplasma->state.ne;
      }
    }
    else
    {
      cooling_adiabatic = mplasma->est.cooling_adiabatic;
    }
  }
  cooling_normalisation += cooling_adiabatic;

  if (mode == KPKT_MODE_ALL)
  {
    cooling_bbtot = mplasma->est.cooling_bbtot;
    cooling_bf_coltot = mplasma->est.cooling_bf_coltot;
  }
  else
  {
    cooling_bbtot = mplasma->derived.cooling_bb_simple_tot;
    cooling_bf_coltot = 0.0;
    /* we don't conaider collisional ionization of simple ions as a cooling process */
  }
  cooling_normalisation += cooling_bbtot + cooling_bf_coltot;


  /* The cooling rates for the recombination and collisional processes are now known.

     Choose which process destroys the k-packet with a random number. */

  destruction_choice = random_number (0.0, 1.0) * cooling_normalisation;

  /* This logic of what follows may not be obvious.  For choosing the basic
   * process, we just look to see if the destruction choice is less than bf, bf+bb, bf+bb+ff
   * etc, but inside the bhe "basic_choices", we iteratively reduce! the destruction
   * choice until we get to one that is less than the cooling associated with a specific
   * tranistion. If we do not find such a transition, within for example the bf if statement
   * we drop all the way down to the Error at the end.
   */


  if (destruction_choice < mplasma->est.cooling_bftot)
  {                             //destruction by bf

    /* JM 1503 -- we used to loop over ntop_phot here,
       but we should really loop over the tabulated Verner Xsections too
       see #86, #141 */
    for (i = 0; i < nphot_total; i++)
    {
      if (destruction_choice < mplasma->est.cooling_bf[i])
      {

        if (i > nphot_total - 1)
        {
          Error ("kpkt (matom.c): trying to destroy k-packet in unknown process. Abort.\n");
          *escape = TRUE;
          p->istat = P_ERROR_MATOM;
          return (0);
        }


        *nres = i + NLINES + 1;
        *escape = TRUE;

        p->freq = matom_select_bf_freq (xwind, i);

        /* if the cross-section corresponds to a simple ion (macro_info == FALSE)
           or if we are treating all ions as simple, then adopt the total emissivity
           approach to choosing photon weights - this means we
           multipy down the photon weight by a factor nu/(nu-nu_0)
           and we force a kpkt to be created */
        if (modes.use_upweighting_of_simple_macro_atoms)
        {
          if (phot_top[i].macro_info == FALSE || geo.macro_simple == TRUE)
          {
            upweight_factor = xplasma->state.recomb_simple_upweight[i];
            p->w *= upweight_factor;

            /* record the amount of energy being extracted from the simple ion ionization pool */
            xplasma->derived.bf_simple_ionpool_out += p->w - (p->w / upweight_factor);
          }
        }

        return (0);
      }
      else
      {
        destruction_choice = destruction_choice - mplasma->est.cooling_bf[i];
      }
    }
  }
  else if (destruction_choice < (mplasma->est.cooling_bftot + cooling_bbtot))
  {
    /* a collisional destruction has occurred and so, if the line is  associated with
       a macro atom, it  must be excited.
     */

    destruction_choice = destruction_choice - mplasma->est.cooling_bftot;
    for (i = 0; i < nlines; i++)
    {
      /* this is a bit inelegant, but whether we want to consider the contribution
         here depends on the mode and type of line */
      if (mode == KPKT_MODE_ALL || line[i].macro_info == FALSE || geo.macro_simple == TRUE)
      {
        cooling_bb_use = mplasma->est.cooling_bb[i];
      }
      else
      {
        cooling_bb_use = 0.0;
      }

      if (destruction_choice < cooling_bb_use && cooling_bb_use != 0.0)
      {
        *nres = line[i].where_in_list;
        if (line[i].macro_info == TRUE && geo.macro_simple == FALSE)
        {
          *escape = FALSE;
        }
        else
        {
          *escape = TRUE;
          p->freq = line[i].freq;
        }
        return (0);
      }
      else
      {
        destruction_choice = destruction_choice - cooling_bb_use;
      }
    }
  }

  else if (destruction_choice < (mplasma->est.cooling_bftot + cooling_bbtot + mplasma->est.cooling_ff))
  {
    /* If reached this point, it is a FF destruction event */
    /* consult issues #187, #492 regarding free-free */
    *escape = TRUE;
    *nres = NRES_FF;
    p->freq = one_ff (xplasma, freqmin, freqmax);
    return (0);
  }
  else if (destruction_choice < (mplasma->est.cooling_bftot + cooling_bbtot + mplasma->est.cooling_ff + mplasma->est.cooling_ff_lofreq))
  {
/*this is ff at a frequency that is so low frequency that it is not worth tracking further */
    *escape = TRUE;
    *nres = NRES_FF;
    p->istat = P_LOFREQ_FF;
    return (0);
  }


  else if (destruction_choice <
           (mplasma->est.cooling_bftot + cooling_bbtot + mplasma->est.cooling_ff + mplasma->est.cooling_ff_lofreq + cooling_adiabatic))
  {
/* It is a k-packat that is destroyed by adiabatic cooling */

    if (geo.adiabatic == 0 || mode == KPKT_MODE_CONTINUUM)
    {
      Error ("kpkt: Destroying kpkt by adiabatic cooling even though it is turned off.\n");
    }
    *escape = TRUE;
    *nres = NRES_FF;
    p->istat = P_ADIABATIC;

    return (0);
  }


  else
  {
    /* It is a k-packed destroyed by collisional ionization in a macro atom. */
    destruction_choice =
      destruction_choice - mplasma->est.cooling_bftot - cooling_bbtot - mplasma->est.cooling_ff - mplasma->est.cooling_ff_lofreq -
      cooling_adiabatic;

    for (i = 0; i < nphot_total; i++)
    {
      if (destruction_choice < mplasma->est.cooling_bf_col[i])
      {
        if (i > nphot_total - 1)
        {
          Error ("kpkt: trying to destroy k-packet in unknown process. Abort.\n");
          *escape = TRUE;
          p->istat = P_ERROR_MATOM;
          return (0);
        }

        *nres = i + NLINES + 1;
        *escape = FALSE;

        return (0);
      }
      else
      {
        destruction_choice = destruction_choice - mplasma->est.cooling_bf_col[i];
      }
    }
  }



  Error ("kpkt: Failed to select a destruction process in kpkt. Abort.\n");
  Error
    ("kpkt: choice %8.4e norm %8.4e cooling_bftot %g, cooling_bbtot %g, cooling_ff %g, cooling_ff_lofreq %g, cooling_bf_coltot %g cooling_adiabatic %g cooling_adiabatic %g\n",
     destruction_choice, cooling_normalisation, mplasma->est.cooling_bftot, cooling_bbtot, mplasma->est.cooling_ff,
     mplasma->est.cooling_ff_lofreq, cooling_bf_coltot, mplasma->est.cooling_adiabatic, cooling_adiabatic);

  *escape = TRUE;
  p->istat = P_ERROR_MATOM;

  return (0);

}



/************************************************************
 **
 * @brief routine for dealing with bound-bound "simple ions" within the hybrid macro-atom framework
 *
 * @param [in,out] PhotPtr p   the packet at the point of activation
 * @param [in] int nres    the process which activates the Macro Atom
 * @param [out] int escape  TRUE if de-activation is
 *                           via an r-packet or FALSE if a k-packet.
 * @return The routine alwas returns 0.
 * However, as noted above the escape indicates
 * whether the deactivation is via an r- or k-packet.
 * If the deactivation is via an r-packet the frequency of the photon is
 * set to that of the line
 *
 * @details
 *
 * fake_matom_bb is the macro atom routine that deals with line events involving
 * simple ions (i.e. ions for which a full macro atom treatment is not employed.
 * When this routine is called a simple line has absorbed a packet. This routine
 * creates a fake two-level macro atom and determines whether the packet energy
 * is simply re-emitted in the line or is thermalised. If it is thermalised it
 * turns into a k-packet and the appropriate routine is called.
 *
************************************************************/

int
fake_matom_bb (PhotPtr p, int *nres, int *escape)
{
  double kprb, rprb;
  WindPtr one;
  PlasmaPtr xplasma;
  struct lines *line_ptr;
  double electron_temperature;
  double normalisation;
  double choice;

  one = &wmain[p->grid];
  xplasma = &plasmamain[one->nplasma];

  line_ptr = lin_ptr[*nres];
  electron_temperature = xplasma->state.t_e;

  /* Upon calling we know that the upper level of our fake two level macro
     atom is excited. Since it's only two-levels there are no jumping probabilities
     and the only decision that needs to be made is whether it de-excites by an r-packet
     in the line or a collisionally produced k-packet. To decide this we compute the
     two emission probabilities. rprb is the probability of the r-packet and kprb
     for the k-packet. */

  /* Since both k- and r-packet emission have the same delta energy associated with them then
     we can ignore the energy factor that is otherwise needed in the de-activation probabilities.
     The upper level population is also factored out.

     The relative rates are then given by
     Einstein-A * escape probability  for radiative deactivations  and
     Collision de-excitation rate coeff. * electron density for where
     The extra factor of (1 - exp (-h nu / k T)) is explained in KSL's notes on Sirocco. It appears because we are
     using the "scattering fraction" formalism for simple ions.

   */

  rprb = a21 (line_ptr) * p_escape (line_ptr, xplasma);

  kprb = q21 (line_ptr, electron_temperature) * xplasma->state.ne * (1. - exp (-H_OVER_K * line_ptr->freq / electron_temperature));

  normalisation = kprb + rprb;

  kprb = kprb / normalisation;

  rprb = rprb / normalisation;

  /* Now just use a random number to decide what happens.

     If "choice" is less than rprb then we have chosen a radiative decay - for this fake macro atom there
     is only one line so there's nothing to do - the energy is re-radiated in the line and that's it. We
     don't need to change nres. Otherwise we've chosen a collisional destruction and so we need to get
     kpkt to deal with the resulting k-packet and give us a new value of nres following thermalisation and
     re-emission of the energy. */

  choice = random_number (0.0, 1.0);

  if (choice < rprb)
  {
    *escape = TRUE;
    p->freq = line_ptr->freq;
  }
  else
  {
    *escape = FALSE;
  }

  return (0);

}


/************************************************************
 **
 * @brief calculate the frequency with for a bound-free transition for "simple ions"
 * within the hybrid macro-atom framework
 *
 * @param [in,out] PhotPtr p   the packet at the point of activation
 * @param [in]     int nres    the process which activates the Macro Atom
 * @param [out]   int escape  in principle this tells us whether de-activation is
 *                             via an r-packet or a k-packet. For this routine at the
 *                             moment only r-packets are possible so it always returns
 *                             escape = TRUE
 *
 * @return The routine always returns 0, and addition escape is always set to TRUE.
 *
 * @details
 * fake_matom_bf is the macro atom routine that deals with photoionisation
 * events involving simple ions (i.e. ions for which a full macro atom treatment
 * is not employed). When this routine is called a photoionisation has absorbed a packet.
 * The idea of this routine is to deal with the subsequenct
 * by creating a fake two-level atom. However, in the absense of collisional recombination
 * (or something similar) there's only one thing that can happen - radiative
 * recombination. Therefore there's no need to do anything here unless
 * collisional recombination is to be introduced at some point in the
 * future. All this routine does for now is choose a new frequency for the
 * emitted r-packet.
 *
************************************************************/

int
fake_matom_bf (PhotPtr p, int *nres, int *escape)
{
  WindPtr one;
//OLD  PlasmaPtr xplasma;

  one = &wmain[p->grid];
//OLD  xplasma = &plasmamain[one->nplasma];

  *escape = TRUE;

  p->freq = matom_select_bf_freq (one, *nres - NLINES - 1);


  return (0);

}




/**********************************************************/
/**
 * @brief   calculate the band limited frequency emitted by a macro-atom during spectral cycles.
 *
 * @param [in]     WindPtr w   the ptr to the structure defining the wind
 * @param [in,out]  PhotPtr p   the packet at the point of activation and deactivation
 * @param [in]     int upper   the upper level that we deactivate from
 * @param [in,out]  int nres    the process by which deactivation occurs
 * @param [in]     double freq_min  the lower limit to the desired frequency
 * @param [in]    double freq_max   the upper limit to the desired frequency
 * @return 0
 *
 * @details  emit matom is called to generate a photon frequency emitted in the
 * wind during spectral cycles in the matom case.  The upper level is known and
 * so the question is how to create a photon in the desired frequency range. For
 * BB transitions this is straightfoward but for recombination this is less clear.
 * and so there it is necessary to test for the correct frequency.
 *
 * ### Notes
 *
 * emit_matom is a scaled down version of matom which deals with the emission due
 * to deactivating macro atoms in the detailed spectrum part of the calculation.
 *
***********************************************************/

int
emit_matom (WindPtr w, PhotPtr p, int *nres, int upper, double freq_min, double freq_max)
{
  struct lines *line_ptr;
  struct topbase_phot *cont_ptr;
  int uplvl;
  double eprbs[NBBJUMPS + NBFJUMPS];
  double penorm;
  double threshold, run_tot;
  double sp_rec_rate;
  int n, m;
  int nbbd, nbfd;
//OLD  double t_e, ne;
  double ne;
  double bb_cont;
  WindPtr one;
  PlasmaPtr xplasma;
  MacroPtr mplasma;


  one = &w[p->grid];
  xplasma = &plasmamain[one->nplasma];
  mplasma = &macromain[one->nplasma];

//OLD  t_e = xplasma->state.t_e;
  ne = xplasma->state.ne;

  uplvl = upper;


  /*  The excited configuration is known. Compute all the probabilities of deactivation
     /jumping from this configuration. Then choose one. */


  nbbd = xconfig[uplvl].n_bbd_jump;     //number of bb downward jumps
  nbfd = xconfig[uplvl].n_bfd_jump;     //number of bf downared jumps


  // Set the emission probabilities and normalization factor to 0

  m = 0;
  penorm = 0.0;
  for (n = 0; n < nbbd + nbfd; n++)
  {
    eprbs[n] = 0;
  }

  /* bb */

  /* First consider BB transitions

     Since we are only interested in making an r-packet here we can (a) ignore collisional
     deactivation and (b) ignore lines outside the frequency range of interest.

   */

  for (n = 0; n < nbbd; n++)
  {
    line_ptr = &line[xconfig[uplvl].bbd_jump[n]];
    if ((line_ptr->freq > freq_min) && (line_ptr->freq < freq_max))     // correct range
    {
      bb_cont = (a21 (line_ptr) * p_escape (line_ptr, xplasma));

      eprbs[m] = bb_cont * (xconfig[uplvl].ex - xconfig[line[xconfig[uplvl].bbd_jump[n]].nconfigl].ex); //energy difference
      penorm += eprbs[m];
    }
    m++;
  }

  /* Now deal with BF
     If the edge is above the frequency range we are interested in then we need not consider this
     bf process.
   */

  for (n = 0; n < nbfd; n++)
  {
    cont_ptr = &phot_top[xconfig[uplvl].bfd_jump[n]];

    if (cont_ptr->freq[0] < freq_max)
    {
      sp_rec_rate = mplasma->est.recomb_sp[xconfig[uplvl].bfd_indx_first + n];

      eprbs[m] = sp_rec_rate * ne * (xconfig[uplvl].ex - xconfig[phot_top[xconfig[uplvl].bfd_jump[n]].nlev].ex);        //energy difference
      penorm += eprbs[m];
    }
    m++;
  }

  /* Probabilities of emission (e) are now known. The normalisation factor pnorm has
     also been recorded. the integer m now gives the total number of possibilities too.
     now select what happens next. Start by choosing the random threshold value at which the
     event will occur. */


  double xfreq;
  xfreq = -1;                   /* need to be below 0 to allow for case when freq_min = 0 */
  while (xfreq < freq_min || xfreq > freq_max)
  {
    threshold = random_number (0.0, 1.0);

    run_tot = 0;
    n = 0;
    threshold = threshold * penorm;     //normalise to total emission prob.
    while (run_tot < threshold)
    {
      run_tot += eprbs[n];
      n++;
    }
    n = n - 1;

    /* n now identifies the jump that occurs - now set nres for the return value. */

    if (n < nbbd)
    {
      line_ptr = &line[xconfig[uplvl].bbd_jump[n]];
      *nres = line[xconfig[uplvl].bbd_jump[n]].where_in_list;
      p->freq = xfreq = line[xconfig[uplvl].bbd_jump[n]].freq;
    }
    else if (n < (nbbd + nbfd))
    {
      *nres = xconfig[uplvl].bfd_jump[n - nbbd] + NLINES + 1;
      xfreq = matom_select_bf_freq (one, xconfig[uplvl].bfd_jump[n - nbbd]);
      p->freq = xfreq;
    }

    else
    {
      Error ("Trying to emit from Macro Atom but no available route (emit_matom). Abort.");
      Exit (0);
    }
  }

  return (0);
}
