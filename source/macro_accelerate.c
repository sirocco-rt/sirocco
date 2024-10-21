#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#include "atomic.h"
#include "sirocco.h"

/**********************************************************/
/**
 * @brief calculate the matrix of probabilities for the accelerated macro-atom scheme
 *
 * @param [in] PlasmaPtr  xplasma
 * @param [in,out] double **matom_matrix
 *        the 2D matrix array we will populate with normalised probabilities
 *
 * @details
 * given an activation state, this routine calculates the probability that a packet
 * will deactivate from a given state. Let's suppose that the \f$(i,j)\f$-th element of matrix
 * \f$Q\f$ contains the **jumping** probability from state i to state j, and the diagonals
 * \f$(i,i)\f$ of matrix \f$R\f$ contain the emission probabilities from each state, then it
 * can be shown that (see short notes from Vogl, or
 * <a href="Ergon et al. 2018">https://ui.adsabs.harvard.edu/abs/2018A%26A...620A.156E</a>)
 * the quantity we want is then \f$B = N R\f$, where \f$N = (I - Q)^{-1}\f$, where \f$I\f$
 * is the identity matrix. This routine does this calculation.
 *
 **********************************************************/

void
calc_matom_matrix (xplasma, matom_matrix)
     PlasmaPtr xplasma;
     double **matom_matrix;
{
  MacroPtr mplasma;
  double t_e, ne;
  int nbbd, nbbu, nbfd, nbfu;
  int uplvl, target_level, escape_dummy;
  double Qcont;
  struct lines *line_ptr;
  struct auger *auger_ptr;
  struct topbase_phot *cont_ptr;
  double rad_rate, coll_rate;
  int n, i, nn, mm, iauger, nauger;
  double Qcont_kpkt, bb_cont, sp_rec_rate, bf_cont, lower_density, density_ratio;
  double kpacket_to_rpacket_rate, norm, Rcont, auger_rate;
  int matrix_error;
  double *a_data, *a_inverse;
  mplasma = &macromain[xplasma->nplasma];       //telling us where in the matom structure we are
  struct photon pdummy;

  t_e = xplasma->t_e;           //electron temperature
  ne = xplasma->ne;             //electron number density

  /* allocate arrays for matrices and normalsiations
     we want to include all macro-atom levels + 1 kpacket levels */
  int nrows = nlevels_macro + 1;
  double **R_matrix = (double **) calloc (sizeof (double *), nrows);
  double **Q_matrix = (double **) calloc (sizeof (double *), nrows);
  double *Q_norm = (double *) calloc (sizeof (double), nrows);

  for (i = 0; i < nrows; i++)
  {
    R_matrix[i] = (double *) calloc (sizeof (double), nrows);
    Q_matrix[i] = (double *) calloc (sizeof (double), nrows);
  }

  /* initialise everything to zero */
  norm = 0.0;
  for (uplvl = 0; uplvl < nrows; uplvl++)
  {
    Q_norm[uplvl] = 0.0;
    for (target_level = 0; target_level < nrows; target_level++)
    {
      Q_matrix[uplvl][target_level] = 0.0;
      R_matrix[uplvl][target_level] = 0.0;
      matom_matrix[uplvl][target_level] = 0.0;
    }
  }

  /* loop over all macro-atom levels and populate the rate matrix */
  for (uplvl = 0; uplvl < nlevels_macro; uplvl++)
  {
    nbbd = xconfig[uplvl].n_bbd_jump;   //store these for easy access -- number of bb downward jumps
    nbbu = xconfig[uplvl].n_bbu_jump;   // number of bb upward jump from this configuration
    nbfd = xconfig[uplvl].n_bfd_jump;   // number of bf downward jumps from this transition
    nbfu = xconfig[uplvl].n_bfu_jump;   // number of bf upward jumps from this transiion
    nauger = xconfig[uplvl].nauger;     /* number of Auger jumps */
    iauger = xconfig[uplvl].iauger;

    /* bound-bound */
    for (n = 0; n < nbbd; n++)
    {

      line_ptr = &line[xconfig[uplvl].bbd_jump[n]];

      rad_rate = a21 (line_ptr) * p_escape (line_ptr, xplasma);
      coll_rate = q21 (line_ptr, t_e);  // this is multiplied by ne below

      bb_cont = rad_rate + (coll_rate * ne);

      target_level = line_ptr->nconfigl;

      //internal jump to another macro atom level
      Q_matrix[uplvl][target_level] += Qcont = bb_cont * xconfig[target_level].ex;      //energy of lower state

      //jump to the k-packet pool (we used to call this "deactivation")
      Q_matrix[uplvl][nlevels_macro] += Qcont_kpkt = (coll_rate * ne) * (xconfig[uplvl].ex - xconfig[target_level].ex); //energy of lower state

      //deactivation back to r-packet
      R_matrix[uplvl][uplvl] += Rcont = rad_rate * (xconfig[uplvl].ex - xconfig[target_level].ex);      //energy difference

      Q_norm[uplvl] += Qcont + Qcont_kpkt + Rcont;
    }

    /* Auger ionization */
    if (xconfig[uplvl].iauger >= 0)
    {
      auger_ptr = &auger_macro[iauger];

      for (n = 0; n < nauger; n++)
      {
        target_level = auger_ptr->nconfig_target[n];
        auger_rate = auger_ptr->Avalue_auger * auger_ptr->branching_ratio[n];

        //internal jump to another macro atom level
        Q_matrix[uplvl][target_level] += Qcont = auger_rate * xconfig[target_level].ex; //energy of lower state

        //jump to the k-packet pool (we used to call this "deactivation")
        Q_matrix[uplvl][nlevels_macro] += Qcont_kpkt = auger_rate * (xconfig[uplvl].ex - xconfig[target_level].ex);     //energy of lower state

        //deactivation back to r-packet isn't possible for the Auger process
        Q_norm[uplvl] += Qcont + Qcont_kpkt;
      }
    }

    /* bound-free */
    for (n = 0; n < nbfd; n++)
    {

      cont_ptr = &phot_top[xconfig[uplvl].bfd_jump[n]]; //pointer to continuum

      sp_rec_rate = mplasma->recomb_sp[xconfig[uplvl].bfd_indx_first + n];      //need this twice so store it
      bf_cont = (sp_rec_rate + q_recomb (cont_ptr, t_e) * ne) * ne;

      target_level = phot_top[xconfig[uplvl].bfd_jump[n]].nlev;

      if (bf_cont > 0.0)
      {

        //internal jump to another macro atom level
        Q_matrix[uplvl][target_level] += Qcont = bf_cont * xconfig[target_level].ex;    //energy of lower state

        //jump to the k-packet pool (we used to call this "deactivation")
        Q_matrix[uplvl][nlevels_macro] += Qcont_kpkt = q_recomb (cont_ptr, t_e) * ne * ne * (xconfig[uplvl].ex - xconfig[target_level].ex);     //energy difference

        //deactivation back to r-packet
        R_matrix[uplvl][uplvl] += Rcont = ne * sp_rec_rate * (xconfig[uplvl].ex - xconfig[target_level].ex);    //energy difference

        Q_norm[uplvl] += Qcont + Qcont_kpkt + Rcont;
      }
    }

    /* Now upwards jumps. */

    /* bound-bound */
    for (n = 0; n < nbbu; n++)
    {
      line_ptr = &line[xconfig[uplvl].bbu_jump[n]];
      rad_rate = (b12 (line_ptr) * mplasma->jbar_old[xconfig[uplvl].bbu_indx_first + n]);

      coll_rate = q12 (line_ptr, t_e);  // this is multiplied by ne below

      target_level = line[xconfig[uplvl].bbu_jump[n]].nconfigu;
      Q_matrix[uplvl][target_level] += Qcont = ((rad_rate) + (coll_rate * ne)) * xconfig[uplvl].ex;     //energy of lower state

      Q_norm[uplvl] += Qcont;
    }

    /* bound-free */
    for (n = 0; n < nbfu; n++)
    {
      /* For bf ionization the jump probability is just gamma * energy
         gamma is the photoionisation rate. Stimulated recombination also included. */
      cont_ptr = &phot_top[xconfig[uplvl].bfu_jump[n]]; //pointer to continuum

      /* first let us take care of the situation where the lower level is zero or close to zero */
      lower_density = den_config (xplasma, cont_ptr->nlev);
      if (lower_density >= DENSITY_PHOT_MIN)
      {
        density_ratio = den_config (xplasma, cont_ptr->uplev) / lower_density;
      }
      else
        density_ratio = 0.0;

      target_level = phot_top[xconfig[uplvl].bfu_jump[n]].uplev;
      Qcont = (mplasma->gamma_old[xconfig[uplvl].bfu_indx_first + n] - (mplasma->alpha_st_old[xconfig[uplvl].bfu_indx_first + n] * xplasma->ne * density_ratio) + (q_ioniz (cont_ptr, t_e) * ne)) * xconfig[uplvl].ex;  //energy of lower state

      /* this error condition can happen in unconverged hot cells where T_R >> T_E.
         for the moment we set to 0 and hope spontaneous recombiantion takes care of things */
      /* note that we check and report this in check_stimulated_recomb() in estimators.c once a cycle */
      if (Qcont < 0.)           //test (can be deleted eventually SS)
      {
        //Error ("Negative probability (matom, 6). Abort?\n");
        Qcont = 0.0;

      }
      Q_matrix[uplvl][target_level] += Qcont;
      Q_norm[uplvl] += Qcont;
    }


  }

  /* Now need to do k-packet processes */
  escape_dummy = 0;
  init_dummy_phot (&pdummy);
  fill_kpkt_rates (xplasma, &escape_dummy, &pdummy);

  /* Cooling due to collisional transitions in lines and collision ionization [for macro atoms] constitute internal transitions from the k-packet pool to macro atom states. */
  kpacket_to_rpacket_rate = 0.0;        // keep track of rate for kpacket_to_rpacket channel

  for (i = 0; i < nlines; i++)
  {
    if (line[i].macro_info == 1 && geo.macro_simple == 0)       //line is for a macro atom
    {
      target_level = line[i].nconfigu;
      Q_matrix[nlevels_macro][target_level] += Qcont = mplasma->cooling_bb[i];
      Q_norm[nlevels_macro] += Qcont;
    }
    else
    {
      /* the idea here is that if it is a simple line then it *must* create an r-packet eventually, so this is essentially
         a k->r transition */
      kpacket_to_rpacket_rate += mplasma->cooling_bb[i];
    }
  }

  for (i = 0; i < nphot_total; i++)
  {
    if (phot_top[i].macro_info == 1 && geo.macro_simple == 0)   //part of macro atom
    {
      target_level = phot_top[i].uplev;
      Q_matrix[nlevels_macro][target_level] += Qcont = mplasma->cooling_bf_col[i];
      Q_norm[nlevels_macro] += Qcont;

    }
    else
    {
      /* XXX - ask stuart about this! */
      kpacket_to_rpacket_rate += mplasma->cooling_bf_col[i];
    }
  }

  /* Cooling due to other processes will always contribute to k-packet -> r-packet channel */

  kpacket_to_rpacket_rate += mplasma->cooling_bftot;
  kpacket_to_rpacket_rate += mplasma->cooling_adiabatic;
  kpacket_to_rpacket_rate += mplasma->cooling_ff + mplasma->cooling_ff_lofreq;
  R_matrix[nlevels_macro][nlevels_macro] += Rcont = kpacket_to_rpacket_rate;
  Q_norm[nlevels_macro] += Rcont;

  /* end kpacket */

  /* now in one step, we multiply by the identity matrix and normalise the probabilities
     this means that what is now stored in Q_matrix is no longer Q, but N=(I - Q) using Vogl
     notation. We check that Q_norm is 0, because some states (ground states) can have 0
     jumping probabilities and so zero normalisation too */
  for (uplvl = 0; uplvl < nrows; uplvl++)
  {
    for (target_level = 0; target_level < nrows; target_level++)
    {
      if (Q_norm[uplvl] > 0.0)
      {
        Q_matrix[uplvl][target_level] = -1. * Q_matrix[uplvl][target_level] / Q_norm[uplvl];
      }
    }
    Q_matrix[uplvl][uplvl] = 1.0;

    if (Q_norm[uplvl] > 0.0)
    {
      R_matrix[uplvl][uplvl] = R_matrix[uplvl][uplvl] / Q_norm[uplvl];
    }
  }

  /* Check normalisation of the matrix. the diagonals of R, minus N which is now stored in
     Q_matrix, should be zero */
  for (uplvl = 0; uplvl < nrows; uplvl++)
  {
    norm = R_matrix[uplvl][uplvl];
    for (target_level = 0; target_level < nlevels_macro + 1; target_level++)
    {
      norm -= Q_matrix[uplvl][target_level];
    }

    /* throw an error if this normalisation is not zero */
    /* note that the ground state is a special case here (improve error check) */
    if ((fabs (norm) > 1e-14 && uplvl != ion[xconfig[uplvl].nion].first_nlte_level) || sane_check (norm))
      Error ("calc_matom_matrix: matom accelerator matrix has bad normalisation for level %d: %8.4e\n", norm, uplvl);
  }

  /* This next line produces an array of the correct size to hold the rate matrix */
  a_data = (double *) calloc (nrows * nrows, sizeof (double));
  a_inverse = (double *) calloc (nrows * nrows, sizeof (double));

  /* We now copy our rate matrix into the prepared matrix */
  for (mm = 0; mm < nrows; mm++)
  {
    for (nn = 0; nn < nrows; nn++)
    {
      a_data[mm * nrows + nn] = Q_matrix[mm][nn];       /* row-major */
    }
  }

  matrix_error = invert_matrix (a_data, a_inverse, nrows);

  if (matrix_error != EXIT_SUCCESS)
  {
    Error ("calc_matom_matrix: error %d whilst inverting Q_matrix in plasma cell %d\n", matrix_error, xplasma->nplasma);
  }

  free (a_data);

  /* We now copy our rate matrix into the prepared matrix */
  for (mm = 0; mm < nrows; mm++)
  {
    for (nn = 0; nn < nrows; nn++)
    {
      /* in Christian Vogl's notation this is doing his equation 3: B = (N * R)
         where N is the inverse matrix we have just calculated. */
      /* the reason this matrix multiplication is so simple here is because R is a diagonal matrix */
      matom_matrix[mm][nn] = a_inverse[mm * nrows + nn] * R_matrix[nn][nn];
    }
  }

  free (a_inverse);

  /* need to free each calloc-ed row of the matrixes */
  for (i = 0; i < nrows; i++)
  {
    free (R_matrix[i]);
    free (Q_matrix[i]);
  }

  free (R_matrix);
  free (Q_matrix);
  free (Q_norm);
}


/**********************************************************/
/**
 * @brief calculate the cooling rates for the conversion of k-packets.
 *
 * @param [in] PlasmaPtr  xplasma
 * @param [in,out] double **matom_matrix
 *        the 2D matrix array we will populate with normalised probabilities
 *
 *
 * @details
 *
 *
 **********************************************************/

int
fill_kpkt_rates (xplasma, escape, p)
     PlasmaPtr xplasma;
     int *escape;
     PhotPtr p;
{

  int i;
  int ulvl;
  double cooling_bf[nphot_total];
  double cooling_bf_col[nphot_total];   //collisional cooling in bf transitions
  double cooling_bb[NLINES];
  double cooling_adiabatic;
  struct topbase_phot *cont_ptr;
  struct lines *line_ptr;
  double cooling_normalisation;
  double electron_temperature;
  double cooling_bbtot, cooling_bftot, cooling_bf_coltot;
  double lower_density, upper_density;
  double cooling_ff;
  WindPtr one;

  MacroPtr mplasma;

  double coll_rate, rad_rate;
  double freqmin, freqmax;

  mplasma = &macromain[xplasma->nplasma];
  one = &wmain[xplasma->nwind];

  electron_temperature = xplasma->t_e;

  /* JM 1511 -- Fix for issue 187. We need band limits for free free packet
     generation (see call to one_ff below) */
  freqmin = xband.f1[0];
  freqmax = ALPHA_FF * xplasma->t_e / H_OVER_K;

  /* ksl This is a Bandaid for when the temperatures are very low */
  /* in this case cooling_ff should be low compared to cooling_ff_lofreq anyway */
  if (freqmax < 1.1 * freqmin)
  {
    freqmax = 1.1 * freqmin;
  }

  /* If the kpkt destruction rates for this cell are not known they are calculated here.  This happens
   * every time the wind is updated */

  if (mplasma->kpkt_rates_known != TRUE)
  {
    cooling_normalisation = 0.0;
    cooling_bftot = 0.0;
    cooling_bbtot = 0.0;
    cooling_ff = 0.0;
    cooling_bf_coltot = 0.0;
    mplasma->cooling_bb_simple_tot = 0.0;

    /* Start of BF calculation */
    /* JM 1503 -- we used to loop over ntop_phot here,
       but we should really loop over the tabulated Verner Xsections too
       see #86, #141 */
    for (i = 0; i < nphot_total; i++)
    {
      cont_ptr = &phot_top[i];
      ulvl = cont_ptr->uplev;

      if (cont_ptr->macro_info == TRUE && geo.macro_simple == FALSE)
      {
        upper_density = den_config (xplasma, ulvl);
        cooling_bf[i] = mplasma->cooling_bf[i] =
          upper_density * PLANCK * cont_ptr->freq[0] * (mplasma->recomb_sp_e[xconfig[ulvl].bfd_indx_first + cont_ptr->down_index]);
      }
      else
      {
        upper_density = xplasma->density[cont_ptr->nion + 1];

        cooling_bf[i] = mplasma->cooling_bf[i] = upper_density * PLANCK * cont_ptr->freq[0] * (xplasma->recomb_simple[i]);
      }

      if (cooling_bf[i] < 0)
      {
        Error ("kpkt: phot %d bf cooling rate negative. Density was %g\n", p->np, upper_density);
        Error ("alpha_sp(cont_ptr, xplasma,2) %g \n", alpha_sp (cont_ptr, xplasma, 2));
        Error ("i, ulvl, nphot_total, nion %d %d %d %d\n", i, ulvl, nphot_total, cont_ptr->nion);
        Error ("nlev, z, istate %d %d %d \n", cont_ptr->nlev, cont_ptr->z, cont_ptr->istate);
        Error ("freq[0] %g\n", cont_ptr->freq[0]);
        cooling_bf[i] = mplasma->cooling_bf[i] = 0.0;
      }
      else
      {
        cooling_bftot += cooling_bf[i];
      }

      cooling_normalisation += cooling_bf[i];

      if (cont_ptr->macro_info == TRUE && geo.macro_simple == FALSE)
      {
        /* Include collisional ionization as a cooling term in macro atoms, but not simple atoms.  */

        lower_density = den_config (xplasma, cont_ptr->nlev);
        cooling_bf_col[i] = mplasma->cooling_bf_col[i] =
          lower_density * PLANCK * cont_ptr->freq[0] * q_ioniz (cont_ptr, electron_temperature);

        cooling_bf_coltot += cooling_bf_col[i];

        cooling_normalisation += cooling_bf_col[i];

      }

    }

    /* End of BF calculation and beginning of BB calculation.  Note that for macro atoms
       the upper level density is stored but for simple atoms it must be calculated. */

    for (i = 0; i < nlines; i++)
    {
      line_ptr = &line[i];
      if (line_ptr->macro_info == TRUE && geo.macro_simple == FALSE)
      {
        cooling_bb[i] = mplasma->cooling_bb[i] =
          den_config (xplasma, line_ptr->nconfigl) * q12 (line_ptr, electron_temperature) * line_ptr->freq * PLANCK;

      }
      else
      {
        two_level_atom (line_ptr, xplasma, &lower_density, &upper_density);

        coll_rate = q21 (line_ptr, electron_temperature) * (1. - exp (-H_OVER_K * line_ptr->freq / electron_temperature));

        cooling_bb[i] =
          (lower_density * line_ptr->gu / line_ptr->gl -
           upper_density) * coll_rate / (exp (H_OVER_K * line_ptr->freq / electron_temperature) - 1.) * line_ptr->freq * PLANCK;

        rad_rate = a21 (line_ptr) * p_escape (line_ptr, xplasma);

        /* Now multiply by the scattering probability - i.e. we are only going to consider bb cooling when
           the photon actually escapes - we don't to waste time by exciting a two-level macro atom only so that
           it makes another k-packet for us! (SS May 04) */

        cooling_bb[i] *= rad_rate / (rad_rate + (coll_rate * xplasma->ne));
        mplasma->cooling_bb[i] = cooling_bb[i];
        mplasma->cooling_bb_simple_tot += cooling_bb[i];
      }

      if (cooling_bb[i] < 0)
      {
        cooling_bb[i] = mplasma->cooling_bb[i] = 0.0;
      }
      else
      {
        cooling_bbtot += cooling_bb[i];
      }
      cooling_normalisation += cooling_bb[i];
    }

    /* end of BB calculation  */


    if (one->inwind >= 0)
    {
      cooling_ff = mplasma->cooling_ff = total_free (xplasma, xplasma->t_e, freqmin, freqmax) / xplasma->vol / xplasma->ne;
      cooling_ff += mplasma->cooling_ff_lofreq = total_free (xplasma, xplasma->t_e, 0.0, freqmin) / xplasma->vol / xplasma->ne;
    }
    else
    {
      /* This should never happen */
      cooling_ff = mplasma->cooling_ff = mplasma->cooling_ff_lofreq = 0.0;
      Error ("kpkt: np %d A scattering event in cell %d with vol = 0???\n", p->np, one->nwind);
      *escape = TRUE;
      p->istat = P_ERROR_MATOM;
      return (0);
    }


    if (cooling_ff < 0)
    {
      Error ("kpkt: np %d ff cooling rate negative. Abort.\n", p->np);
      *escape = TRUE;
      p->istat = P_ERROR_MATOM;
      return (0);
    }
    else
    {
      cooling_normalisation += cooling_ff;
    }


    /* JM -- 1310 -- we now want to add adiabatic cooling as another way of destroying kpkts
       this should have already been calculated and stored in the plasma structure. Note that
       adiabatic cooling does not depend on type of macro atom excited */

    /* note the units here- we divide the total luminosity of the cell by volume and ne to give cooling rate */

    cooling_adiabatic = xplasma->cool_adiabatic / xplasma->vol / xplasma->ne;   // JM 1411 - changed to use filled volume


    if (geo.adiabatic == 0 && cooling_adiabatic > 0.0)
    {
      Error ("kpkt: Adiabatic cooling turned off, but non zero in cell %d", xplasma->nplasma);
    }


    /* If the next error occurs see issue #70.  It is ghought to be fixed, with the possiple
       exception of cells that are partially in the wind.
     */
    if (cooling_adiabatic < 0)
    {
      Error ("kpkt: Photon %d. Adiabatic cooling negative! Major problem if inwind (%d) == 0\n", p->np, one->inwind);
      Error ("kpkt: Photon %d. Setting adiabatic kpkt destruction probability to zero for this matom.\n", p->np);
      cooling_adiabatic = 0.0;
    }


    cooling_normalisation += cooling_adiabatic;

    mplasma->cooling_bbtot = cooling_bbtot;
    mplasma->cooling_bftot = cooling_bftot;
    mplasma->cooling_bf_coltot = cooling_bf_coltot;
    mplasma->cooling_adiabatic = cooling_adiabatic;
    mplasma->cooling_normalisation = cooling_normalisation;
    mplasma->kpkt_rates_known = TRUE;

  }

  return (0);
}




/**********************************************************/
/**
 * @brief a routine which calculates what fraction of a level
 *         emissivity comes out in a given frequency range
 *
 * @param [in]      PlasmaPtr   xplasma       Plasma cell in question
 * @param [in]      int         upper         macro-atom level
 * @param [in]      double      freq_min, freq    Frequency range requested
 *                                            (e.g. spectral cycle freq range)
 *
 * @details similar to routines like emit_matom, except that we calculate the fraction
 * of emission in a band rather than calculating frequencies for photons. Used to be
 * done via a less efficient rejection method.
 *
***********************************************************/

double
f_matom_emit_accelerate (xplasma, upper, freq_min, freq_max)
     PlasmaPtr xplasma;
     int upper;
     double freq_min, freq_max;
{
  struct lines *line_ptr;
  struct topbase_phot *cont_ptr;
  int uplvl;
  double eprbs[NBBJUMPS + NBFJUMPS], eprbs_band[NBBJUMPS + NBFJUMPS];
  double penorm, penorm_band;
  double sp_rec_rate;
  int n, m;
  int nbbd, nbfd;
//OLD  double t_e, ne;
  double ne;
  double bb_cont;
  double flast, fthresh, bf_int_full, bf_int_inrange;

//OLD  t_e = xplasma->t_e;           //electron temperature
  ne = xplasma->ne;             //electron number density

  /* The first step is to identify the configuration that has been excited. */

  uplvl = upper;

  /* Unlike the proper matom routine, we don't want any jumps, only deactivations here. */


  /*  The excited configuration is now known. Now compute all the probabilities of deactivation
     /jumping from this configuration. Then choose one. */


  nbbd = xconfig[uplvl].n_bbd_jump;     //store these for easy access -- number of bb downward jumps
  nbfd = xconfig[uplvl].n_bfd_jump;     // number of bf downared jumps from this transition


  // Start by setting everything to 0

  m = 0;                        //m counts the total number of possible ways to leave the level
  penorm = 0.0;                 //stores the total emission probability
  penorm_band = 0.0;

  for (n = 0; n < nbbd + nbfd; n++)
  {
    eprbs[n] = 0;               //stores the individual emission probabilities SS
    eprbs_band[n] = 0;
  }
  /* Finished zeroing. */


  /* bb downwards jumps */

  for (n = 0; n < nbbd; n++)
  {
    line_ptr = &line[xconfig[uplvl].bbd_jump[n]];
    /* Since we are only interested in making an r-packet here we can (a) ignore collisional
       deactivation and (b) ignore lines outside the frequency range of interest. */
    bb_cont = (a21 (line_ptr) * p_escape (line_ptr, xplasma));

    eprbs[n] = bb_cont * (xconfig[uplvl].ex - xconfig[line[xconfig[uplvl].bbd_jump[n]].nconfigl].ex);   //energy difference

    if (eprbs[n] < 0.)          //test (can be deleted eventually SS)
    {
      //Error ("Negative probability (matom, 2) for number %d at %e. Abort.\n", n, eprbs[n]);
      //Exit (0);
      Error ("Negative probability (matom, 2) for number %d at %e. Setting to 0.\n", n, eprbs[n]);
      eprbs[n] = 0.0;
    }

    penorm += eprbs[n];
    if ((line_ptr->freq > freq_min) && (line_ptr->freq < freq_max))     // correct range
    {
      eprbs_band[m] = eprbs[n];
      penorm_band += eprbs_band[m];
      m++;
    }
  }

  /* bf downwards jumps */
  for (n = 0; n < nbfd; n++)
  {
    cont_ptr = &phot_top[xconfig[uplvl].bfd_jump[n]];   //pointer to continuum

    /* If the edge is above the frequency range we are interested in then we need not consider this
       bf process. */
    sp_rec_rate = alpha_sp (cont_ptr, xplasma, 0);
    eprbs[n + nbbd] = sp_rec_rate * ne * (xconfig[uplvl].ex - xconfig[phot_top[xconfig[uplvl].bfd_jump[n]].nlev].ex);   //energy difference

    if (eprbs[n + nbbd] < 0.)   //test (can be deleted eventually SS)
    {
      //Error ("Negative probability (matom, 4) for %d at %e. Abort.\n", n, eprbs[n + nbbd]);
      //Exit (0);
      Error ("Negative probability (matom, 4) for %d at %e. Setting to 0\n", n, eprbs[n + nbbd]);
      eprbs[n + nbbd] = 0.0;

    }
    penorm += eprbs[n + nbbd];
    if (cont_ptr->freq[0] < freq_max && cont_ptr->freq[cont_ptr->np - 1] > freq_min)    //means that it may contribute
    {
      eprbs_band[m] = eprbs[n + nbbd];
      fthresh = cont_ptr->freq[0];      //first frequency in list
      flast = cont_ptr->freq[cont_ptr->np - 1]; //last frequency in list
      bf_int_full = scaled_alpha_sp_integral_band_limited (cont_ptr, xplasma, 0, fthresh, flast);

      /* the limits differ depending on whether the band covers all or part of the cross-section.
         four possibilities */
      if (fthresh < freq_min && flast > freq_max)
      {
        bf_int_inrange = scaled_alpha_sp_integral_band_limited (cont_ptr, xplasma, 0, freq_min, freq_max);
      }
      else if (fthresh < freq_min && flast < freq_max)
      {
        bf_int_inrange = scaled_alpha_sp_integral_band_limited (cont_ptr, xplasma, 0, freq_min, flast);
      }
      else if (fthresh > freq_min && flast > freq_max)
      {
        bf_int_inrange = scaled_alpha_sp_integral_band_limited (cont_ptr, xplasma, 0, fthresh, freq_max);
      }
      else if (fthresh > freq_min && flast < freq_max)
      {
        bf_int_inrange = scaled_alpha_sp_integral_band_limited (cont_ptr, xplasma, 0, fthresh, flast);
      }
      else
      {
        Error ("Something wrong here: f_matom_emit_accelerate broke the law!");
        bf_int_inrange = 0;
      }
      penorm_band += eprbs_band[m] * bf_int_inrange / bf_int_full;
      m++;
    }
  }

  /* catch the case where penorm == 0 for e.g. ground states with zero emission probability */
  if (penorm > 0)
  {
    return (penorm_band / penorm);
  }
  else
  {
    return (0.0);
  }
}

/* The frequency and the value of nres have been set correctly. All done. */






/**********************************************************/
/**
 * @brief calculate what fraction of the thermal continuum emission comes out in the required band
 *
 * This routine uses the various cooling rates to work out which k->r processes contribute
 * in the frequency range requested (between freq_min and freq_max). This involves doing a series
 * of band-limited integrals, for the bound-free "alpha" estimators, for each cross-section.
 * It also calls functions like total_free() to work out the fraction of free-free emission
 * that emerges in the frequency range.
 *
 *
 *
 * @param [in]      PlasmaPtr   xplasma       Plasma cell in question
 *
 * @param [in]      double      freq_min, freq    Frequency range requested
 *                                            (e.g. spectral cycle freq range)
 *
 * @return penorm_band / penorm   total fraction of k->r emission that emerges in the band
************************************************************/

double
f_kpkt_emit_accelerate (xplasma, freq_min, freq_max)
     PlasmaPtr xplasma;
     double freq_min, freq_max;
{

  int i;
  int escape_dummy;
  struct topbase_phot *cont_ptr;
  MacroPtr mplasma;
  struct photon pdummy;
  double ff_freq_min, ff_freq_max;
  double eprbs, eprbs_band, penorm, penorm_band;
  double flast, fthresh, bf_int_full, bf_int_inrange;
  double total_ff_lofreq, total_ff;

  penorm = 0.0;
  penorm_band = 0.0;

  if (check_plasma (xplasma, "f_kpkt_emit_accelerate"))
  {
    Error ("f_kpkt_emit_accelrate: Calculating emissivity for the dummy plasma cell\n");
    return (0);
  }
  mplasma = &macromain[xplasma->nplasma];

//OLD  electron_temperature = xplasma->t_e;

  /* JM 1511 -- Fix for issue 187. We need band limits for free free packet
     generation (see call to one_ff below) */
  ff_freq_min = xband.f1[0];
  ff_freq_max = ALPHA_FF * xplasma->t_e / H_OVER_K;

  /* ksl This is a Bandaid for when the temperatures are very low */
  /* in this case cooling_ff should be low compared to cooling_ff_lofreq anyway */
  if (ff_freq_max < 1.1 * ff_freq_min)
  {
    ff_freq_max = 1.1 * ff_freq_min;
  }


  /* Now need to do k-packet processes. This subroutine calculates the k-packet
     cooling rates and stores them in mplasma->cooling. Dummy variables are needed
     because this routine is also used in the main kpkt routine */
  escape_dummy = 0;
  init_dummy_phot (&pdummy);
  fill_kpkt_rates (xplasma, &escape_dummy, &pdummy);

  for (i = 0; i < nphot_total; i++)
  {
    cont_ptr = &phot_top[i];    //pointer to continuum

    /* If the edge is above the frequency range we are interested in then we need not consider this
       bf process. */

    eprbs = mplasma->cooling_bf[i];
    penorm += eprbs;
    if (cont_ptr->freq[0] < freq_max && cont_ptr->freq[cont_ptr->np - 1] > freq_min)    //means that it may contribute
    {
      eprbs_band = mplasma->cooling_bf[i];
      fthresh = cont_ptr->freq[0];      //first frequency in list
      flast = cont_ptr->freq[cont_ptr->np - 1]; //last frequency in list
      bf_int_full = scaled_alpha_sp_integral_band_limited (cont_ptr, xplasma, 0, fthresh, flast);
      if (fthresh < freq_min && flast > freq_max)
      {
        bf_int_inrange = scaled_alpha_sp_integral_band_limited (cont_ptr, xplasma, 0, freq_min, freq_max);
      }
      else if (fthresh < freq_min && flast < freq_max)
      {
        bf_int_inrange = scaled_alpha_sp_integral_band_limited (cont_ptr, xplasma, 0, freq_min, flast);
      }
      else if (fthresh > freq_min && flast > freq_max)
      {
        bf_int_inrange = scaled_alpha_sp_integral_band_limited (cont_ptr, xplasma, 0, fthresh, freq_max);
      }
      else if (fthresh > freq_min && flast < freq_max)
      {
        bf_int_inrange = scaled_alpha_sp_integral_band_limited (cont_ptr, xplasma, 0, fthresh, flast);
      }
      else
      {
        Error ("Something wrong here: f_matom_emit_accelerate broke the law!");
        bf_int_inrange = 0;
      }
      penorm_band += eprbs_band * bf_int_inrange / bf_int_full;
    }
  }

  for (i = 0; i < nlines; i++)
  {
    if (line[i].macro_info == 1 && geo.macro_simple == 0)       //line is for a macro atom
    {
      eprbs = 0.0;              //these are not deactivations in this approach any more but jumps to macro atom levels
    }
    else                        //line is not for a macro atom - use simple method
    {
      penorm += eprbs = mplasma->cooling_bb[i];
      if ((line[i].freq > freq_min) && (line[i].freq < freq_max))       // correct range
      {
        penorm_band += eprbs_band = eprbs;
      }
    }
  }



  /* consult issues #187, #492 regarding free-free */
  penorm += eprbs = mplasma->cooling_ff + mplasma->cooling_ff_lofreq;

  total_ff_lofreq = total_free (xplasma, xplasma->t_e, 0, ff_freq_min);
  total_ff = total_free (xplasma, xplasma->t_e, ff_freq_min, ff_freq_max);

  /*
   * Do not increment penorm_band when the total free-free luminosity is zero
   */

  if (freq_min > ff_freq_min)
  {
    if (total_ff > 0)
      penorm_band += total_free (xplasma, xplasma->t_e, freq_min, freq_max) / total_ff * mplasma->cooling_ff;
  }
  else if (freq_max > ff_freq_min)
  {
    if (total_ff > 0)
      penorm_band += total_free (xplasma, xplasma->t_e, ff_freq_min, freq_max) / total_ff * mplasma->cooling_ff;
    if (total_ff_lofreq > 0)
      penorm_band += total_free (xplasma, xplasma->t_e, freq_min, ff_freq_min) / total_ff_lofreq * mplasma->cooling_ff_lofreq;
  }
  else
  {
    if (total_ff_lofreq > 0)
      penorm_band += total_free (xplasma, xplasma->t_e, freq_min, freq_max) / total_ff_lofreq * mplasma->cooling_ff_lofreq;
  }

  penorm += eprbs = mplasma->cooling_adiabatic;

  for (i = 0; i < nphot_total; i++)
  {
    if (phot_top[i].macro_info == 0 || geo.macro_simple == 1)
    {
      penorm += eprbs = mplasma->cooling_bf_col[i];
    }
  }

  if (penorm > 0)
  {
    return (penorm_band / penorm);
  }
  else
  {
    return (0.0);
  }

}

/**********************************************************/
/**
 * @brief  Choose a deactivation process using the matrix scheme
 *         for macro-atom transition probabilities
 *
 * @param[in] PlasmaPtr xplasma        The plasma cell in question
 * @param[in] int       uplvl          The level the macro-atom was activated with
 *
 * @return  int   j  the level the macro-atom will deactivate from
 *
 * @details
 *
 **********************************************************/

int
matom_deactivation_from_matrix (xplasma, uplvl)
     PlasmaPtr xplasma;
     int uplvl;
{
  double z, total;
  int j, i;
  int nrows = nlevels_macro + 1;
  double **matom_matrix;
  MacroPtr mplasma;

  mplasma = &macromain[xplasma->nplasma];

  if (mplasma->matrix_rates_known == FALSE)
  {
    if (mplasma->store_matom_matrix == FALSE)
    {
      /* we aren't storing the macro-atom matrix, so we need to allocate and calculate it */
      matom_matrix = (double **) calloc (sizeof (double *), nrows);
      for (i = 0; i < nrows; i++)
      {
        matom_matrix[i] = (double *) calloc (sizeof (double), nrows);
      }
    }
    else
    {
      matom_matrix = mplasma->matom_matrix;
    }

    calc_matom_matrix (xplasma, matom_matrix);

    /* if we are storing the matrix, flag that we know the rates now */
    if (mplasma->store_matom_matrix == TRUE)
    {
      mplasma->matrix_rates_known = TRUE;
    }
  }
//OLD  else if (mplasma->store_matom_matrix == TRUE)
  else
  {
    matom_matrix = mplasma->matom_matrix;
  }

  /* Now use the B matrix to calculate the outgoing state from activating state "uplvl" */
  /* we draw a random number and sample from the column in the matrix corresponding to uplvl */
  z = random_number (0.0, 1.0);
  j = 0;
  total = 0.0;
  while (total < z)
  {
    total += matom_matrix[uplvl][j];
    j++;
  }

  /* This if statement is added to prevent case where z is essentially 0. */
  if (j > 0)
  {
    j = j - 1;
  }

  if (mplasma->store_matom_matrix == FALSE)
  {
    /* need to free each calloc-ed row of the matrixes */
    for (i = 0; i < nrows; i++)
    {
      free (matom_matrix[i]);
    }
    free (matom_matrix);
  }

  return (j);
}

/**********************************************************/
/**
 * @brief calculate all the macro-atom matrices in advance
 * @return  0
 *
 * @details calculate all the macro-atom B matrices in advance
 * of the ionization cycles, and communicate them between
 * parallel threads if necessary. Populates matom_matrix in
 * macromain.
 **********************************************************/

int
calc_all_matom_matrices (void)
{
  int ndo, my_nmin, my_nmax, n;
  struct timeval timer_t0;
  char message[LINELENGTH];
  MacroPtr mplasma;
  PlasmaPtr xplasma;
#ifdef MPI_ON
  ndo = get_parallel_nrange (rank_global, NPLASMA, np_mpi_global, &my_nmin, &my_nmax);
#else
  my_nmin = 0;
  my_nmax = NPLASMA;
  ndo = NPLASMA;
#endif

  timer_t0 = init_timer_t0 ();

  for (n = my_nmin; n < my_nmax; n++)
  {
    xplasma = &plasmamain[n];
    mplasma = &macromain[n];

    if (mplasma->store_matom_matrix == TRUE)
    {
      calc_matom_matrix (xplasma, mplasma->matom_matrix);
    }
  }

  /* print the time taken for this thread to complete */
  sprintf (message, "calc_all_matom_matrices: thread %d calculated %d matrices in", rank_global, ndo);
  print_timer_duration (message, timer_t0);

  /* this deals with communicating the matrices between threads (does nothing in serial mode) */
  broadcast_macro_atom_state_matrix (my_nmin, my_nmax, ndo);

  /* flag the matrix rates as known */
  for (n = 0; n < NPLASMA; n++)
  {
    macromain[n].matrix_rates_known = TRUE;
  }

  return (0);
}
