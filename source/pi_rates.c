
/***********************************************************/
/** @file  pi_rates.c
 * @author nsh
 * @date   Aug, 2014
 *
 * @brief  Routines associated with calculating photionization rates
 *
 * ???
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

struct topbase_phot *xtop;      //Topbase description of a photoionization x-section - this is the only type we use - tabulated.

double qromb_temp;              //Temperature used in integrations - has to be an external variable so qromb can use it

//Model parameters for integrations, also passed externally so qromb can use them

double xpl_alpha, xpl_w, xpl_logw;
double xexp_temp, xexp_w;


/**********************************************************/
/**
 * @brief      calculate the photoionization rate coefficient for an ion in a given cell
 *
 * @param [in] int  nion   The ion being ionized (or the index to an inner shell cross section)
 * @param [in] PlasmaPtr  xplasma   The cell for which rates are calulated
 * @param [in] int  mode Use either a piecewise model (1), a blackbody model (2), or the fine cell histogram (3) for J_nu
 * @param [in] type Compute either an outer shell (1) or inner shell (2) rate
 * @return     The photioinization rate coefficient for the ion.
 *
 * @details
 *  calc_pi_rate  calculates the photoionization rate coefficient for ion
 * 	nion, based upon the mean intensity stored in cell xplasma
 * 	The mode tells the subroutine wether we are modelling the
 * 	mean intensity as a dilute blackbody (mode2) or as a series
 * 	of power laws and exponentials (mode1). The type allows for the calculation
 *      of inner shell rates - in this case nion is the index into an inner
 * 	shell cross section record - this links to the relevant ion internally.
 *      the reason for the difference is that there is exactly one outer rate
 * 	per ion, but there can be many inner shell rates. Most of the information
 * 	needed for the calculations is stored in the xplasma structure (e.g. temperature
 * 	and spectral model)
 *
 * ### Notes ###
 * This was created in Summer 2014 in preparation for matrix ionization solver. Previously, this code
 * was contained in two subroutines bb_correct_2 and pl_correct_2.The functionality of these two
 * have been combined into one - hence the requirement for the mode parameter. It was further extended
 * to deal with inner shell rates - hence the type parameter
 *
 **********************************************************/

double
calc_pi_rate (int nion, PlasmaPtr xplasma, int mode, int type)
{
  int j;
  double pi_rate;
  int ntmin, nvmin;
  double fthresh, fmax, fmaxtemp;
  double f1, f2;
  double exp_qromb, pl_qromb;


  exp_qromb = 1e-4;             /*These are the two tolerance values for the rhomberg integrals */
  pl_qromb = 1e-4;

  ntmin = nvmin = -1;           /* Initialize these to an unreasonable number. We dont use them all the time */

  if (type == 1)                //We are computing a normal outer shell rate
  {
    if (-1 < nion && nion < nions)      //Get cross section for this specific ion_number
    {
      ntmin = ion[nion].ntop_ground;    /*We only ever use the ground state cross sections. This is for topbase */
      nvmin = ion[nion].nxphot; /* this gets any verner cross section */
    }
    else
    {
      Error ("calc_pi_rate: %d is unacceptable value of nion\n", nion);
      Exit (0);
      return (1.0);
    }
    if (ion[nion].phot_info > 0)        //topbase or hybrid VFKY (GS)+TB excited
    {
      xtop = &phot_top[ntmin];
    }
    else if (ion[nion].phot_info == 0)  // verner
    {
      xtop = &phot_top[nvmin];  //If we don't have a topbase or hybrid cross section we fall back to VFKY
    }
    else
    {
      Error ("calc_pi_rate: No photoionization xsection for ion %d (element %d, ion state %d)\n", nion, ion[nion].z, ion[nion].istate);
      Exit (0);
    }
  }
  else if (type == 2)           //We are computing an inner shell rate - nion refers to an actual cross section.
  {
    if (-1 < nion && nion < n_inner_tot)        //We have a reasonable value for nion_in
    {
      xtop = &inner_cross[nion];        //Set the cross sections for the integral to the inner shell cross section
    }
    else
    {
      Error ("calc_pi_rate: No inner shell xsection for record %d (there are %d)\n", nion, n_inner_tot);
      Exit (0);
    }
  }
  else
  {
    Error ("calc_pi_rate: unknown mode %i\n", type);
  }



/* At this stage, xtop points to either an outer or inner shell photoionization cross section
	We now prepare to explicitly integrate J x csec over modelled bands */

  fthresh = xtop->freq[0];      //The first frequency for which we have a cross section
  fmax = xtop->freq[xtop->np - 1];      //The last frequency
  pi_rate = 0;                  //Initialise the pi rate - it will be computed through several integrals



  if (mode == 1)                //Modelled version of J
  {
    for (j = 0; j < xplasma->state.nbands; j++) //We loop over all the bands in this cell
    {
      xpl_alpha = xplasma->state.pl_alpha[j];   //set the various model parameters to those for this model
      xpl_logw = xplasma->state.pl_log_w[j];
      xexp_temp = xplasma->state.exp_temp[j];
      xexp_w = xplasma->state.exp_w[j];
      if (xplasma->state.spec_mod_type[j] != SPEC_MOD_FAIL)     //Only bother doing the integrals if we have a model in this band
      {
        f1 = xplasma->state.fmin_mod[j];        //NSH 131114 - Set the low frequency limit to the lowest frequency that the model applies to
        f2 = xplasma->state.fmax_mod[j];        //NSH 131114 - Set the high frequency limit to the highest frequency that the model applies to
        if (f1 < fthresh && fthresh < f2 && f1 < fmax && fmax < f2)     //Case 1-
        {
          if (xplasma->state.spec_mod_type[j] == SPEC_MOD_PL)
          {
            pi_rate += num_int (tb_logpow, fthresh, fmax, pl_qromb);

          }
          else
          {
            pi_rate += num_int (tb_exp, fthresh, fmax, exp_qromb);
          }
        }
        else if (f1 < fthresh && fthresh < f2 && f2 < fmax)     //case 2
        {
          if (xplasma->state.spec_mod_type[j] == SPEC_MOD_PL)
          {
            pi_rate += num_int (tb_logpow, fthresh, f2, pl_qromb);
          }
          else
          {
            pi_rate += num_int (tb_exp, fthresh, f2, exp_qromb);
          }
        }
        else if (f1 >= fthresh && f1 < fmax && fmax < f2)       //case 3 (f1 == fthresh when a band edge is the threshold)
        {
          if (xplasma->state.spec_mod_type[j] == SPEC_MOD_PL)
          {
            pi_rate += num_int (tb_logpow, f1, fmax, pl_qromb);
          }
          else
          {
            pi_rate += num_int (tb_exp, f1, fmax, exp_qromb);
          }
        }
        else if (f1 >= fthresh && f2 < fmax)    // case 4
        {
          if (xplasma->state.spec_mod_type[j] == SPEC_MOD_PL)
          {
            pi_rate += num_int (tb_logpow, f1, f2, pl_qromb);
          }
          else
          {
            pi_rate += num_int (tb_exp, f1, f2, exp_qromb);
          }
        }
        else                    //case 5 - should only be the case where the band is outside the range for the integral.
        {
          pi_rate += 0;         // Add nothing - bit of a null statement, but makes the code look nice.
        }
      }                         //End of loop to only integrate in this band if there is power
    }

  }
  else if (mode == 2)           //blackbody mode
  {
    fmaxtemp = xtop->freq[xtop->np - 1];        //Set the maximum frequency temporarily to the maximum cross section frequency
    fmax = check_freq_max (fmaxtemp, xplasma->state.t_r);       /*Check that the requested maximum frequency is sensible - if it is way
                                                                   off the end of the wien tail then the integration can fail - reset if necessary. */
    if (fthresh > fmax)         //The threshold for PI is above the maximum frequency of the radiation
    {
      pi_rate = 0.0;
    }
    else                        //We are OK - do the integral
    {
      qromb_temp = xplasma->state.t_r;
      pi_rate = xplasma->state.w * num_int (tb_planck, fthresh, fmax, 1.e-4);
    }
  }

  else if (mode == 3)           //Fine cell histogram (band model cell_spectrum)
  {
    pi_rate = pi_rate_hist_integral (xplasma, type);
  }

  pi_rate = (4 * PI * pi_rate) / PLANCK;        //We multiply by 4pi and divide by photon energy - the division by nu is done in the integrands




  return (pi_rate);
}


/**********************************************************/
/**
 * @brief      The integrand for working out the PI rate in a BB radiation field
 *
 * @param [in] double  freq   the frequency
 * @param [in] void  params   An extra (unused) variable to make it paletable for the gsl integrator
 * @return     value of sigma B_nu/nu at nu
 *
 * @details
 * tb_planck is the function a\nuB\nu and is called by Qromb in order to integrate over the frequency range where the
 *    ionisation cross section is significant. This is the function for ions with a topbase cross section NSH
 *
 * ### Notes ###
 *
 **********************************************************/

double
tb_planck (double freq, void *params)
{
  double answer, bbe;
  bbe = exp ((PLANCK * freq) / (BOLTZMANN * qromb_temp));
  answer = (2. * PLANCK * pow (freq, 3.)) / (pow (VLIGHT, 2));
  answer *= (1 / (bbe - 1));
//      answer*=weight;
  answer *= sigma_phot (xtop, freq);

  answer /= freq;

  return (answer);
}


/**********************************************************/
/**
 * @brief      The integrand for working out the PI rate in a radiation field modelled by a power law
 *
 * @param [in] double  freq   the frequency
 * @param [in] void  params   An extra (unused) variable to make it paletable for the gsl integrator
 * @return     value of sigma J_nu/nu at nu
 *
 * @details
 * This is the integrand for computing the PI rate coefficient when j_nu is described by a power
 * law. The returned value is J_nu x cross section / nu. The power law is computed in log space to
 * allow a larger range of possible values.
 *
 * ### Notes ###
 *
 **********************************************************/

double
tb_logpow (double freq, void *params)
{
  double answer;
  answer = pow (10, xpl_logw + (xpl_alpha - 1.0) * log10 (freq));       //NB - the alpha-1.0 appears because we divide by nu
  answer *= sigma_phot (xtop, freq);    // and finally multiply by the cross section.

  return (answer);
}


/**********************************************************/
/**
 * @brief      The integrand for working out the PI rate in a radiation field modelled by an exponential
 *
 * @param [in] double  freq   the frequency
 * @param [in] void  params   An extra (unused) variable to make it paletable for the gsl integrator
 * @return     value of sigma J_nu/nu at nu
 *
 * @details
 * This is the integrand for computing the PI rate coefficient when j_nu is described by an exponential
 * law. The returned value is J_nu x cross section / nu.
 *
 * ### Notes ###
 *
 **********************************************************/

double
tb_exp (double freq, void *params)
{
  double answer;

  answer = xexp_w * exp ((-1.0 * PLANCK * freq) / (BOLTZMANN * xexp_temp));
  answer *= sigma_phot (xtop, freq);    // and finally multiply by the cross section.

  answer /= freq;               //then finally finally divide by the frequency
  return (answer);
}



/* Cache of the integrals of sigma/nu over each cell spectrum bin, one array per
   cross section, computed when first needed.  The bins do not change during a run. */

static double **hist_sigma_outer = NULL;
static double **hist_sigma_inner = NULL;


/**********************************************************/
/**
 * @brief      The integral of sigma/nu over each cell spectrum bin, for the
 * cross section xtop
 *
 * @param [in] int  type   1 for an outer-shell (phot_top), 2 for an inner-shell (inner_cross) cross section
 * @return     A pointer to the geo.nbins_in_cell_spec integrals
 *
 * @details
 * Each integral runs over the part of the bin between the threshold and the
 * last frequency of the cross section, and is done with Simpson's rule in
 * ln nu on 8 intervals.  The bins have a boundary at every threshold when
 * this is used (cell_spec_threshold_grid), so sigma is smooth within each bin.
 **********************************************************/

static double *
hist_sigma_integrals (int type)
{
  double **cache;
  double *s;
  double lo, hi, h, x, f, sum;
  int index;
  int ncache;
  int i, k;

  if (type == 1)
  {
    if (hist_sigma_outer == NULL)
    {
      hist_sigma_outer = calloc (nphot_total, sizeof (double *));
    }
    cache = hist_sigma_outer;
    index = xtop - phot_top;
    ncache = nphot_total;
  }
  else
  {
    if (hist_sigma_inner == NULL)
    {
      hist_sigma_inner = calloc (n_inner_tot, sizeof (double *));
    }
    cache = hist_sigma_inner;
    index = xtop - inner_cross;
    ncache = n_inner_tot;
  }

  if (index < 0 || index >= ncache)
  {
    Error ("hist_sigma_integrals: cross section index %d out of range (%d)\n", index, ncache);
    Exit (1);
  }

  if (cache[index] != NULL)
  {
    return (cache[index]);
  }

  s = calloc (geo.nbins_in_cell_spec, sizeof (double));
  for (i = 0; i < geo.nbins_in_cell_spec; i++)
  {
    lo = geo.cell_freq[i];
    hi = geo.cell_freq[i + 1];
    if (lo < xtop->freq[0])
    {
      lo = xtop->freq[0];
    }
    if (hi > xtop->freq[xtop->np - 1])
    {
      hi = xtop->freq[xtop->np - 1];
    }
    if (hi <= lo)
    {
      s[i] = 0.0;
      continue;
    }

    /* integral of sigma dnu/nu = integral of sigma d(ln nu) */
    h = (log (hi) - log (lo)) / 8.;
    sum = 0.0;
    for (k = 0; k <= 8; k++)
    {
      x = log (lo) + k * h;
      f = sigma_phot (xtop, exp (x));
      if (k == 0 || k == 8)
      {
        sum += f;
      }
      else if (k % 2 == 1)
      {
        sum += 4. * f;
      }
      else
      {
        sum += 2. * f;
      }
    }
    s[i] = sum * h / 3.;
  }

  cache[index] = s;
  return (s);
}


/**********************************************************/
/**
 * @brief      The integrals of sigma/nu over each cell spectrum bin for
 * the cross section calc_pi_rate would use
 *
 * @param [in] int  n   The ion (type 1) or the inner-shell cross section (type 2)
 * @param [in] int  type   1 for the outer-shell (ground state), 2 for an inner-shell cross section
 * @return     A pointer to the geo.nbins_in_cell_spec integrals, or NULL if
 * the ion has no photoionization cross section
 *
 * @details
 * Used by windsave2fits, so that the photoionization rate of any description
 * of J_nu can be computed outside Sirocco: the rate is
 * 4 pi / h * sum over bins of J_nu (bin) * integral.  geo.cell_freq must be set.
 **********************************************************/

double *
pi_rate_bin_integrals (int n, int type)
{
  if (type == 1)
  {
    if (n < 0 || n >= nions)
    {
      return (NULL);
    }
    if (ion[n].phot_info > 0)
    {
      xtop = &phot_top[ion[n].ntop_ground];
    }
    else if (ion[n].phot_info == 0)
    {
      xtop = &phot_top[ion[n].nxphot];
    }
    else
    {
      return (NULL);
    }
  }
  else
  {
    if (n < 0 || n >= n_inner_tot)
    {
      return (NULL);
    }
    xtop = &inner_cross[n];
  }
  return (hist_sigma_integrals (type));
}


/**********************************************************/
/**
 * @brief      The integral of J_nu sigma / nu over frequency, with J_nu
 * from the fine cell histogram
 *
 * @param [in] PlasmaPtr  xplasma   The cell
 * @param [in] int  type   1 for an outer-shell, 2 for an inner-shell cross section (xtop)
 * @return     The integral; calc_pi_rate multiplies it by 4 pi / h
 *
 * @details
 * J_nu is taken as constant in each bin of the normalised cell spectrum
 * (cell_spec_flux).  Used by calc_pi_rate mode 3, for band model
 * cell_spectrum and for the development table written with -pi_rate_table.
 **********************************************************/

double
pi_rate_hist_integral (PlasmaPtr xplasma, int type)
{
  double *s;
  double sum;
  int i;

  s = hist_sigma_integrals (type);

  sum = 0.0;
  for (i = 0; i < geo.nbins_in_cell_spec; i++)
  {
    if (s[i] > 0)
    {
      sum += xplasma->est.cell_spec_flux[i] * s[i];
    }
  }
  return (sum);
}


/**********************************************************/
/**
 * @brief      The direct-estimator photoionization rate of an ion
 *
 * @param [in] int  nion   The ion
 * @param [in] PlasmaPtr  xplasma   The cell
 * @return     The rate per ion (s^-1), or 0 if the ion has no density
 *
 * @details
 * The number of photoionizations counted during transport (est.ioniz)
 * divided by the density of the ion, as in matrix_est.  Whether
 * excited levels are counted is set by -direct_ground_only.
 **********************************************************/

double
pi_rate_direct (int nion, PlasmaPtr xplasma)
{
  if (xplasma->state.density[nion] > 0)
  {
    return (xplasma->est.ioniz[nion] / xplasma->state.density[nion]);
  }
  return (0.0);
}


/**********************************************************/
/**
 * @brief      The photoionization rate of an ion in matrix_pow and matrix_multishot
 *
 * @param [in] int  nion   The ion
 * @param [in] PlasmaPtr  xplasma   The cell
 * @return     The rate per ion (s^-1)
 *
 * @details
 * Normally the band models (calc_pi_rate mode 1); with band model cell_spectrum the
 * fine cell histogram (calc_pi_rate mode 3).
 **********************************************************/

double
pi_rate_matrix_pow (int nion, PlasmaPtr xplasma)
{
  if (geo.band_model == BAND_MODEL_CELL_SPECTRUM)
  {
    return (calc_pi_rate (nion, xplasma, 3, 1));
  }
  return (calc_pi_rate (nion, xplasma, 1, 1));
}


/**********************************************************/
/**
 * @brief      The inner-shell photoionization rate for one inner-shell
 * cross section in matrix_pow and matrix_multishot
 *
 * @param [in] int  n   The inner-shell cross section (index into inner_cross)
 * @param [in] PlasmaPtr  xplasma   The cell
 * @return     The rate per ion (s^-1)
 *
 * @details
 * As pi_rate_matrix_pow, but for inner shells.
 **********************************************************/

double
inner_rate_matrix_pow (int n, PlasmaPtr xplasma)
{
  if (geo.band_model == BAND_MODEL_CELL_SPECTRUM)
  {
    return (calc_pi_rate (n, xplasma, 3, 2));
  }
  return (calc_pi_rate (n, xplasma, 1, 2));
}


/**********************************************************/
/**
 * @brief      Write the band, histogram and direct photoionization
 * rates of every ion in a cell
 *
 * @param [in] PlasmaPtr  xplasma   The cell
 * @return     Always 0
 *
 * @details
 * Development (-pi_rate_table).  Called once per cell per ionization cycle,
 * after the band models have been fitted and before the ionization is
 * calculated, so all three rates come from the same photons.  Each MPI
 * rank writes the cells it updates to root.pi_rates.<rank>.txt.
 *
 * Columns: the cycle, the cell (nplasma), the ion, z, the ionization state,
 * the ion density, and the outer-shell photoionization rates (s^-1) from
 * the band models, the fine histogram and the direct estimators.
 * Ions in their highest ionization state are left out.
 **********************************************************/

int
pi_rate_table_write (PlasmaPtr xplasma)
{
  static FILE *fptr = NULL;
  char filename[2 * LINELENGTH];
  int nion;
  double r_bands;
  double r_hist;
  double r_direct;

  if (fptr == NULL)
  {
    snprintf (filename, sizeof (filename), "%s.pi_rates.%d.txt", files.root, rank_global);
    fptr = fopen (filename, "w");
    if (fptr == NULL)
    {
      Error ("pi_rate_table_write: could not open %s\n", filename);
      return (0);
    }
    fprintf (fptr, "# Photoionization rates (s^-1) from the band models, the fine cell histogram and the direct estimators\n");
    fprintf (fptr, "# direct_ground_only %d; histogram %d bins\n", modes.direct_ground_only, geo.nbins_in_cell_spec);
    fprintf (fptr, "cycle nplasma nion z istate density rate_bands rate_hist rate_direct\n");
  }

  for (nion = 0; nion < nions; nion++)
  {
    if (ion[nion].istate == ele[ion[nion].nelem].istate_max)
    {
      continue;
    }
    if (ion[nion].phot_info < 0)
    {
      continue;
    }
    r_bands = calc_pi_rate (nion, xplasma, 1, 1);
    r_hist = calc_pi_rate (nion, xplasma, 3, 1);
    r_direct = pi_rate_direct (nion, xplasma);
    fprintf (fptr, "%d %d %d %d %d %.4e %.6e %.6e %.6e\n", geo.wcycle, xplasma->nplasma, nion, ion[nion].z, ion[nion].istate,
             xplasma->state.density[nion], r_bands, r_hist, r_direct);
  }
  fflush (fptr);

  return (0);
}
