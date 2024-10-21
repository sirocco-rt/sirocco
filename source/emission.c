
/***********************************************************/
/** @file  emission.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  Subroutines having to do with radiation from the wind,
 * including those that generate photons for the wind
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

/**********************************************************/
/**
 * @brief      calculate the energy radiated by the wind
 * wind between freqencies f1 and f2 either using a the co-moving
 * frame time step or an observer frame time step
 *
 * @param [in] double  f1         The minimum frequency
 * @param [in] double  f2         The maximum frequency
 * @param [in] int     mode       A switch indicating whether energy
 *  is to be calculated in the observer (MODE_OBSERVER_TIME) or
 *  co-moving frame (MODE_CMF_TIME).
 * @return     The luminosity of the entire wind
 *
 * @details
 *
 * The name of the routine is something of a misnomer.  The
 * routine calculates amount of energy radiated per unit time,
 * from the emissivity of the plasma x its volume in the co-moving
 * frame.  It does not account for the fact that these photons
 * will be Doppler shifted as a result of the motion of the wind
 * cells.
 *
 * In MODE_CMF_TIME, the time step for all cells is 1 second.
 *
 * In MODE_OBSERVER_FRAME_TIME, the time step is reduced by gamma
 * to accoun for time dilation in the observer frame.
 *
 * The routine simply calls total_emission for each wind cell with
 * a positive volume in the wind
 *
 * The routine also populates several luminosity related variables in
 * geo, which give the luminosity for separate processes, e.g free-free
 * and free-bound emission.
 *
 * ### Notes ###
 * @bug The do loop might be simpler if made over the plasma cells
 * instead of the wind, but one should be careful of the dummy cell
 *
 *
 * CK20180801:
 *
 *           in non-macro atom mode, the only continuum process treated as scattering is
 *           electron scattering, and this is assigned nres = -1. The only valid values
 *           of nres in non-macro-atom mode are therefore nres = -1 and 0 <= nres <= nlines-1
 *           (with the lattter range covering the lines).
 *
 *           in macro atom mode, nres = -1 indicates electron scattering,
 *           nres = -2 indicates ff, and nres > NLINES indicates bound-free.
 * 	     [nres == NLINES is never used. Note also that NLINES is the *max* number of lines, whereas nlines
 *	     is the *actual* number of lines. So, actually, it's not just nres = NLINES that's never used, but
 *	     the entire range of nlines <= nres <= NLINES]
 *
 **********************************************************/

double
wind_luminosity (double f1, double f2, int mode)
{
  int n_plasma;
  int n_start;
  int n_stop;
  double total_lum;
  double lum_lines;
  double lum_rad_recomb;
  double lum_free_free;
  double gamma_factor;
  int n_cells_rank;

#ifdef MPI_ON
  n_cells_rank = get_parallel_nrange (rank_global, NPLASMA, np_mpi_global, &n_start, &n_stop);
#else
  n_start = 0;
  n_stop = NPLASMA;
  n_cells_rank = 0;
#endif

  /* Each rank will find the total emission for a subset of the wind grid */
  for (n_plasma = n_start; n_plasma < n_stop; ++n_plasma)
  {
    if (wmain[plasmamain[n_plasma].nwind].inwind < W_ALL_INWIND)
    {
      Error ("wind_luminosty: Trying to calculate luminosity for a wind cell %d that has plasma cell %d but is not in the wind\n",
             plasmamain[n_plasma].nwind, n_plasma);
    }

    total_emission (&plasmamain[n_plasma], f1, f2);
  }

  broadcast_wind_luminosity (n_start, n_stop, n_cells_rank);

  total_lum = 0.0;
  lum_lines = 0.0;
  lum_free_free = 0.0;
  lum_rad_recomb = 0.0;

  for (n_plasma = 0; n_plasma < NPLASMA; n_plasma++)
  {

    if (mode == MODE_OBSERVER_FRAME_TIME)
    {
      gamma_factor = 1.0 / plasmamain[n_plasma].xgamma; /* this is dt_cmf */
    }
    else
    {
      gamma_factor = 1.0;
    }

    lum_lines += plasmamain[n_plasma].lum_lines * gamma_factor;
    lum_rad_recomb += plasmamain[n_plasma].lum_rr * gamma_factor;
    lum_free_free += plasmamain[n_plasma].lum_ff * gamma_factor;
  }

  total_lum = lum_lines + lum_rad_recomb + lum_free_free;

  if (mode == MODE_CMF_TIME)
  {
    geo.lum_lines = lum_lines;
    geo.lum_rr = lum_rad_recomb;
    geo.lum_ff = lum_free_free;
  }

  return (total_lum);
}


/**********************************************************/
/**
 * @brief      Calculate the band-limited emission of a single cell
 *
 * @param [in] PlasmaPtr  xplasma   The plasma cell of interest
 * @param [in] double  f1   The minimum frequency for the calculation
 * @param [in] double  f2   The maximum frequency for the calculation
 * @return
 * It returns the total luminosity (within frequency limits)
 *
 * The routine also stores
 * the luminosity due
 * to various emission processes, e.g ff, fb, lines, compton into
 * varius variables in thea associated Plasma cell
 *
 * @details
 *
 * ### Notes ###
 * Total emission gives the total enery loss due to photons.  It does
 * not include other cooling sources, e. g. adiabatic expansion.
 *
 * The name total emission is a misnomer.  The returns a
 * band limited luminosity.  This is because the routine can be used
 * to establish the number of photons to be emitted by the wind.
 *
 * Comment:  Compton cooling is not included here because this
 * does not result in new photons.
 *
 *
 **********************************************************/

double
total_emission (xplasma, f1, f2)
     PlasmaPtr xplasma;
     double f1, f2;
{
  double t_e;

  t_e = xplasma->t_e;


  if (f2 < f1)
  {
    xplasma->lum_tot = xplasma->lum_lines = xplasma->lum_ff = xplasma->lum_rr = 0;
  }
  else
  {
    if (geo.rt_mode == RT_MODE_MACRO)
    {
      xplasma->lum_rr = (total_fb_matoms (xplasma, t_e, f1, f2) + total_fb (xplasma, t_e, f1, f2, FB_FULL, OUTER_SHELL));       //outer shellrecombinations

      /*
       *The first term here is the fb cooling due to macro ions and the second gives
       *the fb cooling due to simple ions.
       *total_fb has been modified to exclude recombinations treated using macro atoms.
       */
      xplasma->lum_tot = xplasma->cool_rr;
      /* Note: This the fb_matom call makes no use of f1 or f2. They are passed for
       * now in case they should be used in the future. But they could
       * also be removed.
       * (SS)
       */
      xplasma->lum_lines = total_bb_cooling (xplasma, t_e);
      xplasma->lum_tot += xplasma->lum_lines;
      /* total_bb_cooling gives the total cooling rate due to bb transisions whether they
         are macro atoms or simple ions. */
      xplasma->lum_ff = total_free (xplasma, t_e, f1, f2);
      xplasma->lum_tot += xplasma->lum_ff;


    }
    else                        //default (non-macro atoms) (SS)
    {
      xplasma->lum_tot = xplasma->lum_lines = total_line_emission (xplasma, f1, f2);
      xplasma->lum_tot += xplasma->lum_ff = total_free (xplasma, t_e, f1, f2);
      /* We compute the radiative recombination luminosity - this is not the same as the rr cooling rate and
         so is stored in a separate variable */
      xplasma->lum_tot += xplasma->lum_rr = total_fb (xplasma, t_e, f1, f2, FB_FULL, OUTER_SHELL);      //outer shell recombinations
    }
  }


  return (xplasma->lum_tot);
}




/**********************************************************/
/**
 * @brief      generates
 * 	photons within the wind between two freqencies and stores them in the photon array
 *
 * @param [out] PhotPtr  p   The entire photon stucture
 * @param [in] double  weight   weight of each photon to be generated
 * @param [in] double  freqmin   The minimum frequency
 * @param [in] double  freqmax   The maximum frequency
 * @param [in] int  photstart   The place in the photon stucture where the first photon will be stored
 * @param [in] int  nphot   The number of photon bundles to generate
 * @return     Normally returns the number of photons created.
 *
 * When the routine is completed new  photons will exist from p[istart]
 * to p[istart+return value].
 *
 * @details
 *
 * The routine first generates a random number which is used to determine
 * in which wind cell  should be generated, and then  determines the
 * type of photon to generate.  Once this is done the routine cycles through
 * the PlasmaCells generatating all of the photons for each cell at once.
 *
 *
 * ### Notes ###
 * This logic was adopted for speed related reasons.
 *
 * First the routine determines how many photons should be generated in each (plasma) cell
 * of each type
 *
 * Then it generates the photons on a cell by cell basis.
 *
 * If photo_gen_wind tries to create more photons than exist in the photon structure the
 * program will stop (rather than continue incorrectly or start blasting away memory.)
 *
 **********************************************************/

#define FREE_FREE  0
#define FREE_BOUND 1
#define BOUND_BOUND 2

int
photo_gen_wind (p, weight, freqmin, freqmax, photstart, nphot)
     PhotPtr p;
     double weight;
     double freqmin, freqmax;
     int photstart, nphot;
{
  int nn, np;
  int kkk;
  int photstop;
  double xlum, xlumsum, lum, dt_cmf;
  int icell, icell_old;
  int nplasma = 0;
  int nnscat;
  int ptype[NPLASMA][3];        //Store for the types of photons to be generated in each cell, ff first, fb next, line third
  dt_cmf = 0.0;

  for (nplasma = 0; nplasma < NPLASMA; nplasma++)
  {
    for (nn = 0; nn < 3; nn++)
      ptype[nplasma][nn] = 0;
  }

  limit_lines (freqmin, freqmax);

  photstop = photstart + nphot;
  Log_silent ("photo_gen_wind creates nphot %5d photons from %5d to %5d \n", nphot, photstart, photstop);

  for (kkk = photstart; kkk < photstop; kkk++)
  {
    p[kkk].nres = -1;
    p[kkk].nnscat = 1;

    /* Locate the wind_cell in which the photon bundle originates, allowing for the
       fact that we want to create the correct number of photons using observer
       an obsever frame time step. geo.f_wind is the energy generated by the wind
       in the desired band in the observer frame, wherease plaama_main[].lum_tot is
       was calculated in the CMF. */

    xlum = random_number (0.0, 1.0) * geo.f_wind;

    xlumsum = 0;
    nplasma = 0;

    while (xlumsum < xlum)
    {
      dt_cmf = 1.0 / plasmamain[nplasma].xgamma;
      xlumsum += plasmamain[nplasma].lum_tot * dt_cmf;
      nplasma++;
    }
    nplasma--;

    if (nplasma == NPLASMA)
    {
      Error ("photo_gen_wind: Generating photon in a cell that should not be in wind (%.2e > %.2e) geo.f_wind %.2e\n", xlumsum, xlum,
             geo.f_wind);
    }

    /* At this point we know the cell in which the photon will be generated */

    plasmamain[nplasma].nrad += 1;

    /*Determine the type of photon this photon will be and increment ptype, which stores the total number of
     * each photon type to be made in each cell. We don't need to account for time
     * dilation heere, because all processes are in the same cell*/

    lum = plasmamain[nplasma].lum_tot;
    xlum = lum * random_number (0.0, 1.0);
    xlumsum = 0;

    if ((xlumsum += plasmamain[nplasma].lum_ff) > xlum)
    {
      ptype[nplasma][FREE_FREE]++;
    }
    else if ((xlumsum += plasmamain[nplasma].lum_rr) > xlum)
    {
      ptype[nplasma][FREE_BOUND]++;
    }
    else
    {
      ptype[nplasma][BOUND_BOUND]++;
    }
  }



/* Now generate the photons looping over the Plasma cells */

  photstop = photstart;

  icell_old = (-1);
  for (nplasma = 0; nplasma < NPLASMA; nplasma++)
  {

    photstart = photstop;
    photstop = photstart + ptype[nplasma][FREE_FREE] + ptype[nplasma][FREE_BOUND] + ptype[nplasma][BOUND_BOUND];

    icell = plasmamain[nplasma].nplasma;

    for (np = photstart; np < photstop; np++)
    {

      if (np < photstart + ptype[nplasma][FREE_FREE])
      {
        p[np].freq = one_ff (&plasmamain[nplasma], freqmin, freqmax);
        if (p[np].freq <= 0.0)
        {
          Error_silent
            ("photo_gen_wind: On return from one_ff: icell %d vol %g t_e %g\n", nplasma, plasmamain[nplasma].vol, plasmamain[nplasma].t_e);
          p[np].freq = 0.0;
        }
      }
      else if (np < photstart + ptype[nplasma][FREE_FREE] + ptype[nplasma][FREE_BOUND])
      {
        p[np].freq = one_fb (&plasmamain[nplasma], freqmin, freqmax);
      }
      else
      {

        /* fill the lin_ptr->pow array. This must be done because it is not stored for all cells.
           The if statement is intended to prevent recalculating the power if more than one
           line photon is generated from this cell in this cycle.
         */
        if (icell != icell_old)
        {
          lum_lines (&plasmamain[nplasma], nline_min, nline_max);
          icell_old = icell;
        }

        p[np].freq = one_line (&plasmamain[nplasma], &p[np].nres);
        if (p[np].freq == 0)
        {
          Error ("photo_gen_wind: one_line returned 0 for freq %g %g\n", freqmin, freqmax);
        }
      }

      p[np].w = weight;
      get_random_location (plasmamain[nplasma].nwind, p[np].x);
      p[np].grid = plasmamain[nplasma].nwind;

      nnscat = 1;

      /* Select a direction for the photon, depending on the scattering mode and/or
         the type of photon that was generated
       */

      if (p[np].nres < 0 || geo.scatter_mode == SCATTER_MODE_ISOTROPIC)
      {
        randvec (p[np].lmn, 1.0);
      }
      else if (geo.scatter_mode == SCATTER_MODE_THERMAL)
      {
        randwind_thermal_trapping (&p[np], &nnscat);
      }
      p[np].nnscat = nnscat;

      /* Photons are generated in the local frame and so must be Doppler shifted into
         the observer frame.
       */

      p[np].istat = P_INWIND;
      p[np].frame = F_LOCAL;
      p[np].tau = p[np].nscat = p[np].nrscat = p[np].nmacro = 0;
      p[np].origin = PTYPE_WIND;

      if (local_to_observer_frame (&p[np], &p[np]))
      {
        Error ("photo_gen_wind: frame tranformation error\n");
      }

      /* Extra processing for reveberation calculations */
      switch (geo.reverb)
      {
      case REV_WIND:
      case REV_MATOM:
        wind_paths_gen_phot (&wmain[plasmamain[nplasma].nwind], &p[np]);
        break;
      case REV_PHOTON:
        simple_paths_gen_phot (&p[np]);
        break;
      case REV_NONE:
      default:
        break;
      }

    }
  }


  return (nphot);
}




/**********************************************************/
/**
 * @brief      gets the frequency of a
 * single collisionally excited line photon in a particular cell
 * of the wind.
 *
 * @param [in] WindPtr  one   The wind cell
 * @param [out] int *  nres   The number associated with the transition
 * @return     The frequency of the transition
 *
 * @details
 * The routine generates a random number and uses this to select
 * the transition that was excited
 *
 * ### Notes ###
 *
 *
 **********************************************************/

double
one_line (xplasma, nres)
     PlasmaPtr xplasma;
     int *nres;
{
  double xlum, xlumsum;
  int m;
  /* Put in a bunch of checks */
  if (xplasma->lum_lines <= 0)
  {
    Error ("one_line: requesting a line when line lum is 0\n");
    return (0);
  }
  if (nline_min == nline_max)
  {
    Error ("one_line: no lines %d %d\n", nline_min, nline_max);
    return (0);
  }


  xlum = xplasma->lum_lines * random_number (0.0, 1.0);

  xlumsum = 0;
  m = nline_min;
  while (xlumsum < xlum && m < nline_max)
  {
    xlumsum += lin_ptr[m]->pow;
    m++;
  }
  m--;
  *nres = m;
  return (lin_ptr[m]->freq);
}


/**********************************************************/
/** Next section deals with bremsstrahlung radiation
 * 4*PI* 8/3* sqrt(2*PI/3)*e**6/m**2/c**3 sqrt(m/kT) or
 * 4*PI times the normal constant to dL_nu/dnu
 **********************************************************/
#define BREMS_CONSTANT 6.85e-38


/**********************************************************/
/**
 * @brief      calculates the band-limited ff luminosity of a cell.
 *
 * @param [in] PlasmaPtr  xplasma   A plasma cell
 * @param [in] double  t_e   The electron temperature to use
 * @param [in] double  f1   The minimum frequency
 * @param [in] double  f2   The maximum frequency
 * @return     The free-free luminosity between f1 and f2
 *
 * @details
 *
 * ### Notes ###
 * The program uses an integral formula rather than integrating on
 * the fly which is not ideal but saves time
 *
 * The gaunt factor used depends on whether data for this is read in
 * from the atomic data files.  If it is, then the gaunt factor determined
 * using data from Sutherland (1998).  If not, the gaunt factor is set
 * to 1.
 *
 * Total free free emission in this function is calculated CMF.
 *
 **********************************************************/

double
total_free (xplasma, t_e, f1, f2)
     PlasmaPtr xplasma;
     double t_e;
     double f1, f2;
{
  double g_ff_h, g_ff_he;
  double gaunt;
  double x, sum;
  double gsqrd;                 /*The scaled inverse temperature  - used to compute the gaunt factor */
  int nion;
  if (f2 < f1)
  {
    Error ("total_free: band limited ff  emissivity requested by f1 %g > f2 %g\n", f1, f2);
    return (0.0);
  }

  if (ALPHA_FF * xplasma->t_e / H_OVER_K < f1)
  {
    return (0.0);
  }

  if (ALPHA_FF * xplasma->t_e / H_OVER_K < f2)
  {
    f2 = ALPHA_FF * xplasma->t_e / H_OVER_K;
  }

  if (t_e < TMIN)
  {
    return (0.0);
  }


  if (gaunt_n_gsqrd == 0)       //fallback soln if atomic data files do not include gaunt factor.
  {
    g_ff_h = g_ff_he = 1.0;
    if (nelements > 1)
    {
      x = BREMS_CONSTANT * xplasma->ne * (xplasma->density[1] * g_ff_h + 4. * xplasma->density[4] * g_ff_he) / H_OVER_K;
    }
    else
    {
      x = BREMS_CONSTANT * xplasma->ne * (xplasma->density[1] * g_ff_h) / H_OVER_K;
    }
  }
  else
  {
    sum = 0.0;
    for (nion = 0; nion < nions; nion++)
    {
      if (ion[nion].istate != 1)        //neutral ions do not contribute
      {
        gsqrd = ((ion[nion].istate - 1) * (ion[nion].istate - 1) * RYD2ERGS) / (BOLTZMANN * t_e);
        gaunt = gaunt_ff (gsqrd);
        sum += xplasma->density[nion] * (ion[nion].istate - 1) * (ion[nion].istate - 1) * gaunt;
      }
    }
    x = BREMS_CONSTANT * xplasma->ne * (sum) / H_OVER_K;
  }

  /* JM 1604 -- The reason why this is proportional to t_e**1/2,
     rather than t_e**(-1/2) as in equation 40 of LK02 is because
     one gets an extra factor of (k*t_e/h) when one does the integral */
  x *= sqrt (t_e) * xplasma->vol;
  x *= (exp (-H_OVER_K * f1 / t_e) - exp (-H_OVER_K * f2 / t_e));
  return (x);
}



/**********************************************************/
/**
 * @brief      calculate f_nu for free free emisssion
 *
 * @param [in] PlasmaPtr  one   A plasma  cell
 * @param [in] double  t_e   The temperature of the plasma
 * @param [in] double  freq   The frequency at which f_nu is to be calculated
 * @return     f_nu for the specific freqency and temperature requested and densities
 * contained in the cell indicated
 *
 * @details
 *
 * ### Notes ###
 *
 * Within sirocco, this routine is accessed through one_ff
 *
 * Most of the computation in this routine arises from calculating
 * a the gaunt factor, so this version of the code checks to
 * see if that can be avoided.
 **********************************************************/

double ff_constant = 0;
int ff_nplasma = -100;
double ff_t_e = -100.;

double
ff (xplasma, t_e, freq)
     PlasmaPtr xplasma;
     double t_e, freq;
{
  double g_ff_h, g_ff_he;
  double fnu;
  double gsqrd, gaunt, sum;
  int nplasma;
  int nion;

  nplasma = xplasma->nplasma;


  if (t_e < TMIN)
    return (0.0);

  if (nplasma != ff_nplasma || t_e != ff_t_e)
  {
    ff_nplasma = nplasma;
    ff_t_e = t_e;

    if (gaunt_n_gsqrd == 0)     //case  where gaunt factors have not been provided in atomic data files
    {
      g_ff_h = g_ff_he = 1.0;
      if (nelements > 1)
      {
        fnu = BREMS_CONSTANT * xplasma->ne * (xplasma->density[1] * g_ff_h + 4. * xplasma->density[4] * g_ff_he);
      }
      else
      {
        fnu = BREMS_CONSTANT * xplasma->ne * (xplasma->density[1] * g_ff_h);
      }
    }
    else
    {
      sum = 0.0;
      for (nion = 0; nion < nions; nion++)
      {
        if (ion[nion].istate != 1)      //The neutral ion does not contribute
        {
          gsqrd = ((ion[nion].istate - 1) * (ion[nion].istate - 1) * RYD2ERGS) / (BOLTZMANN * t_e);
          gaunt = gaunt_ff (gsqrd);
          sum += xplasma->density[nion] * (ion[nion].istate - 1) * (ion[nion].istate - 1) * gaunt;
        }
        else
        {
          sum += 0.0;
        }
      }
      fnu = BREMS_CONSTANT * xplasma->ne * (sum);
    }

    ff_constant = fnu * xplasma->vol;
  }


  fnu = ff_constant * exp (-H_OVER_K * freq / t_e) / sqrt (t_e);
  return (fnu);
}



/// We initialise the arrays that will contain the unscaled PDF
double ff_x[ARRAY_PDF], ff_y[ARRAY_PDF];

/// Old values
double one_ff_f1, one_ff_f2, one_ff_te;


/**********************************************************/
/**
 * @brief      randomly generate the frequency of a
 * 	ff photon within the frequency interval f1 and f2
 *
 * @param [in] PlasmaPtr  xplasma   A specific plasma cell
 * @param [in] double  f1   The minimum frequency
 * @param [in] double  f2   The maximum frequency
 * @return     A randomly geneted frequncy for a ff photon given
 * the conditions specified
 *
 * @details
 *
 * ### Notes ###
 * one is the windcell where the photon will be created.  It is needed
 * only for the temperature.  The code would be simplified
 * if simply the temperature were transmitted.
 *
 * The routine creeates a cdf given the conditions specified and then
 * samples this.  When the routine is re-entered it checks if the
 * conditions are the same as previously.  Given this it is desirable
 * to generate all of the photons in a single cell in one go.
 *
 **********************************************************/

double
one_ff (xplasma, f1, f2)
     PlasmaPtr xplasma;         /* a single cell */
     double f1, f2;             /* freqmin and freqmax */
{
  double freq, dfreq;
  int n;
  int echeck;

  if (f2 < f1)
  {
    Error ("one_ff: Bad inputs f2 %g < f1 %g returning 0.0  t_e %g\n", f2, f1, xplasma->t_e);
    return (-1.0);
  }

  /* Check to see if we have already generated a pdf */

  if (xplasma->t_e != one_ff_te || f1 != one_ff_f1 || f2 != one_ff_f2)
  {                             /* Generate a new pdf */
    dfreq = (f2 - f1) / (ARRAY_PDF - 1);
    for (n = 0; n < ARRAY_PDF - 1; n++)
    {
      ff_x[n] = f1 + dfreq * n;
      ff_y[n] = ff (xplasma, xplasma->t_e, ff_x[n]);
    }

    ff_x[ARRAY_PDF - 1] = f2;
    ff_y[ARRAY_PDF - 1] = ff (xplasma, xplasma->t_e, ff_x[ARRAY_PDF - 1]);
    if ((echeck = cdf_gen_from_array (&cdf_ff, ff_x, ff_y, ARRAY_PDF, f1, f2)) != 0)
    {
      Error
        ("one_ff: cdf_gen_from_array error %d : nplasma %d f1 %g f2 %g te %g ne %g nh %g vol %g\n",
         echeck, xplasma->nplasma, f1, f2, xplasma->t_e, xplasma->ne, xplasma->density[1], xplasma->vol);
      Exit (0);
    }
    one_ff_te = xplasma->t_e;
    one_ff_f1 = f1;
    one_ff_f2 = f2;
  }
  freq = cdf_get_rand (&cdf_ff);
  return (freq);
}




/**********************************************************/
/**
 * @brief      computes the frequency averaged gaunt factor for ff emissiona
 * 	       based on scaled value of the temperature from Sutherland (1988).
 *
 * @param [in] double  gsquared   The variable on which
 * Sutherland's ff gaunt factors are interpolated.
 * @return     An estimate of the Gaunt factor (for a particular ion)
 *
 * @details
 * It interpolates simply interpolates between the tabulated
 * factors from Table 3 of Sutherland (which is a spline fit
 * to the calculated gaunt factors at various values of gsquared)
 *
 * ### Notes ###
 * The reference for this is Sutherland, R.~S. 1998, MNRAS, 300, 321
 *
 * gsquared is the scaled inverse temperature experienced by an ion,
 * Z**2/kT(Ry).
 *
 * See issue #750 for more explanation
 **********************************************************/

double
gaunt_ff (gsquared)
     double gsquared;
{
  int i, index;
  double gaunt;
  double log_g2;
  double delta;                 //The log difference between our G2 and the one in the table

  delta = 0.0;
  index = 0;
  log_g2 = log10 (gsquared);    //The data is in log format
  if (log_g2 < gaunt_total[0].log_gsqrd || log_g2 > gaunt_total[gaunt_n_gsqrd - 1].log_gsqrd)
  {
    return (1.0);
  }

  i = 0;
  while (gaunt_total[i].log_gsqrd < log_g2)
  {
    i++;
  }

  index = i - 1;
  delta = log_g2 - gaunt_total[index].log_gsqrd;

  gaunt = gaunt_total[index].gff + delta * (gaunt_total[index].s1 + delta * (gaunt_total[index].s2 + gaunt_total[index].s3));
  return (gaunt);
}
