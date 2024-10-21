/** ************************************************************************* */
/**
 * @file     py_optical_depth_transport.c
 * @author   Edward Parkinson
 * @date     April 2019
 *
 * @brief    Functions related to the transport of photons and the main
 *           algorithm functions.
 *
 * ************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "atomic.h"
#include "sirocco.h"
#include "sirocco_optd.h"
static const double MAXDIFF = VCHECK / VLIGHT;  // For linear velocity requirement for photon transport

/* ************************************************************************* */
/**
 * @brief  Calculate the total optical depth a photon experiences across the
 *         cell of distance SMAX_FRAC * smax.
 * @param[in]  photon  The photon packet
 * @param[in,out]  *c_column_density  The column density the photon has moved
 *                                    through
 * @param[in,out]  *c_optical_depth  The optical depth experienced by the photon
 *
 * @return p_istat  The current photon status or EXIT_FAILURE on failure.
 *
 * @details
 *
 * This function is concerned with finding the opacity of the photon's current
 * cell, the distance the photon can move in the cell and hence it increments
 * the optical depth tau a photon has experienced as it moves through the wind.
 *
 * ************************************************************************** */

int
integrate_tau_across_cell (PhotPtr photon, double *c_column_density, double *c_optical_depth)
{
  int p_istat;
  int n_domain, n_plasma;
  double kappa_total, density;
  double smax, diff;
  double freq_inner, freq_outer, mean_freq;
  WindPtr c_wind_cell;
  PlasmaPtr c_plasma_cell;
  struct photon p_start, p_stop, p_now;

  photon->grid = where_in_grid (wmain[photon->grid].ndom, photon->x);
  if (photon->grid < 0)
  {
    printf ("integrate_tau_across_cell: photon is not in a grid cell\n");
    return EXIT_FAILURE;
  }

  /*
   * Determine where in the plasma grid the cell is. This is required so a
   * column density can be calculated. Also calculate the maximum distance the
   * photon can traverse across the cell
   */

  c_wind_cell = &wmain[photon->grid];
  n_domain = c_wind_cell->ndom;
  n_plasma = wmain[photon->grid].nplasma;
  c_plasma_cell = &plasmamain[n_plasma];

  if (COLUMN_MODE == COLUMN_MODE_RHO)
  {
    density = c_plasma_cell->rho;
  }
  else
  {
    density = c_plasma_cell->density[COLUMN_MODE_ION_NUMBER];
  }

  smax = smax_in_cell (photon) * SMAX_FRAC;
  if (smax < 0)
  {
    errormsg ("smax %e < 0 in cell %d\n", smax, photon->grid);
    return EXIT_FAILURE;
  }

  /*
   * Transform the photon into the CMF and create a photon in the CMF at the
   * end of the path of length smax
   */

  observer_to_local_frame (photon, &p_start);
  stuff_phot (photon, &p_stop);
  move_phot (&p_stop, smax);
  observer_to_local_frame (&p_stop, &p_stop);

  /* At this point p_start and p_stop are in the local frame
   * at the and p_stop is at the maximum distance it can
   * travel. We want to check that the frequency shift is
   * not too great along the path that a linear approximation
   * to the change in frequency is not reasonable
   */

  while (smax > DFUDGE)
  {
    stuff_phot (photon, &p_now);
    move_phot (&p_now, smax * 0.5);
    observer_to_local_frame (&p_now, &p_now);
    diff = fabs (p_now.freq - 0.5 * (p_start.freq + p_stop.freq)) / p_start.freq;
    if (diff < MAXDIFF)
      break;
    stuff_phot (&p_now, &p_stop);
    smax *= 0.5;
  }

  freq_inner = p_start.freq;
  freq_outer = p_stop.freq;
  mean_freq = 0.5 * (freq_inner + freq_outer);

  /*
   * Now we can finally calculate the opacity due to all the continuum
   * processes. In macro-atom mode, we need to calculate the continuum opacity
   * using kappa_bf and kappa_ff using the macro treatment. For simple mode, we
   * can simply use radiation which **SHOULD** return the continuum opacity
   * as well, plus something from induced Compton heating. In either cases,
   * we still then need to add the optical depth from electron scattering at
   * the end.
   */

  kappa_total = 0;

  if (RUN_MODE != RUN_MODE_ES_PHOTOSPHERE)
  {
    if (geo.rt_mode == RT_MODE_2LEVEL)
    {
      kappa_total += radiation (photon, smax);
    }
    else                        // macro atom case
    {
      if (c_wind_cell->vol > 0)
      {
        kappa_total += kappa_bf (c_plasma_cell, freq_inner, 0);
        kappa_total += kappa_ff (c_plasma_cell, freq_inner);
      }
    }
  }

  if (RUN_MODE != RUN_MODE_NO_ES_OPACITY)
  {
    kappa_total += klein_nishina (mean_freq) * c_plasma_cell->ne * zdom[n_domain].fill;
  }

  /*
   * Increment the optical depth and column density variables and move the
   * photon to the edge of the cell
   */

  photon->nscat++;

  *c_column_density += smax * density;
  *c_optical_depth += smax * kappa_total;
  move_phot (photon, smax);
  p_istat = photon->istat;

  return p_istat;
}

/* ************************************************************************* */
/**
 * @brief           Extract the optical depth the photon packet porig must
 *                  travel through to reach the observer.
 *
 * @param[in]  photon  The photon packet to extract
 * @param[out]  *c_column_density  The column depth of the extracted photon angle
 * @param[out]  *c_optical_depth  The optical depth from photon origin to
 *                                the observer
 *
 * @return  EXIT_SUCCESS or EXIT_FAILURE
 *
 * @details
 *
 * The direction of the observer is set in porig.lmn and should be set prior
 * to passing a photon to this function. If any values of lmn are negative, the
 * photon will translate in the negative direction. However, have no fear as this
 * is normal and is fine due to the assumed symmetry of models in Python.
 *
 * ************************************************************************** */

int
integrate_tau_across_wind (PhotPtr photon, double *c_column_density, double *c_optical_depth)
{
  int err;
  int n_dom, where;
  enum istat_enum p_istat;
  const int max_translate_in_space = 10;
  int n_in_space;
  double norm[3];
  struct photon p_extract;

  p_istat = P_INWIND;           // assume photon is in wind for initialisation reasons
  stuff_phot (photon, &p_extract);

  n_in_space = 0;
  while (p_istat == P_INWIND)
  {
    where = where_in_wind (p_extract.x, &n_dom);

    if (where < 0)
    {
      translate_in_space (&p_extract);
      if (++n_in_space > max_translate_in_space)
      {
        errormsg ("something has gone wrong as this photon has translated in space %d times\n", n_in_space);
        return EXIT_FAILURE;
      }
    }
    else if ((p_extract.grid = where_in_grid (n_dom, p_extract.x)) >= 0)
    {
      err = integrate_tau_across_cell (&p_extract, c_column_density, c_optical_depth);
      if (err)
        return EXIT_FAILURE;
    }
    else
    {
      errormsg ("photon in unknown location grid stat %i\n", p_extract.grid);
      return EXIT_FAILURE;
    }

    p_istat = walls (&p_extract, photon, norm);

    if (RUN_MODE == RUN_MODE_ES_PHOTOSPHERE)
    {
      if (*c_optical_depth >= TAU_DEPTH)
      {
        p_istat = P_ABSORB;
      }
    }
  }

  /*
   * If we are in RUN_MODE_ES_PHOTOSPHERE, then we shouldn't care about hitting the
   * star or disc, since we are aiming for the origin of the system
   */

  if (RUN_MODE == RUN_MODE_TAU_INTEGRATE)
  {
    if (p_istat == P_HIT_STAR || p_istat == P_HIT_DISK)
    {
      errormsg ("photon hit central source or disk incorrectly istat = %i\n", p_istat);
      return EXIT_FAILURE;
    }
  }
  else
  {
    stuff_phot (&p_extract, photon);
    if (p_istat == P_HIT_DISK)
    {
      errormsg ("the photon hit the disk whilst in RUN_MODE_ES_PHOTOSPHERE when it should hit the central source\n");
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
