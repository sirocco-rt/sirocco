/** ************************************************************************* */
/**
 * @file     py_optical_depth_util.c
 * @author   Edward Parkinson
 * @date     May 2021
 *
 * @brief    Functions which aren't related to the transport of photons, or
 *           creation of the results.
 *
 * ************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "atomic.h"
#include "sirocco.h"
#include "sirocco_optd.h"

/* ************************************************************************* */
/**
 * @brief  Initialize the inclination angles for a 2D outward run.
 *
 * @param[out]  int n_angles  The number of angles initialized
 *
 * @return  SightLines_t *inclinations  The initialize inclinations structure
 *
 * @details
 *
 * If xxpsec has been read in, i.e. there were spectral cycles run the model,
 * then the same inclination angles in xxspec are used. Otherwise, a set of
 * default angles will be used instead, with a default phase of 0.5. The phase
 * should not matter, as it works for, i.e., 0.5 and 1 and gives the same
 * results as of writing this.
 *
 * ************************************************************************** */

SightLines_t *
outward_initialize_2d_model_angles (int *n_angles, double *input_inclinations)
{
  int i;
  int len;
  const double default_phase = 0.5;
  const double default_angles[] = { 0.0, 10.0, 30.0, 45.0, 60.0, 75.0, 85.0, 90.0 };
  const int n_default_angles = sizeof default_angles / sizeof default_angles[0];

  /*
   * Count the number of inclination angles provided by the -i command. Since
   * the array is initialized as -1, this assumes any values > -1 is a valid
   * inclination angle.
   */

  int n_input = 0;
  for (i = 0; i < MAX_CUSTOM_ANGLES; ++i)
  {
    if (input_inclinations[i] > -1)
    {
      n_input++;
    }
  }

  /*
   * Use the angles specified for by the user for spectrum generation, this
   * requires for xxspec to be initialised. Otherwise use the pre-defined angles
   * above.
   */

  SightLines_t *inclinations = NULL;

  if (n_input > 0)
  {
    *n_angles = n_input;
    inclinations = calloc (n_input, sizeof *inclinations);
    if (inclinations == NULL)
    {
      errormsg ("cannot allocate %lu bytes for observers array\n", n_input * sizeof *inclinations);
      exit (EXIT_FAILURE);
    }

    for (i = 0; i < n_input; i++)
    {
      len = snprintf (inclinations[i].name, NAMELEN, "A%02.0fP%04.2f", input_inclinations[i], default_phase);
      if (len < 0)
      {
        errormsg ("there was an error writing the name to the sight lines array\n");
        exit (EXIT_FAILURE);
      }

      inclinations[i].lmn[0] = sin (input_inclinations[i] / RADIAN) * cos (-default_phase * 360.0 / RADIAN);
      inclinations[i].lmn[1] = sin (input_inclinations[i] / RADIAN) * sin (-default_phase * 360.0 / RADIAN);
      inclinations[i].lmn[2] = cos (input_inclinations[i] / RADIAN);
      inclinations[i].angle = input_inclinations[i];
    }
  }
  else if (xxspec != NULL && geo.nangles > 0)
  {
    *n_angles = geo.nangles;
    inclinations = calloc (geo.nangles, sizeof *inclinations);
    if (inclinations == NULL)
    {
      errormsg ("cannot allocate %lu bytes for observers array\n", geo.nangles * sizeof *inclinations);
      exit (EXIT_FAILURE);
    }

    for (i = MSPEC; i < MSPEC + geo.nangles; i++)
    {
      strcpy (inclinations[i - MSPEC].name, xxspec[i].name);
      stuff_v (xxspec[i].lmn, inclinations[i - MSPEC].lmn);
      inclinations[i].angle = -1;       // todo: implement way to get angle xxspec
    }
  }
  else
  {
    printf ("\nNo spec.save file has been found, using a default set of inclination angles\n\n");

    *n_angles = n_default_angles;
    inclinations = calloc (n_default_angles, sizeof *inclinations);
    if (inclinations == NULL)
    {
      errormsg ("cannot allocate %lu bytes for observers array\n", n_default_angles * sizeof *inclinations);
      exit (EXIT_FAILURE);
    }

    for (i = 0; i < n_default_angles; i++)
    {
      len = snprintf (inclinations[i].name, NAMELEN, "A%02.0fP%04.2f", default_angles[i], default_phase);
      if (len < 0)
      {
        errormsg ("there was an error writing the name to the sight lines array\n");
        exit (EXIT_FAILURE);
      }

      inclinations[i].lmn[0] = sin (default_angles[i] / RADIAN) * cos (-default_phase * 360.0 / RADIAN);
      inclinations[i].lmn[1] = sin (default_angles[i] / RADIAN) * sin (-default_phase * 360.0 / RADIAN);
      inclinations[i].lmn[2] = cos (default_angles[i] / RADIAN);
      inclinations[i].angle = default_angles[i];
    }
  }

  return inclinations;
}

/* ************************************************************************* */
/**
 * @brief  Initialize the inclination angles for a 1D outward run.
 *
 * @param[out]  int n_angles  The number of angles initialized
 *
 * @return  SightLines_t *inclinations  The initialize inclinations structure
 *
 * @details
 *
 * Since this is for a 1D model, it does not matter which inclination angle
 * we use. We have opted to used an angle of 45 degrees, as this is also the
 * angle we use to define a 1D spherical grid. It should not matter either
 * way.
 *
 * ************************************************************************** */

SightLines_t *
outward_initialize_1d_model_angles (int *n_angles)
{
  int len;
  const int n_default_angles = 1;
  const double default_angle = 45.0;
  const double default_phase = 0.5;

  *n_angles = n_default_angles;
  SightLines_t *inclinations = calloc (n_default_angles, sizeof (SightLines_t));

  if (inclinations == NULL)
  {
    errormsg ("unable to allocate %ld bytes for observers array\n", sizeof *inclinations);
    exit (EXIT_FAILURE);
  }

  len = snprintf (inclinations[0].name, NAMELEN, "A%02.0fP%04.2f", default_angle, default_phase);
  if (len < 0)
  {
    errormsg ("there was an error writing the name to the sight lines array\n");
    exit (EXIT_FAILURE);
  }

  inclinations[0].lmn[1] = sin (default_angle / RADIAN) * sin (-default_phase * 360.0 / RADIAN);
  inclinations[0].lmn[0] = sin (default_angle / RADIAN) * cos (-default_phase * 360.0 / RADIAN);
  inclinations[0].lmn[2] = cos (default_angle / RADIAN);
  inclinations[0].angle = default_angle;

  return inclinations;
}

/* ************************************************************************* */
/**
 * @brief  Initialize the inclination angles to find the photosphere for a
 *         2d model.
 *
 * @param[out]  int n_angles  The number of angles initialized
 *
 * @return  SightLines_t *inclinations  The initialize inclinations structure
 *
 * @details
 *
 * The same function is called for both 1D and 2D models. This creates extra
 * work for 1D model, but as the algorithm takes very little time to run, it
 * does not matter.
 *
 * 500 inclination angles are defined, to very finely resolve the photosphere
 * surface. This is fixed for now. I think 500 is probably far too many, but
 * it takes absolutely no time to run. The results from 500 angles probably
 * need smoothing if the grid is coarse.
 *
 * ************************************************************************** */

SightLines_t *
photosphere_2d_initialize_angles (int *n_angles)
{
  int i;
  double default_angle;
  const double default_phase = 1.0;
  const int n_default_angles = 500;
  const double d_theta = 90.0 / (double) n_default_angles;

  *n_angles = n_default_angles;
  SightLines_t *inclinations = calloc (n_default_angles, sizeof (SightLines_t));

  if (inclinations == NULL)
  {
    errormsg ("unable to allocate memory for sight lines array\n");
    exit (EXIT_FAILURE);
  }

  for (i = 0; i < n_default_angles; i++)
  {
    default_angle = i * d_theta;
    inclinations[i].lmn[0] = sin (default_angle / RADIAN) * cos (-default_phase * 360.0 / RADIAN);
    inclinations[i].lmn[1] = sin (default_angle / RADIAN) * sin (-default_phase * 360.0 / RADIAN);
    inclinations[i].lmn[2] = cos (default_angle / RADIAN);
    inclinations[i].angle = default_angle;
  }

  return inclinations;
}

/* ************************************************************************* */
/**
 * @brief  Initialize the inclination angles to find the photosphere.
 *
 * @param[out]  int n_angles  The number of angles initialized
 *
 * @return  SightLines_t *inclinations  The initialize inclinations structure
 *
 * @details
 *
 * The same function is called for both 1D and 2D models. This creates extra
 * work for 1D model, but as the algorithm takes very little time to run, it
 * does not matter.
 *
 * 500 inclination angles are defined, to very finely resolve the photosphere
 * surface. This is fixed for now. I think 500 is probably far too many, but
 * it takes absolutely no time to run. The results from 500 angles probably
 * need smoothing if the grid is coarse.
 *
 * ************************************************************************** */

SightLines_t *
photosphere_1d_initialize_angles (int *n_angles)
{
  const double default_phase = 1.0;
  const int n_default_angles = 1;
  const double default_angle = 45;

  *n_angles = n_default_angles;
  SightLines_t *inclinations = calloc (n_default_angles, sizeof (SightLines_t));

  if (inclinations == NULL)
  {
    errormsg ("unable to allocate memory for sight lines array\n");
    exit (EXIT_FAILURE);
  }

  inclinations[0].lmn[0] = sin (default_angle / RADIAN) * cos (-default_phase * 360.0 / RADIAN);
  inclinations[0].lmn[1] = sin (default_angle / RADIAN) * sin (-default_phase * 360.0 / RADIAN);
  inclinations[0].lmn[2] = cos (default_angle / RADIAN);
  inclinations[0].angle = default_angle;

  return inclinations;
}

/* ************************************************************************* */
/**
 * @brief  Wrapper function for initializing the inclination angles depending
 *         on the run mode of the program.
 *
 * @param[out]  int n_angles  The number of angles initialized
 *
 * @return  SightLines_t *inclinations  The initialize inclinations structure
 *
 * ************************************************************************** */

SightLines_t *
initialize_inclination_angles (int *n_angles, double *input_inclinations)
{
  SightLines_t *inclinations;

  if (RUN_MODE == RUN_MODE_ES_PHOTOSPHERE)
  {
    if (zdom[N_DOMAIN].coord_type == SPHERICAL)
    {
      inclinations = photosphere_1d_initialize_angles (n_angles);
    }
    else
    {
      inclinations = photosphere_2d_initialize_angles (n_angles);
    }
  }
  else
  {
    if (zdom[N_DOMAIN].coord_type == SPHERICAL)
    {
      inclinations = outward_initialize_1d_model_angles (n_angles);
    }
    else
    {
      inclinations = outward_initialize_2d_model_angles (n_angles, input_inclinations);
    }
  }

  return inclinations;
}

/* ************************************************************************* */
/**
 * @brief  Generate a photon packet with a given frequency nu for use with the
 *         optical depth diagnostic routines.
 *
 * @param[in,out]  p_out  The photon packet to initialise
 * @param[in]  freq  The frequency of the photon packet
 *
 * @return  EXIT_SUCCESS or EXIT_FAILURE
 *
 * @details
 *
 * This routine assumes that the user wishes for the photon to be generated from
 * the surface of the central source, taking into account the current direction
 * of the sight-line being extracted. It's currently not possible, without a
 * tedious re-write, to place the photon at the very origin of the grid when
 * there is a central source because of how we check boundaries.
 *
 * Note that photons are initialised with a weight of f_tot as photons are
 * required to have weight, but since functions do not care about the weight of
 * the photon, it is set to something large to make sure it does not get
 * destroyed by accident somehow.
 *
 * ************************************************************************** */

int
create_photon (PhotPtr p_out, double freq, double *lmn)
{
  int i;

  if (freq < 0)
  {
    errormsg ("photon can't be created with negative frequency\n");
    return EXIT_FAILURE;
  }

  p_out->freq = p_out->freq_orig = freq;
  p_out->origin = p_out->origin_orig = PTYPE_DISK;
  p_out->istat = P_INWIND;
  p_out->w = p_out->w_orig = geo.f_tot;
  p_out->tau = 0.0;
  p_out->frame = F_OBSERVER;
  p_out->x[0] = p_out->x[1] = p_out->x[2] = 0.0;
  stuff_v (lmn, p_out->lmn);

  if (RUN_MODE == RUN_MODE_ES_PHOTOSPHERE)
  {
    move_phot (p_out, zdom[N_DOMAIN].rmax - DFUDGE);
    for (i = 0; i < 3; ++i)
    {
      p_out->lmn[i] *= -1.0;    // Make the photon point inwards
    }

  }
  else
  {
    move_phot (p_out, geo.rstar + DFUDGE);
  }

  return EXIT_SUCCESS;
}
