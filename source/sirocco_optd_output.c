/** ************************************************************************* */
/**
 * @file     py_optical_depth_output.c
 * @author   Edward Parkinson
 * @date     May 2021
 *
 * @brief    Functions for writing the output to stdout and to file.
 *
 * ************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"
#include "sirocco_optd.h"

/* ************************************************************************* */
/**
 * @brief  Write a generic header to the provided file pointer.
 *
 * @details
 *
 * The version sirocco, commit hash and date are written to the file pointer
 * provided. This assumes that the file pointer is pointing to the top of the
 * file, as this is supposed to be writing a header.
 *
 * ************************************************************************** */

void
write_generic_file_header (FILE *fp)
{
  char time_string[LINELENGTH];

  get_time (time_string);
  fprintf (fp, "# Python Version %s\n", VERSION);
  fprintf (fp, "# Git commit hash %s\n", GIT_COMMIT_HASH);
  fprintf (fp, "# Date %s\n", time_string);
  fprintf (fp, "#\n");
}

/* ************************************************************************** */
/**
 * @brief  Print the various optical depths for the photoionization edges
 *
 * @param[in]  SightLines_t  *inclinations  The inclination angles used
 * @param[in]  int  n_inclinations          The number of inclinations
 * @param[in]  Edges_t  edges[]             The photoionization edges
 * @param[in]  int  n_edges                 The number of edges evaluated
 * @param[in]  double *optical_depth        The optical depth of the edges
 * @param[in]  double *column_density       The column density of the inclinations
 *
 * @details
 *
 * Prints the different optical depths for each angle and optical depth.
 * Historically, this used to also print to file. But I found that it was mostly
 * useless to do this.
 *
 * ************************************************************************** */

void
print_optical_depths (SightLines_t *inclinations, int n_inclinations, Edges_t edges[], int n_edges, double *optical_depth,
                      double *column_density)
{
  int i, j;
  int len, c_linelen;
  char str[LINELENGTH];
  const int MAX_COL = 120;

  printf ("\nOptical depths along the defined line of sights for domain %i:\n\n", N_DOMAIN);

  for (i = 0; i < n_inclinations; i++)
  {
    if (COLUMN_MODE == COLUMN_MODE_RHO)
    {
      printf ("%-8s: Mass column density     : %3.2e cm^-2\n", inclinations[i].name, column_density[i]);
      printf ("%-8s: Hydrogen column density : %3.2e cm^-2\n", "", column_density[i] * rho2nh);
    }
    else
    {
      printf ("%-8s: %s %i column density    : %3.2e cm^-2\n", inclinations[i].name, ele[ion[COLUMN_MODE_ION_NUMBER].nelem].name,
              ion[COLUMN_MODE_ION_NUMBER].istate, column_density[i]);
    }

    c_linelen = 0;
    for (j = 0; j < n_edges; j++)
    {
      len = snprintf (str, LINELENGTH, "tau_%-9s: %3.2e  ", edges[j].name, optical_depth[i * n_edges + j]);
      if (len < 0)
      {
        errormsg ("error when trying to write to string for output\n");
        exit (EXIT_FAILURE);
      }

      c_linelen += len;
      if (c_linelen > MAX_COL)
      {
        c_linelen = 0;
        printf ("\n");
      }

      printf ("%s", str);
    }
    printf ("\n\n");
  }
}

/* ************************************************************************* */
/**
 * @brief  Write the optical depth spectrum to file.
 *
 * @param[in]  SightLines_t  *inclination  The inclinations angles the optical
 *                                         depth was extracted from.
 * @param[in]  int  n_inclinations         The number of inclination angles
 * @param[in]  double  *tau_spectrum       The optical depth spectrum values
 * @param[in]  double  freq_min            The starting frequency of the
 *                                         spectrum
 * @param[in]  double  dfreq               The frequency spacing of the spectrum
 *
 * @details
 *
 * Simply writes the optical depth spectra to the file named
 * root.tau_spec.
 *
 * ************************************************************************** */

void
write_optical_depth_spectrum (SightLines_t *inclinations, int n_inclinations, double *tau_spectrum, double freq_min, double d_freq)
{
  int i, j;
  double c_wavelength, c_frequency;
  char filename[LINELENGTH];
  FILE *fp;

  int len = snprintf (filename, LINELENGTH, "%s.spec_tau", files.root);
  if (len < 0)
  {
    errormsg ("error when creating filename string\n");
    exit (EXIT_FAILURE);
  }

  fp = fopen (filename, "w");
  if (fp == NULL)
  {
    errormsg ("unable to open %s in write mode\n", filename);
    exit (EXIT_FAILURE);
  }

  write_generic_file_header (fp);
  fprintf (fp, "%-15s %-15s ", "Freq.", "Lambda");
  for (i = 0; i < n_inclinations; i++)
  {
    fprintf (fp, "%-15s ", inclinations[i].name);
  }
  fprintf (fp, "\n");

  c_frequency = log10 (freq_min);
  for (i = 0; i < NUM_FREQUENCY_BINS; i++)
  {
    c_wavelength = VLIGHT / pow (10, c_frequency) / ANGSTROM;
    fprintf (fp, "%-15e %-15e ", pow (10, c_frequency), c_wavelength);

    for (j = 0; j < n_inclinations; j++)
    {
      fprintf (fp, "%-15e ", tau_spectrum[j * NUM_FREQUENCY_BINS + i]);
    }

    fprintf (fp, "\n");
    c_frequency += d_freq;
  }

  if (fclose (fp))
  {
    errormsg ("unable to close %s, output may be unfinished!\n", filename);
  }
}

/* ************************************************************************* */
/**
 * @brief  Write the photosphere location points to file.
 *
 * @param[in]  Positions_t *positions  An array of positions
 * @param[in]  int n_angles            The number of angles which were used
 *                                     in the photosphere calculation.
 *
 * @details
 *
 * This uses the Positions_t type which keeps track of the x, y and z locations
 * of where a photon reached a cumulative optical depth of TAU_DEPTH.
 *
 * ************************************************************************** */

void
write_photosphere_location_to_file (Positions_t *positions, int n_angles)
{
  int i;
  double pos1d[3];
  char filename[LINELENGTH];
  FILE *fp;

  int len = snprintf (filename, LINELENGTH, "%s.photosphere", files.root);
  if (len < 0)
  {
    errormsg ("error when creating filename string\n");
    exit (EXIT_FAILURE);
  }

  fp = fopen (filename, "w");
  if (fp == NULL)
  {
    errormsg ("unable to open %s in write mode\n", filename);
    exit (EXIT_FAILURE);
  }

  write_generic_file_header (fp);
  fprintf (fp, "# Electron scatter photosphere locations for tau_es = %f\n#\n", TAU_DEPTH);

  if (zdom[N_DOMAIN].coord_type != SPHERICAL)
  {
    fprintf (fp, "%-15s %-15s %-15s\n", "x", "y", "z");
  }
  else
  {
    fprintf (fp, "%-15s\n", "r");
  }

  for (i = 0; i < n_angles; i++)
  {
    if (zdom[N_DOMAIN].coord_type != SPHERICAL)
    {
      fprintf (fp, "%-15e %-15e %-15e\n", positions[i].x, positions[i].y, positions[i].z);
    }
    else
    {
      pos1d[0] = positions[0].x;
      pos1d[1] = positions[0].y;
      pos1d[2] = positions[0].z;
      fprintf (fp, "%-15e\n", length (pos1d));
    }
  }

  if (fclose (fp))
  {
    errormsg ("unable to close %s, output may be unfinished!\n", filename);
  }
}
