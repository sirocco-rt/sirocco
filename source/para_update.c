
/***********************************************************/
/** @file  para_update.c
 * @author ksl, jm
 * @date   January, 2018
 *
 * @brief  routines for communicating MC estimators and spectra between MPI threads.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/** 
 * @brief      communicates the MC estimators between tasks
 *
 * @details
 * communicates the MC estimators between tasks relating to 
 * spectral models, heating and cooling and cell diagnostics like IP. 
 * In the case of some variables, the quantities are maxima and minima so the 
 * flag MPI_MAX or MPI_MIN is used in MPI_Reduce. For summed
 * quantities like heating we use MPI_SUM.
 *  
 * This routine should only do anything if the MPI_ON flag was present 
 * in compilation. It communicates all the information
 * required for the spectral model ionization scheme, and 
 * also heating and cooling quantities in cells.
 **********************************************************/

int
communicate_estimators_para ()
{
#ifdef MPI_ON                   // these routines should only be called anyway in parallel but we need these to compile

  int mpi_i, mpi_j;
  double *maxfreqhelper, *maxfreqhelper2;
  /*NSH 131213 the next line introduces new helper arrays for the max and min frequencies in bands */
  double *maxbandfreqhelper, *maxbandfreqhelper2, *minbandfreqhelper, *minbandfreqhelper2;
  double *redhelper, *redhelper2, *qdisk_helper, *qdisk_helper2;
  double *ion_helper, *ion_helper2;
  double *inner_ion_helper, *inner_ion_helper2;
  double *flux_helper, *flux_helper2;

//OLD  int nspec, size_of_commbuffer;
  int size_of_commbuffer;

  int *iredhelper, *iredhelper2, *iqdisk_helper, *iqdisk_helper2;
  // int size_of_helpers;
  int plasma_double_helpers, plasma_int_helpers;

  /* The size of the helper array for doubles. We transmit 10 numbers 
     for each cell, plus three arrays, each of length NXBANDS */

  plasma_double_helpers = (39 + 3 * NXBANDS) * NPLASMA;

  /* The size of the helper array for integers. We transmit 7 numbers 
     for each cell, plus one array of length NXBANDS */
  plasma_int_helpers = (7 + NXBANDS) * NPLASMA;


  maxfreqhelper = calloc (sizeof (double), NPLASMA);
  maxfreqhelper2 = calloc (sizeof (double), NPLASMA);
  /* NSH 131213 - allocate memory for the band limited max and min frequencies */
  maxbandfreqhelper = calloc (sizeof (double), NPLASMA * NXBANDS);
  maxbandfreqhelper2 = calloc (sizeof (double), NPLASMA * NXBANDS);
  minbandfreqhelper = calloc (sizeof (double), NPLASMA * NXBANDS);
  minbandfreqhelper2 = calloc (sizeof (double), NPLASMA * NXBANDS);
  redhelper = calloc (sizeof (double), plasma_double_helpers);
  redhelper2 = calloc (sizeof (double), plasma_double_helpers);

  ion_helper = calloc (sizeof (double), NPLASMA * nions);
  ion_helper2 = calloc (sizeof (double), NPLASMA * nions);
  inner_ion_helper = calloc (sizeof (double), NPLASMA * n_inner_tot);
  inner_ion_helper2 = calloc (sizeof (double), NPLASMA * n_inner_tot);
  /* JM -- added routine to average the qdisk quantities. The 2 is because
     we only have two doubles to worry about (heat and ave_freq) and 
     two integers (nhit and nphot) */
  qdisk_helper = calloc (sizeof (double), NRINGS * 2);
  qdisk_helper2 = calloc (sizeof (double), NRINGS * 2);

  flux_helper = calloc (sizeof (double), NPLASMA * NFLUX_ANGLES * 3);
  flux_helper2 = calloc (sizeof (double), NPLASMA * NFLUX_ANGLES * 3);

  MPI_Barrier (MPI_COMM_WORLD);
  // the following blocks gather all the estimators to the zeroth (Master) thread


  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    maxfreqhelper[mpi_i] = plasmamain[mpi_i].max_freq;
    redhelper[mpi_i] = plasmamain[mpi_i].j / np_mpi_global;
    redhelper[mpi_i + NPLASMA] = plasmamain[mpi_i].ave_freq / np_mpi_global;
    redhelper[mpi_i + 2 * NPLASMA] = plasmamain[mpi_i].cool_tot / np_mpi_global;
    redhelper[mpi_i + 3 * NPLASMA] = plasmamain[mpi_i].heat_tot / np_mpi_global;
    redhelper[mpi_i + 4 * NPLASMA] = plasmamain[mpi_i].heat_lines / np_mpi_global;
    redhelper[mpi_i + 5 * NPLASMA] = plasmamain[mpi_i].heat_ff / np_mpi_global;
    redhelper[mpi_i + 6 * NPLASMA] = plasmamain[mpi_i].heat_comp / np_mpi_global;
    redhelper[mpi_i + 7 * NPLASMA] = plasmamain[mpi_i].heat_ind_comp / np_mpi_global;
    redhelper[mpi_i + 8 * NPLASMA] = plasmamain[mpi_i].heat_photo / np_mpi_global;
    redhelper[mpi_i + 9 * NPLASMA] = plasmamain[mpi_i].ip / np_mpi_global;
    redhelper[mpi_i + 10 * NPLASMA] = plasmamain[mpi_i].j_direct / np_mpi_global;
    redhelper[mpi_i + 11 * NPLASMA] = plasmamain[mpi_i].j_scatt / np_mpi_global;
    redhelper[mpi_i + 12 * NPLASMA] = plasmamain[mpi_i].ip_direct / np_mpi_global;
    redhelper[mpi_i + 13 * NPLASMA] = plasmamain[mpi_i].ip_scatt / np_mpi_global;
    redhelper[mpi_i + 14 * NPLASMA] = plasmamain[mpi_i].heat_auger / np_mpi_global;
    redhelper[mpi_i + 15 * NPLASMA] = plasmamain[mpi_i].rad_force_es[0] / np_mpi_global;
    redhelper[mpi_i + 16 * NPLASMA] = plasmamain[mpi_i].rad_force_es[1] / np_mpi_global;
    redhelper[mpi_i + 17 * NPLASMA] = plasmamain[mpi_i].rad_force_es[2] / np_mpi_global;
    redhelper[mpi_i + 18 * NPLASMA] = plasmamain[mpi_i].rad_force_es[3] / np_mpi_global;
    redhelper[mpi_i + 19 * NPLASMA] = plasmamain[mpi_i].F_vis[0] / np_mpi_global;
    redhelper[mpi_i + 20 * NPLASMA] = plasmamain[mpi_i].F_vis[1] / np_mpi_global;
    redhelper[mpi_i + 21 * NPLASMA] = plasmamain[mpi_i].F_vis[2] / np_mpi_global;
    redhelper[mpi_i + 22 * NPLASMA] = plasmamain[mpi_i].F_vis[3] / np_mpi_global;
    redhelper[mpi_i + 23 * NPLASMA] = plasmamain[mpi_i].F_UV[0] / np_mpi_global;
    redhelper[mpi_i + 24 * NPLASMA] = plasmamain[mpi_i].F_UV[1] / np_mpi_global;
    redhelper[mpi_i + 25 * NPLASMA] = plasmamain[mpi_i].F_UV[2] / np_mpi_global;
    redhelper[mpi_i + 26 * NPLASMA] = plasmamain[mpi_i].F_UV[3] / np_mpi_global;
    redhelper[mpi_i + 27 * NPLASMA] = plasmamain[mpi_i].F_Xray[0] / np_mpi_global;
    redhelper[mpi_i + 28 * NPLASMA] = plasmamain[mpi_i].F_Xray[1] / np_mpi_global;
    redhelper[mpi_i + 29 * NPLASMA] = plasmamain[mpi_i].F_Xray[2] / np_mpi_global;
    redhelper[mpi_i + 30 * NPLASMA] = plasmamain[mpi_i].F_Xray[3] / np_mpi_global;
    redhelper[mpi_i + 31 * NPLASMA] = plasmamain[mpi_i].rad_force_bf[0] / np_mpi_global;
    redhelper[mpi_i + 32 * NPLASMA] = plasmamain[mpi_i].rad_force_bf[1] / np_mpi_global;
    redhelper[mpi_i + 33 * NPLASMA] = plasmamain[mpi_i].rad_force_bf[2] / np_mpi_global;
    redhelper[mpi_i + 34 * NPLASMA] = plasmamain[mpi_i].rad_force_bf[3] / np_mpi_global;
    redhelper[mpi_i + 35 * NPLASMA] = plasmamain[mpi_i].rad_force_ff[0] / np_mpi_global;
    redhelper[mpi_i + 36 * NPLASMA] = plasmamain[mpi_i].rad_force_ff[1] / np_mpi_global;
    redhelper[mpi_i + 37 * NPLASMA] = plasmamain[mpi_i].rad_force_ff[2] / np_mpi_global;
    redhelper[mpi_i + 38 * NPLASMA] = plasmamain[mpi_i].rad_force_ff[3] / np_mpi_global;

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      redhelper[mpi_i + (39 + mpi_j) * NPLASMA] = plasmamain[mpi_i].xj[mpi_j] / np_mpi_global;
      redhelper[mpi_i + (39 + NXBANDS + mpi_j) * NPLASMA] = plasmamain[mpi_i].xave_freq[mpi_j] / np_mpi_global;
      redhelper[mpi_i + (39 + 2 * NXBANDS + mpi_j) * NPLASMA] = plasmamain[mpi_i].xsd_freq[mpi_j] / np_mpi_global;
      /* 131213 NSH populate the band limited min and max frequency arrays */
      maxbandfreqhelper[mpi_i * NXBANDS + mpi_j] = plasmamain[mpi_i].fmax[mpi_j];
      minbandfreqhelper[mpi_i * NXBANDS + mpi_j] = plasmamain[mpi_i].fmin[mpi_j];
    }
    for (mpi_j = 0; mpi_j < nions; mpi_j++)
    {
      ion_helper[mpi_i * nions + mpi_j] = plasmamain[mpi_i].ioniz[mpi_j] / np_mpi_global;
    }
    for (mpi_j = 0; mpi_j < n_inner_tot; mpi_j++)
    {
      inner_ion_helper[mpi_i * n_inner_tot + mpi_j] = plasmamain[mpi_i].inner_ioniz[mpi_j] / np_mpi_global;
    }
    for (mpi_j = 0; mpi_j < NFLUX_ANGLES; mpi_j++)
    {
      flux_helper[mpi_i * (3 * NFLUX_ANGLES) + mpi_j] = plasmamain[mpi_i].F_UV_ang_x[mpi_j] / np_mpi_global;
      flux_helper[mpi_i * (3 * NFLUX_ANGLES) + NFLUX_ANGLES + mpi_j] = plasmamain[mpi_i].F_UV_ang_y[mpi_j] / np_mpi_global;
      flux_helper[mpi_i * (3 * NFLUX_ANGLES) + 2 * NFLUX_ANGLES + mpi_j] = plasmamain[mpi_i].F_UV_ang_z[mpi_j] / np_mpi_global;
    }
  }

  for (mpi_i = 0; mpi_i < NRINGS; mpi_i++)
  {
    qdisk_helper[mpi_i] = qdisk.heat[mpi_i] / np_mpi_global;
    qdisk_helper[mpi_i + NRINGS] = qdisk.ave_freq[mpi_i] / np_mpi_global;
  }


  /* 131213 NSH communicate the min and max band frequencies these use MPI_MIN or MPI_MAX */
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Reduce (minbandfreqhelper, minbandfreqhelper2, NPLASMA * NXBANDS, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (maxbandfreqhelper, maxbandfreqhelper2, NPLASMA * NXBANDS, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (maxfreqhelper, maxfreqhelper2, NPLASMA, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (redhelper, redhelper2, plasma_double_helpers, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  //XXXX Duplicate MPI_Reduce (redhelper, redhelper2, plasma_double_helpers, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (flux_helper, flux_helper2, NPLASMA * 3 * NFLUX_ANGLES, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce (ion_helper, ion_helper2, NPLASMA * nions, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (inner_ion_helper, inner_ion_helper2, NPLASMA * n_inner_tot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /* JM 1607 -- sum up the qdisk values */
  MPI_Reduce (qdisk_helper, qdisk_helper2, 2 * NRINGS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank_global == 0)
  {

    Log_parallel ("Zeroth thread successfully received the normalised estimators. About to broadcast.\n");
  }

  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Bcast (redhelper2, plasma_double_helpers, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (maxfreqhelper2, NPLASMA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  /* 131213 NSH Send out the global min and max band limited frequencies to all threads */
  MPI_Bcast (minbandfreqhelper2, NPLASMA * NXBANDS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (maxbandfreqhelper2, NPLASMA * NXBANDS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (flux_helper2, NPLASMA * 3 * NFLUX_ANGLES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (ion_helper2, NPLASMA * nions, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (inner_ion_helper2, NPLASMA * n_inner_tot, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* JM 1607 -- send out the qdisk values to all threads */
  MPI_Bcast (qdisk_helper2, NRINGS, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    plasmamain[mpi_i].max_freq = maxfreqhelper2[mpi_i];
    plasmamain[mpi_i].j = redhelper2[mpi_i];
    plasmamain[mpi_i].ave_freq = redhelper2[mpi_i + NPLASMA];
    plasmamain[mpi_i].cool_tot = redhelper2[mpi_i + 2 * NPLASMA];
    plasmamain[mpi_i].heat_tot = redhelper2[mpi_i + 3 * NPLASMA];
    plasmamain[mpi_i].heat_lines = redhelper2[mpi_i + 4 * NPLASMA];
    plasmamain[mpi_i].heat_ff = redhelper2[mpi_i + 5 * NPLASMA];
    plasmamain[mpi_i].heat_comp = redhelper2[mpi_i + 6 * NPLASMA];
    plasmamain[mpi_i].heat_ind_comp = redhelper2[mpi_i + 7 * NPLASMA];
    plasmamain[mpi_i].heat_photo = redhelper2[mpi_i + 8 * NPLASMA];
    plasmamain[mpi_i].ip = redhelper2[mpi_i + 9 * NPLASMA];
    plasmamain[mpi_i].j_direct = redhelper2[mpi_i + 10 * NPLASMA];
    plasmamain[mpi_i].j_scatt = redhelper2[mpi_i + 11 * NPLASMA];
    plasmamain[mpi_i].ip_direct = redhelper2[mpi_i + 12 * NPLASMA];
    plasmamain[mpi_i].ip_scatt = redhelper2[mpi_i + 13 * NPLASMA];
    plasmamain[mpi_i].heat_auger = redhelper2[mpi_i + 14 * NPLASMA];
    plasmamain[mpi_i].rad_force_es[0] = redhelper2[mpi_i + 15 * NPLASMA];
    plasmamain[mpi_i].rad_force_es[1] = redhelper2[mpi_i + 16 * NPLASMA];
    plasmamain[mpi_i].rad_force_es[2] = redhelper2[mpi_i + 17 * NPLASMA];
    plasmamain[mpi_i].rad_force_es[3] = redhelper2[mpi_i + 18 * NPLASMA];
    plasmamain[mpi_i].F_vis[0] = redhelper2[mpi_i + 19 * NPLASMA];
    plasmamain[mpi_i].F_vis[1] = redhelper2[mpi_i + 20 * NPLASMA];
    plasmamain[mpi_i].F_vis[2] = redhelper2[mpi_i + 21 * NPLASMA];
    plasmamain[mpi_i].F_vis[3] = redhelper2[mpi_i + 22 * NPLASMA];
    plasmamain[mpi_i].F_UV[0] = redhelper2[mpi_i + 23 * NPLASMA];
    plasmamain[mpi_i].F_UV[1] = redhelper2[mpi_i + 24 * NPLASMA];
    plasmamain[mpi_i].F_UV[2] = redhelper2[mpi_i + 25 * NPLASMA];
    plasmamain[mpi_i].F_UV[3] = redhelper2[mpi_i + 26 * NPLASMA];
    plasmamain[mpi_i].F_Xray[0] = redhelper2[mpi_i + 27 * NPLASMA];
    plasmamain[mpi_i].F_Xray[1] = redhelper2[mpi_i + 28 * NPLASMA];
    plasmamain[mpi_i].F_Xray[2] = redhelper2[mpi_i + 29 * NPLASMA];
    plasmamain[mpi_i].F_Xray[3] = redhelper2[mpi_i + 30 * NPLASMA];
    plasmamain[mpi_i].rad_force_bf[0] = redhelper2[mpi_i + 31 * NPLASMA];
    plasmamain[mpi_i].rad_force_bf[1] = redhelper2[mpi_i + 32 * NPLASMA];
    plasmamain[mpi_i].rad_force_bf[2] = redhelper2[mpi_i + 33 * NPLASMA];
    plasmamain[mpi_i].rad_force_bf[3] = redhelper2[mpi_i + 34 * NPLASMA];
    plasmamain[mpi_i].rad_force_ff[0] = redhelper2[mpi_i + 35 * NPLASMA];
    plasmamain[mpi_i].rad_force_ff[1] = redhelper2[mpi_i + 36 * NPLASMA];
    plasmamain[mpi_i].rad_force_ff[2] = redhelper2[mpi_i + 37 * NPLASMA];
    plasmamain[mpi_i].rad_force_ff[3] = redhelper2[mpi_i + 38 * NPLASMA];

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      plasmamain[mpi_i].xj[mpi_j] = redhelper2[mpi_i + (39 + mpi_j) * NPLASMA];
      plasmamain[mpi_i].xave_freq[mpi_j] = redhelper2[mpi_i + (39 + NXBANDS + mpi_j) * NPLASMA];
      plasmamain[mpi_i].xsd_freq[mpi_j] = redhelper2[mpi_i + (39 + NXBANDS * 2 + mpi_j) * NPLASMA];

      /* 131213 NSH And unpack the min and max banded frequencies to the plasma array */
      plasmamain[mpi_i].fmax[mpi_j] = maxbandfreqhelper2[mpi_i * NXBANDS + mpi_j];
      plasmamain[mpi_i].fmin[mpi_j] = minbandfreqhelper2[mpi_i * NXBANDS + mpi_j];
    }
    for (mpi_j = 0; mpi_j < nions; mpi_j++)
    {
      plasmamain[mpi_i].ioniz[mpi_j] = ion_helper2[mpi_i * nions + mpi_j];
    }
    for (mpi_j = 0; mpi_j < n_inner_tot; mpi_j++)
    {
      plasmamain[mpi_i].inner_ioniz[mpi_j] = inner_ion_helper2[mpi_i * n_inner_tot + mpi_j];
    }
    for (mpi_j = 0; mpi_j < NFLUX_ANGLES; mpi_j++)
    {
      plasmamain[mpi_i].F_UV_ang_x[mpi_j] = flux_helper2[mpi_i * (3 * NFLUX_ANGLES) + mpi_j];
      plasmamain[mpi_i].F_UV_ang_y[mpi_j] = flux_helper2[mpi_i * (3 * NFLUX_ANGLES) + NFLUX_ANGLES + mpi_j];
      plasmamain[mpi_i].F_UV_ang_z[mpi_j] = flux_helper2[mpi_i * (3 * NFLUX_ANGLES) + 2 * NFLUX_ANGLES + mpi_j];
    }
  }


  for (mpi_i = 0; mpi_i < NRINGS; mpi_i++)
  {
    qdisk.heat[mpi_i] = qdisk_helper2[mpi_i];
    qdisk.ave_freq[mpi_i] = qdisk_helper2[mpi_i + NRINGS];
  }
  Log_silent ("Thread %d happy after broadcast.\n", rank_global);

  /* now we've done all the doubles so we can free their helper arrays */
  free (qdisk_helper);
  free (qdisk_helper2);
  free (redhelper);
  free (redhelper2);
  free (maxfreqhelper);
  free (maxfreqhelper2);
  free (maxbandfreqhelper);
  free (maxbandfreqhelper2);
  free (minbandfreqhelper);
  free (minbandfreqhelper2);

  free (ion_helper);
  free (ion_helper2);

  free (inner_ion_helper);
  free (inner_ion_helper2);
  /* allocate the integer helper arrays, set a barrier, then do all the integers. */
  iqdisk_helper = calloc (sizeof (int), NRINGS * 2);
  iqdisk_helper2 = calloc (sizeof (int), NRINGS * 2);
  iredhelper = calloc (sizeof (int), plasma_int_helpers);
  iredhelper2 = calloc (sizeof (int), plasma_int_helpers);
  MPI_Barrier (MPI_COMM_WORLD);

  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    iredhelper[mpi_i] = plasmamain[mpi_i].ntot;
    iredhelper[mpi_i + NPLASMA] = plasmamain[mpi_i].ntot_star;
    iredhelper[mpi_i + 2 * NPLASMA] = plasmamain[mpi_i].ntot_bl;
    iredhelper[mpi_i + 3 * NPLASMA] = plasmamain[mpi_i].ntot_disk;
    iredhelper[mpi_i + 4 * NPLASMA] = plasmamain[mpi_i].ntot_wind;
    iredhelper[mpi_i + 5 * NPLASMA] = plasmamain[mpi_i].ntot_agn;
    iredhelper[mpi_i + 6 * NPLASMA] = plasmamain[mpi_i].nioniz;

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      iredhelper[mpi_i + (7 + mpi_j) * NPLASMA] = plasmamain[mpi_i].nxtot[mpi_j];
    }
  }

  for (mpi_i = 0; mpi_i < NRINGS; mpi_i++)
  {
    iqdisk_helper[mpi_i] = qdisk.nphot[mpi_i];
    iqdisk_helper[mpi_i + NRINGS] = qdisk.nhit[mpi_i];
  }

  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Reduce (iredhelper, iredhelper2, plasma_int_helpers, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (iqdisk_helper, iqdisk_helper2, 2 * NRINGS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank_global == 0)
  {
    Log_parallel ("Zeroth thread successfully received the integer sum. About to broadcast.\n");
  }

  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Bcast (iredhelper2, plasma_int_helpers, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (iqdisk_helper2, NRINGS, MPI_INT, 0, MPI_COMM_WORLD);


  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    plasmamain[mpi_i].ntot = iredhelper2[mpi_i];
    plasmamain[mpi_i].ntot_star = iredhelper2[mpi_i + NPLASMA];
    plasmamain[mpi_i].ntot_bl = iredhelper2[mpi_i + 2 * NPLASMA];
    plasmamain[mpi_i].ntot_disk = iredhelper2[mpi_i + 3 * NPLASMA];
    plasmamain[mpi_i].ntot_wind = iredhelper2[mpi_i + 4 * NPLASMA];
    plasmamain[mpi_i].ntot_agn = iredhelper2[mpi_i + 5 * NPLASMA];
    plasmamain[mpi_i].nioniz = iredhelper2[mpi_i + 6 * NPLASMA];

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      plasmamain[mpi_i].nxtot[mpi_j] = iredhelper2[mpi_i + (7 + mpi_j) * NPLASMA];
    }
  }

  for (mpi_i = 0; mpi_i < NRINGS; mpi_i++)
  {
    qdisk.nphot[mpi_i] = iqdisk_helper2[mpi_i];
    qdisk.nhit[mpi_i] = iqdisk_helper2[mpi_i + NRINGS];
  }

  free (iredhelper);
  free (iredhelper2);
  free (iqdisk_helper);
  free (iqdisk_helper2);


/* Now during ionization cycles, process the cell spectra */

  /* The size of the commbuffers need to be the number of spectra x the length of each */

  if (geo.ioniz_or_extract == CYCLE_IONIZ)
  {
    MPI_Barrier (MPI_COMM_WORLD);
    size_of_commbuffer = NPLASMA * NBINS_IN_CELL_SPEC;
//OLD    nspec = NPLASMA;

    redhelper = calloc (sizeof (double), size_of_commbuffer);
    redhelper2 = calloc (sizeof (double), size_of_commbuffer);

    for (mpi_i = 0; mpi_i < NBINS_IN_CELL_SPEC; mpi_i++)
    {
      for (mpi_j = 0; mpi_j < NPLASMA; mpi_j++)
      {
        redhelper[mpi_i * NPLASMA + mpi_j] = plasmamain[mpi_j].cell_spec_flux[mpi_i] / np_mpi_global;

      }
    }



    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Reduce (redhelper, redhelper2, size_of_commbuffer, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast (redhelper2, size_of_commbuffer, MPI_DOUBLE, 0, MPI_COMM_WORLD);



    for (mpi_i = 0; mpi_i < NBINS_IN_CELL_SPEC; mpi_i++)
    {
      for (mpi_j = 0; mpi_j < NPLASMA; mpi_j++)
      {
        plasmamain[mpi_j].cell_spec_flux[mpi_i] = redhelper2[mpi_i * NPLASMA + mpi_j];

      }
    }
    MPI_Barrier (MPI_COMM_WORLD);



    free (redhelper);
    free (redhelper2);


  }


#endif
  return (0);
}


/**********************************************************/
/**
 * @brief sum up the synthetic and cell spectra between threads.
 *
 * @details
 * sum up the synthetic spectra between threads. Does an
 * MPI_Reduce then an MPI_Bcast for each element of the
 * linear and log spectra arrays (xxspec)
 *
 **********************************************************/


int
gather_spectra_para ()
{
#ifdef MPI_ON
  double *redhelper, *redhelper2;
  int mpi_i, mpi_j;
  int size_of_commbuffer, nspec;

  if (geo.ioniz_or_extract == CYCLE_EXTRACT)
  {
    nspec = MSPEC + geo.nangles;
  }
  else
  {
    nspec = MSPEC;
  }

  size_of_commbuffer = 4 * nspec * NWAVE_MAX;   //we need space for all 4 separate spectra we are normalizing

  redhelper = calloc (sizeof (double), size_of_commbuffer);
  redhelper2 = calloc (sizeof (double), size_of_commbuffer);

  for (mpi_i = 0; mpi_i < NWAVE_MAX; mpi_i++)
  {
    for (mpi_j = 0; mpi_j < nspec; mpi_j++)
    {
      redhelper[mpi_i * nspec + mpi_j] = xxspec[mpi_j].f[mpi_i] / np_mpi_global;
      redhelper[mpi_i * nspec + mpi_j + (NWAVE_MAX * nspec)] = xxspec[mpi_j].lf[mpi_i] / np_mpi_global;
      redhelper[mpi_i * nspec + mpi_j + (2 * NWAVE_MAX * nspec)] = xxspec[mpi_j].f_wind[mpi_i] / np_mpi_global;
      redhelper[mpi_i * nspec + mpi_j + (3 * NWAVE_MAX * nspec)] = xxspec[mpi_j].lf_wind[mpi_i] / np_mpi_global;
    }
  }

  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Reduce (redhelper, redhelper2, size_of_commbuffer, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast (redhelper2, size_of_commbuffer, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for (mpi_i = 0; mpi_i < NWAVE_MAX; mpi_i++)
  {
    for (mpi_j = 0; mpi_j < nspec; mpi_j++)
    {
      xxspec[mpi_j].f[mpi_i] = redhelper2[mpi_i * nspec + mpi_j];
      xxspec[mpi_j].lf[mpi_i] = redhelper2[mpi_i * nspec + mpi_j + (NWAVE_MAX * nspec)];
      xxspec[mpi_j].f_wind[mpi_i] = redhelper2[mpi_i * nspec + mpi_j + (2 * NWAVE_MAX * nspec)];
      xxspec[mpi_j].lf_wind[mpi_i] = redhelper2[mpi_i * nspec + mpi_j + (3 * NWAVE_MAX * nspec)];
    }
  }
  MPI_Barrier (MPI_COMM_WORLD);

  free (redhelper);
  free (redhelper2);
#endif

  return (0);
}





/**********************************************************/
/**
 * @brief
 *
 * @details averages the macro-atom estimators between tasks using MPI_Reduce.
 *   It should only be called if the MPI_ON flag was present
 *   in compilation, and returns 0 immediately if no macro atom levels.
 *   This should probably be improved by working out exactly
 *   what is needed in simple-ion only mode.
 *
 **********************************************************/

int
communicate_matom_estimators_para ()
{
#ifdef MPI_ON                   // these routines should only be called anyway in parallel but we need these to compile

  int n, mpi_i;
  double *gamma_helper, *alpha_helper;
  double *level_helper, *cell_helper, *jbar_helper;
  double *gamma_helper2, *alpha_helper2;
  double *level_helper2, *cell_helper2, *jbar_helper2;
  double *cooling_bf_helper, *cooling_bb_helper;
  double *cooling_bf_helper2, *cooling_bb_helper2;

  if (nlevels_macro == 0 && geo.nmacro == 0)
  {
    /* in this case no space would have been allocated for macro-atom estimators */
    Log ("No need to communicate matom estimators as no macro-atoms!\n");
    return (0);
  }



  /* allocate helper arrays for the estimators we want to communicate */
  /* the sizes of these arrays should match the allocation in calloc_estimators in gridwind.c */
  /* we need two arrays for each set of variables. Note that we stick all estimators of
     the same size in the same helper array */
  /* could put some error conditions here to check memory allocation worked */
  jbar_helper = calloc (sizeof (double), NPLASMA * size_Jbar_est);
  gamma_helper = calloc (sizeof (double), NPLASMA * 4 * size_gamma_est);
  alpha_helper = calloc (sizeof (double), NPLASMA * 2 * size_alpha_est);
  level_helper = calloc (sizeof (double), NPLASMA * nlevels_macro);
  cell_helper = calloc (sizeof (double), 8 * NPLASMA);
  cooling_bf_helper = calloc (sizeof (double), NPLASMA * 2 * nphot_total);
  cooling_bb_helper = calloc (sizeof (double), NPLASMA * nlines);

  jbar_helper2 = calloc (sizeof (double), NPLASMA * size_Jbar_est);
  gamma_helper2 = calloc (sizeof (double), NPLASMA * 4 * size_gamma_est);
  alpha_helper2 = calloc (sizeof (double), NPLASMA * 2 * size_alpha_est);
  level_helper2 = calloc (sizeof (double), NPLASMA * nlevels_macro);
  cell_helper2 = calloc (sizeof (double), 8 * NPLASMA);
  cooling_bf_helper2 = calloc (sizeof (double), NPLASMA * 2 * nphot_total);
  cooling_bb_helper2 = calloc (sizeof (double), NPLASMA * nlines);


  /* set an mpi barrier before we start */
  MPI_Barrier (MPI_COMM_WORLD);



  /* now we loop through each cell and copy the values of our variables
     into our helper arrays */
  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    /* one kpkt_abs quantity per cell */
    cell_helper[mpi_i] = plasmamain[mpi_i].kpkt_abs / np_mpi_global;

    /* each of the cooling sums and normalisations also have one quantity per cell */
    cell_helper[mpi_i + NPLASMA] = macromain[mpi_i].cooling_normalisation / np_mpi_global;
    cell_helper[mpi_i + 2 * NPLASMA] = macromain[mpi_i].cooling_bftot / np_mpi_global;
    cell_helper[mpi_i + 3 * NPLASMA] = macromain[mpi_i].cooling_bf_coltot / np_mpi_global;
    cell_helper[mpi_i + 4 * NPLASMA] = macromain[mpi_i].cooling_bbtot / np_mpi_global;
    cell_helper[mpi_i + 5 * NPLASMA] = macromain[mpi_i].cooling_ff / np_mpi_global;
    cell_helper[mpi_i + 6 * NPLASMA] = macromain[mpi_i].cooling_ff_lofreq / np_mpi_global;
    cell_helper[mpi_i + 7 * NPLASMA] = macromain[mpi_i].cooling_adiabatic / np_mpi_global;



    for (n = 0; n < nlevels_macro; n++)
    {
      level_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].matom_abs[n] / np_mpi_global;
    }

    for (n = 0; n < size_Jbar_est; n++)
    {
      jbar_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].jbar[n] / np_mpi_global;
    }

    for (n = 0; n < size_gamma_est; n++)
    {
      gamma_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].alpha_st[n] / np_mpi_global;
      gamma_helper[mpi_i + ((n + size_gamma_est) * NPLASMA)] = macromain[mpi_i].alpha_st_e[n] / np_mpi_global;
      gamma_helper[mpi_i + ((n + 2 * size_gamma_est) * NPLASMA)] = macromain[mpi_i].gamma[n] / np_mpi_global;
      gamma_helper[mpi_i + ((n + 3 * size_gamma_est) * NPLASMA)] = macromain[mpi_i].gamma_e[n] / np_mpi_global;
    }

    for (n = 0; n < size_alpha_est; n++)
    {
      alpha_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].recomb_sp[n] / np_mpi_global;
      alpha_helper[mpi_i + ((n + size_alpha_est) * NPLASMA)] = macromain[mpi_i].recomb_sp_e[n] / np_mpi_global;
    }

    for (n = 0; n < nphot_total; n++)
    {
      cooling_bf_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].cooling_bf[n] / np_mpi_global;
      cooling_bf_helper[mpi_i + ((n + nphot_total) * NPLASMA)] = macromain[mpi_i].cooling_bf_col[n] / np_mpi_global;
    }

    for (n = 0; n < nlines; n++)
    {
      cooling_bb_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].cooling_bb[n] / np_mpi_global;
    }
  }

  /* because in the above loop we have already divided by number of processes, we can now do a sum
     with MPI_Reduce, passing it MPI_SUM as an argument. This will give us the mean across threads */
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Reduce (cell_helper, cell_helper2, NPLASMA * 8, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (level_helper, level_helper2, NPLASMA * nlevels_macro, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (jbar_helper, jbar_helper2, NPLASMA * size_Jbar_est, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (gamma_helper, gamma_helper2, NPLASMA * 4 * size_gamma_est, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (alpha_helper, alpha_helper2, NPLASMA * 2 * size_alpha_est, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (cooling_bf_helper, cooling_bf_helper2, NPLASMA * 2 * nphot_total, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (cooling_bb_helper, cooling_bb_helper2, NPLASMA * nlines, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


  if (rank_global == 0)
  {
    Log_parallel ("Zeroth thread successfully received the macro-atom estimators. About to broadcast.\n");
  }

  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Bcast (cell_helper2, NPLASMA * 8, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (level_helper2, NPLASMA * nlevels_macro, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (jbar_helper2, NPLASMA * size_Jbar_est, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (gamma_helper2, NPLASMA * 4 * size_gamma_est, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (alpha_helper2, NPLASMA * 2 * size_alpha_est, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (cooling_bf_helper2, NPLASMA * 2 * nphot_total, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (cooling_bb_helper2, NPLASMA * nlines, MPI_DOUBLE, 0, MPI_COMM_WORLD);




  /* We now need to copy these reduced variables to the plasma structure in each thread */


  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    /* one kpkt_abs quantity per cell */
    plasmamain[mpi_i].kpkt_abs = cell_helper2[mpi_i];

    /* each of the cooling sums and normalisations also have one quantity per cell */
    macromain[mpi_i].cooling_normalisation = cell_helper2[mpi_i + NPLASMA];
    macromain[mpi_i].cooling_bftot = cell_helper2[mpi_i + 2 * NPLASMA];
    macromain[mpi_i].cooling_bf_coltot = cell_helper2[mpi_i + 3 * NPLASMA];
    macromain[mpi_i].cooling_bbtot = cell_helper2[mpi_i + 4 * NPLASMA];
    macromain[mpi_i].cooling_ff = cell_helper2[mpi_i + 5 * NPLASMA];
    macromain[mpi_i].cooling_ff_lofreq = cell_helper2[mpi_i + 6 * NPLASMA];
    macromain[mpi_i].cooling_adiabatic = cell_helper2[mpi_i + 7 * NPLASMA];


    for (n = 0; n < nlevels_macro; n++)
    {
      macromain[mpi_i].matom_abs[n] = level_helper2[mpi_i + (n * NPLASMA)];
    }

    for (n = 0; n < size_Jbar_est; n++)
    {
      macromain[mpi_i].jbar[n] = jbar_helper2[mpi_i + (n * NPLASMA)];
    }

    for (n = 0; n < size_gamma_est; n++)
    {
      macromain[mpi_i].alpha_st[n] = gamma_helper2[mpi_i + (n * NPLASMA)];
      macromain[mpi_i].alpha_st_e[n] = gamma_helper2[mpi_i + ((n + size_gamma_est) * NPLASMA)] / np_mpi_global;
      macromain[mpi_i].gamma[n] = gamma_helper2[mpi_i + ((n + 2 * size_gamma_est) * NPLASMA)];
      macromain[mpi_i].gamma_e[n] = gamma_helper2[mpi_i + ((n + 3 * size_gamma_est) * NPLASMA)];
    }

    for (n = 0; n < size_alpha_est; n++)
    {
      macromain[mpi_i].recomb_sp[n] = alpha_helper2[mpi_i + (n * NPLASMA)];
      macromain[mpi_i].recomb_sp_e[n] = alpha_helper2[mpi_i + ((n + size_alpha_est) * NPLASMA)];
    }

    for (n = 0; n < nphot_total; n++)
    {
      macromain[mpi_i].cooling_bf[n] = cooling_bf_helper2[mpi_i + (n * NPLASMA)];
      macromain[mpi_i].cooling_bf_col[n] = cooling_bf_helper2[mpi_i + ((n + nphot_total) * NPLASMA)];
    }

    for (n = 0; n < nlines; n++)
    {
      macromain[mpi_i].cooling_bb[n] = cooling_bb_helper2[mpi_i + (n * NPLASMA)];
    }
  }


  /* at this stage each thread should have the correctly averaged estimators */
  Log_silent ("Thread %d happy after broadcast.\n", rank_global);



  /* set a barrier and free the memory for the helper arrays */

  MPI_Barrier (MPI_COMM_WORLD);

  free (cell_helper);
  free (level_helper);
  free (jbar_helper);
  free (gamma_helper);
  free (alpha_helper);
  free (cooling_bf_helper);
  free (cooling_bb_helper);

  free (cell_helper2);
  free (level_helper2);
  free (jbar_helper2);
  free (gamma_helper2);
  free (alpha_helper2);
  free (cooling_bf_helper2);
  free (cooling_bb_helper2);
#endif


  return (0);
}

/**********************************************************/
/**
 * @brief helper routine for splitting up tasks in MPI
 * @param   [in]      int   rank       processor rank (typically set from rank_global)
 * @param   [in]      int   ntotal     total number of tasks, e.g. NPLASMA
 * @param   [in]       int   nproc      total number of MPI processors
 * @param   [in,out]  int   *my_nmax   pointer to integer value of first task
 * @param   [in,out]  int   *my_nmax   pointer to integer value of final task
 * @return            int   ndo        number of tasks this thread is working on 

 * @details  For a process with ntotal tasks, 
 * this routine calculates which thread will be given each task.
 * typically ntotal is NPLASMA and the thread is splitting up wind cells.
 * The routine deals with remainders by distributing the remainder over the 
 * threads if the cells do not divide evenly by thread
 **********************************************************/

int
get_parallel_nrange (rank, ntotal, nproc, my_nmin, my_nmax)
     int rank;
     int ntotal;
     int nproc;
     int *my_nmin;
     int *my_nmax;
{
  /* divide the cells between the threads */
  int ndo;
  int num_mpi_cells = floor (ntotal / nproc);

  /* the remainder from the above division */
  int num_mpi_extra = ntotal - (nproc * num_mpi_cells);

  /* this section distributes the remainder over the threads if the cells
     do not divide evenly by thread */
  if (rank < num_mpi_extra)
  {
    *my_nmin = rank_global * (num_mpi_cells + 1);
    *my_nmax = (rank_global + 1) * (num_mpi_cells + 1);
  }
  else
  {
    *my_nmin = num_mpi_extra * (num_mpi_cells + 1) + (rank - num_mpi_extra) * (num_mpi_cells);
    *my_nmax = num_mpi_extra * (num_mpi_cells + 1) + (rank - num_mpi_extra + 1) * (num_mpi_cells);
  }
  ndo = *my_nmax - *my_nmin;

  return ndo;
}


/**********************************************************/
/**
 * @brief communicates the macro-atom B matrices between threads
 *
 *
 * @details communicates the macro-atom B matrices between threads 
 * using MPI_Reduce. Only does anything if MPI_ON flag is on, 
 * and should only be called if geo.rt_mode == RT_MODE_MACRO 
 * and nlevels_macro > 0
 *
 **********************************************************/

int
communicate_matom_matrices ()
{
#ifdef MPI_ON
  int size_of_commbuffer, nrows, n_mpi, n_mpi2, num_comm;
  int my_nmax, my_nmin, ndo, n, position, i;
  char *commbuffer;
  ndo = get_parallel_nrange (rank_global, NPLASMA, np_mpi_global, &my_nmin, &my_nmax);
  Log_parallel ("communicate_matom_matrices: Thread %d is communicating matom matrices for cells %d to %d.\n", rank_global, my_nmin,
                my_nmax);

  nrows = nlevels_macro + 1;
  size_of_commbuffer = 8 * ((nrows * nrows) + 2) * (floor (NPLASMA / np_mpi_global) + 1);
  commbuffer = (char *) malloc (size_of_commbuffer * sizeof (char));

  for (n_mpi = 0; n_mpi < np_mpi_global; n_mpi++)
  {
    position = 0;

    if (rank_global == n_mpi)
    {
      MPI_Pack (&ndo, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
      for (n = my_nmin; n < my_nmax; n++)
      {
        MPI_Pack (&n, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);

        /* we only communicate the matrix if it is being stored in this cell */
        if (macromain[n].store_matom_matrix == TRUE)
        {
          for (i = 0; i < nrows; i++)
          {
            MPI_Pack (macromain[n].matom_matrix[i], nrows, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
          }
        }
      }
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Bcast (commbuffer, size_of_commbuffer, MPI_PACKED, n_mpi, MPI_COMM_WORLD);
    MPI_Barrier (MPI_COMM_WORLD);
    Log_silent ("communicate_matom_matrices: MPI task %d survived broadcasting matom_matrix update information.\n", rank_global);

    position = 0;

    if (rank_global != n_mpi)
    {
      MPI_Unpack (commbuffer, size_of_commbuffer, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
      for (n_mpi2 = 0; n_mpi2 < num_comm; n_mpi2++)
      {
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &n, 1, MPI_INT, MPI_COMM_WORLD);

        /* we only communicate the matrix if it is being stored in this cell */
        if (macromain[n].store_matom_matrix == TRUE)
        {
          for (i = 0; i < nrows; i++)
          {
            MPI_Unpack (commbuffer, size_of_commbuffer, &position, macromain[n].matom_matrix[i], nrows, MPI_DOUBLE, MPI_COMM_WORLD);
          }
        }
      }
    }
  }

  free (commbuffer);
#endif
  return (0);
}
