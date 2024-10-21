
/***********************************************************/
/** @file  unit_test.c
 * @author ksl
 * @date   January, 2021
 *
 * @brief  This file contains main and various related routines
 * to test one or more low level routines in Python
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/stat.h>
#include <string.h>

#include <string.h>
#include <math.h>
#include "atomic.h"
#include "sirocco.h"

char inroot[LINELENGTH];

#define COMPTON_TEST 0
#define NUMBER_TEST  1
#define LUM_TEST 0

int
main (argc, argv)
     int argc;
     char *argv[];
{

//  FILE *fptr, *fopen ();

  int my_rank;                  // these two variables are used regardless of parallel mode
  int np_mpi;                   // rank and number of processes, 0 and 1 in non-parallel



#ifdef MPI_ON
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &np_mpi);
#else
  my_rank = 0;
  np_mpi = 1;
#endif

  np_mpi_global = np_mpi;       // Global variable which holds the number of MPI processes
  rank_global = my_rank;        // Global variable which holds the rank of the active MPI process
  Log_set_mpi_rank (my_rank, np_mpi);   // communicates my_rank to kpar

  init_rand (1084515760 + (13 * rank_global));

  printf ("Beginning unit test\n");


  xsignal ("unit_test", "%-20s Initializing variables for %s\n", "NOK", "unit_test");

#if NUMBER_TEST == 1
  int n;
  int nmax = 1000000;
  double x[nmax + 1];
  double q;

  double m, mrepeat;
  double j, jrepeat;

  mrepeat = 1e2;
  jrepeat = 1e2;;

  for (n = 0; n <= nmax; n++)
  {
    x[n] = .1 + random_number (0, 100);
  }

  xsignal ("unit_test", "%-20s at start %s\n", "NOK", "unit_test");


  for (j = 0; j < jrepeat; j++)
  {
    for (m = 0; m < mrepeat; m++)
    {
      for (n = 0; n < nmax; n++)
      {
        q = x[n] + x[n + 1];
      }
    }
  }

  xsignal ("unit_test", "%-20s after addition %s\n", "NOK", "unit_test");


  for (j = 0; j < jrepeat; j++)
  {
    for (m = 0; m < mrepeat; m++)
    {
      for (n = 0; n < nmax; n++)
      {
        q = x[n] * x[n + 1];
      }
    }
  }

  xsignal ("unit_test", "%-20s after multiplication %s\n", "NOK", "unit_test");


  for (j = 0; j < jrepeat; j++)
  {
    for (m = 0; m < mrepeat; m++)
    {
      for (n = 0; n < nmax; n++)
      {
        q = x[n] / x[n + 1];
      }
    }
  }
  xsignal ("unit_test", "%-20s after division %s\n", "NOK", "unit_test");


  for (j = 0; j < jrepeat; j++)
  {
    for (m = 0; m < mrepeat; m++)
    {
      for (n = 0; n < nmax; n++)
      {
        q = exp (-x[n]);
      }
    }
  }
  xsignal ("unit_test", "%-20s after expon %s\n", "NOK", "unit_test");



  for (j = 0; j < jrepeat; j++)
  {
    for (m = 0; m < mrepeat; m++)
    {
      for (n = 0; n < nmax; n++)
      {
        q = log10 (x[n]);
      }
    }
  }
  xsignal ("unit_test", "%-20s after log10 %s\n", "NOK", "unit_test");


  double fibo = 1;
  for (j = 0; j < jrepeat; j++)
  {
    for (m = 0; m < mrepeat; m++)
    {
      for (n = 0; n < nmax; n++)
      {
        fibo += fibo;
      }
    }
  }
  xsignal ("unit_test", "%-20s after fibo %s\n", "NOK", "unit_test");




  printf ("q %e fibo %e\n", q, fibo);



#endif



#if COMPTON_TEST == 1
  int n;
  double velocity[3];
  double t;
  double v;
  double speed = 0;
  double speed2 = 0;
  t = 1.0e4;

  for (n = 0; n < 1000000; n++)
  {
    compton_get_thermal_velocity (t, velocity);
    printf ("%10.3e %10.3e %10.3e\n", velocity[0], velocity[1], velocity[2]);
    v = sqrt (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);
    speed += v;
    speed2 += v * v;

  }

  printf (" Average speed %e  average v**2  %e\n", speed / n, speed2 / n);

  printf (" The expected average is %e\n", 2. / sqrt (3.141592) * sqrt (2. * BOLTZMANN * t / MELEC));
#endif


#if LUM_TEST == 1
/* Create a way to do luminosities in parallel */

  char infile[LINELENGTH];
  double lum_one;
  int zparse ();

  double f1 = 1e12;
  double f2 = 1e18;


  zparse (argc, argv);
  sprintf (infile, "%.150s.wind_save", inroot);
  printf ("Reading %s \n", infile);
  wind_read (infile);


  int number;





  xsignal ("unit_test", "%-20s before old wind_lum for %s\n", "NOK", "unit_test");

  for (number = 0; number < 10; number++)
  {
    lum_one = wind_luminosity (f1, f2, MODE_CMF_TIME);
  }

  if (my_rank == 0)
  {
    printf ("Old version: %e %e %e %e\n", lum_one, geo.lum_lines, geo.lum_rr, geo.lum_ff);
  }



  double par_wind_luminosity ();

  xsignal ("unit_test", "%-20s before new wind luminosity %s\n", "NOK", "unit_test");

  for (number = 0; number < 10; number++)
  {
    lum_one = par_wind_luminosity (f1, f2, MODE_CMF_TIME);
  }

  xsignal ("unit_test", "%-20s after new wind luminosity %s\n", "NOK", "unit_test");

  Log_parallel ("new version: %d %e %e %e %e\n", my_rank, lum_one, geo.lum_lines, geo.lum_rr, geo.lum_ff);
#endif

  printf ("Finished unit test\n");


  return (0);

}


int
zparse (argc, argv)
     int argc;
     char *argv[];
{
  char dummy[LINELENGTH];

  if (argc != 2)
  {
    printf ("usage: unit_test root\n");
    exit (1);
  }


  strcpy (dummy, argv[1]);
  get_root (inroot, dummy);



  return (0);

}



double
par_wind_luminosity (f1, f2, mode)
     double f1, f2;
     int mode;
{
  double lum, lum_lines, lum_rr, lum_ff, factor;
  int nplasma;
  int ndo, my_nmin, my_nmax, n;


  lum = lum_lines = lum_rr = lum_ff = factor = 0.0;



#ifdef MPI_ON
  ndo = get_parallel_nrange (rank_global, NPLASMA, np_mpi_global, &my_nmin, &my_nmax);
  // Log parallel gets each thread as does printf
  Log_parallel ("xxxwind_luminosity: thread %d lum %d to %d %d \n", rank_global, my_nmin, my_nmax, ndo);
  Log ("wind_luminosity: thread %d lum %d to %d %d \n", rank_global, my_nmin, my_nmax, ndo);
#else
  my_nmin = 0;
  my_nmax = NPLASMA;
  ndo = NPLASMA;
#endif


  Log_parallel ("Hello world %d\n", rank_global);


  for (nplasma = my_nmin; nplasma < my_nmax; nplasma++)
  {

    if (wmain[plasmamain[nplasma].nwind].inwind < 0)
    {
      Error ("wind_luminosty: Trying to calculate luminosity for a wind cell %d that has plasma cell %d but is not in the wind\n",
             plasmamain[nplasma].nwind, nplasma);
    }


    total_emission (&plasmamain[nplasma], f1, f2);
  }

  /* So at this point I have calculated the lumiosities in one of the threads for all of the plasma cells,
     but i now need to get this informationm to all of the threads.   Presumably the values are undetermined
     in the cells I have not calculated
   */

#ifdef MPI_ON
  int size_of_commbuffer, n_mpi, n_mpi2, num_comm;
  int position;
  char *commbuffer;

  // We are currently transmtting 1 integer (4) and 1 double, but we need to add 4 bits for the process number

  /* We need to transmit 
     process #  only 4 bits  once

     The remainder need to be transimitted for each cell

     cell number     4
     lum_lines       8
     lum_rr          8
     lum_ff          8

     Total          28
   */



  size_of_commbuffer = 28 * (floor (NPLASMA / np_mpi_global) + 1) + 4;

  commbuffer = (char *) malloc (size_of_commbuffer * sizeof (char));

  Log ("commbuffer size %d  %d\n", size_of_commbuffer, (floor (NPLASMA / np_mpi_global) + 1));


  MPI_Barrier (MPI_COMM_WORLD);

  for (n_mpi = 0; n_mpi < np_mpi_global; n_mpi++)
  {
    Log ("MPI task %d is working on cells %d to max %d (total size %d).\n", rank_global, my_nmin, my_nmax, NPLASMA);

    position = 0;

    if (rank_global == n_mpi)
    {
      // First tansmit ndo, whikch is the number of tasks (elements) this thread is working on (4)

      MPI_Pack (&ndo, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
      // Log ("Position1 %d %d\n", n, position);
      for (n = my_nmin; n < my_nmax; n++)
      {
        // Next  transimit number of the the plasma cell (4) 
        MPI_Pack (&n, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        //   Log ("Position2 %d %d\n", n, position);
        // Now transimit the values we want (8)
        MPI_Pack (&plasmamain[n].lum_lines, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_rr, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_ff, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);

        //    Log ("Position3 %d %d\n", n, position);
      }
    }
    Log ("Luminoisity,MPI task %d broadcasting plasma update information.\n", rank_global);

    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Bcast (commbuffer, size_of_commbuffer, MPI_PACKED, n_mpi, MPI_COMM_WORLD);
    MPI_Barrier (MPI_COMM_WORLD);
    Log ("Luminosity: MPI task %d survived plasma update information.\n", rank_global);

    position = 0;

    if (rank_global != n_mpi)
    {
      MPI_Unpack (commbuffer, size_of_commbuffer, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
      for (n_mpi2 = 0; n_mpi2 < num_comm; n_mpi2++)
      {
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &n, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);

      }
    }
  }
  free (commbuffer);

#endif


  lum = lum_lines = lum_rr = lum_ff = factor = 0.0;

  for (nplasma = 0; nplasma < NPLASMA; nplasma++)
  {

    if (mode == MODE_OBSERVER_FRAME_TIME)
      factor = 1.0 / plasmamain[nplasma].xgamma;        /* this is dt_cmf */
    else if (mode == MODE_CMF_TIME)
      factor = 1.0;

    lum_lines += plasmamain[nplasma].lum_lines * factor;
    lum_rr += plasmamain[nplasma].lum_rr * factor;
    lum_ff += plasmamain[nplasma].lum_ff * factor;
  }

  lum = lum_lines + lum_rr + lum_ff;

  if (mode == MODE_CMF_TIME)
  {
    geo.lum_lines = lum_lines;
    geo.lum_rr = lum_rr;
    geo.lum_ff = lum_ff;
  }

  return (lum);
}
