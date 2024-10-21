
/***********************************************************/
/** @file  sirocco.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  This file contains main and various related routines
 * that are central to the operation of Python
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_errno.h>

#include "atomic.h"
#include "sirocco.h"
#include "models.h"

#define NSPEC	20


/**********************************************************/
/**
 * @brief     The main routine for Python, which supervises the ingest of data defining a model, the initializtion
 * of the model, and actual run of the model, and writing the data to the disk
 *
 *
 * @param [in] int  argc   The number of command line arguments
 * @param [in] char *  argv[]   The command line arguments
 * @return     Always returns 0
 *
 * @details
 *
 * ### Notes ###
 *
 * The main routine of Python is fairly simple conceptually.  One gather the data, one allocates and intializes
 * all of the data structures, one runs the ionization cycles and one run the detalied spectral cycles.
 * The main ruoutin of Python supervises all of these things
 *
 * The main routine is complicated (mainly) due to the multiple ways a model can be run.
 *
 *  * A model can be run from scratch
 *  * A model can be contined from an earlier run (increasing only the number of ionization or spectral cycles)
 *  * A new model can be run beginning with an old windsave file (allowing various parameters of the radiation souces to be set differently
 *  from the old model)
 *
 *  Additionally, there is logic associated with the fact that different types of models require different data.  An AGN for example does
 *  not have (the possibility of) a secondary star.  Most of the inputs from Python come from a parameter file, but Python also has a robust
 *  set of command line switches (See parse_command_line).
 *
 *  Once all the inputs are obtained, the main routine calles routines to allocate data structures needed  to hold the model and to complete
 *  the intialization of both the data structues and various other variables.
 *
 *  The it runs the ionization cycles (See calculate_ionization );
 *
 *  Finally it runs the cycles that calulate the detailed spectra (See make_spectra);
 *
 *
 **********************************************************/

int
main (argc, argv)
     int argc;
     char *argv[];
{
  WindPtr w;

  double freqmin, freqmax;
  int n;
  char values[LINELENGTH], answer[LINELENGTH];
  int get_models ();            // Note: Needed because get_models cannot be included in templates.h
  int dummy_spectype;
  int opar_stat, restart_stat;
  double time_max;
  double lstar;

/* Initialize MPI */

  int my_rank;                  // these two variables are used regardless of parallel mode
  int np_mpi;                   // rank and number of processes, 0 and 1 in non-parallel

#ifdef MPI_ON
  int mpi_err = MPI_Init (&argc, &argv);
  if (mpi_err != EXIT_SUCCESS)
  {
    Error ("Failed to initialise MPI\n");
    Exit (mpi_err);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &np_mpi);
#else
  my_rank = 0;
  np_mpi = 1;
#endif

#ifdef CUDA_ON
  int cuda_err = cusolver_create ();
  if (cuda_err != EXIT_SUCCESS)
  {
    Error ("Failed to initialise CUDA/cuSolver\n");
    Exit (cuda_err);
  }
#endif

  /* Turn off gsl error handling so that the code does not abort on error */
  gsl_set_error_handler_off ();

  np_mpi_global = np_mpi;       // Global variable which holds the number of MPI processes
  rank_global = my_rank;        // Global variable which holds the rank of the active MPI process
  Log_set_mpi_rank (my_rank, np_mpi);   // communicates my_rank to kpar

/* This completes the initialisation of mpi */

  opar_stat = 0;                /* Initialize opar_stat to indicate that if we do not open a rdpar file,
                                   the assumption is that we are reading from the command line */
  restart_stat = FALSE;         /* Assume initially that these is a new run from scratch, and not
                                   a restart */
  time_max = 13.8e9 * 3.2e7;    /* The maximum time the program will run without stopping.  This
                                   is initially set to the lifetime of the universe */

  set_max_time (files.root, time_max);

  rel_mode = REL_MODE_FULL;
  run_xtest = FALSE;
  NWAVE_MAX = (int) NWAVE_IONIZ;

  /* Set the verbosity level for logging.  To get more info raise the verbosity level to a higher number. To
     get less set the verbosity to a lower level. The verbosity can be reset from the comamnd line */

  verbosity = 3;
  Log_set_verbosity (verbosity);

  /* initialise the advanced mode flags (all set to 0) which is a structure in sirocco.h */

  init_advanced_modes ();

  /* Parse the command line. Get the root. create files.diagfolder + diagfiles */


  restart_stat = parse_command_line (argc, argv);

  /* If the restart flag has been set, we check to see if a windsave file exists.  If it doues we will
     we will restart from that point.  If the windsave file does not exist we will start from scratch */

  init_log_and_windsave (restart_stat);

  Log ("Thread %d starting.\n", my_rank);

  /* Start logging of errors and comments */

  Log ("!!Python Version %s \n", VERSION);      //54f -- ksl -- Now read from version.h
  Log ("!!Git commit hash %s\n", GIT_COMMIT_HASH);

  /* warn the user if there are uncommited changes */

  int git_diff_status = GIT_DIFF_STATUS;
  if (git_diff_status > 0)
    Log ("!!Git: This version was compiled with %i files with uncommitted changes.\n", git_diff_status);

  Log ("!!Python is running with %d processors\n", np_mpi_global);
  Log ("This is MPI task number %d (a total of %d tasks are running).\n", rank_global, np_mpi_global);

  Debug ("Debug statements are on. To turn off use lower verbosity (< 5).\n");

  xsignal (files.root, "%-20s Initializing variables for %s\n", "NOK", files.root);

  opar_stat = setup_created_files ();

/* Allocate the domain structure */

  zdom = (DomainPtr) calloc (sizeof (domain_dummy), MAX_DOM);

  /* BEGIN GATHERING INPUT DATA */

  /* Describe the basic calculation in terms of the number of iterations which will
     be used to calculate the wind parameters and the number of iterations and wavelength
     range which will be used for the final spectrom.  Also describe the observer's views
     of the system */


  if (restart_stat == TRUE)     /* We want to continue a run. This is generally used
                                   because we had to limit to runtime of sirocco or we decided
                                   we needed more ionization or spectral cycles */
  {
    Log ("Continuing a previous run of %s \n", files.root);
    strcpy (files.old_windsave, files.root);
    strcat (files.old_windsave, ".wind_save");

    if (wind_read (files.old_windsave) < 0)
    {
      Error ("sirocco: Unable to open %s\n", files.old_windsave);
      Exit (0);
    }
    w = wmain;

    geo.run_type = RUN_TYPE_RESTART;

    xsignal (files.root, "%-20s Read %s\n", "COMMENT", files.old_windsave);

    if (geo.model_count > 0)    //We have previously used models - we need to read them in again
    {
      for (n = 0; n < geo.model_count; n++)
      {
        get_models (geo.model_list[n], 2, &dummy_spectype);
      }
    }
    if (geo.disk_tprofile == DISK_TPROFILE_READIN)      //We also need to re-read in any previously used disk temperature profile
    {
      rdstr ("Disk.T_profile_file", files.tprofile);
      read_non_standard_disk_profile (files.tprofile);
    }
    if (geo.pcycle > 0)
    {
      spec_read (files.specsave);
      xsignal (files.root, "%-20s Read %s\n", "COMMENT", files.specsave);
    }
  }

  else if (restart_stat == FALSE)       /* We are starting a new run, which is the normal mode of operation */
  {

    /* First,  establish the overall system type.  System type should be a physical system,
     * to make things easier for the user.  So really want system types to be something like
     * CV, YSO, AGN so that defaults can be set.  This is now issue #420
     */

    geo.system_type = SYSTEM_TYPE_STAR;
    geo.run_type = RUN_TYPE_NEW;
    geo.mono_freq = 0;


    strcpy (answer, "star");
    sprintf (values, "%d,%d,%d,%d,%d", SYSTEM_TYPE_STAR, SYSTEM_TYPE_CV, SYSTEM_TYPE_BH, SYSTEM_TYPE_AGN, SYSTEM_TYPE_PREVIOUS);
    geo.system_type = rdchoice ("System_type(star,cv,bh,agn,previous)", values, answer);


    if (geo.system_type == SYSTEM_TYPE_PREVIOUS)
    {

      /* This option is for the confusing case where we want to start with a previous wind
         model,(presumably because that run produced a wind close to the one we are looking for,
         but we are going to change some parameters that do not affect the wind geometry,
         We will write use new filenames for the results, so all of the previous work is still saved,
         Note that wind_read also reads the atomic data file that was used to create the previous run of the data.
       */

      strcpy (files.old_windsave, "earlier.run");
      rdstr ("Wind.old_windfile(root_only)", files.old_windsave);
      strcat (files.old_windsave, ".wind_save");


      Log ("Starting a new run from scratch starting using a previous windfile\n");


      if (wind_read (files.old_windsave) < 0)
      {
        Error ("sirocco: Unable to open %s\n", files.old_windsave);
        Exit (0);
      }

      geo.run_type = RUN_TYPE_PREVIOUS;

      w = wmain;
      geo.wcycle = 0;
      geo.pcycle = 0;
      geo.model_count = 0;
    }


    if (geo.run_type == RUN_TYPE_NEW || geo.run_type == RUN_TYPE_PREVIOUS)
    {
      /* This option is the most common one, where we are starting to define a completely new system.
       */

      if (geo.run_type == RUN_TYPE_NEW)
      {
        init_geo ();
      }

      /* get_stellar_params gets information like mstar, rstar, tstar etc.
         it returns the luminosity of the star */

      lstar = get_stellar_params ();

      /* Describe the disk */

      get_disk_params ();


      /* describe the boundary layer / agn components to the spectrum if they exist.
         So that initial condiditions for the bl and agn are initialized sensibly this has
         to come after the disk is defined.
       */

      get_bl_and_agn_params (lstar);

      /* At this point we check whether we have any sources of radiation and exit if we do not */

      if (!geo.star_radiation && !geo.disk_radiation && !geo.bl_radiation && !geo.wind_radiation && !geo.agn_radiation)
      {
        Error ("sirocco: No radiation sources so nothing to do but quit!\n");
        Exit (0);
      }

      /* Describe the wind (or more correctly the various domains).
       */

      rdpar_comment ("Parameters describing the various winds or coronae in the system");

      if (geo.run_type == RUN_TYPE_NEW)
      {
        geo.ndomain = 1;
        rdint ("Wind.number_of_components", &geo.ndomain);

        if (geo.ndomain > MAX_DOM)
        {
          Error ("Maximum number of wind components allowed is %d\n", MAX_DOM);
          Exit (EXIT_FAILURE);
        }

        for (n = 0; n < geo.ndomain; n++)
        {
          get_domain_params (n);
        }
      }

    }
  }

  /* zdom only temporarily needs to be MAX_DOM. Now that we know the number of domains, we'll reallocate
   * the domain to make it smaller */
  zdom = realloc (zdom, sizeof (domain_dummy) * geo.ndomain);
  if (zdom == NULL)
  {
    Error ("sirocco: unable to re-allocate space for domain structure from %d domains to %d domains\n", MAX_DOM, geo.ndomain);
    Exit (EXIT_FAILURE);
  }

/* Get the remainder of the input data.  Note that the next few lines are read from the input file whether or not the windsave file was read in,
   because these are things one would like to be able to change even if we have read in an old windsave file.  init_photons reads in
   the numbers of ionization and spectral cycles and then instatiates PhotPtr.  It is possilbe that this should be moved to else where in
   the flow of reading in data.
 */

  rdpar_comment ("Parameters associated with photon number, cycles,ionization and radiative transfer options");
  init_photons ();
  init_ionization ();

  /* Note: ksl - At this point, SYSTEM_TYPE_PREVIOUS refers both to a restart and to a situation where
   * one is starting from an early wind file as implemented this is quite restrictive about what one
   * can change in the previous case.   */

  if (geo.run_type == RUN_TYPE_NEW)
  {

    /* Getting the atomic data has been moved here as when it was deeply nested in init_ionization, it made it very
     * difficult to write unit tests or other tests. The key issue is that init_ionization is initialising too
     * much all at once, which makes it very difficult to initialise or modify specific things for testing or debug
     * purposes */
    rdstr ("Atomic_data", geo.atomic_filename);
    setup_atomic_data (geo.atomic_filename);

    /* Describe the wind, by calling get_wind_params one or more times
       and then gets params by calling e.g. get_sv_wind_params() */

    for (n = 0; n < geo.ndomain; n++)
    {
      rdpar_comment ("Parameters for Domain %d", n);
      get_wind_params (n);
    }

  }                             // End of block to define a model for the first time
  else if (modes.zeus_connect == 1)     /* We are in rad-hydro mode, we want the new density and temperature */
  {
    /* Hydro takes the wind domain number as an argument in the current domains setup */
    Log ("We are going to read in the density and temperature from a zeus file\n");
    get_hydro (geo.hydro_domain_number);        //This line just populates the hydro structures
  }


  /* Calculate parameters associated with the binary star system */

  if (geo.binary)
    binary_basics ();

  /* Check that the parameters which have been supplied for the star, disk and boundary layer will
     allow generation of photons where that is appropriate */

  if (geo.tstar <= 0.0)
    geo.star_radiation = 0;
  if (geo.disk_mdot <= 0.0 && geo.disk_tprofile == DISK_TPROFILE_STANDARD)
    geo.disk_radiation = 0;
  if (geo.t_bl <= 0.0 || geo.lum_bl <= 0.0)
    geo.bl_radiation = 0;

  /* If the disk radius is <0, assume no disk was intended. */

  if (geo.disk_rad_max <= 0.0)
  {
    geo.disk_type = DISK_NONE;
    geo.disk_radiation = 0;
  }

  if (geo.star_radiation)
    Log ("There is a star which radiates\n");
  else
    Log ("The star in the system does not radiate\n");

  if (!geo.disk_type)
    Log ("There is no disk in the system \n");
  else if (!geo.disk_radiation)
    Log ("The disk exists, but only absorbs photons\n");
  else
    Log ("There is a disk which radiates and absorbs\n");

  if (geo.bl_radiation || (geo.agn_radiation && geo.system_type != SYSTEM_TYPE_AGN))
    Log ("There is a boundary layer which radiates\n");
  else
    Log ("There is no boundary layer\n");

  if (geo.agn_radiation)
    Log ("There is a BH  which radiates\n");
  else
    Log ("There is no BH \n");

  if (geo.rt_mode == RT_MODE_MACRO)
  {
    if (nlevels_macro == 0)
    {
      Error ("THIS IS A MACROATOM CALCULATION WITH NO MACROLEVELS\n");
    }
    else
    {
      Log ("This is a macro-atom calculation\n");
    }
  }
  else
  {
    Log ("This is a simple atom calculation\n");
  }

  /* Describe the spectra which will be extracted and the way it will be extracted */


  if (geo.pcycles > 0 && geo.pcycle == 0)
  {

    rdpar_comment ("Parameters defining the spectra seen by observers\n");

    if (geo.star_radiation)
    {
      geo.star_spectype = geo.star_ion_spectype;
      get_spectype (geo.star_radiation, "Central_object.rad_type_in_final_spectrum(bb,models,uniform,mono)", &geo.star_spectype);
    }

    if (geo.disk_radiation)
    {
      geo.disk_spectype = geo.disk_ion_spectype;
      get_spectype (geo.disk_radiation, "Disk.rad_type_in_final_spectrum(bb,models,uniform,mono,mod_bb)", &geo.disk_spectype);
    }

    if (geo.bl_radiation)
    {
      geo.bl_spectype = geo.bl_ion_spectype;
      get_spectype (geo.bl_radiation, "Boundary_layer.rad_type_in_final_spectrum(bb,models,uniform)", &geo.bl_spectype);
    }

    if (geo.agn_radiation)
    {
      // This block will run for both AGN, *and* some versions of a boundary layer.
      // Even though we're setting the same params, we need to change the wording based on the system, unfortunately.
      geo.agn_spectype = geo.agn_ion_spectype;

      // If there is 'AGN radiation' that genuinely *is* AGN radiation (and not a star boundary layer
      if (geo.system_type == SYSTEM_TYPE_AGN || geo.system_type == SYSTEM_TYPE_BH)
      {
        get_spectype (geo.agn_radiation, "Central_object.rad_type_in_final_spectrum(bb,models,power,cloudy,brems,mono)", &geo.agn_spectype);
      }
      else
      {
        get_spectype (geo.agn_radiation, "Boundary_layer.rad_type_in_final_spectrum(power)", &geo.agn_spectype);
      }


      if (geo.agn_spectype >= 0 && comp[geo.agn_spectype].nmods != 1)
      {
        Error ("sirocco: When using models with an AGN, there should be exactly 1 model, we have %i for spectrum cycles\n",
               comp[geo.agn_ion_spectype].nmods);
        exit (0);
      }
    }
    init_observers ();
  }

  geo.matom_radiation = 0;      //initialise for ionization cycles - don't use pre-computed emissivities for macro-atom levels/ k-packets.
  get_standard_care_factors ();
  get_meta_params ();

/* Establish limits on the frequency intervals to be used by the ionization cycles and
 * the fraquency bands for stratified sampling. These bands are alos used as the spectral
 * intervals for creating crude spectra in each of the cells
 *
 * This section of inputs might logically go earlier in the code, but
 * was put here so it would add on to existing .pf files.  It would be reasonble to consider moving
 * it to a more logical location
 */

  rdpar_comment ("Other parameters");

  bands_init (-1, &xband);
  freqmin = xband.f1[0];
  freqmax = xband.f2[xband.nbands - 1];

  if (modes.iadvanced)
  {
    /* Do we require extra diagnostics or not */
    strcpy (answer, "no");
    modes.diag_on_off = rdchoice ("@Diag.extra(yes,no)", "1,0", answer);
    if (modes.diag_on_off)
    {
      get_extra_diagnostics ();
    }
  }


  /* Wrap up and save all the inputs */


  if (strncmp (files.root, "mod", 3) == 0)
  {
    cpar ("mod.pf");
  }

  else if (opar_stat == 1)
  {
    cpar (files.input);
  }
  else
  {
    cpar (files.new_pf);
  }

  /* At this point, all inputs have been obtained at this point and the inputs have been copied to "mod.pf" or "sirocco.pf"
   * If we have used, the -i flag, we quit; otherwise we continue on to run the model */
  if (modes.quit_after_inputs)
  {
    Log ("This was was run with the -i or --dry-run flag set, so quitting now inputs have been gathered.\n");
    error_summary ("dry run.");
    exit (0);
  }


  if (rdpar_check ())
  {
    Log ("Some of the input have not been updated for the current version of Python.  Please correct and rerun\n");
    exit (0);
  }


  /* INPUTS ARE FINALLY COMPLETE */


  Log ("There are %d domains\n", geo.ndomain);
  for (n = 0; n < geo.ndomain; n++)
  {
    Log ("%20s type: %3d  ndim: %3d mdim: %3d ndim2: %4d\n", zdom[n].name, zdom[n].wind_type, zdom[n].ndim, zdom[n].mdim, zdom[n].ndim2);
  }


  /* Now define the wind cones generically. modifies the global windcone structure */
  setup_windcone ();

  /* initialize the random number generator */
  /* By default, the random number generator start with fixed seeds (differnt
   * for each processor, but this can be changed using a command line
   * switch.
   *
   * An exception is when we are in zeus mode, where it would be inappropriate
   * to use the same phtons in ezch cycle.  There we initiate the seeds unsing
   * the clock
   */
  if (modes.rand_seed_usetime)
  {
    n = (unsigned int) clock () * (rank_global + 1);
    init_rand (n);

  }
  else
  {
    init_rand (1084515760 + (13 * rank_global));
  }

  /* DFUDGE is a distance that assures we can "push through" boundaries.  setup_dfudge
     sets the push through distance depending on the size of the system.
   */

  DFUDGE = setup_dfudge ();

  /* Next line finally defines the wind if this is the initial time this model is being run */

  if (geo.run_type == RUN_TYPE_NEW)
  {
    define_wind ();
  }

  Log ("DFUDGE set to %e based on geo.rmax\n", DFUDGE);

  if (modes.zeus_connect == 1)  //We have restarted, but are in zeus connect mode, so we want to update density, temp and velocities
  {

    /* Hydro takes the wind domain number as an argument in the current domains setup */
    hydro_restart (geo.hydro_domain_number);
  }

  /* this routine checks, somewhat crudely, if the grid is well enough resolved */
  check_grid ();

  w = wmain;
  if (modes.extra_diagnostics)
  {
    init_extra_diagnostics ();
  }


/* The next section sets up two represenations of the disk structure
 *
 * disk_init is called primarily to get
 * a defined set of annular rings which are kept throughout the
 * ionization calculation.
 * disk_init calculates the flux from the disk in the energy range set by
 * freqmin and freqmax, and uses is this to identify the position of the
 * rings in the disk, so that each ring contributes the same amount to
 * the flux
 *
 * A second structure qdisk is needed
 * because in the process of generating photons in various bands
 * the annular rings are changed
 *
 */


  disk_init (geo.disk_rad_min, geo.disk_rad_max, geo.mstar, geo.disk_mdot, freqmin, freqmax, 0, &geo.f_disk);
  qdisk_init (geo.disk_rad_min, geo.disk_rad_max, geo.mstar, geo.disk_mdot);
  xsignal (files.root, "%-20s Finished initialization for %s\n", "NOK", files.root);
  check_time (files.root);

  /* allow the user to quit after the wind has been defined */
  if (modes.quit_after_wind_defined)
  {
    wind_save (files.windsave);
    Log ("This was was run with the ---grid-only flag set, so quitting now wind has been defined.\n");
    error_summary ("wind definition only (--grid-only).");
#ifdef MPI_ON
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize ();
#endif
    return EXIT_SUCCESS;
  }


  /* Allow for the possibility of running a special diagnostic mode in
     a stand alone routine xtest. This will happen with the command line
     option -xtest.  */
  if (run_xtest)
  {
    xtest ();
  }


  if (geo.wcycles == 0 && geo.pcycles == 0)
  {
    wind_save (files.windsave);
    Log ("Both ionization and spectral cycles are set to 0; Saving windfile but then exiting\n");
    exit (1);                   //There is really nothing to do!
  }

/* XXXX - THE CALCULATION OF THE IONIZATION OF THE WIND */

  geo.ioniz_or_extract = CYCLE_IONIZ;

  calculate_ionization (restart_stat);

/* XXXX - END OF CYCLE TO CALCULATE THE IONIZATION OF THE WIND */
  Log (" Completed wind creation.  The elapsed TIME was %f\n", timer ());
  /* Evaluate wind paths for last iteration */
  if (geo.reverb == REV_WIND || geo.reverb == REV_MATOM)
  {                             //If this is a mode in which we keep wind arrays, update them
    wind_paths_evaluate (w, my_rank);
  }

/* XXXX - THE CALCULATION OF A DETAILED SPECTRUM IN A SPECIFIC REGION OF WAVELENGTH SPACE */

  freqmax = VLIGHT / (geo.swavemin * 1.e-8);
  freqmin = VLIGHT / (geo.swavemax * 1.e-8);

  /* Perform the initilizations required to handle macro-atoms during the detailed
     calculation of the spectrum.

     Next lines turns off macro atom estimators and other portions of the code that are
     unnecessary during spectrum cycles.  */

  geo.ioniz_or_extract = CYCLE_EXTRACT;


/* Next step speeds up extraction stage */

  if (!modes.keep_photoabs)
  {
    DENSITY_PHOT_MIN = -1.0;    // Do not calculated photoabsorption in detailed spectrum
  }

  /*Switch on k-packet/macro atom emissivities  SS June 04 */

  if (geo.rt_mode == RT_MODE_MACRO)
  {
    geo.matom_radiation = 1;
  }

  /* Finished initializations required for macro-atom approach */

  /* Calculate and store which bf processess need to be considered in each cell
   * Note that this is not macro-specific but is just to speed the program up.
   */

  kbf_need (freqmin, freqmax);

  /* XXXX - Execute  CYCLES TO CREATE THE DETAILED SPECTRUM */

  make_spectra (restart_stat);


/* Finally done */

#ifdef MPI_ON
  char dummy[LINELENGTH];
  sprintf (dummy, "End of program, Thread %d only", rank_global);       // added so we make clear these are just errors for thread ngit status
  error_summary (dummy);        // Summarize the errors that were recorded by the program
  Log ("Run py_error.py for full error report.\n");
  MPI_Finalize ();
#else
  error_summary ("End of program");     // Summarize the errors that were recorded by the program
#endif

#ifdef CUDA_ON
  cusolver_destroy ();
#endif

  xsignal (files.root, "%-20s %s\n", "COMPLETE", files.root);
  Log ("\nBrief Run Summary\nAt program completion, the elapsed TIME was %f\n", timer ());
  Log ("There were %d of %d ionization cycles and %d of %d spectral cycles run\n", geo.wcycle, geo.wcycles, geo.pcycle, geo.pcycles);
  if (geo.rt_mode == RT_MODE_MACRO)
  {
    if (nlevels_macro == 0)
    {
      Log ("THIS WAS A MACROATOM CALCULATION WITH NO MACROLEVELS. (Use for diagnostics only)\n");
    }
    else
    {
      Log ("This was a macro-atom calculation\n");
    }
  }
  else
  {
    Log ("This was a simple atom calculation\n");
  }

  Log ("Convergence statistics for the wind after the ionization calculation:\n");
  check_convergence ();
  Log ("Information about luminosities and apparent fluxes due to various portions of the system:\n");
  phot_status ();

  clean_on_exit ();

  return (0);
}
