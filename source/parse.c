
/***********************************************************/
/** @file  parse.c
 * @author ksl
 * @date   February, 2018
 *
 * @brief  Routines associated with parsing the command line
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/** 
 * @brief      parses the command line options
 *
 * @param [in]  int  argc   the number of command line arguments
 * @param [in]  char *  argv[]   The command line arguments
 * @return      restart_stat   1 if restarting a previous model,
 * 0 in all other cases.
 *
 * Python has a fairly rich set of command line options, which
 * are parsed by this routine
 *
 * The routine also creates the diag folder, which is where most
 * of the log files are written 
 *
 * ###Notes###
 *
 * The general purpose of each of the command line options
 * should be fairly obvious from reading the code.
 *
 * If changes to the command line interface are made they should
 * be described in the routine help 
 *
 * Although this routine uses the standard Log and Error commands
 * the diag files have not been created yet and so this information
 * is really simply written to the terminal.  
 *
 **********************************************************/

int
parse_command_line (argc, argv)
     int argc;
     char *argv[];
{
  int restart_stat, verbosity, max_errors, i;
  int j = 0;
  char dummy[LINELENGTH];
  int mkdir ();
  double time_max;
  char *fgets_rc;
  double x;

  restart_stat = 0;

  if (argc == 1)
  {
    printf ("Parameter file name (e.g. my_model.pf, or just my_model):");
    fgets_rc = fgets (dummy, LINELENGTH, stdin);
    if (!fgets_rc)
    {
      Error ("Input rootname is NULL or invalid\n");
      exit (1);
    }
    get_root (files.root, dummy);
    strcpy (files.diag, files.root);
    strcat (files.diag, ".diag");
  }
  else
  {

    for (i = 1; i < argc; i++)
    {
      if (strcmp (argv[i], "-h") == 0)
      {
        help ();
      }
      else if (strcmp (argv[i], "-r") == 0)
      {
        restart_stat = 1;
        j = i;
      }
      else if (strcmp (argv[i], "-t") == 0)
      {
        if (sscanf (argv[i + 1], "%lf", &time_max) != 1)
        {
          Error ("sirocco: Expected time after -t switch\n");
          exit (1);
        }
        set_max_time (files.root, time_max);
        i++;
        j = i;
        Log ("Program will stop after time %f\n", time_max);

      }
      else if (strcmp (argv[i], "-v") == 0)
      {
        if (sscanf (argv[i + 1], "%d", &verbosity) != 1)
        {
          Error ("sirocco: Expected verbosity after -v switch\n");
          exit (0);
        }
        Log_set_verbosity (verbosity);
        i++;
        j = i;
        Log ("Verbosity level set to %d\n", verbosity);

      }
      else if (strcmp (argv[i], "-e") == 0)
      {
        if (sscanf (argv[i + 1], "%lf", &x) != 1)
        {
          Error ("sirocco: Expected max errors after -e switch\n");
          exit (1);
        }
        max_errors = x;
        Log_quit_after_n_errors (max_errors);
        i++;
        j = i;
        Log ("Setting the maximum number of errors of a type before quitting to %d\n", max_errors);

      }
      else if (strcmp (argv[i], "-e_write") == 0)
      {
        if (sscanf (argv[i + 1], "%lf", &x) != 1)
        {
          Error ("sirocco: Expected max errors after -e switch\n");
          exit (1);
        }
        max_errors = x;
        Log_print_max (max_errors);
        i++;
        j = i;
        Log ("Setting the maximum number of errors of a type to print out to  %d\n", max_errors);

      }
      else if (strcmp (argv[i], "-d") == 0)
      {
        modes.iadvanced = 1;
        Log ("Enabling advanced/diagnostic inputs (@ commands)\n");
        j = i;
      }
      else if (strcmp (argv[i], "-gamma") == 0)
      {
        rel_mode = REL_MODE_FULL;
        Log ("Using special relativity for Doppler shifts, etc., and including co-moving frame effects.\n");
        j = i;
      }
      else if (strcmp (argv[i], "-sr_doppler_only") == 0)
      {
        rel_mode = REL_MODE_SR_FREQ;
        Log ("Using full special relativity only for Doppler shifts, etc., and not considering co-moving frame effects.");
        j = i;
      }
      else if (strcmp (argv[i], "-nonrel") == 0)
      {
        rel_mode = REL_MODE_LINEAR;
        Log ("Using only old approach with linear Doppler shifts, etc. and not considering co-moving frame effects.\n");
        j = i;
      }
      else if (strcmp (argv[i], "-no-matrix-storage") == 0)
      {
        modes.store_matom_matrix = FALSE;
        Log ("Not storing the macro-atom matrix (on-the-fly method) if Matom.ransition_mode is matrix.\n");
        j = i;
      }

      else if (strcmp (argv[i], "--version") == 0)
      {
        /* give information about the sirocco version, such as commit hash */
        Log ("Python Version %s \n", VERSION);  //54f -- ksl -- Now read from version.h
        Log ("Built from git commit hash %s\n", GIT_COMMIT_HASH);
        /* warn the user if there are uncommited changes */
        int git_diff_status = GIT_DIFF_STATUS;
        if (git_diff_status > 0)
          Log ("This version was compiled with %i files with uncommitted changes.\n", git_diff_status);
        exit (1);
      }

      else if (strcmp (argv[i], "-xtest") == 0)
      {
        run_xtest = TRUE;
        Log ("Run xstest, usually instead of normal Python.\n");
        j = i;
      }
      else if (strcmp (argv[i], "-include_partial_cells") == 0)
      {
        modes.partial_cells = PC_INCLUDE;
        Log ("Cells partially in the wind will be included.\n");
        j = i;
      }
      else if (strcmp (argv[i], "-ignore_partial_cells") == 0)
      {
        modes.partial_cells = PC_ZERO_DEN;
        Log ("Cells partially in the wind will be ingnored.\n");
        j = i;
      }
      else if (strcmp (argv[i], "-f") == 0)
      {
        modes.fixed_temp = 1;
        Log ("Invoking fixed temperature mode\n");
        j = i;
      }
      else if (strcmp (argv[i], "--rseed") == 0)
      {
        modes.rand_seed_usetime = 1;
        j = i;
        Log ("Using a random seed in random number generator\n");
      }
      else if (strcmp (argv[i], "--rng") == 0)
      {
        modes.save_rng = 1;
        modes.load_rng = 1;
        j = i;
        Log ("Using a persistent RNG state\n");
      }
      else if (strcmp (argv[i], "-z") == 0)
      {
        modes.zeus_connect = 1;
        Log ("Setting zeus_connect to %i\n", modes.zeus_connect);
        j = i;
      }
      else if (strcmp (argv[i], "-i") == 0)
      {
        modes.quit_after_inputs = 1;
        j = i;
      }
      else if (strcmp (argv[i], "-p") == 0)
      {
        Log ("Logarithmic photon stepping enabled\n");
        modes.photon_speedup = 1;
        if (sscanf (argv[i + 1], "%lf", &PHOT_RANGE) == 1)
        {
          i++;
        }
        else
        {
          PHOT_RANGE = 1.;
        }
        j = i;
      }
      else if (strcmp (argv[i], "--dry-run") == 0)
      {
        modes.quit_after_inputs = 1;
        j = i;
      }

      else if (strcmp (argv[i], "--grid-only") == 0)
      {
        modes.quit_after_wind_defined = 1;
        j = i;
      }

      else if (strcmp (argv[i], "--version") == 0)
      {
        /* give information about the pyhon version, such as commit hash */
        Log ("Python Version %s \n", VERSION);  //54f -- ksl -- Now read from version.h
        Log ("Built from git commit hash %s\n", GIT_COMMIT_HASH);
        /* warn the user if there are uncommited changes */
        int git_diff_status = GIT_DIFF_STATUS;
        if (git_diff_status > 0)
          Log ("This version was compiled with %i files with uncommitted changes.\n", git_diff_status);
        exit (1);
      }

      else if (strncmp (argv[i], "-", 1) == 0)
      {
        Error ("sirocco: Unknown switch %s\n", argv[i]);
        help ();
      }
    }

    /* The last command line variable is always the .pf file */

    if (j + 1 == argc)
    {
      Error ("All of the command line has been consumed without specifying a parameter file name, so exiting\n");
      exit (1);
    }

    strcpy (dummy, argv[argc - 1]);
    get_root (files.root, dummy);

    /* This completes the parsing of the command line */

    /* Create a subdirectory to store diaganostic files */

    sprintf (files.diagfolder, "diag_%.100s/", files.root);
    mkdir (files.diagfolder, 0777);

    sprintf (dummy, "_%02d.diag", rank_global);

    sprintf (files.diag, "%.50s/%.50s%.50s", files.diagfolder, files.root, dummy);

    /* Set up the directory structure for storing the rng state */

    if (modes.save_rng)
    {
      init_rng_directory (files.root, rank_global);
    }
  }

  if (restart_stat)
  {
    Log ("\n*****Restarting %s *****\n", files.root);
    Log
      ("WARNING: With  the execption of the number of ionization and spectral cycles\n any changes in the parameter files will be ignored\n\n");
  }

  return (restart_stat);
}




/**********************************************************/
/** 
 * @brief      print out some basic help concering command line options for the program
 *
 * @return     0 in all cases
 *
 * This simply prints out a hardwired multi-line string to the command line.
 *
 * ###Notes###
 *
 * An improved version of this routine would simply read and print and external 
 * file.
 *
 * If one chooses to add to what is printed out then one needs to be careful 
 * concerning the ends of lines.
 *
 **********************************************************/

void
help ()
{
#ifdef MPI_ON
  if (rank_global == 0)
  {
#endif
    char *some_help;

    some_help = "\
\n\
This program simulates radiative transfer in a (biconical) CV, YSO, quasar, TDE or (spherical) stellar wind \n\
\n\
Usage:  py [-h] [-r] [-t time_max] [-v n] [--dry-run] [-i] [--version] [--rseed] [-p n_steps] [-nonrel] [-sr_doppler_only] xxx  or simply py \n\
\n\
where xxx is the rootname or full name of a parameter file, e. g. test.pf \n\
\n\
and the switches have the following meanings \n\
\n\
 -h                     Print this help message and stop \n\
 -r                     Restart a run of the progarm reading the file xxx.windsave \n\
 -t time_max            Limit the total time to approximately time_max seconds.  Note that the program checks \n\
                        for this limit somewhat infrequently, usually at the ends of cycles, because it \n\
                        is attempting to save the program outputs so that the program can be restarted with \n\
                        -r if that is desired. \n\
 -v n                   Increase or decrease the amount of print out.  The default is 5.  Larger numbers increase  \n\
                        the amount printed; smaller numbers decrease it.   \n\
 --dry-run              Create a new .pf file and stop \n\
 -i                     Same as --dry-run \n\
 --grid-only            Define the wind grid and save to wind_save file, then stop \n\
 --version              Print out sirocco version, commit hash and if there were files with uncommitted \n\
                        changes and stop \n\
 --rseed                Set the random number seed to be time-based, rather than fixed. \n\
 --rng                  Save or load the RNG state to file, to allow persistent RNG states between restarts\n\
\n\
Other switches exist but these are not intended for the general user.\n\
These are largely diagnostic or for special cases. These include\n\
 -d                     Enable advanced/diagnostic inputs (normally for debugging purposes) \n\
                        Python will then query the user for information about what to do with a series of \n\
                        inputs beginning with @ \n\
 -e                     Change the maximum number of errors of one type (by default 100,000) before the program will quit\n\
 -e_write               Change the maximum number of errors of one type (by default 100) to print out before recording errors silently\n\
 -f                     Invoke a fixed temperature mode, used for runs with Zeus or Plutu \n\
 -z                     Invoke a special mode for that causes Python to start with a run from Zeus or Plutu\n\
 -p [range]             Vary the number of photons in ionization cycles logarthmically building up to the final value\n\
                        Range is in powers of 10, the difference beween the number of photons in the first cycle \n\
                        compared to the last. If range is missing, range is assumed to be 1, in which case the  \n\
                        number of photons will in the first cycle will be one order of magniude less than in the last cycle \n\
 -nonrel                Use Python in its old non-relativistic configuration, with linear Doppler shifts, etc., and where co-moving frame\n\
                        effects are not taken into account.\n\
 -sr_doppler_only       Use Python with full special relativity for Doppler shifts, etc., but do not include any co-moving frame\n\
                        effects.\n\
 -ignore_partial_cells  Ignore wind cells that are only partially filled by the wind (This is now the default)  \n\
 -include_partial_cells Include wind cells that are only partially filled by the wind   \n\
 -no-matrix-storage     Do not store macro-atom transition matrices if using the macro-atom line transfer and the matrix matom_transition_mode.\n\
\n\
 -xtest                 Instead of running sirocco, call the routine xtest so that one can diagnose issues associted with the \n\
                        setup.  This is only useful to devlopers \n\
If one simply types py or pyZZ where ZZ is the version number, one is queried for a name \n\
of the parameter file and inputs will be requested from the command line. \n\
\n\
\n\
";                              // End of string to provide one with help

    printf ("%s\n", some_help);
#ifdef MPI_ON
  }
#endif

  Exit (0);                     // Note that here we simply do want to exit, not use Exit
}
