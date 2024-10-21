/***********************************************************/
/** @file  setup_files.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  A series of routines mostly associated with opening
 * files needed by Python
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
 * @brief      init_files
 *
 * @param [in] int  restart_stat   A flag indicating whether 
 * this run is a restart or a new run.  TRUE for a restart,
 * FALSE for a new run.
 *
 * @return     Always returns 0                   
 *
 * @details
 *
 * This routine open the log file and signal files. 
 *
 * On restarts, it checks whether there is an appropriate
 * windsave file to use.  
 *
 *
 * ### Notes ###
 *
 * @bug This routine was refactored out of sirocco.c in an attempt
 * to make that routine simpler, but it is not obvious that 
 * what done here actually makes sense. For example, although
 * the routine checks if the windsave file exists, and sets
 * restart_stat to 0 if it does not.  This information is
 * not transmitted back to any other portion of the program.
 *
 **********************************************************/

int
init_log_and_windsave (restart_stat)
     int restart_stat;
{
  FILE *fopen (), *qptr;

  if (restart_stat == FALSE)
  {                             // Then we are simply running from a new model
    xsignal_rm (files.root);
    xsignal (files.root, "%-20s %s \n", "START", files.root);
    Log_init (files.diag);
  }
  else
  {
    /* Note that although we check that we can open the windsave file, it is not read here.   */

    strcpy (files.windsave, files.root);
    strcat (files.windsave, ".wind_save");
    qptr = fopen (files.windsave, "r");

    if (qptr != NULL)
    {
      /* Then the file does exist and we can restart */
      fclose (qptr);
      xsignal (files.root, "%-20s %s\n", "RESTART", files.root);
      Log_append (files.diag);
    }
    else
    {
      /* It does not exist and so we start from scratch */
      xsignal_rm (files.root);
      xsignal (files.root, "%-20s %s \n", "START", files.root);
      Log_init (files.diag);
    }
  }

  return (0);
}




/**********************************************************/
/** 
 * @brief      Use the root name of the .pf file to create
 * names for the various output files written by Python
 *
 * @return     Always returns 0
 *
 * @details
 * 
 * This names of the files sirocco uses are stored in he filenames
 * structure.  This routine uses the root of the .pf file to
 * create names for all of the files that Python will write and
 * stores them in the files structure
 *
 * The routine also opens the parameter file if it exists
 *
 * ### Notes ###
 *
 * The routine establishes the names of files, but does not
 * open them.
 * 
 * 
 *
 **********************************************************/

int
setup_created_files ()
{
  int opar_stat;

  opar_stat = 0;                /* Initialize opar_stat to indicate that if we do not open a rdpar file, 
                                   the assumption is that we are reading from the command line */


  if (strncmp (files.root, "stdin", 5) == 0 || strncmp (files.root, "rdpar", 5) == 0 || files.root[0] == ' ' || strlen (files.root) == 0)
  {
    strcpy (files.root, "mod");
    Log ("Proceeding in interactive mode\n Output files will have rootname mod\n");
  }

  else
  {
    strcpy (files.input, files.root);
    strcat (files.input, ".pf");

    if ((opar_stat = opar (files.input)) == 2)
    {
      Log ("Reading data from file %s\n", files.input);
    }
    else
    {
      Log ("Creating a new parameter file %s\n", files.input);
    }

  }


  /* Now create the names of all the files which will be written.  Note that some files
     have the same root as the input file, while others have a generic name of sirocco.
     This is intended so that files which you really want to keep have unique names, while
     those which are for short-term diagnostics are overwritten.   */

  strcpy (files.lwspec, files.root);    //generated photon in log space

  strcpy (files.lwspec_wind, files.root);

  strcpy (files.spec, files.root);
  strcpy (files.lspec, files.root);

  strcpy (files.spec_wind, files.root);
  strcpy (files.lspec_wind, files.root);

  strcpy (files.new_pf, files.root);
  strcat (files.new_pf, ".out.pf");


  strcpy (files.windrad, "sirocco");
  strcpy (files.windsave, files.root);
  strcpy (files.specsave, files.root);

  /* save sirocco.phot and disk.diag files under diag_root folder */
  strcpy (files.phot, files.diagfolder);
  strcat (files.phot, "sirocco");

  sprintf (files.disk, "%.100s/%.100s", files.diagfolder, files.root);

  strcat (files.lwspec, ".log_spec_tot");

  strcat (files.lwspec_wind, ".log_spec_tot_wind");


  strcat (files.spec, ".spec");
  strcat (files.lspec, ".log_spec");

  strcat (files.spec_wind, ".spec_wind");
  strcat (files.lspec_wind, ".log_spec_wind");


  strcat (files.windrad, ".wind_rad");
  strcat (files.windsave, ".wind_save");
  strcat (files.specsave, ".spec_save");
  strcat (files.phot, ".phot");
  strcat (files.disk, ".disk.diag");



  return (opar_stat);
}
