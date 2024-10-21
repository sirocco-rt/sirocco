
/***********************************************************/
/** @file  windsave2table.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  A standalone routine for writing a standard set of data
 * from a windsavefile to ascii tables all readable with as astropy.io.ascii
 *
 * This routine is run from the command line, as follows
 *
 * windsave2table  windsave_root
 *
 * where windsave_root is the root name for a sirocco run, or more precisely
 * the rootname of a windsave file, as the .pf file is not read.
 *
 * The routine reads the windsavefile and then writes out a selected 
 * set of variables into a variety of number of files, each of which
 * are readable as astropy tables, where each row corresponds to one
 * element of the wind.   
 *
 * All of the files begin
 * with the rootname and the domain number.
 * If Python run that created the windsave file
 * has multiple domains a sepearate output file is created for each
 * domain
 *
 * All of the files, begin with
 *
 *        x        z    i    j inwind 
 *
 * where x and z correspond to the position of the cell, i and j correspond to 
 * the cell number and inwind says whether the cell was in the wind, so that it is
 * fairly straightforward to create plotting routines for various parameters
 * that are contained in the remaining columns.  
 *
 * The files include a so-called masterfile which records basic parameters
 * like n_e, velocity,  rho, t_e, t_r, and ionization fractions
 * for a few key ions in each cell when the windsave file was written
 *
 * ### Notes ###
 *
 * Whereas swind is intended to be run interactively, windsave2table is
 * entirely hardwired so that it produces a standard set of output
 * files.  To change the outputs one has to modify the routine
 *
 * This file just contains the driving routine.  All of the 
 * real work is carried out in windsave2table_sub.c  Indeed so 
 * little of the real work is done here that it might be sensible
 * to move some portion of that code here.
 *
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
 * @brief Parse the command line arguments given to windsave2table
 *
 * @param[in] int argc        The number of arguments in the command line
 * @param[in] char *argv[]    The command line arguments
 * @param[out] char root[]    The rootname of the Python simulation
 *
 * @return void
 *
 * @details
 *
 * Parse arguments provided in the command line whilst invoking windsave2table.
 *
 * The allowed switches include 
 *
 *  --version    Print out information about the version number
 *
 *  -d     Write out ion densisties, rather than ion fractions in the cell
 *  -s     Write out the number of scatters per unit volume  by an ion in a cell, instead of the 
 *         ion fraction
 *  -x     windcell Writes out the detailed spectra in a specific windcell 
 *  -xall  Writes out the detiled windcell spectra for all of the cells that are acutally in 
 *         the wind
 *
 * The switches only affect the ion tables not the master table
 * This was originally implemented to enable somebody to query which version of
 * Python windsave2table was compiled with. Works in a similar fashion to how
 * the version information is stored and viewed in Python.
 * 
 * 
 *
 **********************************************************/

char windsave2table_help[] = "Usage: windsave2table [-r or -s] [-a] [-x wincell_no] [-xall] [-h] [--version] rootname \n\
-d             Return densities instead of ion fraction in ion tables \n\
-s             Return number of scatters per unit volume of an ion instead if ion fractions \n\
-a             Print additional tables with more information about ions  \n\
-edge          Include edge cells in the various tables  \n\
--version      Return version info and quit \n\
-x windcell    In addition to the normal tables, print out the detailed spectra in a specific windcell\n\
-xall          In addition to the normal tables, print out a large file containing all of the detailed cell spectra\n\
               for those cells that are in the wind\n\
-h             get this help message and quit\n\
";

void
parse_arguments (int argc, char *argv[], char root[], int *ion_switch, int *spec_switch, int *edge_switch)
{
  int i;
  char *fget_rc;
  char input[LINELENGTH];

  *ion_switch = 0;
  *spec_switch = -1;
  *edge_switch = FALSE;


  if (argc == 1)
  {
    printf ("Root for wind file :");
    fget_rc = fgets (input, LINELENGTH, stdin);
    if (!fget_rc)
    {
      printf ("No root file provided, exiting\n\n");
      printf ("%s", windsave2table_help);

      exit (0);
    }
    get_root (root, input);
  }
  else
  {
    for (i = 1; i < argc; i++)
    {
      if (!strcmp (argv[i], "--version"))
      {
        printf ("Python Version %s\n", VERSION);
        printf ("windsave2table built from git commit hash %s\n", GIT_COMMIT_HASH);
        if (GIT_DIFF_STATUS)
          printf ("This version was compiled with %i files with uncommitted changes.\n", GIT_DIFF_STATUS);
        exit (0);
      }
      else if (!strncmp (argv[i], "-xall", 5))
      {
        *spec_switch = -2;

      }
      else if (!strncmp (argv[i], "-x", 2))
      {
        if (sscanf (argv[i + 1], "%d", spec_switch) != 1)
        {
          Error ("sirocco: wind cell number  after -x switch\n");
          exit (1);
        }

        printf ("Cell spectrum in wind cell %d will be printed\n", *spec_switch);
        i++;
      }
      else if (!strncmp (argv[i], "-d", 2))
      {
        *ion_switch = 1;
        printf ("Ion outputs will be densities");
      }
      else if (!strncmp (argv[i], "-s", 2))
      {
        *ion_switch = 2;
        printf ("Ion outputs will be the number of scatters for this ion in a cell");
      }
      else if (!strncmp (argv[i], "-a", 2))
      {
        *ion_switch = 99;
        printf ("Various files detailing information about each ion in a cell will be created\n");
      }
      else if (!strncmp (argv[i], "-edge", 5))
      {
        *edge_switch = TRUE;
        printf ("Files will include edge cells\n");
      }
      else if (!strncmp (argv[i], "-h", 2))
      {
        printf ("%s", windsave2table_help);
        exit (0);
      }
      else if (!strncmp (argv[i], "-", 1))
      {
        printf ("Unknown switch %s\nExiting\n\n", argv[i]);
        printf ("%s", windsave2table_help);
        exit (0);
      }
      else
      {
        strcpy (input, argv[argc - 1]);
        get_root (root, input);
      }
    }
  }
}




/**********************************************************/
/** 
 * @brief      windsave2table writes key variables in a windsave file 
 * to an astropy table calculated Python.  This is the  main routine.
 *
 * @param [in] int  argc   The number of argments in the command line
 * @param [in] char *  argv[]   The command line
 * @return     Always returns 0, unless the windsave file is not
 * found in which case the routine will issue and error before
 * exiting.  
 *
 * @details
 * argc and argv[] are the standard variables provided to main 
 * in a c-program.  Only the first command line argument is 
 * parsed and this should be rootname of the windsave file
 *
 * ### Notes ###
 *
 * This routine is a supervisory routine. The real work 
 * is in do_windsave2table
 *
 * The routine does not read the .pf file.  It does read  
 * the windsave file and the associated atomic data file
 * before calling do_windsave2table.
 *
 **********************************************************/

int
main (argc, argv)
     int argc;
     char *argv[];
{
  char root[LINELENGTH], xroot[LINELENGTH];
  char outputfile[LINELENGTH];
  char windsavefile[LINELENGTH];
  char parameter_file[LINELENGTH];
  int ion_switch;
  int spec_switch;
  int edge_switch;
  int ndom;


  strcpy (parameter_file, "NONE");

  /* Next command stops Debug statements printing out in swind */
  Log_set_verbosity (3);

  /*
   * EP: added some extra argument parsing for windsave2table - specifically
   * because I was having some trouble with knowing when windsave2table was
   * last compiled and on what commit this was
   */

  parse_arguments (argc, argv, root, &ion_switch, &spec_switch, &edge_switch);

  printf ("Reading data from file %s\n", root);

  /* Now create the names of all the files which will be written */

  strcpy (windsavefile, root);
  strcpy (outputfile, root);

  strcat (windsavefile, ".wind_save");
  strcat (outputfile, ".txt");


/* Read in the wind file */

  zdom = calloc (MAX_DOM, sizeof (domain_dummy));
  if (zdom == NULL)
  {
    printf ("Unable to allocate memory for domain\n");
    return EXIT_FAILURE;
  }

  if (wind_read (windsavefile) < 0)
  {
    Error ("swind: Could not open %s", windsavefile);
    exit (0);
  }

  printf ("Read wind_file %s\n", windsavefile);
  printf ("Read Atomic data from %s\n", geo.atomic_filename);

  do_windsave2table (root, ion_switch, edge_switch);

  if (spec_switch == -2)
  {
    for (ndom = 0; ndom < geo.ndomain; ndom++)
    {
      if (geo.ndomain > 1)
      {
        sprintf (xroot, "%.200s.%d", root, ndom);
      }
      else
      {
        sprintf (xroot, "%s", root);
      }

      create_big_detailed_spec_table (ndom, xroot);
    }
  }

  else if (spec_switch >= 0)
  {
    printf ("Cell spectrum in wind cell %d will be printed\n", spec_switch);
    create_detailed_cell_spec_table (spec_switch, root);
  }


  return (0);
}
