
/***********************************************************/
/** @file  swind.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  swind is a program which can be used to display various parameters of a wind
as calculated by sirocco.  This is the main routine.
 *
 *
 * Arguments:
 *
 * swind [-h] [-s] [-p parameter_file] [root]
 *
 * where
 * * -h 	prints out a short help file and exits (see help routine below)
 * * -s	causes certain parameters in the windsave file to be printed out
 * as individual ascii files after which the program exits
 * * -d	In cases where a windsavefile was made for each ionization
 * cycle, this prints out ascii files for "certain" parameters
 * for each of the ionization cycles, as well as the ascii files
 * for the final ionization cycle (if windsave files for each
 * ionization cycle were not created then -s and -d are
 * equivalent).
 * * -p parameter file
 * Instead of reading the choices from the command line read them
 * from a parameter file
 * root	optional root name of wind_save file.  If this is not given,
 * the user is queried for this
 *
 *
 * Description:
 *
 * Py_wind simply reads and then displays portions of the wind file created by sirocco.
 * It can select various parameters from the wind and it can write them to files so that
 * the variables can be plotted.
 *
 * The normal mode of running swind is to run it interactively.  As you run it interactively
 * the variables you select are displayed ont the screen.  The variables that you display
 * can also be written to files (depending on the answer to the question Make_files)
 *
 * The commands that were executed in the interactive will be stored in swind.pf (if you end
 * with a "q", and not an EOF response to the choice question)  EOF terminates the program
 * at that point before the command file is written to swind.pf
 *
 * The swind.pf file is useful if you want to run the exact same set of commands on another
 * windfile. You should rename swind.pf to something_else.pf and run swind on that
 * data set using the -p something_else.pf option
 *
 * The command_line switches -d and -s are intended to produce a standard set of output files sufficient
 * for many purposes.
 *
 * ###Notes####
 *
 * The files that are produced  can contain either the original gridding which was used by sirocco, in which
 * case the file prefix will be "x.", or "z", in which case it will be regridded to a
 * linear array.  This option is intended so one can create a contour plot more easily.
 *
 * The structure of swind is intended to make it easy to add additional variables to print to the
 * scen or a file.  Basically to add a new option, one needs to update the list of possibilities,
 * add a new choice statement in one_choice, and then write the routine that grabs a new variable (or
 * calculates a value from several), and add this to either swind_subs or swind_macro_subs.
 *
 * Because this was intended as a diagnositic routine, some of the more esoteric choices may seem
 * a bit odd today.
 *
 * IMPORTANT: swind has not been adapted to handle multiple domains.  Use windsave2table instead.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

//char *choice_options;


/**********************************************************/
/** 111125 - ksl - Replaced print statements giving choices with a string. The point is to be able to include
 * the choices in the help string  Note carefullly the format if you revise this
 * lines should end with \n\  to make the string continue.  Do not leave any trailing spaces after the last
 * \ to avoid warnings
 *
 **********************************************************/
char *choice_options = "\n\
    1=onefile summary 2=all ions in a given cell\n\
 	 n=ne,  R=rho,  v=vel,         i=ion info, j=ave_tau, f=ave_freq, p=nphot, S=sim_alpha\n\
   r=t_r, t=t_e,  w=rad_weight,  s=vol,     l=lum,     C=cooling/heating,  b=adiabatic cooling\n\
   a=abs,         L=line lum,    g=photo,   h=recomb,  k=tau H,     m=F_rad, x=total, y=mod_te,\n\
   o=overview,    e=everything,  P=Partial emission meas,  I=Ionisation parameter\n\
   W=wind_region, D=dvds_ave,    X=position summary, M=macro atom info *=shock heating \n\
   d=convergence status  E=convergence_all_info   B=PlasmaPtr  J=Radiation density\n\
   H=All Heating and Cooling mechanisms in one shot  O=Spectral model parameters S=Spectral models\n\
   z=Zoom,u=unZoom,Z=switch to/from raw and yz projected modes, F=Create files, A=Change file write defaults\n\
   #=Wind grid    N=new.windfile q=quit (preferred over EOF) &=Coll Strengths Q=switch domain\n";


/**********************************************************/
/**
 * @brief      swind is a program which can be used to display various parameters of a wind
 * 		as calculated by sirocco.  This is the  main routine
 *
 * @param [in] int  argc   The number of command line arguments
 * @param [in] char *  argv[]   The command line arguments
 * @return     Always returns 0
 *
 * @details
 * Py_wind simply reads and then displays portions of the wind file created by sirocco.
 * 	It can select various parameters from the wind and it can write them to files so that
 * 	the variables can be plotted.
 *
 * 	The normal mode of running swind is to run it interactively.  As you run it interactively
 * 	the variables you select are displayed ont the screen.  The variables that you display
 * 	can also be written to files (depending on the answer to the question Make_files)
 *
 * 	The commands that were executed in the interactive will be stored in swind.pf (if you end
 * 	with a "q", and not an EOF response to the choice question)  EOF terminates the program
 * 	at that point before the command file is written to swind.pf
 *
 * 	The swind.pf file is useful if you want to run the exact same set of commands on another
 * 	windfile. You should rename swind.pf to something_else.pf and run swind on that
 * 	data set using the -p something_else.pf option
 *
 * 	The command_line switches -d and -s are intended to produce a standard set of output files sufficient
 * 	for many purposes.
 *
 * 	The main routine gets the input data, and repeately presents the user with a set of choices (in interactive
 * 	mode) until the users dicides she has had enough.  The case statements that are used to select what to do based
 * 	on user inputs are (now) contained in the subroutine one_choice
 *
 * ### Notes ###
 * The files that are produced  can contain either the original gridding which was used by sirocco, in which
 * case the file prefix will be "x.", or "z", in which case it will be regridded to a
 * linear array.  This option is intended so one can create a contour plot more easily.
 *
 **********************************************************/

int
main (argc, argv)
     int argc;
     char *argv[];
{


  int i;
  int ochoice;
  char c;

  char root[LINELENGTH], input[LINELENGTH], wspecfile[LINELENGTH], specfile[LINELENGTH];
  char windradfile[LINELENGTH], windsavefile[LINELENGTH];
  char parameter_file[LINELENGTH];
  char photfile[LINELENGTH];
  char *fgets_rc;
  double freq;
  int interactive;


  // swind uses rdpar, but only in an interactive mode. As a result
  // there is no associated .pf file

  interactive = 1;              /* Default to the standard operating mofe for swind */
  strcpy (parameter_file, "NONE");

  /* Next command stops Debug statements printing out in swind */
  Log_set_verbosity (3);

  /* Parse the command line.  It is Important that swind_help
   * below is updated if this changes
   */

  if (argc == 1)
  {
    printf ("Root for wind file :");
    fgets_rc = fgets (input, LINELENGTH, stdin);
    if (!fgets_rc)
    {
      Error ("Input rootname is NULL or EOF\n");
      exit (1);                 // Exit if NULL returned
    }
    get_root (root, input);
  }
  else
  {
    for (i = 1; i < argc; i++)
    {
      if (strcmp (argv[i], "-h") == 0)
      {
        swind_help ();
      }
      else if (strcmp (argv[i], "-d") == 0)
      {
        interactive = -1;
      }

      else if (strcmp (argv[i], "-s") == 0)
      {
        interactive = 0;
      }
      else if (strcmp (argv[i], "-p") == 0)
      {
        interactive = 0;
        i = i + 1;
        strcpy (parameter_file, argv[i]);
      }
      else if (strncmp (argv[i], "-", 1) == 0)
      {
        Error ("swind: unknown switch %s\n", argv[i]);
        swind_help ();
      }
    }

    strcpy (input, argv[argc - 1]);
    get_root (root, input);
  }


  printf ("Reading data from file %s\n", root);

  /* Now create the names of all the files which will be written */

  strcpy (wspecfile, root);
  strcpy (specfile, root);
  strcpy (windradfile, root);
  strcpy (windsavefile, root);
  strcpy (photfile, root);

  strcat (wspecfile, ".spec_tot");
  strcat (specfile, ".spec");
  strcat (windradfile, ".wind_rad");
  strcat (windsavefile, ".wind_save");
  strcat (photfile, ".phot");


  /* Initialize other variables here */

  swind_project = 1;          // The default is to try to project onto a yz plane

/* Read in the wind file */

/* Note that wind_read allocates the space for the WindPtr array.  The
reason that this is done in wind_read is because wind_read also reads
the geo structure.  The geo struc contains the dimensions of the wind
array and so until it is read the space for the WindPtr structure cannot
be allocated.  Howver, one cannot assign w in a subroutine.  Therefore
the old call, namely wind_read(w,windsavefile, will not work.  The statement
w=wmain is a way around this since wmain is external and therefore can
be assigned.  Then w can be set to this ptr and all is restored. The
use of w is endemic in the program. and it is always called through main.
I did not change this now.  Though it could be done.  02apr ksl */

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

/* aaa is used to store variable for writing to files for the purpose of plotting*/
  aaa = calloc (sizeof (freq), NDIM2);

  printf ("Read wind_file %s\n", windsavefile);

  get_atomic_data (geo.atomic_filename);

  printf ("Read Atomic data from %s\n", geo.atomic_filename);


  printf ("There were %d domains.  Starting with domain 0\n", geo.ndomain);
  current_domain = 0;
  printf ("JM COORD TYPE %d\n\n", zdom[current_domain].coord_type);
  /*Set the current domain to zero */
/* Produce a standard set of output files and exit*/
  if (interactive == 0 && strcmp (parameter_file, "NONE") == 0)
  {
    zoom (1);                   /* This affects the logfile */
    ochoice = 1;
    complete_file_summary (wmain, root, ochoice);
    exit (0);
  }
  else if (interactive == -1)
  {
    /* In cases, where the windsave file was written out for
     * eaach ionization cycle Write the sumary ascii files
     * for each of the ionization cycles as
     * as well as the final cycle */
    zoom (1);                   /* This affects the logfile */
    ochoice = 1;
    complete_file_summary (wmain, root, ochoice);
    i = 0;
    strcpy (root, "");
    sprintf (root, "sirocco%02d", i);
    strcpy (windsavefile, "");
    sprintf (windsavefile, "sirocco%02d.wind_save", i);
    while (wind_read (windsavefile) > 0)
    {
      Log ("Trying %s %s\n", windsavefile, root);
      complete_file_summary (wmain, root, ochoice);
      strcpy (root, "");
      sprintf (root, "sirocco%02d", i);
      strcpy (windsavefile, "");
      sprintf (windsavefile, "sirocco%02d.wind_save", i);
      i++;
    }
    exit (0);
  }




  if (strcmp (parameter_file, "NONE") != 0)
  {
    zoom (1);                   /* This affects the logfile */
    ochoice = 1;
    opar (parameter_file);
  }


/* Choices */
  ochoice = 0;
  rdint ("Make_files(0=no,1=original,2=regrid_to_linear)", &ochoice);
  c = 'i';
  zoom (1);



/* 111126 - Consolodated choice statements into a string that is an external variable
 * so that it can be printed out whenever necessary.
 */

  printf ("%s\n", choice_options);
  printf ("Model %s   :\n", root);
  rdchar ("Choice", &c);

  while (c != EOF)
  {
    if (c == 'Q')
    {
      printf ("There were %d domains.\n", geo.ndomain);
      rdint ("Domain_to_examine", &current_domain);
      if (current_domain < 0 || current_domain >= geo.ndomain)
      {
        printf ("Unkown Domain (forcing current_domain to 0\n");
        current_domain = 0;
      }
      else
      {
        printf ("Swithching to domain %d\n", current_domain);
      }
      zoom (1);                 /* Unzoom */
      rdchar ("Choice", &c);
    }

    one_choice (c, root, ochoice);
    printf ("%s\n", choice_options);
    rdchar ("Choice", &c);
  }


  return (0);
}




/**********************************************************/
/**
 * @brief      Process a request to display one variable
 *
 * @param [in] char  choice   A character indicating the
 * variable to be displayed
 * @param [in] char *  root   The rootname of the windsave file
 * @param [in] int  ochoice   An integer indicating whether the data
 * are only to be displayed, only to be written to a file, or both.
 * @return     Always returns 0
 *
 * @details
 * The routine simply contains a large case statement for
 * interpreting what to do based on the variable choice
 *
 * ### Notes ###
 *
 * In some cases, additional imputs are required.
 *
 **********************************************************/

int
one_choice (choice, root, ochoice)
     char choice;
     char *root;
     int ochoice;
{
  double lambda, freq;
  int n, istate, iswitch;
  char windsavefile[LINELENGTH];

  iswitch = 0;

  /* JM 1312 --initialise variables to avoid compilation warnings */
  istate = 0;
  n = 0;


  switch (choice)
  {
  case 'a':                    /* Energy absorbed */
    abs_summary (wmain, root, ochoice);
    break;
  case 'A':                    // Change the file defaults
    rdint ("Make_files(0=no,1=original,2=regrid_to_linear)", &ochoice);
    break;
  case 'b':                    /*Adiabatic cooling */
    adiabatic_cooling_summary (wmain, root, ochoice);
    break;
  case 'B':
    plasma_cell (wmain, root, ochoice);
    break;
  case 'c':
    flux_summary (wmain, root, ochoice);
    break;
  case 'C':                    /*the ratio cooling to heating */
    coolheat_summary (wmain, root, ochoice);
    break;
  case 'd':
    convergence_summary (wmain, root, ochoice);
    break;
  case 'D':                    /* dvds summary */
    dvds_summary (wmain, root, ochoice);
    break;
  case 'E':
    convergence_all (wmain, root, ochoice);
    break;
  case 'e':                    /* print out everything about an element */
    wind_element (wmain);
    break;
  case 'f':                    /* Electron summary */
    freq_summary (wmain, root, ochoice);
    break;
  case 'F':                    /* Complete file summary */
    complete_file_summary (wmain, root, ochoice);
    break;
  case 'g':                    /*n photo */
    photo_summary (wmain, root, ochoice);
    break;
//  case 'G':    Removed relevant data from sirocco so option removed May 18         /* inner shell summary */
//    inner_shell_summary (wmain, root, ochoice);
//    break;
  case 'h':                    /*n photo */
    Log ("Don't get discouraged.  This takes a little while!");
    recomb_summary (wmain, root, ochoice);
    break;
  case 'H':                    /* heating and cooling mechanisms breakdown */
    heatcool_summary (wmain, root, ochoice);
    break;
  case 'i':                    /* Allow user to display information about ions in the wind */

    rdint ("Ion_info_type(0=fraction,1=density,2=scatters,3=abs", &iswitch);

    n = 6;
    istate = 4;

    while (rdint ("element(0=return)", &n) != EOF)
    {
      if (n <= 0)
        break;
      rdint ("ion", &istate);
      ion_summary (wmain, n, istate, iswitch, root, ochoice);   // 0 implies ion fractions
    }
    break;
  case 'I':
    IP_summary (wmain, root, ochoice);
    break;

  case 'j':                    /* Calculate the average tau at the center of a cell */
    n = 6;
    istate = 4;
    lambda = 1550;

    rddoub ("wavelength", &lambda);
    freq = VLIGHT / (lambda * 1.e-8);

    while (rdint ("element(0=return)", &n) != EOF)
    {
      if (n <= 0)
        break;
      rdint ("ion", &istate);
      tau_ave_summary (wmain, n, istate, freq, root, ochoice);
    }
    break;
  case 'J':                    /* radiation density in cell */
    J_summary (wmain, root, ochoice);
    break;
  case 'k':                    /* tau at H edge */
    tau_h_summary (wmain, root, ochoice);
    break;
  case 'K':                    /* cell J split by direct photons and scattered photons */
    J_scat_summary (wmain, root, ochoice);
    break;
  case 'l':                    /* Lum of shell */
    lum_summary (wmain, root, ochoice);
    break;
  case 'L':                    /*Line emission */
    line_summary (wmain, root, ochoice);
    break;
  case 'm':                    /* Radiation force */
    mo_summary (wmain, root, ochoice);
    break;
  case 'M':
    macro_summary (wmain, root, ochoice);
    break;
  case 'n':                    /* Electron summary */
    electron_summary (wmain, root, ochoice);
    break;
  case 'N':                    /* Read a different wind save file */
    rdstr ("New.rootname", root);
    strcpy (windsavefile, root);
    strcat (windsavefile, ".wind_save");
    if (wind_read (windsavefile) < 0)
    {
      Error ("one_choice: Could not read %s", windsavefile);
    }


/* aaa is used to store variable for writing to files for the purpose of plotting*/
    if (aaa != NULL)
    {
      free (aaa);
    }
    aaa = calloc (sizeof (freq), NDIM2);

    printf ("Read wind_file %s\n", windsavefile);

    get_atomic_data (geo.atomic_filename);

    printf ("Read Atomic data from %s\n", geo.atomic_filename);

    break;
  case 'o':                    /* overview */
    overview (wmain, root);
    break;
  case 'O':                    /* spectral model parameters */
    model_bands (wmain, root, ochoice);
    break;
  case 'p':                    /* nphot summary */
    nphot_summary (wmain, root, ochoice);
    break;
  case 'P':                    /* Allow user to display information about the wind */

    n = 6;
    istate = 4;

    while (rdint ("element(0=return)", &n) != EOF)
    {
      if (n <= 0)
        break;
      rdint ("ion", &istate);
      partial_measure_summary (wmain, n, istate, root, ochoice);
    }
    break;
  case 'r':                    /* Temp summary */
    temp_rad (wmain, root, ochoice);
    break;
  case 'R':                    /* Rho summary */
    rho_summary (wmain, root, ochoice);
    break;
  case 's':                    /* Volume summary */
    vol_summary (wmain, root, ochoice);
    break;
  case 'S':
    alpha_summary (wmain, root, ochoice);
    break;
  case 't':                    /* Temp summary */
    temp_summary (wmain, root, ochoice);
    break;
  case 'T':
    thompson (wmain, root, ochoice);
    break;
  case 'v':                    /* Velocity summary */
    velocity_summary (wmain, root, ochoice);
    break;
  case 'V':                    /* Split of scatters in the cell between electron and resonant */
    nscat_split (wmain, root, ochoice);
    break;
  case 'w':                    /* inten weight summary */
    weight_summary (wmain, root, ochoice);
    break;
  case 'W':                    /*Show regions in the wind */
    wind_reg_summary (wmain, root, ochoice);
    break;
  case 'x':                    /*Total emission */
    total_emission_summary (root, ochoice);
    break;
  case 'X':                    /* Position summary */
    position_summary (wmain);
    break;
  case 'y':                    /* Recalculate temperatures */
    modify_te (wmain, root, ochoice);
    break;
  case 'Y':                    /* Split of photons from different sources */
    phot_split (wmain, root, ochoice);
    break;
  case 'z':                    /* inspect a specific region */
    zoom (0);
    break;
  case 'Z':                    /* Switch between raw and projected mode */
    if (swind_project == 0)
    {
      swind_project = 1;
      Log ("Switching to raw display, whatever the coordinate system");
    }
    else if (swind_project == 1)
    {
      swind_project = 0;
      Log ("Switching to projected y z display");
    }
    break;

    /* JM -- typing 1 gives you a summary of everything in one file with
       astropy.io.ascii compliant headers */

  case '1':
    complete_physical_summary (wmain, root, ochoice);   //
    break;


  case 'u':                    /* Go back to full image */
    zoom (1);
    break;
  case '2':
    complete_ion_summary (wmain, root, ochoice);        //
    break;
  case '#':
    grid_summary (wmain, root, ochoice);        //
    break;
  case '&':
    collision_summary (wmain, root, ochoice);   //
    break;
  case '*':
    shock_heating_summary (wmain, root, ochoice);
    break;

  case 'q':                    /* quit */
    /* Write out a parameterfile that gives all of the commands used in this run */
    cpar ("swind.pf");
    exit (0);
    break;

  }

  return (0);
}




/**********************************************************/
/**
 * @brief      print out some helpful information to the screen
 * 	and exits
 *
 * @return     N/A since the program exits once the help is
 * written
 *
 * @details
 * This routine simply writes information to the screen.
 * It should be updated whenever the program changes in
 * a major way.
 *
 * It is called if swind is invoked with -h, or if there
 * if it is given a switch that it does not understand.
 *
 * ### Notes ###
 * Unfortunately unlike python 
 * c is not self documenting
 *
 **********************************************************/

void
swind_help ()
{

  char *some_help;

/* Beginning of a string to provide some help for swind */
  some_help = "\
\n\
This program reads a wind save file created by sirocco and examine the wind structure\\
\n\
	Usage: swind [-h] [-s] [-p parameter_file] [root] \n\
\n\
	where\n\
		-h 	prints out a short help file and exits (see help routine below) \n\
		-s	causes certain parameters to be printed out in individual files \n\
			after which the program exits \n\
		-d	searches for the diagnostic wind save files and prints information \n\
			to ascii files \n\
		-p	Gets information from a parameter file \n\
		root	optional root name of wind_save file.  If this is not given, \n\
			the user is queried for this \n\
\n\
";

  printf ("%s\n", some_help);

  printf ("Choices are indicated below\n%s\n", choice_options);

  exit (0);
}
