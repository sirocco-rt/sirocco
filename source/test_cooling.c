
/***********************************************************/
/** @file  inspect_wind.c
 * @author ksl
 * @date   October, 2021
 *
 * @brief  Routines to inspect variables in a wind structure 
 *
 *###Notes###
 * This is intended just as a diagnostic routine 
 * so that one can print out whatever variables in
 * a windstrucutre one wants in order to diagnose
 * a problem.  It was written so that we could inspect
 * some of the macro atom variables in paralell mode
 * in diagnosing issue #898 and #910, but anyon 
 * should change it so that other problems might 
 * be addressed.
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"



char inroot[LINELENGTH], outroot[LINELENGTH], model_file[LINELENGTH], folder[LINELENGTH];
int model_flag, ksl_flag, cmf2obs_flag, obs2cmf_flag;

double line_matom_lum_single (double lum[], PlasmaPtr xplasma, int uplvl);
int line_matom_lum (int uplvl);
int create_matom_level_map ();

/**********************************************************/
/**
 * @brief      parses the command line options
 *
 * @param [in]  int  argc   the number of command line arguments
 * @param [in]  char *  argv[]   The command line arguments
 *
 *
 * ###Notes###
 *
 * The general purpose of each of the command line options
 * should be fairly obvious from reading the code.
 *
 *
 * Although this routine uses the standard Log and Error commands
 * the diag files have not been created yet and so this information
 * is really simply written to the terminal.
 *
 **********************************************************/

int
xparse_command_line (argc, argv)
     int argc;
     char *argv[];
{
  int j = 0;
  int i;
  char dummy[LINELENGTH];
  int mkdir ();
  char *fgets_rc;


  sprintf (outroot, "%s", "new");

  model_flag = ksl_flag = obs2cmf_flag = cmf2obs_flag = 0;

  if (argc == 1)
  {
    printf ("Parameter file name (e.g. my_model.pf, or just my_model):");
    fgets_rc = fgets (dummy, LINELENGTH, stdin);
    if (!fgets_rc)
    {
      printf ("Input rootname is NULL or invalid\n");
      exit (1);
    }
    get_root (inroot, dummy);
  }
  else
  {

    for (i = 1; i < argc; i++)
    {
      if (strcmp (argv[i], "-out_root") == 0)
      {
        if (sscanf (argv[i + 1], "%s", dummy) != 1)
        {
          printf ("sirocco: Expected out_root after -out_root switch\n");
          exit (0);
        }

        get_root (outroot, dummy);
        i++;
        j = i;

      }
      if (strcmp (argv[i], "-model_file") == 0)
      {
        if (sscanf (argv[i + 1], "%s", dummy) != 1)
        {
          printf ("sirocco: Expected a model file containing density, velocity and temperature after -model_file switch\n");
          exit (0);
        }
        get_root (model_file, dummy);
        i++;
        j = i;
        printf ("got a model file %s\n", model_file);
        model_flag = 1;
      }
      else if (strcmp (argv[i], "-ksl") == 0)
      {
        printf ("Carrying out a simple hard wired ion modification\n");
        ksl_flag = 1;
      }
      else if (strcmp (argv[i], "--dry-run") == 0)
      {
        modes.quit_after_inputs = 1;
        j = i;
      }
      else if (strcmp (argv[i], "-cmf") == 0)
      {
        obs2cmf_flag = 1;
      }
      else if (strcmp (argv[i], "-obs") == 0)
      {
        cmf2obs_flag = 1;
      }
      else if (strncmp (argv[i], "-", 1) == 0)
      {
        printf ("sirocco: Unknown switch %s\n", argv[i]);
        exit (0);
      }
    }

    /* The last command line variable is always the windsave file */

    if (j + 1 == argc)
    {
      printf ("All of the command line has been consumed without specifying a parameter file name, so exiting\n");
      exit (1);
    }
    strcpy (dummy, argv[argc - 1]);
    get_root (inroot, dummy);

  }

  return (0);
}



/* An externall pointer reference used by zero_emit.  */
PlasmaPtr xxxplasma;



double
xcalc_te (PlasmaPtr xplasma, double tmin, double tmax)
{
  double z1, z2;
  int ierr = FALSE;


  /* we assign a plasma pointer here to a fixed structure because
   * we need to call zbrent and we cannot pass the xplasma ptr directly
   */

  xxxplasma = xplasma;

  xxxplasma->heat_tot += xxxplasma->heat_ch_ex;

  xplasma->t_e = tmin;
  z1 = zero_emit (tmin);
  xplasma->t_e = tmax;
  z2 = zero_emit (tmax);

  /* The way this works is that if we have a situation where the cooling
   * at tmax and tmin brackets the heating, then we use zbrent to improve
   * the estimated temperature, but if not we chose the best direction
   */

  if ((z1 * z2 < 0.0))
  {                             // Then the interval is bracketed
    xplasma->t_e = zero_find (zero_emit2, tmin, tmax, 50., &ierr);
    if (ierr)
    {
      Error ("calc_te: zero_find failed to find a temperature\n");
    }

  }
  else if (fabs (z1) < fabs (z2))
  {
    xplasma->t_e = tmin;
  }
  else
  {
    xplasma->t_e = tmax;
  }
  /* With the new temperature in place for the cell, get the correct value of heat_tot.
     SS June  04 */

  /* At this point we know the temperature that balances heating and cooling
   * within the constraints set by tmin and tmax.
   */


  /* Update heat_tot and heat_lines for macro_bb_heating at the new temperature. 
   * We subtract the current value and then compute at the new temperature and
   * add this back */

  xplasma->heat_tot -= xplasma->heat_lines_macro;
  xplasma->heat_lines -= xplasma->heat_lines_macro;

  xplasma->heat_lines_macro = macro_bb_heating (xplasma, xplasma->t_e);

  xplasma->heat_tot += xplasma->heat_lines_macro;
  xplasma->heat_lines += xplasma->heat_lines_macro;

  /* Similaryly for macro_atom_bf_heating */

  xplasma->heat_tot -= xplasma->heat_photo_macro;
  xplasma->heat_photo -= xplasma->heat_photo_macro;

  xplasma->heat_photo_macro = macro_bf_heating (xplasma, xplasma->t_e);

  xplasma->heat_tot += xplasma->heat_photo_macro;
  xplasma->heat_photo += xplasma->heat_photo_macro;


  return (xplasma->t_e);

}

/**********************************************************/
/**
 * @brief      the main routine which carries out the effort
 *
 * @param [in]  int  argc   the number of command line arguments
 * @param [in]  char *  argv[]   The command line arguments
 *
 *
 * ###Notes###
 *
 * This routine oversees the effort.  The basic steps are
 *
 * - parse the command line to get the names of files
 * - read the old windsave file
 * - read the densities from in this case H
 * - modify the densities
 * - write out the new windsave file
 *
 *
 **********************************************************/


int
main (argc, argv)
     int argc;
     char *argv[];
{

  char infile[LINELENGTH], outfile[LINELENGTH];
  int n, i;
  FILE *fptr, *fopen ();
  int ii, jj, ndom, nnwind;
  int mkdir ();


  xparse_command_line (argc, argv);

  sprintf (infile, "%.150s.wind_save", inroot);
  sprintf (outfile, "%.150s.txt", inroot);

  printf ("Reading %s and writing to %s\n", infile, outfile);

  zdom = calloc (MAX_DOM, sizeof (domain_dummy));
  if (zdom == NULL)
  {
    printf ("Unable to allocate memory for domain\n");
    return EXIT_FAILURE;
  }

  wind_read (infile);

  if (nlevels_macro == 0)
  {
    printf ("Currently this routine only looks at macro atom values, and this is a simple atom file\n");
    exit (0);
  }



  fptr = fopen (outfile, "w");

  fprintf (fptr, "# Results for %s\n", infile);


  printf ("te %.3e\n", plasmamain[0].t_e);
  fprintf (fptr, "te %.3e\n", plasmamain[0].t_e);

  double t, t_new;
  t = plasmamain[0].t_e;

  t_new = xcalc_te (&plasmamain[0], 0.7 * t, 1.3 * t);
  printf ("te_new %.3e\n", t_new);
  fprintf (fptr, "te_new %.3e\n", t_new);

  fclose (fptr);

  exit (0);

}
