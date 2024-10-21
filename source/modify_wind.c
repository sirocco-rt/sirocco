
/***********************************************************/
/** @file  modify_wind.c
 * @author ksl
 * @date   June, 2020
 *
 * @brief  Routines to modfify a wind structure to for example
 * change densities of certain ions
 *
 *###Notes###
 * At present this is just a prototype of a more general 
 * purpose routine to paint sections of a windsavefile
 * with new values.  At present the code just supports
 * changing densities, although it would be fairly 
 * straighfoward to extend to propeties like T_e
 *
 * At present the changes have to be hardcoded into 
 * main, and the program lacks an interface to for
 * example files of the tupe produced by windsave2table
 * 
 * The intent is to make the structure of this program
 * like the inverse of windsave2table.
 *
 * The initial commit that KSL write is 1485d5f6669704dde55e437f5c0cb061ec202653
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

char inroot[LINELENGTH], outroot[LINELENGTH], model_file[LINELENGTH];
int model_flag, ksl_flag, cmf2obs_flag, obs2cmf_flag;

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

  double *den;
  char name[LINELENGTH];        /* file name extension */
  char infile[LINELENGTH], outfile[LINELENGTH];
  int put_ion ();
  int apply_model ();
  int frame_transform ();
  int ndom;
  int i;

  ndom = 0;

  Log_set_verbosity (3);
  xparse_command_line (argc, argv);


  sprintf (infile, "%.150s.wind_save", inroot);
  sprintf (outfile, "%.150s.wind_save", outroot);

  printf ("Reading %s and writing to %s\n", infile, outfile);

  zdom = calloc (MAX_DOM, sizeof (domain_dummy));
  if (zdom == NULL)
  {
    printf ("Unable to allocate memory for domain\n");
    return EXIT_FAILURE;
  }

  wind_read (infile);

  if (model_flag)
  {
    apply_model (ndom, model_file);
  }
  if (obs2cmf_flag || cmf2obs_flag)
  {
    printf ("BOOM - we are going to be converting from one frame to another\n");
    frame_transform (ndom);
  }
  if (ksl_flag)
  {
    den = get_ion (0, 1, 1, 1, name);

    /* Having gotten the densities, modify them */

    printf ("%d\n", zdom[0].ndim2);

    for (i = 0; i < zdom[0].ndim2; i++)
    {
      printf ("%e\n", den[i]);
      den[i] *= 1e-6;
    }

    /* Now paint the densities into the appropriatle
       places in the plasma structure
     */
    put_ion (0, 1, 1, den);
  }


  printf ("outputting to %s\n", outfile);


  wind_save (outfile);


  printf ("gotcha %s\n", files.root);

  exit (0);

}


/**********************************************************/
/**
 * @brief      update a specific ion with new densities
 *
 * @param [in]  int ndom   the domain number
 * @param [in]  int element the element(z) for which 
 * @param [in]  int istate  the ion number
 * @param [in]  double den  an array with the new densities
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
put_ion (ndom, element, istate, den)
     int ndom, element, istate;
     double *den;
{
  int i, n;
  int nion;
  int nstart, ndim2;
  int nplasma;



  for (i = 0; i < zdom[0].ndim2; i++)
  {
//    printf ("%e\n", den[i]);
  }
  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;


  /* Find the ion */

  nion = 0;
  while (nion < nions && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;

  printf ("Found %d\n", nion);


  for (n = 0; n < ndim2; n++)
  {
    nplasma = wmain[nstart + n].nplasma;
    plasmamain[nplasma].density[nion] = den[n];
  }

  return (0);

}


int
apply_model (ndom, filename)
     int ndom;
     char *filename;
{
  int ndim, mdim;
//OLD  int nstart, n, nion, nplasma;
  int n, nion, nplasma;

  printf ("We have been given a model file - we will be using this for new densities in domain 0\n");
  ndim = zdom[ndom].ndim;
  mdim = zdom[ndom].mdim;
  printf ("Current dimensions are %i %i\n", ndim, mdim);
  import_wind2 (ndom, model_file);
  printf ("Model dimensions are %i %i\n", imported_model[ndom].ndim, imported_model[ndom].mdim);
  if (ndim == imported_model[ndom].ndim && mdim == imported_model[ndom].mdim)
  {
    printf ("The model dimensions match the current file - proceeding\n");
    if (zdom[ndom].coord_type == SPHERICAL)
    {
      printf ("We have a spherical model\n");
//OLD      nstart = zdom[ndom].nstart;
      for (n = 0; n < ndim; n++)
      {
        wmain[n].v[0] = imported_model[ndom].v_r[n];    //we need a value for v_r for all cells including ghosts
        if (wmain[n].inwind > -1)
        {
          nplasma = wmain[n].nplasma;
          for (nion = 0; nion < nions; nion++)  //Change the absolute number densities, fractions remain the same
          {
            plasmamain[nplasma].density[nion] =
              plasmamain[nplasma].density[nion] * (imported_model[ndom].mass_rho[n] / plasmamain[nplasma].rho);
          }
          plasmamain[nplasma].rho = imported_model[ndom].mass_rho[n];
          if (imported_model[ndom].init_temperature == FALSE)
          {
            plasmamain[nplasma].t_e = imported_model[ndom].t_e[n];
            plasmamain[nplasma].t_r = imported_model[ndom].t_r[n];
          }
        }
      }
    }
    else if (zdom[ndom].coord_type == RTHETA)
    {
      printf ("We have an r-theta model %i %i\n", ndim, mdim);
      for (n = 0; n < imported_model[ndom].ncell; n++)  //Loop over all the cells in the model
      {
        wmain[n].v[0] = imported_model[ndom].v_x[n];    //we need a value for v_r for all cells including ghosts
        wmain[n].v[1] = imported_model[ndom].v_y[n];    //we need a value for v_r for all cells including ghosts
        wmain[n].v[2] = imported_model[ndom].v_z[n];    //we need a value for v_r for all cells including ghosts
        if (wmain[n].inwind > -1)
        {
          nplasma = wmain[n].nplasma;
          for (nion = 0; nion < nions; nion++)  //Change the absolute number densities, fractions remain the same
          {
            plasmamain[nplasma].density[nion] =
              plasmamain[nplasma].density[nion] * (imported_model[ndom].mass_rho[n] / plasmamain[nplasma].rho);
          }
          plasmamain[nplasma].rho = imported_model[ndom].mass_rho[n];
          if (imported_model[ndom].init_temperature == FALSE)
          {
            plasmamain[nplasma].t_e = imported_model[ndom].t_e[n];
            plasmamain[nplasma].t_r = imported_model[ndom].t_r[n];
          }
        }
      }
    }
  }
  else
  {
    printf ("The model doesnt match the current windsave aborting\n");
    exit (0);
  }
  return (0);
}


int
frame_transform (ndom)
     int ndom;
{
  int n, nion;
  double factor;                //This will either be gamma or 1/gamma
  int nplasma;

  printf ("BOOM got here - frame=%i\n", geo.frame);

  if (cmf2obs_flag)
  {
    printf ("We are going from cmf to observer frame\n");
    if (geo.frame == CMF_FRAME)
    {
      printf ("BOOM IN CMF frame - converting to OBS\n");
      geo.frame = OBS_FRAME;
    }
    else if (geo.frame == OBS_FRAME)
    {
      printf ("BOOM IN OBS frame already - exiting\n");
      exit (0);
    }
  }
  else if (obs2cmf_flag)
  {
    printf ("We are going from observer to cmf frame\n");
    if (geo.frame == CMF_FRAME)
    {
      printf ("BOOM IN CMF frame already - exiting\n");
      exit (0);
    }
    else if (geo.frame == OBS_FRAME)
    {
      printf ("BOOM IN OBS frame - converting to CMF\n");
      geo.frame = CMF_FRAME;
    }

  }
  printf ("domain %i %i\n", zdom[ndom].ndim, zdom[ndom].mdim);

  //First deal with wind quantities
  printf ("looping over %i %i\n", zdom[ndom].nstart, zdom[ndom].nstop);
  for (n = zdom[ndom].nstart; n < zdom[ndom].nstop; n++)
  {
    if (cmf2obs_flag)
      factor = 1. / wmain[n].xgamma_cen;
//OLD    else if (obs2cmf_flag)
    else
      factor = wmain[n].xgamma_cen;
    wmain[n].vol *= factor;
    if (wmain[n].inwind == W_ALL_INWIND || wmain[n].inwind == W_PART_INWIND)
    {
      printf ("This cell is in the wind - transforming plasma variables\n");
      nplasma = wmain[n].nplasma;
      plasmamain[nplasma].vol *= factor;
      plasmamain[nplasma].rho /= factor;
      plasmamain[nplasma].ne /= factor;
      for (nion = 0; nion < nions; nion++)
      {
        plasmamain[nplasma].density[nion] /= factor;
      }

    }
  }



  return (0);
}
