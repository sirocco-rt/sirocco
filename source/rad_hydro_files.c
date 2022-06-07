
/***********************************************************/
/** @file  rad_hydro_files.c
 * @author nsh
 * @date   April, 2020
 *
 * @brief  A standalone routine for writing a set of files required for
 * rad_hydro runs from a windsavefile 
 *
 * This routine is run from the command line, as follows
 *
 * rad_hydro_files  windsave_root
 *
 * where windsave_root is the root name for a python run, or more precisely
 * the rootname of a windsave file, as the .pf file is not read.
 *
 * The routine reads the windsavefile and then writes out a set of files
 * that are required to calculate the force multiplier, the continuum
 * driving force and the heating ad cooling rates.
 *


 * ### Notes ###
 *
 * Whereas py_wind is intended to be run interactively, rad_hydro_files is
 * entirely hardwired so that it produces a standard set of output
 * files.  To change the outputs one has to modify the routine
 *
 * This file just contains the driving routine.  All of the 
 * real work is carried out in rad_hydro_files_sub.c  Indeed so 
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
#include "python.h"


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
 *  -d   Write out ion densisties, rather than ion fractions in the cell
 *  -s   Write out the number of scatters per unit volume  by an ion in a cell, instead of the 
 *       ion fraction
 *
 * The switches only affect the ion tables not the master table
 * This was originally implemented to enable somebody to query which version of
 * Python windsave2table was compiled with. Works in a similar fashion to how
 * the version information is stored and viewed in Python.
 * 
 * 
 *
 **********************************************************/

char rad_hydro_files_help[] = "Usage: rad_hydro_files rootname \n\
-h get this help message and quit\n\
";

void
xparse_arguments (int argc, char *argv[], char root[], int *ion_switch)
{
  char *fget_rc;
  char input[LINELENGTH];

  *ion_switch = 0;


  if (argc == 1)
  {
    printf ("Root for wind file :");
    fget_rc = fgets (input, LINELENGTH, stdin);
    if (!fget_rc)
    {
      printf ("No root file provided, exiting\n\n");
      printf ("%s", rad_hydro_files_help);

      exit (0);
    }
    get_root (root, input);
  }
  else
  {
    strcpy (input, argv[argc - 1]);
    get_root (root, input);
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
  char root[LINELENGTH];
  char outputfile[LINELENGTH];
  char windsavefile[LINELENGTH];
  char parameter_file[LINELENGTH];
  int ion_switch, nwind, nplasma;
  int i, j, ii, domain;
  double vol, kappa_es, lum_sum, cool_sum;
  double t_opt, t_UV, t_Xray, v_th, fhat[3];    /*This is the dimensionless optical depth parameter computed for communication to rad-hydro. */

  struct photon ptest;          //We need a test photon structure in order to compute t

  FILE *fptr, *fptr2, *fptr3, *fptr4, *fptr5, *fptr6, *fptr7, *fptr8, *fptr9, *fptr10, *fopen ();       /*This is the file to communicate with zeus */
  domain = geo.hydro_domain_number;


  strcpy (parameter_file, "NONE");

  /* Next command stops Debug statements printing out in py_wind */
  Log_set_verbosity (3);

  /*
   * EP: added some extra argument parsing for windsave2table - specifically
   * because I was having some trouble with knowing when windsave2table was
   * last compiled and on what commit this was
   */

  xparse_arguments (argc, argv, root, &ion_switch);

  printf ("Reading data from file %s\n", root);
  /* Now create the names of all the files which will be written */
  strcpy (windsavefile, root);
  strcpy (outputfile, root);

  strcat (windsavefile, ".wind_save");
  strcat (outputfile, ".txt");


/* Read in the wind file */

  if (wind_read (windsavefile) < 0)
  {
    Error ("py_wind: Could not open %s", windsavefile);
    exit (0);
  }


  printf ("Read wind_file %s\n", windsavefile);

  get_atomic_data (geo.atomic_filename);

  printf ("Read Atomic data from %s\n", geo.atomic_filename);


  cool_sum = wind_cooling ();   /*We call wind_cooling here to obtain an up to date set of cooling rates */
  lum_sum = wind_luminosity (0.0, VERY_BIG, MODE_CMF_TIME);     /*and we also call wind_luminosity to get the luminosities */

  fptr = fopen ("py_heatcool.dat", "w");
  fptr2 = fopen ("py_driving.dat", "w");
  fptr3 = fopen ("py_ion_data.dat", "w");
  fptr4 = fopen ("py_spec_data.dat", "w");
  fptr5 = fopen ("py_pcon_data.dat", "w");
  fptr6 = fopen ("py_debug_data.dat", "w");
  fptr7 = fopen ("py_fluxes.dat", "w");
  fptr8 = fopen ("directional_flux_x.dat", "w");
  fptr9 = fopen ("directional_flux_y.dat", "w");
  fptr10 = fopen ("directional_flux_z.dat", "w");

//  fptr11 = fopen ("idomega.dat", "w");


  if (zdom[domain].coord_type == SPHERICAL || zdom[domain].coord_type == RTHETA)
    fprintf (fptr, "i j rcen thetacen vol temp xi ne heat_xray heat_comp heat_lines heat_ff cool_comp cool_lines cool_ff rho n_h\n");
  else if (zdom[domain].coord_type == CYLIND)
    fprintf (fptr, "i j rcen zcen vol temp xi ne heat_xray heat_comp heat_lines heat_ff cool_comp cool_lines cool_ff rho n_h\n");



  fprintf (fptr3, "nions %i\n", nions);
  for (i = 0; i < nions; i++)
  {
    fprintf (fptr3, "ion %i %s %i %i\n", i, ele[ion[i].nelem].name, ion[i].z, ion[i].istate);
  }
  fprintf (fptr3, "nplasma %i\n", NPLASMA);


  fprintf (fptr4, "model %i\n", geo.ioniz_mode);
  if (geo.ioniz_mode == IONMODE_MATRIX_SPECTRALMODEL)
  {
    fprintf (fptr4, "nbands %i\n", geo.nxfreq);
    fprintf (fptr4, "nplasma %i\n", NPLASMA);
    for (i = 0; i < geo.nxfreq + 1; i++)
      fprintf (fptr4, "%e ", geo.xfreq[i]);     //hard wired band edges
    fprintf (fptr4, "\n ");
  }
  else if (geo.ioniz_mode == IONMODE_MATRIX_BB)
    fprintf (fptr4, "nplasma %i\n", NPLASMA);

  fprintf (fptr5, "nplasma %i\n", NPLASMA);

  if (zdom[domain].coord_type == SPHERICAL)
    fprintf (fptr6, "i j rcen thetacen v_th dvdr J\n");
  else if (zdom[domain].coord_type == CYLIND)
    fprintf (fptr6, "i j rcen zcen v_th dvz_dz J\n");
  else if (zdom[domain].coord_type == RTHETA)
    fprintf (fptr6, "i j rcen thetacen v_th dvr_dr J\n");

  printf ("Set up files\n");

  if (geo.hydro_domain_number > 0)
  {
    domain = geo.hydro_domain_number;
  }
  else
  {
    domain = 0;
  }

  if (zdom[domain].coord_type == SPHERICAL)
  {
    fprintf (fptr2, "i j rcen thetacen F_vis_r F_UV_r F_Xray_r es_f_r bf_f_r\n");       //directional flux by band
    fprintf (fptr7, "i j rcen thetacen F_vis_r F_UV_r F_Xray_r \n");    //directional flux by band
  }
  else if (zdom[domain].coord_type == CYLIND)
  {
    fprintf (fptr2, "i j rcen zcen vol rho ne F_vis_x F_vis_y F_vis_z F_vis_mod F_UV_x F_UV_y F_UV_z F_UV_mod F_Xray_x F_Xray_y F_Xray_z F_Xray_mod F_vis_x2 F_vis_y2 F_vis_z2 F_vis_mod2 F_UV_x2 F_UV_y2 F_UV_z2 F_UV_mod2 F_Xray_x2 F_Xray_y2 F_Xray_z2 F_Xray_mod2 es_f_x es_f_y es_f_z es_f_mod bf_f_x bf_f_y bf_f_z bf_f_mod\n");  //directional flux by band
    fprintf (fptr7, "i j rcen zcen F_vis_x F_vis_y F_vis_z F_vis_mod F_UV_x F_UV_y F_UV_z F_UV_mod F_Xray_x F_Xray_y F_Xray_z F_Xray_mod\n");   //directional flux by band
  }
  else if (zdom[domain].coord_type == RTHETA)
  {
    fprintf (fptr2, "i j rcen thetacen vol rho ne F_vis_x F_vis_y F_vis_z F_vis_mod F_UV_x F_UV_y F_UV_z F_UV_mod F_Xray_x F_Xray_y F_Xray_z F_Xray_mod es_f_x es_f_y es_f_z es_f_mod bf_f_x bf_f_y bf_f_z bf_f_mod\n");        //directional flux by band
    fprintf (fptr7, "i j rcen thetacen F_vis_x F_vis_y F_vis_z F_vis_mod F_UV_x F_UV_y F_UV_z F_UV_mod F_Xray_x F_Xray_y F_Xray_z F_Xray_mod\n");       //directional flux by band
  }

  fprintf (fptr8, "NANGLES %i\n", NFLUX_ANGLES);
  fprintf (fptr9, "NANGLES %i\n", NFLUX_ANGLES);
  fprintf (fptr10, "NANGLES %i\n", NFLUX_ANGLES);


  for (nwind = zdom[domain].nstart; nwind < zdom[domain].nstop; nwind++)
  {
    if (wmain[nwind].vol > 0.0)
    {
      nplasma = wmain[nwind].nplasma;
      wind_n_to_ij (domain, plasmamain[nplasma].nwind, &i, &j);

      if (zdom[domain].coord_type == SPHERICAL)
        i = i - 1;              //There is an extra radial 'ghost zone' in spherical coords in python, we need to make our i,j agree with zeus
      vol = wmain[plasmamain[nplasma].nwind].vol;
      if (zdom[domain].coord_type == SPHERICAL || zdom[domain].coord_type == RTHETA)
        fprintf (fptr, "%d %d %e %e %e ", i, j, wmain[plasmamain[nplasma].nwind].rcen, wmain[plasmamain[nplasma].nwind].thetacen / RADIAN, vol);        //output geometric things
      else if (zdom[domain].coord_type == CYLIND)
        fprintf (fptr, "%d %d %e %e %e ", i, j, wmain[plasmamain[nplasma].nwind].xcen[0], wmain[plasmamain[nplasma].nwind].xcen[2], vol);       //output geometric things

      fprintf (fptr, "%e %e %e ", plasmamain[nplasma].t_e, plasmamain[nplasma].xi, plasmamain[nplasma].ne);     //output temp, xi and ne to ease plotting of heating rates
      fprintf (fptr, "%e ", (plasmamain[nplasma].heat_photo + plasmamain[nplasma].heat_auger) / vol);   //Xray heating - or photoionization
      fprintf (fptr, "%e ", (plasmamain[nplasma].heat_comp) / vol);     //Compton heating
      fprintf (fptr, "%e ", (plasmamain[nplasma].heat_lines) / vol);    //Line heating 28/10/15 - not currently used in zeus
      fprintf (fptr, "%e ", (plasmamain[nplasma].heat_ff) / vol);       //FF heating 28/10/15 - not currently used in zeus
      fprintf (fptr, "%e ", (plasmamain[nplasma].cool_comp) / vol);     //Compton cooling
      fprintf (fptr, "%e ", (plasmamain[nplasma].lum_lines + plasmamain[nplasma].cool_rr + plasmamain[nplasma].cool_dr) / vol); //Line cooling must include all recombination cooling
      fprintf (fptr, "%e ", (plasmamain[nplasma].lum_ff) / vol);        //ff cooling
      fprintf (fptr, "%e ", plasmamain[nplasma].rho);   //density
      fprintf (fptr, "%e \n", plasmamain[nplasma].rho * rho2nh);        //hydrogen number density

      fprintf (fptr8, "%d %d %e %e ", i, j, wmain[plasmamain[nplasma].nwind].rcen, wmain[plasmamain[nplasma].nwind].thetacen / RADIAN); //output geometric things
      fprintf (fptr9, "%d %d %e %e ", i, j, wmain[plasmamain[nplasma].nwind].rcen, wmain[plasmamain[nplasma].nwind].thetacen / RADIAN); //output geometric things
      fprintf (fptr10, "%d %d %e %e ", i, j, wmain[plasmamain[nplasma].nwind].rcen, wmain[plasmamain[nplasma].nwind].thetacen / RADIAN);        //output geometric things
//      fprintf (fptr11, "%d %d %e %e %e ", i, j, wmain[plasmamain[nplasma].nwind].rcen, wmain[plasmamain[nplasma].nwind].thetacen / RADIAN, wmain[plasmamain[nplasma].nwind].theta / RADIAN);    //output geometric things
      for (ii = 0; ii < NFLUX_ANGLES; ii++)
      {
        fprintf (fptr8, "%e ", plasmamain[nplasma].F_UV_ang_x_persist[ii]);
        fprintf (fptr9, "%e ", plasmamain[nplasma].F_UV_ang_y_persist[ii]);
        fprintf (fptr10, "%e ", plasmamain[nplasma].F_UV_ang_z_persist[ii]);

//        fprintf (fptr11, "%e ", plasmamain[nplasma].idomega[ii]);
      }
      fprintf (fptr8, "\n");
      fprintf (fptr9, "\n");
      fprintf (fptr10, "\n");
//      fprintf (fptr11, "\n");



      if (zdom[domain].coord_type == SPHERICAL || zdom[domain].coord_type == RTHETA)
      {
        fprintf (fptr2, "%d %d %e %e %e ", i, j, wmain[plasmamain[nplasma].nwind].rcen, wmain[plasmamain[nplasma].nwind].thetacen / RADIAN, vol);       //output geometric things
        fprintf (fptr5, "%d %d %e %e ", i, j, wmain[plasmamain[nplasma].nwind].rcen, wmain[plasmamain[nplasma].nwind].thetacen / RADIAN);       //output geometric things
        fprintf (fptr2, "%e ", plasmamain[nplasma].rho);        //density
        fprintf (fptr2, "%e ", plasmamain[nplasma].ne);
        fprintf (fptr7, "%d %d %e %e ", i, j, wmain[plasmamain[nplasma].nwind].rcen, wmain[plasmamain[nplasma].nwind].thetacen / RADIAN);       //output geometric things
      }
      else if (zdom[domain].coord_type == CYLIND)
      {
        fprintf (fptr2, "%d %d %e %e %e ", i, j, wmain[plasmamain[nplasma].nwind].xcen[0], wmain[plasmamain[nplasma].nwind].xcen[2], vol);      //output geometric things
        fprintf (fptr5, "%d %d %e %e ", i, j, wmain[plasmamain[nplasma].nwind].xcen[0], wmain[plasmamain[nplasma].nwind].xcen[2]);      //output geometric things
        fprintf (fptr2, "%e ", plasmamain[nplasma].rho);        //density
        fprintf (fptr2, "%e ", plasmamain[nplasma].ne);
        fprintf (fptr7, "%d %d %e %e ", i, j, wmain[plasmamain[nplasma].nwind].xcen[0], wmain[plasmamain[nplasma].nwind].xcen[2]);      //output geometric things
      }
      if (zdom[domain].coord_type == SPHERICAL)
      {
        fprintf (fptr2, "%e ", plasmamain[nplasma].F_vis[0]);   //directional flux by band
        fprintf (fptr2, "%e ", plasmamain[nplasma].F_UV[0]);    //directional flux by band
        fprintf (fptr2, "%e ", plasmamain[nplasma].F_Xray[0]);  //directional flux by band
        fprintf (fptr2, "%e ", plasmamain[nplasma].rad_force_es[0]);    //electron scattering radiation force in the w(x) direction
        fprintf (fptr2, "%e\n", plasmamain[nplasma].rad_force_bf_persist[0]);   //bound free scattering radiation force in the w(x) direction    
        fprintf (fptr7, "%e ", plasmamain[nplasma].F_vis[0]);   //directional flux by band
        fprintf (fptr7, "%e ", plasmamain[nplasma].F_UV[0]);    //directional flux by band
        fprintf (fptr7, "%e ", plasmamain[nplasma].F_Xray[0]);  //directional flux by band      
      }
      else
      {
        {
          fprintf (fptr2, "%e %e %e %e ", plasmamain[nplasma].F_vis[0], plasmamain[nplasma].F_vis[1], plasmamain[nplasma].F_vis[2], plasmamain[nplasma].F_vis[3]);      //directional flux by band
          fprintf (fptr2, "%e %e %e %e ", plasmamain[nplasma].F_UV[0], plasmamain[nplasma].F_UV[1], plasmamain[nplasma].F_UV[2], plasmamain[nplasma].F_UV[3]);  //directional flux by band
          fprintf (fptr2, "%e %e %e %e ", plasmamain[nplasma].F_Xray[0], plasmamain[nplasma].F_Xray[1], plasmamain[nplasma].F_Xray[2], plasmamain[nplasma].F_Xray[3]);  //directional flux by band
          fprintf (fptr2, "%e ", plasmamain[nplasma].rad_force_es[0]);  //electron scattering radiation force in the w(x) direction
          fprintf (fptr2, "%e ", plasmamain[nplasma].rad_force_es[1]);  //electron scattering radiation force in the phi(rotational) directionz direction
          fprintf (fptr2, "%e ", plasmamain[nplasma].rad_force_es[2]);  //electron scattering radiation force in the z direction
          fprintf (fptr2, "%e ", plasmamain[nplasma].rad_force_es[3]);  //sum of magnitude of electron scattering radiation force
          fprintf (fptr2, "%e ", plasmamain[nplasma].rad_force_bf_persist[0]);  //bound free scattering radiation force in the w(x) direction
          fprintf (fptr2, "%e ", plasmamain[nplasma].rad_force_bf_persist[1]);  //bound free scattering radiation force in the phi(rotational) direction
          fprintf (fptr2, "%e ", plasmamain[nplasma].rad_force_bf_persist[2]);  //bound free scattering radiation force in the z direction
          fprintf (fptr2, "%e \n", plasmamain[nplasma].rad_force_bf_persist[3]);        //sum of magnitude of bound free scattering radiation force 
        }
        fprintf (fptr7, "%e %e %e %e ", plasmamain[nplasma].F_vis_persistent[0], plasmamain[nplasma].F_vis_persistent[1], plasmamain[nplasma].F_vis_persistent[2], plasmamain[nplasma].F_vis_persistent[3]);    //directional flux by band
        fprintf (fptr7, "%e %e %e %e ", plasmamain[nplasma].F_UV_persistent[0], plasmamain[nplasma].F_UV_persistent[1], plasmamain[nplasma].F_UV_persistent[2], plasmamain[nplasma].F_UV_persistent[3]);        //directional flux by band
        fprintf (fptr7, "%e %e %e %e\n ", plasmamain[nplasma].F_Xray_persistent[0], plasmamain[nplasma].F_Xray_persistent[1], plasmamain[nplasma].F_Xray_persistent[2], plasmamain[nplasma].F_Xray_persistent[3]);      //directional flux by band
      }
      fprintf (fptr3, "%d %d ", i, j);  //output geometric things               
      for (ii = 0; ii < nions; ii++)
        fprintf (fptr3, "%e ", plasmamain[nplasma].density[ii]);
      fprintf (fptr3, "\n");
      fprintf (fptr4, "%d %d ", i, j);  //output geometric things 
      if (geo.ioniz_mode == IONMODE_MATRIX_SPECTRALMODEL)
      {
        for (ii = 0; ii < geo.nxfreq; ii++)
          fprintf (fptr4, "%e %e %i %e %e %e %e ",
                   plasmamain[nplasma].fmin_mod[ii], plasmamain[nplasma].fmax_mod[ii], plasmamain[nplasma].spec_mod_type[ii],
                   plasmamain[nplasma].pl_log_w[ii], plasmamain[nplasma].pl_alpha[ii], plasmamain[nplasma].exp_w[ii],
                   plasmamain[nplasma].exp_temp[ii]);
      }
      else if (geo.ioniz_mode == IONMODE_MATRIX_BB)
        fprintf (fptr4, "%e %e ", plasmamain[nplasma].t_r, plasmamain[nplasma].w);
      fprintf (fptr4, "\n ");


      //We need to compute the g factor for this cell and output it.


      v_th = pow ((2. * BOLTZMANN * plasmamain[nplasma].t_e / MPROT), 0.5);     //We need the thermal velocity for hydrogen
//      v_th = 4.2e5;
      stuff_v (wmain[plasmamain[nplasma].nwind].xcen, ptest.x); //place our test photon at the centre of the cell
      ptest.grid = nwind;       //We need our test photon to know where it is 
      kappa_es = THOMPSON * plasmamain[nplasma].ne / plasmamain[nplasma].rho;
      kappa_es = THOMPSON / MPROT;

      //First for the optcial band (up to 4000AA)     
      if (length (plasmamain[nplasma].F_vis) > 0.0)     //Only makes sense if flux in this band is non-zero
      {
        if (zdom[domain].coord_type == SPHERICAL)       //We have to do something special here - because flux is r, theta, phi in sphericals
        {
          fhat[0] = sqrt (length (plasmamain[nplasma].F_vis_persistent));
          fhat[1] = 0.0;
          fhat[2] = sqrt (length (plasmamain[nplasma].F_vis_persistent));
        }
        else
        {
          stuff_v (plasmamain[nplasma].F_vis_persistent, fhat);
        }
        if (renorm (fhat, 1.) == -1)    //A unit vector in the direction of the flux - this can be treated as the lmn vector of a pretend photon
        {
          t_opt = 0.0;
        }
        else
        {
          stuff_v (fhat, ptest.lmn);    //place our test photon at the centre of the cell            
          t_opt = kappa_es * plasmamain[nplasma].rho * v_th / fabs (dvwind_ds_cmf (&ptest));
        }
      }
      else
        t_opt = 0.0;            //Essentually a flag that there is no way of computing t (and hence M) in this cell.

      //Now for the UV band (up to 4000AA->100AA)                                             
      if (length (plasmamain[nplasma].F_UV) > 0.0)      //Only makes sense if flux in this band is non-zero
      {
        if (zdom[domain].coord_type == SPHERICAL)       //We have to do something special here - because flux is r, theta, phi in sphericals
        {
          fhat[0] = sqrt (length (plasmamain[nplasma].F_UV_persistent));
          fhat[1] = 0.0;
          fhat[2] = sqrt (length (plasmamain[nplasma].F_UV_persistent));
        }
        else
        {
          stuff_v (plasmamain[nplasma].F_UV_persistent, fhat);
        }
        if (renorm (fhat, 1.) == -1)    //A unit vector in the direction of the flux - this can be treated as the lmn vector of a pretend photon
        {
          t_UV = 0.0;
        }
        else
        {
          stuff_v (fhat, ptest.lmn);    //place our test photon at the centre of the cell            
          t_UV = kappa_es * plasmamain[nplasma].rho * v_th / fabs (dvwind_ds_cmf (&ptest));
        }
      }
      else
        t_UV = 0.0;             //Essentually a flag that there is no way of computing t (and hence M) in this cell.


      //And finally for the Xray band (up to 100AA and up)
      if (length (plasmamain[nplasma].F_Xray) > 0.0)    //Only makes sense if flux in this band is non-zero
      {
        if (zdom[domain].coord_type == SPHERICAL)       //We have to do something special here - because flux is r, theta, phi in sphericals
        {
          fhat[0] = sqrt (length (plasmamain[nplasma].F_Xray_persistent));
          fhat[1] = 0.0;
          fhat[2] = sqrt (length (plasmamain[nplasma].F_Xray_persistent));
        }
        else
        {
          stuff_v (plasmamain[nplasma].F_Xray_persistent, fhat);
        }
        if (renorm (fhat, 1.) == -1)    //A unit vector in the direction of the flux - this can be treated as the lmn vector of a pretend photon
        {
          t_Xray = 0.0;
        }
        else
        {
          stuff_v (fhat, ptest.lmn);    //place our test photon at the centre of the cell            
          t_Xray = kappa_es * plasmamain[nplasma].rho * v_th / fabs (dvwind_ds_cmf (&ptest));
        }
      }
      else
        t_Xray = 0.0;           //Essentually a flag that there is no way of computing t (and hence M) in this cell.                

      fprintf (fptr5, " %e %e %e %e %e %e %e\n", plasmamain[nplasma].t_e, plasmamain[nplasma].rho,
               plasmamain[nplasma].rho * rho2nh, plasmamain[nplasma].ne, t_opt, t_UV, t_Xray);

      fprintf (fptr6, "%d %d %e %e %e %e %e\n", i, j, wmain[plasmamain[nplasma].nwind].rcen, wmain[plasmamain[nplasma].nwind].thetacen / RADIAN, v_th, fabs (dvwind_ds_cmf (&ptest)), plasmamain[nplasma].j);   //output geometric things
    }
  }
  fclose (fptr);
  fclose (fptr2);
  fclose (fptr3);
  fclose (fptr4);
  fclose (fptr5);
  fclose (fptr6);
  fclose (fptr7);
  fclose (fptr8);
  fclose (fptr9);
  fclose (fptr10);

  exit (0);
}
