
/***********************************************************/
/** @file  windsave2table_sub.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  Subroutines for windsave2table
 *
 * These are various routines used by windsave2table which reads
 * a windsave file and writes various variables of that windsave
 * file to ascii files which can be read with astropy.io ascii
 * as tables.
 *
 * Unlike swind, windsave2table is hardwired and to change what
 * is written one must actually modify the routines themselves.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

int xedge = FALSE;


/**********************************************************/
/**
 * @brief      The main routine associated with windsave2table, this routine calls
 * various other routines which write individual files
 *
 * @param [in] char *  root   The rootname of the windsave file
 * @param [in] int     ion_switch Choose what type of ion data to print
 * @param [in] int    edge_switch If TRUE, include edge cells          
 * @return     Always returns 0
 *
 * @details
 *
 *
 *
 * ### Notes ###
 *
 * The routine cycles through the various domains and calles subroutines
 * that write the individual files for each domain.
 *
 *
 **********************************************************/

int
do_windsave2table (root, ion_switch, edge_switch)
     char *root;
     int ion_switch;
     int edge_switch;
{
  int ndom, i;
  char rootname[LINELENGTH];
  int all[7] = { 0, 4, 5, 6, 7, 8, 9 };

  xedge = edge_switch;          //if TRUE include edge cells


  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {

    if (geo.ndomain > 1)
    {
      sprintf (rootname, "%s.%d", root, ndom);
    }
    else
    {
      sprintf (rootname, "%s", root);
    }

    create_master_table (ndom, rootname);
    create_heat_table (ndom, rootname);
    create_convergence_table (ndom, rootname);
    create_velocity_gradient_table (ndom, rootname);
    create_spec_table (ndom, rootname);

    if (ion_switch != 99)
    {
      create_ion_table (ndom, rootname, 1, ion_switch);
      create_ion_table (ndom, rootname, 2, ion_switch);
      create_ion_table (ndom, rootname, 6, ion_switch);
      create_ion_table (ndom, rootname, 7, ion_switch);
      create_ion_table (ndom, rootname, 8, ion_switch);
      create_ion_table (ndom, rootname, 11, ion_switch);
      create_ion_table (ndom, rootname, 14, ion_switch);
      create_ion_table (ndom, rootname, 20, ion_switch);
      create_ion_table (ndom, rootname, 26, ion_switch);
    }
    else
    {
      for (i = 0; i < 7; i++)
      {
        create_ion_table (ndom, rootname, 1, all[i]);
        create_ion_table (ndom, rootname, 2, all[i]);
        create_ion_table (ndom, rootname, 6, all[i]);
      }

    }


  }
  return (0);
}






/**********************************************************/
/**
 * @brief      writes specific  variables of a windsaave
 * that are intended to be of general interest to
 * file which has the format of an astropy table
 *
 * 	It is intended to be easily modifible.
 *
 * @param [in] int  ndom   A domain number
 * @param [in] char  rootname   The rootname of the master file
 * @return   Always returns 0
 *
 * @details
 *
 * The master table is contains basic information for each
 * cell in the wind, such as the electron density, the density,
 * the ionization parameter, and the radiative and electron temperature
 *
 * The routine takes data directly from wmain, and then calls
 * get_one or get_ion multiple times to fet info from the Plasma
 * structure.
 *
 * It then writes the data to an ascii file which can be read as
 * an  astropy table
 *
 * ### Notes ###
 * To add a variable one just needs to define the column_name
 * and send the appropriate call to either get_one or get_ion.
 *
 * There is some duplicated code in the routine that pertains
 * to whether one is dealing with a spherecial or a 2d coordinate
 * system.  It should be possible to delete this
 *
 **********************************************************/

int
create_master_table (ndom, rootname)
     int ndom;
     char rootname[];
{
  char filename[132];
  double *c[50], *converge;
  char column_name[50][20];
  char one_line[1024], start[1024], one_value[20];
  char name[132];               /* file name extension */


  int i, ii, jj;
  int nstart, ndim2;
  int n, ncols;
  FILE *fptr;

  strcpy (filename, rootname);
  strcat (filename, ".master.txt");


  fptr = fopen (filename, "w");

  /* Get the variables that one needs */

  c[0] = get_one (ndom, "vol");
  strcpy (column_name[0], "vol");

  c[1] = get_one (ndom, "rho");
  strcpy (column_name[1], "rho");

  c[2] = get_one (ndom, "ne");
  strcpy (column_name[2], "ne");

  c[3] = get_one (ndom, "t_e");
  strcpy (column_name[3], "t_e");

  c[4] = get_one (ndom, "t_r");
  strcpy (column_name[4], "t_r");

  c[5] = get_ion (ndom, 1, 1, 0, name);
  strcpy (column_name[5], "h1");

  c[6] = get_ion (ndom, 2, 2, 0, name);
  strcpy (column_name[6], "he2");

  c[7] = get_ion (ndom, 6, 4, 0, name);
  strcpy (column_name[7], "c4");

  c[8] = get_ion (ndom, 7, 5, 0, name);
  strcpy (column_name[8], "n5");

  c[9] = get_ion (ndom, 8, 6, 0, name);
  strcpy (column_name[9], "o6");

  c[10] = get_one (ndom, "dmo_dt_x");
  strcpy (column_name[10], "dmo_dt_x");


  c[11] = get_one (ndom, "dmo_dt_y");
  strcpy (column_name[11], "dmo_dt_y");

  c[12] = get_one (ndom, "dmo_dt_z");
  strcpy (column_name[12], "dmo_dt_z");

  c[13] = get_one (ndom, "ip");
  strcpy (column_name[13], "ip");

  c[14] = get_one (ndom, "xi");
  strcpy (column_name[14], "xi");

  c[15] = get_one (ndom, "ntot");
  strcpy (column_name[15], "ntot");

  c[16] = get_one (ndom, "nrad");
  strcpy (column_name[16], "nrad");

  c[17] = get_one (ndom, "nioniz");
  strcpy (column_name[17], "nioniz");


  /* This should be the maxium number above +1 */
  ncols = 18;


  converge = get_one (ndom, "converge");

  /* At this point oll of the data has been collected */


  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;


  if (zdom[ndom].coord_type == SPHERICAL)
  {


    /*
     * First assemble the header line
     */

    sprintf (start, "%9s %9s %4s %6s %6s %9s %9s %9s ", "r", "rcen", "i", "inwind", "converge", "v_x", "v_y", "v_z");
    strcpy (one_line, start);
    n = 0;
    while (n < ncols)
    {
      sprintf (one_value, "%9.9s ", column_name[n]);
      strcat (one_line, one_value);

      n++;
    }
    fprintf (fptr, "%s\n", one_line);


    /* Now assemble the lines of the table */

    for (i = 0; i < ndim2; i++)
    {
      //This line is different from the two d case
      sprintf (start, "%9.3e %9.3e %4d %6d %8.0f %9.2e %9.2e %9.2e ",
               wmain[nstart + i].r, wmain[nstart + i].rcen, i, wmain[nstart + i].inwind,
               converge[i], wmain[nstart + i].v[0], wmain[nstart + i].v[1], wmain[nstart + i].v[2]);
      strcpy (one_line, start);
      n = 0;
      while (n < ncols)
      {
        sprintf (one_value, "%9.2e ", c[n][i]);
        strcat (one_line, one_value);
        n++;
      }
      fprintf (fptr, "%s\n", one_line);
    }
  }
  else if (zdom[ndom].coord_type == CYLIND)
  {

    /* First assemble the header line */

    sprintf (start, "%8s %8s %8s %8s %4s %4s %6s %8s %9s %9s %9s ", "x", "z", "xcen", "zcen", "i", "j", "inwind", "converge", "v_x", "v_y",
             "v_z");
    strcpy (one_line, start);
    n = 0;
    while (n < ncols)
    {
      sprintf (one_value, "%9.9s ", column_name[n]);
      strcat (one_line, one_value);

      n++;
    }
    fprintf (fptr, "%s\n", one_line);


    /* Now assemble the lines of the table */

    for (i = 0; i < ndim2; i++)
    {
      wind_n_to_ij (ndom, nstart + i, &ii, &jj);
      sprintf (start,
               "%8.2e %8.2e %8.2e %8.2e %4d %4d %6d %8.0f %9.2e %9.2e %9.2e ",
               wmain[nstart + i].x[0], wmain[nstart + i].x[2], wmain[nstart + i].xcen[0], wmain[nstart + i].xcen[2], ii,
               jj, wmain[nstart + i].inwind, converge[i], wmain[nstart + i].v[0], wmain[nstart + i].v[1], wmain[nstart + i].v[2]);
      strcpy (one_line, start);
      n = 0;
      while (n < ncols)
      {
        sprintf (one_value, "%9.2e ", c[n][i]);
        strcat (one_line, one_value);
        n++;
      }
      fprintf (fptr, "%s\n", one_line);
    }
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {

    /* First assemble the header line */

    sprintf (start, "%8s %8s %8s %9s %8s %8s %8s %8s %4s %4s %6s %8s %9s %9s %9s ", "r", "theta", "r_cen", "theta_cen", "x", "z", "xcen",
             "zcen", "i", "j", "inwind", "converge", "v_x", "v_y", "v_z");
    strcpy (one_line, start);
    n = 0;
    while (n < ncols)
    {
      sprintf (one_value, "%9.9s ", column_name[n]);
      strcat (one_line, one_value);

      n++;
    }
    fprintf (fptr, "%s\n", one_line);


    /* Now assemble the lines of the table */

    for (i = 0; i < ndim2; i++)
    {
      wind_n_to_ij (ndom, nstart + i, &ii, &jj);
      sprintf (start,
               "%8.2e %8.2e %8.2e %9.2e %8.2e %8.2e %8.2e %8.2e %4d %4d %6d %8.0f %9.2e %9.2e %9.2e ",
               wmain[nstart + i].r, wmain[nstart + i].theta, wmain[nstart + i].rcen, wmain[nstart + i].thetacen,
               wmain[nstart + i].x[0], wmain[nstart + i].x[2], wmain[nstart + i].xcen[0], wmain[nstart + i].xcen[2], ii,
               jj, wmain[nstart + i].inwind, converge[i], wmain[nstart + i].v[0], wmain[nstart + i].v[1], wmain[nstart + i].v[2]);
      strcpy (one_line, start);
      n = 0;
      while (n < ncols)
      {
        sprintf (one_value, "%9.2e ", c[n][i]);
        strcat (one_line, one_value);
        n++;
      }
      fprintf (fptr, "%s\n", one_line);
    }
  }
  else
  {
    printf ("Error: Cannot print out files for coordinate system type %d\n", zdom[ndom].coord_type);
  }


  fclose (fptr);

  return (0);
}




/**********************************************************/
/**
 * @brief      writes selected variables related to heating and cooling
 * processes to an ascii file which can be read as an astropy table
 *
 *
 * @param [in] int  ndom   The domain of interest
 * @param [in, out] char  rootname[]   The rootname of the windsave file
 * @return     Always returns 0
 *
 * @details
 *
 * ### Notes ###
 *
 * To add a variable one just needs to define the column_name
 * and send the appropriate call to either get_one or get_ion.
 *
 *
 **********************************************************/

int
create_heat_table (ndom, rootname)
     int ndom;
     char rootname[];
{
  char filename[132];
  double *c[50], *converge;
  char column_name[50][22];
  char one_line[1024], start[1024], one_value[20];


  int i, ii, jj;
  int nstart, ndim2;
  int n, ncols;
  FILE *fptr;

  strcpy (filename, rootname);
  strcat (filename, ".heat.txt");


  fptr = fopen (filename, "w");

  /* Get the variables that one needs */

  c[0] = get_one (ndom, "vol");
  strcpy (column_name[0], "vol");

  c[1] = get_one (ndom, "rho");
  strcpy (column_name[1], "rho");

  c[2] = get_one (ndom, "ne");
  strcpy (column_name[2], "ne");

  c[3] = get_one (ndom, "t_e");
  strcpy (column_name[3], "t_e");

  c[4] = get_one (ndom, "t_r");
  strcpy (column_name[4], "t_r");

  c[5] = get_one (ndom, "w");
  strcpy (column_name[5], "w");

  c[6] = get_one (ndom, "ave_freq");
  strcpy (column_name[6], "ave_freq");

  c[7] = get_one (ndom, "J");
  strcpy (column_name[7], "J");

  c[8] = get_one (ndom, "J_direct");
  strcpy (column_name[8], "J_direct");

  c[9] = get_one (ndom, "J_scatt");
  strcpy (column_name[9], "J_scatt");

  c[10] = get_one (ndom, "lum_tot");
  strcpy (column_name[10], "lum_tot");

  c[11] = get_one (ndom, "heat_tot");
  strcpy (column_name[11], "heat_tot");

  c[12] = get_one (ndom, "heat_comp");
  strcpy (column_name[12], "heat_comp");

  c[13] = get_one (ndom, "heat_lines");
  strcpy (column_name[13], "heat_lines");

  c[14] = get_one (ndom, "heat_ff");
  strcpy (column_name[14], "heat_ff");

  c[15] = get_one (ndom, "heat_photo");
  strcpy (column_name[15], "heat_photo");

  c[16] = get_one (ndom, "heat_auger");
  strcpy (column_name[16], "heat_auger");

  c[17] = get_one (ndom, "cool_tot");
  strcpy (column_name[17], "cool_tot");

  c[18] = get_one (ndom, "cool_comp");
  strcpy (column_name[18], "cool_comp");

  c[19] = get_one (ndom, "lum_lines");
  strcpy (column_name[19], "lum_lines");

  c[20] = get_one (ndom, "cool_dr");
  strcpy (column_name[20], "cool_dr");

  c[21] = get_one (ndom, "lum_ff");
  strcpy (column_name[21], "lum_ff");

  c[22] = get_one (ndom, "lum_rr");
  strcpy (column_name[22], "lum_rr");


  c[23] = get_one (ndom, "cool_rr");
  strcpy (column_name[23], "cool_rr");

  c[24] = get_one (ndom, "cool_adiab");
  strcpy (column_name[24], "cool_adiab");

  c[25] = get_one (ndom, "heat_shock");
  strcpy (column_name[25], "heat_shock");

  c[26] = get_one (ndom, "heat_lines_macro");
  strcpy (column_name[26], "ht_ln_macro");

  c[27] = get_one (ndom, "heat_photo_macro");
  strcpy (column_name[27], "ht_ph_macro");

  /* This should be the maximum number above +1 */
  ncols = 28;


  converge = get_one (ndom, "converge");

  /* At this point all of the data has been collected */


  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;


  if (zdom[ndom].coord_type == SPHERICAL)
  {


    /*
     * First assemble the header line
     */

    sprintf (start, "%9s %4s %6s %6s %9s %9s %9s ", "r", "i", "inwind", "converge", "v_x", "v_y", "v_z");
    strcpy (one_line, start);
    n = 0;
    while (n < ncols)
    {
      sprintf (one_value, "%9.9s ", column_name[n]);
      strcat (one_line, one_value);

      n++;
    }
    fprintf (fptr, "%s\n", one_line);


    /* Now assemble the lines of the table */

    for (i = 0; i < ndim2; i++)
    {
      //This line is different from the two d case
      sprintf (start, "%9.3e %4d %6d %8.0f %9.2e %9.2e %9.2e ",
               wmain[nstart + i].r, i, wmain[nstart + i].inwind,
               converge[i], wmain[nstart + i].v[0], wmain[nstart + i].v[1], wmain[nstart + i].v[2]);
      strcpy (one_line, start);
      n = 0;
      while (n < ncols)
      {
        sprintf (one_value, "%9.2e ", c[n][i]);
        strcat (one_line, one_value);
        n++;
      }
      fprintf (fptr, "%s\n", one_line);
    }
  }
  else
  {

    /* First assemble the header line */

    sprintf (start, "%8s %8s %4s %4s %6s %8s %9s %9s %9s ", "x", "z", "i", "j", "inwind", "converge", "v_x", "v_y", "v_z");
    strcpy (one_line, start);
    n = 0;
    while (n < ncols)
    {
      sprintf (one_value, "%9.9s ", column_name[n]);
      strcat (one_line, one_value);

      n++;
    }
    fprintf (fptr, "%s\n", one_line);


    /* Now assemble the lines of the table */

    for (i = 0; i < ndim2; i++)
    {
      wind_n_to_ij (ndom, nstart + i, &ii, &jj);
      sprintf (start,
               "%8.2e %8.2e %4d %4d %6d %8.0f %9.2e %9.2e %9.2e ",
               wmain[nstart + i].xcen[0], wmain[nstart + i].xcen[2], ii,
               jj, wmain[nstart + i].inwind, converge[i], wmain[nstart + i].v[0], wmain[nstart + i].v[1], wmain[nstart + i].v[2]);
      strcpy (one_line, start);
      n = 0;
      while (n < ncols)
      {
        sprintf (one_value, "%9.2e ", c[n][i]);
        strcat (one_line, one_value);
        n++;
      }
      fprintf (fptr, "%s\n", one_line);
    }
  }

  fclose (fptr);

  return (0);
}




/**********************************************************/
/**
 * @brief      writes selected variables related to issues about
 * convergence to an ascii file which can be read as an astropy table
 *
 *
 * @param [in] int  ndom   The domain of interest
 * @param [in, out] char  rootname[]   The rootname of the windsave file
 * @return     Always returns 0
 *
 * @details
 *
 * ### Notes ###
 *
 * To add a variable one just needs to define the column_name
 * and send the appropriate call to either get_one or get_ion.
 *
 *
 **********************************************************/

int
create_convergence_table (ndom, rootname)
     int ndom;
     char rootname[];
{
  char filename[132];
  double *c[50], *converge;
  char column_name[50][20];
  char one_line[1024], start[1024], one_value[20];


  int i, ii, jj;
  int nstart, ndim2;
  int n, ncols;
  FILE *fptr;

  strcpy (filename, rootname);
  strcat (filename, ".converge.txt");


  fptr = fopen (filename, "w");

  /* Get the variables that one needs */

  c[0] = get_one (ndom, "vol");
  strcpy (column_name[0], "vol");

  c[1] = get_one (ndom, "rho");
  strcpy (column_name[1], "rho");

  c[2] = get_one (ndom, "ne");
  strcpy (column_name[2], "ne");

  c[3] = get_one (ndom, "t_e");
  strcpy (column_name[3], "t_e");

  c[4] = get_one (ndom, "t_e_old");
  strcpy (column_name[4], "t_e_old");

  c[5] = get_one (ndom, "dt_e");
  strcpy (column_name[5], "dt_e");

  c[6] = get_one (ndom, "dt_e_old");
  strcpy (column_name[6], "dt_e_old");

  c[7] = get_one (ndom, "t_r");
  strcpy (column_name[7], "t_r");

  c[8] = get_one (ndom, "t_r_old");
  strcpy (column_name[8], "t_r_old");

  c[9] = get_one (ndom, "w");
  strcpy (column_name[9], "w");

  c[10] = get_one (ndom, "heat_tot");
  strcpy (column_name[10], "heat_tot");

  c[11] = get_one (ndom, "heat_tot_old");
  strcpy (column_name[11], "heat_tot_old");

  c[12] = get_one (ndom, "cool_tot");
  strcpy (column_name[12], "cool_tot");

  c[13] = get_one (ndom, "ntot");
  strcpy (column_name[13], "ntot");

  c[14] = get_one (ndom, "ip");
  strcpy (column_name[14], "ip");

  c[15] = get_one (ndom, "nioniz");
  strcpy (column_name[15], "nioniz");

  c[16] = get_one (ndom, "gain");
  strcpy (column_name[16], "gain");

  c[17] = get_one (ndom, "macro_bf_in");
  strcpy (column_name[17], "macro_bf_in");

  c[18] = get_one (ndom, "macro_bf_out");
  strcpy (column_name[18], "macro_bf_out");

  /* This should be the maxium number above +1 */
  ncols = 19;


  converge = get_one (ndom, "converge");

  /* At this point all of the data has been collected */


  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;


  if (zdom[ndom].coord_type == SPHERICAL)
  {


    /*
     * First assemble the header line
     */

    sprintf (start, "%9s %4s %6s %6s %9s %9s %9s ", "r", "i", "inwind", "converge", "v_x", "v_y", "v_z");
    strcpy (one_line, start);
    n = 0;
    while (n < ncols)
    {
      sprintf (one_value, "%9.9s ", column_name[n]);
      strcat (one_line, one_value);

      n++;
    }
    fprintf (fptr, "%s\n", one_line);


    /* Now assemble the lines of the table */

    for (i = 0; i < ndim2; i++)
    {
      //This line is different from the two d case
      sprintf (start, "%9.3e %4d %6d %8.0f %9.2e %9.2e %9.2e ",
               wmain[nstart + i].r, i, wmain[nstart + i].inwind,
               converge[i], wmain[nstart + i].v[0], wmain[nstart + i].v[1], wmain[nstart + i].v[2]);
      strcpy (one_line, start);
      n = 0;
      while (n < ncols)
      {
        sprintf (one_value, "%9.2e ", c[n][i]);
        strcat (one_line, one_value);
        n++;
      }
      fprintf (fptr, "%s\n", one_line);
    }
  }
  else
  {

    /* First assemble the header line */

    sprintf (start, "%8s %8s %4s %4s %6s %8s %9s %9s %9s ", "x", "z", "i", "j", "inwind", "converge", "v_x", "v_y", "v_z");
    strcpy (one_line, start);
    n = 0;
    while (n < ncols)
    {
      sprintf (one_value, "%9.9s ", column_name[n]);
      strcat (one_line, one_value);

      n++;
    }
    fprintf (fptr, "%s\n", one_line);


    /* Now assemble the lines of the table */

    for (i = 0; i < ndim2; i++)
    {
      wind_n_to_ij (ndom, nstart + i, &ii, &jj);
      sprintf (start,
               "%8.2e %8.2e %4d %4d %6d %8.0f %9.2e %9.2e %9.2e ",
               wmain[nstart + i].xcen[0], wmain[nstart + i].xcen[2], ii,
               jj, wmain[nstart + i].inwind, converge[i], wmain[nstart + i].v[0], wmain[nstart + i].v[1], wmain[nstart + i].v[2]);
      strcpy (one_line, start);
      n = 0;
      while (n < ncols)
      {
        sprintf (one_value, "%9.2e ", c[n][i]);
        strcat (one_line, one_value);
        n++;
      }
      fprintf (fptr, "%s\n", one_line);
    }
  }

  fclose (fptr);

  return (0);
}



/**********************************************************/
/**
 * @brief      writes selected variables related to issues about
 * velocity gradients to an ascii file which can be read as an astropy table
 *
 *
 * @param [in] int  ndom   The domain of interest
 * @param [in, out] char  rootname[]   The rootname of the windsave file
 * @return     Always returns 0
 *
 * @details
 *
 * ### Notes ###
 *
 * To add a variable one just needs to define the column_name
 * and send the appropriate call to either get_one or get_ion.
 *
 *
 **********************************************************/

int
create_velocity_gradient_table (ndom, rootname)
     int ndom;
     char rootname[];
{
  char filename[132];
  double *c[50], *converge;
  char column_name[50][20];
  char one_line[1024], start[1024], one_value[20];


  int i, ii, jj;
  int nstart, ndim2;
  int n, ncols;
  FILE *fptr;

  strcpy (filename, rootname);
  strcat (filename, ".gradient.txt");


  fptr = fopen (filename, "w");

  /* Get the variables that one needs */

  c[0] = get_one (ndom, "dv_x_dx");
  strcpy (column_name[0], "dv_x_dx");

  c[1] = get_one (ndom, "dv_y_dx");
  strcpy (column_name[1], "dv_y_dx");

  c[2] = get_one (ndom, "dv_z_dx");
  strcpy (column_name[2], "dv_z_dx");

  c[3] = get_one (ndom, "dv_x_dy");
  strcpy (column_name[3], "dv_x_dy");

  c[4] = get_one (ndom, "dv_y_dy");
  strcpy (column_name[4], "dv_y_dy");

  c[5] = get_one (ndom, "dv_z_dy");
  strcpy (column_name[5], "dv_z_dy");

  c[6] = get_one (ndom, "dv_x_dz");
  strcpy (column_name[6], "dv_x_dz");

  c[7] = get_one (ndom, "dv_y_dz");
  strcpy (column_name[7], "dv_y_dz");

  c[8] = get_one (ndom, "dv_z_dz");
  strcpy (column_name[8], "dv_z_dz");

  c[9] = get_one (ndom, "div_v");
  strcpy (column_name[9], "div_v");

  c[10] = get_one (ndom, "dvds_max");
  strcpy (column_name[10], "dvds_max");

  c[11] = get_one (ndom, "gamma");
  strcpy (column_name[11], "gamma");

  c[12] = get_one (ndom, "dfudge");
  strcpy (column_name[12], "dfudge");

  /* This should be the maxium number above +1 */
  ncols = 13;



  /* At this point all of the data has been collected */


  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;

  converge = get_one (ndom, "converge");


  if (zdom[ndom].coord_type == SPHERICAL)
  {


    /*
     * First assemble the header line
     */

    sprintf (start, "%9s %4s %6s %6s %9s %9s %9s ", "r", "i", "inwind", "converge", "v_x", "v_y", "v_z");
    strcpy (one_line, start);
    n = 0;
    while (n < ncols)
    {
      sprintf (one_value, "%9.9s ", column_name[n]);
      strcat (one_line, one_value);

      n++;
    }
    fprintf (fptr, "%s\n", one_line);


    /* Now assemble the lines of the table */

    for (i = 0; i < ndim2; i++)
    {
      //This line is different from the two d case
      sprintf (start, "%9.3e %4d %6d %8.0f %9.2e %9.2e %9.2e ",
               wmain[nstart + i].r, i, wmain[nstart + i].inwind,
               converge[i], wmain[nstart + i].v[0], wmain[nstart + i].v[1], wmain[nstart + i].v[2]);
      strcpy (one_line, start);
      n = 0;
      while (n < ncols)
      {
        sprintf (one_value, "%9.2e ", c[n][i]);
        strcat (one_line, one_value);
        n++;
      }
      fprintf (fptr, "%s\n", one_line);
    }
  }
  else
  {

    /* First assemble the header line */

    sprintf (start, "%8s %8s %4s %4s %6s %8s %9s %9s %9s ", "x", "z", "i", "j", "inwind", "converge", "v_x", "v_y", "v_z");
    strcpy (one_line, start);
    n = 0;
    while (n < ncols)
    {
      sprintf (one_value, "%9.9s ", column_name[n]);
      strcat (one_line, one_value);

      n++;
    }
    fprintf (fptr, "%s\n", one_line);


    /* Now assemble the lines of the table */

    for (i = 0; i < ndim2; i++)
    {
      wind_n_to_ij (ndom, nstart + i, &ii, &jj);
      sprintf (start,
               "%8.2e %8.2e %4d %4d %6d %8.0f %9.2e %9.2e %9.2e ",
               wmain[nstart + i].xcen[0], wmain[nstart + i].xcen[2], ii,
               jj, wmain[nstart + i].inwind, converge[i], wmain[nstart + i].v[0], wmain[nstart + i].v[1], wmain[nstart + i].v[2]);
      strcpy (one_line, start);
      n = 0;
      while (n < ncols)
      {
        sprintf (one_value, "%9.2e ", c[n][i]);
        strcat (one_line, one_value);
        n++;
      }
      fprintf (fptr, "%s\n", one_line);
    }
  }

  fclose (fptr);

  return (0);
}



/**********************************************************/
/**
 * @brief      makes an astropy table containing the relative abundances
 *  of a given element as a function of the position in the grid
 *
 * @param [in] int  ndom   The domain number
 * @param [in] char  rootname[]   The rootname of the windsave file
 * @param [in] int  iz   atomic number of the element
 * @return     0 on completion
 *
 * @details
 *
 * ### Notes ###
 *
 * The routine calls get_ion mulitiple times, once for
 * each ion of an atom
 *
 **********************************************************/

int
create_ion_table (ndom, rootname, iz, ion_switch)
     int ndom;
     char rootname[];
     int iz;
     int ion_switch;
{
  char filename[132];
  double *c[100];
  int first_ion, number_ions;
  char element_name[20];
  int istate[50];
  char one_line[1024], start[1024], one_value[20];
  int nstart, ndim2;
  char name[132];               /* Name of extension for a
                                 * particular ion table */


  int i, ii, jj, n;
  FILE *fptr;


  i = 0;
  while (i < nelements)
  {
    if (ele[i].z == iz)
    {
      break;
    }
    i++;
  }


  first_ion = ele[i].firstion;
  number_ions = ele[i].nions;
  strcpy (element_name, ele[i].name);


  i = 0;
  while (i < number_ions)
  {
    istate[i] = ion[first_ion + i].istate;

    c[i] = get_ion (ndom, iz, istate[i], ion_switch, name);
    i++;
  }


  /* Open the output file */

  sprintf (filename, "%.100s.%.10s.%.10s.txt", rootname, element_name, name);

  fptr = fopen (filename, "w");



  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;



  if (zdom[ndom].coord_type == SPHERICAL)
  {


    /*
     * First assemble the header line
     */

    sprintf (start, "%8s %4s %6s ", "r", "i", "inwind");

    strcpy (one_line, start);
    n = 0;
    while (n < number_ions)
    {
      sprintf (one_value, "     i%02d ", istate[n]);
      strcat (one_line, one_value);

      n++;
    }


    fprintf (fptr, "%s\n", one_line);


    /* Now assemble the lines of the table */

    for (i = 0; i < ndim2; i++)
    {
      //This line is different from the two d case
      sprintf (start, "%8.2e %4d %6d ", wmain[nstart + i].r, i, wmain[nstart + i].inwind);
      strcpy (one_line, start);

      n = 0;
      while (n < number_ions)
      {
        sprintf (one_value, "%8.2e ", c[n][i]);
        strcat (one_line, one_value);
        n++;
      }


      fprintf (fptr, "%s\n", one_line);
    }
  }
  else
  {


    /* First assemble the header line */

    sprintf (start, "%8s %8s %4s %4s %6s ", "x", "z", "i", "j", "inwind");
    strcpy (one_line, start);
    n = 0;
    while (n < number_ions)
    {
      sprintf (one_value, "     i%02d ", istate[n]);
      strcat (one_line, one_value);

      n++;
    }

    fprintf (fptr, "%s\n", one_line);

    /* Now assemble the lines of the table */

    for (i = 0; i < ndim2; i++)
    {
      wind_n_to_ij (ndom, nstart + i, &ii, &jj);
      sprintf (start, "%8.2e %8.2e %4d %4d %6d ", wmain[nstart + i].xcen[0], wmain[nstart + i].xcen[2], ii, jj, wmain[nstart + i].inwind);
      strcpy (one_line, start);
      n = 0;
      while (n < number_ions)
      {
        sprintf (one_value, "%8.2e ", c[n][i]);
        strcat (one_line, one_value);
        n++;
      }
      fprintf (fptr, "%s\n", one_line);
    }
  }

  fclose (fptr);
  return (0);

}



/**********************************************************/
/**
 * @brief      Get get density, etc for one particular ion
 *
 * @param [in] int  ndom   the domain number
 * @param [in] int  element   the element number
 * @param [in] int  istate   the ionization state
 * @param [in] int  iswitch   a switch controlling exactly what is returned for that ion
 * @return     Normally returns an array with values associated with what is requested
 *    	This will return an array with all zeros if there is no such ion
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

double *
get_ion (ndom, element, istate, iswitch, name)
     int ndom, element, istate, iswitch;
     char *name;
{
  int nion, nelem;
  int n;
  int nplasma;
  double *x;
  int nstart, ndim2;
  double nh;


  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;

  x = (double *) calloc (sizeof (double), ndim2);

  /* Find the ion */

  nion = 0;
  while (nion < nions && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;
  if (nion == nions)
  {
    Log ("Error--element %d ion %d not found in define_wind\n", element, istate);
    return (x);
  }
  nelem = 0;
  while (nelem < nelements && ele[nelem].z != element)
    nelem++;


  /* Now populate the array */

  for (n = 0; n < ndim2; n++)
  {
    x[n] = 0;
    nplasma = wmain[nstart + n].nplasma;
    if (wmain[nstart + n].inwind >= 0 && plasmamain[nplasma].rho > 0.0)
    {
      if (iswitch == 0)
      {
        x[n] = plasmamain[nplasma].density[nion];
        nh = rho2nh * plasmamain[nplasma].rho;
        x[n] /= (nh * ele[nelem].abun);
        strcpy (name, "frac");
      }
      else if (iswitch == 1)
      {
        x[n] = plasmamain[nplasma].density[nion];
        strcpy (name, "den");
      }
      else if (iswitch == 2)
      {
        x[n] = (double) plasmamain[nplasma].scatters[nion] / plasmamain[nplasma].vol;
        strcpy (name, "scat");
      }
      else if (iswitch == 3)
      {
        x[n] = plasmamain[nplasma].xscatters[nion];
        strcpy (name, "ion_frac");
      }
      else if (iswitch == 4)
      {
        x[n] = plasmamain[nplasma].ioniz[nion];
        strcpy (name, "ioniz");
      }
      else if (iswitch == 5)
      {
        x[n] = plasmamain[nplasma].recomb[nion];
        strcpy (name, "recomb");
      }
      else if (iswitch == 6)
      {
        x[n] = plasmamain[nplasma].heat_ion[nion];
        strcpy (name, "heat");
      }
      else if (iswitch == 7)
      {
        x[n] = plasmamain[nplasma].cool_rr_ion[nion];
        strcpy (name, "cool_rr");
      }
      else if (iswitch == 8)
      {
        x[n] = plasmamain[nplasma].lum_rr_ion[nion];
        strcpy (name, "lum_rr");
      }
      else if (iswitch == 9)
      {
        x[n] = plasmamain[nplasma].cool_dr_ion[nion];
        strcpy (name, "cool_dr");
      }
      else
      {
        Error ("get_ion : Unknown switch %d \n", iswitch);
        exit (0);
      }
    }
  }

  return (x);
}



/**********************************************************/
/**
 * @brief      Get a simple variable from the PlasmaPtr array
 *
 * @param [in] int  ndom   The domain in question
 * @param [in] char  variable_name[]   The name of the variable
 * @return     The values in the plasma pointer for this variable. A double
 * 	will be returned even if the PlasmaPtr varible is an integer
 *
 *
 * @details
 *
 * A simple variable is a variable that is just a number, not an array
 *
 * The routine performes a simple tranlation of the character name
 * to a variable in the PlasmaPtr.
 *
 * ### Notes ###
 * Normally returns non-zero values only if a cell is in the wind
 * but this can be changed if external variable xedge is TRUE.
 *
 * Only selected variables are returned, but new variables are easy
 * to add using the template of the other variables
 *
 **********************************************************/

double *
get_one (ndom, variable_name)
     int ndom;
     char variable_name[];
{
  int n;
  int nplasma;
  double *x;
  int ndim2;
  int nstart;

  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;

  x = (double *) calloc (sizeof (double), ndim2);

  for (n = 0; n < ndim2; n++)
  {
    x[n] = 0;
    if (wmain[n + nstart].inwind >= 0 || xedge)
    {
      nplasma = wmain[n + nstart].nplasma;


      if (strcmp (variable_name, "ne") == 0)
      {
        x[n] = plasmamain[nplasma].ne;
      }
      else if (strcmp (variable_name, "rho") == 0)
      {
        x[n] = plasmamain[nplasma].rho;
      }
      else if (strcmp (variable_name, "vol") == 0)
      {
        x[n] = plasmamain[nplasma].vol;
      }
      else if (strcmp (variable_name, "t_e") == 0)
      {
        x[n] = plasmamain[nplasma].t_e;
      }
      else if (strcmp (variable_name, "t_r") == 0)
      {
        x[n] = plasmamain[nplasma].t_r;
      }
      else if (strcmp (variable_name, "t_e_old") == 0)
      {
        x[n] = plasmamain[nplasma].t_e_old;
      }
      else if (strcmp (variable_name, "t_r_old") == 0)
      {
        x[n] = plasmamain[nplasma].t_r_old;
      }
      else if (strcmp (variable_name, "dt_e") == 0)
      {
        x[n] = plasmamain[nplasma].dt_e;
      }
      else if (strcmp (variable_name, "dt_e_old") == 0)
      {
        x[n] = plasmamain[nplasma].dt_e_old;
      }
      else if (strcmp (variable_name, "J") == 0)
      {
        x[n] = plasmamain[nplasma].j;
      }
      else if (strcmp (variable_name, "J_direct") == 0)
      {
        x[n] = plasmamain[nplasma].j_direct;
      }
      else if (strcmp (variable_name, "J_scatt") == 0)
      {
        x[n] = plasmamain[nplasma].j_scatt;
      }
      else if (strcmp (variable_name, "ave_freq") == 0)
      {
        x[n] = plasmamain[nplasma].ave_freq;
      }
      else if (strcmp (variable_name, "converge") == 0)
      {
        x[n] = plasmamain[nplasma].converge_whole;
      }
      else if (strcmp (variable_name, "dmo_dt_x") == 0)
      {
        x[n] = plasmamain[nplasma].dmo_dt[0];
      }
      else if (strcmp (variable_name, "dmo_dt_y") == 0)
      {
        x[n] = plasmamain[nplasma].dmo_dt[1];
      }
      else if (strcmp (variable_name, "dmo_dt_z") == 0)
      {
        x[n] = plasmamain[nplasma].dmo_dt[2];
      }
      else if (strcmp (variable_name, "ntot") == 0)
      {
        x[n] = plasmamain[nplasma].ntot;
      }
      else if (strcmp (variable_name, "ip") == 0)
      {
        x[n] = plasmamain[nplasma].ip;
      }
      else if (strcmp (variable_name, "xi") == 0)
      {
        x[n] = plasmamain[nplasma].xi;
      }
      else if (strcmp (variable_name, "heat_tot") == 0)
      {
        x[n] = plasmamain[nplasma].heat_tot;
      }
      else if (strcmp (variable_name, "heat_tot_old") == 0)
      {
        x[n] = plasmamain[nplasma].heat_tot_old;
      }
      else if (strcmp (variable_name, "heat_comp") == 0)
      {
        x[n] = plasmamain[nplasma].heat_comp;
      }
      else if (strcmp (variable_name, "heat_lines") == 0)
      {
        x[n] = plasmamain[nplasma].heat_lines;
      }
      else if (strcmp (variable_name, "heat_ff") == 0)
      {
        x[n] = plasmamain[nplasma].heat_ff;
      }
      else if (strcmp (variable_name, "heat_photo") == 0)
      {
        x[n] = plasmamain[nplasma].heat_photo;
      }
      else if (strcmp (variable_name, "heat_auger") == 0)
      {
        x[n] = plasmamain[nplasma].heat_auger;
      }
      else if (strcmp (variable_name, "cool_comp") == 0)
      {
        x[n] = plasmamain[nplasma].cool_comp;
      }
      else if (strcmp (variable_name, "lum_tot") == 0)
      {
        x[n] = plasmamain[nplasma].lum_tot;
      }
      else if (strcmp (variable_name, "lum_lines") == 0)
      {
        x[n] = plasmamain[nplasma].lum_lines;
      }
      else if (strcmp (variable_name, "lum_ff") == 0)
      {
        x[n] = plasmamain[nplasma].lum_ff;
      }
      else if (strcmp (variable_name, "lum_rr") == 0)
      {
        x[n] = plasmamain[nplasma].lum_rr;
      }
      else if (strcmp (variable_name, "cool_rr") == 0)
      {
        x[n] = plasmamain[nplasma].cool_rr;
      }
      else if (strcmp (variable_name, "cool_dr") == 0)
      {
        x[n] = plasmamain[nplasma].cool_dr;
      }
      else if (strcmp (variable_name, "cool_tot") == 0)
      {
        x[n] = plasmamain[nplasma].cool_tot;
      }
      else if (strcmp (variable_name, "w") == 0)
      {
        x[n] = plasmamain[nplasma].w;
      }
      else if (strcmp (variable_name, "nrad") == 0)
      {
        x[n] = plasmamain[nplasma].nrad;
      }
      else if (strcmp (variable_name, "nioniz") == 0)
      {
        x[n] = plasmamain[nplasma].nioniz;
      }
      else if (strcmp (variable_name, "heat_shock") == 0)
      {
        x[n] = plasmamain[nplasma].heat_shock;
      }
      else if (strcmp (variable_name, "cool_adiab") == 0)
      {
        x[n] = plasmamain[nplasma].cool_adiabatic;
      }
      else if (strcmp (variable_name, "heat_lines_macro") == 0)
      {
        x[n] = plasmamain[nplasma].heat_lines_macro;
      }
      else if (strcmp (variable_name, "heat_photo_macro") == 0)
      {
        x[n] = plasmamain[nplasma].heat_photo_macro;
      }
      else if (strcmp (variable_name, "gain") == 0)
      {
        x[n] = plasmamain[nplasma].gain;
      }
      else if (strcmp (variable_name, "macro_bf_in") == 0)
      {
        x[n] = plasmamain[nplasma].bf_simple_ionpool_in;
      }
      else if (strcmp (variable_name, "macro_bf_out") == 0)
      {
        x[n] = plasmamain[nplasma].bf_simple_ionpool_out;
      }
      else if (strcmp (variable_name, "dv_x_dx") == 0)
      {
        x[n] = wmain[n].v_grad[0][0];
      }
      else if (strcmp (variable_name, "dv_x_dy") == 0)
      {
        x[n] = wmain[n].v_grad[0][1];
      }
      else if (strcmp (variable_name, "dv_x_dz") == 0)
      {
        x[n] = wmain[n].v_grad[0][2];
      }
      else if (strcmp (variable_name, "dv_y_dx") == 0)
      {
        x[n] = wmain[n].v_grad[1][0];
      }
      else if (strcmp (variable_name, "dv_y_dy") == 0)
      {
        x[n] = wmain[n].v_grad[1][1];
      }
      else if (strcmp (variable_name, "dv_y_dz") == 0)
      {
        x[n] = wmain[n].v_grad[1][2];
      }
      else if (strcmp (variable_name, "dv_z_dx") == 0)
      {
        x[n] = wmain[n].v_grad[2][0];
      }
      else if (strcmp (variable_name, "dv_z_dy") == 0)
      {
        x[n] = wmain[n].v_grad[2][1];
      }
      else if (strcmp (variable_name, "dv_z_dz") == 0)
      {
        x[n] = wmain[n].v_grad[2][2];
      }
      else if (strcmp (variable_name, "dvds_max") == 0)
      {
        x[n] = wmain[n].dvds_max;
      }
      else if (strcmp (variable_name, "div_v") == 0)
      {
        x[n] = wmain[n].div_v;
      }
      else if (strcmp (variable_name, "gamma") == 0)
      {
        x[n] = wmain[n].xgamma;
      }
      else if (strcmp (variable_name, "dfudge") == 0)
      {
        x[n] = wmain[n].dfudge;
      }
      else
      {
        Error ("get_one: Unknown variable %s\n", variable_name);
      }
    }
  }

  return (x);

}

/**********************************************************/
/**
 * @brief      Get a simple array element from the PlasmaPtr array
 *
 * @param [in] int  ndom   The domain in question
 * @param [in] char  variable_name[]   The name of the variable
 * @param [in] int  array_dim  The number of element in the array
 * @param [out] *double  xval An array containing the requested values
 * @return  Always returns 0   
 *  The values in the plasma pointer for this variable. A double
 * 	will be returned even if the PlasmaPtr variable is an integer
 *      The results are a 1-d array, where the various array elements
 *      have effectively been concatenated
 *
 *
 * @details
 *
 * An  array element here means a one element of a 1d array
 *
 * The routine performes a simple translation of the character name
 * to a variable in the PlasmaPtr.
 * 
 * The values in the plasma pointer for this variable. A double
 * will be returned even if the PlasmaPtr variable is an integer
 * The results are a 1-d array, where the various array elements
 * have effectively been concatenated
 *
 * ### Notes ###
 * This was written to get information about crude fits to the 
 * spectra in a cell, where one wants to create a single astropy
 * table that contains all, for example, photons in each band. One
 * would only use it in a situation, where making a separate column
 * for each band would create a table that is too wide, for normal use
 * or when one does not know before hand, how many elements are in the
 * array.

 * Only selected variables are returned, but new variables are easy
 * to add using the template of the other variables
 *
 **********************************************************/

int
get_one_array_element (ndom, variable_name, array_dim, xval)
     int ndom;
     char variable_name[];
     int array_dim;
     double xval[];
{
  int j, n;
  int nplasma;
  int ndim2;
  int nstart;
  int m;

  /* In order for a routine to return an array, one must
     declare it as static */


  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;


  m = 0;
  for (j = 0; j < array_dim; j++)
  {
    for (n = 0; n < ndim2; n++)
    {
      xval[m] = 0;
      if (wmain[n + nstart].inwind >= 0)
      {
        nplasma = wmain[n + nstart].nplasma;


        if (strcmp (variable_name, "xj") == 0)
        {
          xval[m] = plasmamain[nplasma].xj[j];
        }
        else if (strcmp (variable_name, "xave_freq") == 0)
        {
          xval[m] = plasmamain[nplasma].xave_freq[j];
        }
        else if (strcmp (variable_name, "fmin") == 0)
        {
          xval[m] = plasmamain[nplasma].fmin[j];
        }
        else if (strcmp (variable_name, "fmax") == 0)
        {
          xval[m] = plasmamain[nplasma].fmax[j];
        }
        else if (strcmp (variable_name, "fmin_mod") == 0)
        {
          xval[m] = plasmamain[nplasma].fmin_mod[j];
        }
        else if (strcmp (variable_name, "fmax_mod") == 0)
        {
          xval[m] = plasmamain[nplasma].fmax_mod[j];
        }
        else if (strcmp (variable_name, "xsd_freq") == 0)
        {
          xval[m] = plasmamain[nplasma].xsd_freq[j];
        }
        else if (strcmp (variable_name, "nxtot") == 0)
        {
          xval[m] = plasmamain[nplasma].nxtot[j];
        }
        else if (strcmp (variable_name, "spec_mod_type") == 0)
        {
          xval[m] = plasmamain[nplasma].spec_mod_type[j];
        }
        else if (strcmp (variable_name, "pl_alpha") == 0)
        {
          xval[m] = plasmamain[nplasma].pl_alpha[j];
        }
        else if (strcmp (variable_name, "pl_log_w") == 0)
        {
          xval[m] = plasmamain[nplasma].pl_log_w[j];
        }
        else if (strcmp (variable_name, "exp_temp") == 0)
        {
          xval[m] = plasmamain[nplasma].exp_temp[j];
        }
        else if (strcmp (variable_name, "exp_w") == 0)
        {
          xval[m] = plasmamain[nplasma].exp_w[j];
        }
        else
        {
          Error ("get_on_array_elemente: Unknown variable %s\n", variable_name);
        }
      }
      m++;
    }
  }

  return (0);

}



/**********************************************************/
/**
 * @brief      writes specific  variables associated with
 * model spectrua of a windsaave
 * which are intended to be of general interest to
 * file which has the format of an astropy table
 *
 * 	It is intended to be easily modifible.
 *
 * @param [in] int  ndom   A domain number
 * @param [in] char  rootname   The rootname of the master file
 * @return   Always returns 0
 *
 * @details
 *
 * The master table is contains basic information for each
 * cell in the wind, such as the electron density, the density,
 * the ionization parameter, and the radiative and electron temperature
 *
 * The routine takes data directly from wmain, and then calls
 * get_one or get_ion multiple times to fet info from the Plasma
 * structure.
 *
 * It then writes the data to an ascii file which can be read as
 * an  astropy table
 *
 * ### Notes ###
 * To add a variable one just needs to define the column_name
 * and send the appropriate call to either get_one or get_ion.
 *
 * There is some duplicated code in the routine that pertains
 * to whether one is dealing with a spherecial or a 2d coordinate
 * system.  It should be possible to delete this
 *
 **********************************************************/

int
create_spec_table (ndom, rootname)
     int ndom;
     char rootname[];
{
  char filename[132];
  double *c[50], *converge;
  char column_name[50][20];
  char one_line[1024], start[1024], one_value[20];

  int nxfreq;


  int i, ii, jj;
  int nstart, ndim2;
  int n, ncols;
  int nx;
  int j;
  FILE *fptr;

  nxfreq = geo.nxfreq;
  ndim2 = zdom[ndom].ndim2;

  strcpy (filename, rootname);
  strcat (filename, ".spec.txt");


  fptr = fopen (filename, "w");

  /* Get the variables that one needs */

  for (i = 0; i < 50; i++)
  {
    c[i] = (double *) calloc (sizeof (double), ndim2 * nxfreq);
  }

  get_one_array_element (ndom, "nxtot", nxfreq, c[0]);
  strcpy (column_name[0], "nxtot");

  get_one_array_element (ndom, "xj", nxfreq, c[1]);
  strcpy (column_name[1], "xj");

  get_one_array_element (ndom, "xave_freq", nxfreq, c[2]);
  strcpy (column_name[2], "xave_freq");

  get_one_array_element (ndom, "xsd_freq", nxfreq, c[3]);
  strcpy (column_name[3], "xsd_freq");

  get_one_array_element (ndom, "fmin", nxfreq, c[4]);
  strcpy (column_name[4], "fmin");

  get_one_array_element (ndom, "fmax", nxfreq, c[5]);
  strcpy (column_name[5], "fmax");

  get_one_array_element (ndom, "fmin_mod", nxfreq, c[6]);
  strcpy (column_name[6], "fmin_mod");

  get_one_array_element (ndom, "fmax_mod", nxfreq, c[7]);
  strcpy (column_name[7], "fmax_mod");

  get_one_array_element (ndom, "spec_mod_type", nxfreq, c[8]);
  strcpy (column_name[8], "spec_mod_type");

  get_one_array_element (ndom, "pl_alpha", nxfreq, c[9]);
  strcpy (column_name[9], "pl_alpha");

  get_one_array_element (ndom, "pl_log_w", nxfreq, c[10]);
  strcpy (column_name[10], "pl_log_w");

  get_one_array_element (ndom, "exp_temp", nxfreq, c[11]);
  strcpy (column_name[11], "exp_temp");

  get_one_array_element (ndom, "exp_w", nxfreq, c[12]);
  strcpy (column_name[12], "exp_w");

  /* This should be the maximum number above +1 */
  ncols = 13;


  converge = get_one (ndom, "converge");

  /* At this point oll of the data has been collected */


  nstart = zdom[ndom].nstart;


  if (zdom[ndom].coord_type == SPHERICAL)
  {


    /*
     * First assemble the header line
     */

    sprintf (start, "%9s %9s %4s %6s %6s %9s  ", "r", "rcen", "i", "inwind", "converge", "nband");
    strcpy (one_line, start);
    n = 0;
    while (n < ncols)
    {
      sprintf (one_value, "%9.9s ", column_name[n]);
      strcat (one_line, one_value);

      n++;
    }
    fprintf (fptr, "%s\n", one_line);


    /* Now assemble the lines of the table */

    j = 0;
    for (nx = 0; nx < nxfreq; nx++)
    {

      for (i = 0; i < ndim2; i++)
      {
        //This line is different from the two d case
        sprintf (start, "%9.3e %9.3e %4d %6d %8.0f %6d ",
                 wmain[nstart + i].r, wmain[nstart + i].rcen, i, wmain[nstart + i].inwind, converge[i], nx);
        strcpy (one_line, start);
        n = 0;
        while (n < ncols)
        {
          sprintf (one_value, "%9.2e ", c[n][j]);
          strcat (one_line, one_value);
          n++;
        }
        fprintf (fptr, "%s\n", one_line);
        j++;
      }
    }
  }
  else if (zdom[ndom].coord_type == CYLIND)
  {

    /* First assemble the header line */

    sprintf (start, "%8s %8s %8s %8s %4s %4s %6s %8s %6s ", "x", "z", "xcen", "zcen", "i", "j", "inwind", "converge", "nband");
    strcpy (one_line, start);
    n = 0;
    while (n < ncols)
    {
      sprintf (one_value, "%9.9s ", column_name[n]);
      strcat (one_line, one_value);

      n++;
    }
    fprintf (fptr, "%s\n", one_line);


    /* Now assemble the lines of the table */

    j = 0;
    for (nx = 0; nx < nxfreq; nx++)
    {
      for (i = 0; i < ndim2; i++)
      {
        wind_n_to_ij (ndom, nstart + i, &ii, &jj);
        sprintf (start,
                 "%8.2e %8.2e %8.2e %8.2e %4d %4d %6d %8.0f %6d  ",
                 wmain[nstart + i].x[0], wmain[nstart + i].x[2], wmain[nstart + i].xcen[0], wmain[nstart + i].xcen[2], ii,
                 jj, wmain[nstart + i].inwind, converge[i], nx);
        strcpy (one_line, start);
        n = 0;
        while (n < ncols)
        {
          sprintf (one_value, "%9.2e ", c[n][j]);
          strcat (one_line, one_value);
          n++;
        }
        fprintf (fptr, "%s\n", one_line);
        j++;
      }
    }
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {

    /* First assemble the header line */

    sprintf (start, "%8s %8s %8s %9s %8s %8s %8s %8s %4s %4s %6s %8s %6s ", "r", "theta", "r_cen", "theta_cen", "x", "z", "xcen",
             "zcen", "i", "j", "inwind", "converge", "nband");
    strcpy (one_line, start);
    n = 0;
    while (n < ncols)
    {
      sprintf (one_value, "%9.9s ", column_name[n]);
      strcat (one_line, one_value);

      n++;
    }
    fprintf (fptr, "%s\n", one_line);


    /* Now assemble the lines of the table */

    j = 0;
    for (nx = 0; nx < nxfreq; nx++)
    {
      for (i = 0; i < ndim2; i++)
      {
        wind_n_to_ij (ndom, nstart + i, &ii, &jj);
        sprintf (start,
                 "%8.2e %8.2e %8.2e %9.2e %8.2e %8.2e %8.2e %8.2e %4d %4d %6d %8.0f %6d ",
                 wmain[nstart + i].r, wmain[nstart + i].theta, wmain[nstart + i].rcen, wmain[nstart + i].thetacen,
                 wmain[nstart + i].x[0], wmain[nstart + i].x[2], wmain[nstart + i].xcen[0], wmain[nstart + i].xcen[2], ii,
                 jj, wmain[nstart + i].inwind, converge[i], nx);
        strcpy (one_line, start);
        n = 0;
        while (n < ncols)
        {
          sprintf (one_value, "%9.2e ", c[n][j]);
          strcat (one_line, one_value);
          n++;
        }
        fprintf (fptr, "%s\n", one_line);
        j++;
      }
    }
  }
  else
  {
    printf ("Error: Cannot print out files for coordinate system type %d\n", zdom[ndom].coord_type);
  }


  fclose (fptr);

  return (0);
}


/**********************************************************/
/**
 * @brief      write the detailed spectrum in a cell
 * to a file
 * @param [in] int  ncell   The number of a cell in the wind
 * @param [in] char  rootname   The rootname of the master file
 * @return   Always returns 0
 *
 * @details
 *
 *
 **********************************************************/

int
create_detailed_cell_spec_table (ncell, rootname)
     int ncell;
     char rootname[];
{
  FILE *fptr;
  char filename[132];

  double freq[NBINS_IN_CELL_SPEC];
  double flux[NBINS_IN_CELL_SPEC];
  int i, nplasma;

  printf ("%e %e\n", geo.cell_log_freq_min, geo.cell_delta_lfreq);

  sprintf (filename, "%s.xspec.%d.txt", rootname, ncell);

  for (i = 0; i < NBINS_IN_CELL_SPEC; i++)
  {
    freq[i] = pow (10., geo.cell_log_freq_min + i * geo.cell_delta_lfreq);
  }

  nplasma = wmain[ncell].nplasma;

  for (i = 0; i < NBINS_IN_CELL_SPEC; i++)
  {
    flux[i] = plasmamain[nplasma].cell_spec_flux[i];
  }


  fptr = fopen (filename, "w");

  fprintf (fptr, "Freq.          Flux\n");

  for (i = 0; i < NBINS_IN_CELL_SPEC; i++)
  {
    fprintf (fptr, "%10.3e  %10.3e\n", freq[i], flux[i]);
  }


  fclose (fptr);

  return (0);
}

#define MAX_COLUMNS 10000
#define MAX_IN_TABLE 1000


/**********************************************************/
/**
 * @brief      write the detailed spectra in cells in a particular
 * domain to a file or files
 * @param [in] int  ncell   The domain number 
 * @param [in] char  rootname   The rootname of the master file
 * @return   Always returns 0
 *
 * @details
 *
 * The routine simply reads data in stored in the cell_spec_flux
 * array of the Plasma structure
 *
 * Notes:
 *
 * The column names for the files represent the wind cell number
 * (i,j) for for cylindrical and rtheta coordinates, (i) for 
 * spherical coordinates.  
 *
 * Multiple files are written out if there are so many cells
 * in the wind that the length of the lines seems impossibly
 * large.  Whether this is a real limit or not is debatable.
 *
 * If there are multiple domains, the domain number needs
 * to be incorporated into the rootname.
 *
 **********************************************************/

int
create_big_detailed_spec_table (ndom, rootname)
     int ndom;
     char *rootname;
{
  char column_name[MAX_COLUMNS][20];
  int nplasma[MAX_COLUMNS], ii, jj, ncols;
  int i, n, nwind_start, nwind_stop;

  FILE *fptr;
  char filename[132];
  double freq[NBINS_IN_CELL_SPEC];
  int nstart, nstop;

  /* Identify the range of wind cells for this domain */

  nwind_start = zdom[ndom].nstart;
  nwind_stop = zdom[ndom].nstop;

  ncols = 0;
  for (n = nwind_start; n < nwind_stop; n++)
  {
    if (wmain[n].inwind >= 0)
    {
      nplasma[ncols] = wmain[n].nplasma;
      if (zdom[ndom].coord_type == SPHERICAL)
      {
        sprintf (column_name[ncols], "F%03d      ", n - nwind_start);
      }
      else
      {
        wind_n_to_ij (ndom, n, &ii, &jj);
        sprintf (column_name[ncols], "F%03d_%03d  ", ii, jj);
      }

      ncols++;
      if (ncols == MAX_COLUMNS)

      {
        printf ("There are more than %d cells in the wind, increase MAX_COLUMNS to get the rest\n", MAX_COLUMNS);
        break;
      }

    }
  }

  /* Now we have the column names and we know which fluxes to print out, and we
     can write out the header and the data */





  /* Calculate the frequencies */
  for (i = 0; i < NBINS_IN_CELL_SPEC; i++)
  {
    freq[i] = pow (10., geo.cell_log_freq_min + i * geo.cell_delta_lfreq);
  }

  nstart = 0;
  nstop = nstart + MAX_IN_TABLE;
  while (nstart < ncols)
  {
    if (nstop > ncols)
    {
      nstop = ncols;
    }

    if (ncols < MAX_IN_TABLE)
    {
      sprintf (filename, "%s.xspec.all.txt", rootname);
    }
    else
    {
      sprintf (filename, "%s.xspec.%02d.txt", rootname, nstart / MAX_IN_TABLE);
    }

    fptr = fopen (filename, "w");

    fprintf (fptr, "Freq.       ");
    for (n = nstart; n < nstop; n++)
    {
      fprintf (fptr, "%-10s ", column_name[n]);
    }
    fprintf (fptr, "\n");

    for (i = 0; i < NBINS_IN_CELL_SPEC; i++)
    {
      fprintf (fptr, "%10.3e ", freq[i]);

      for (n = nstart; n < nstop; n++)
      {
        fprintf (fptr, "%10.3e ", plasmamain[nplasma[n]].cell_spec_flux[i]);
      }

      fprintf (fptr, "\n");
    }


    fclose (fptr);

    nstart += MAX_IN_TABLE;
    nstop = nstart + MAX_IN_TABLE;
  }

  return (0);






}
