
/***********************************************************/
/** @file  windsave.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  Routines to save and read in the structues the constitute
 * a model and and spectra which have been genrerated
 *
 * 
 * The first two routines in this file write and read the wind structure.  		
 * The second two routines do the same thing for the spectrum structure
 * 
 * ### Notes ###
 *
 * The files here are all written out as binary files.  They are
 * used for restars, and also by routines like swind and windsave2talbe
 * which inspect what is happening in the wind.
 *
 * There are separate ascii_writing 
 * routines for writing the spectra out for plotting.)
 * 
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/** 
 * @brief      Save all of the strutures associated with the 
 * wind to a file
 *
 * @param [in] char  filename[]   The name of the file to write to
 * @return     The number of successful writes
 *
 * @details
 *
 * ### Notes ###
 *
 * For the most part, adding a variable to the structures geo,
 * or plasma, does not require changes to this routine, unless
 * new variable length arrays are involved.
 *
 **********************************************************/

int
wind_save (filename)
     char filename[];
{
  FILE *fptr;
  char header[LINELENGTH];
  int ndom;
  int m;
  int n;

  if ((fptr = fopen (filename, "w")) == NULL)
  {
    Error ("wind_save: Unable to open %s\n", filename);
    Exit (0);
  }

  sprintf (header, "Version %s\n", VERSION);
  n = fwrite (header, sizeof (header), 1, fptr);
  n += fwrite (&geo, sizeof (geo), 1, fptr);

  n += fwrite (zdom, sizeof (domain_dummy), geo.ndomain, fptr);
  for (ndom = 0; ndom < geo.ndomain; ++ndom)
  {
    n += fwrite (zdom[ndom].wind_x, sizeof (double), zdom[ndom].ndim, fptr);
    n += fwrite (zdom[ndom].wind_z, sizeof (double), zdom[ndom].mdim, fptr);
    n += fwrite (zdom[ndom].wind_midx, sizeof (double), zdom[ndom].ndim, fptr);
    n += fwrite (zdom[ndom].wind_midz, sizeof (double), zdom[ndom].mdim, fptr);

    if (zdom[ndom].coord_type == CYLVAR)
    {
      n += fwrite (zdom[ndom].wind_z_var, sizeof (double), zdom[ndom].ndim * zdom[ndom].mdim, fptr);
      n += fwrite (zdom[ndom].wind_midz_var, sizeof (double), zdom[ndom].ndim * zdom[ndom].mdim, fptr);
    }
  }

  n += fwrite (wmain, sizeof (wind_dummy), NDIM2, fptr);
  n += fwrite (&disk, sizeof (disk), 1, fptr);
  n += fwrite (&qdisk, sizeof (disk), 1, fptr);
  n += fwrite (plasmamain, sizeof (plasma_dummy), NPLASMA, fptr);

/* Write out the variable length arrays
in the plasma structure */

  for (m = 0; m < NPLASMA; m++)
  {
    n += fwrite (plasmamain[m].density, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].partition, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].ioniz, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].recomb, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].inner_recomb, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].scatters, sizeof (int), nions, fptr);
    n += fwrite (plasmamain[m].xscatters, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].heat_ion, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].cool_rr_ion, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].cool_dr_ion, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].lum_rr_ion, sizeof (double), nions, fptr);
    n += fwrite (plasmamain[m].levden, sizeof (double), nlte_levels, fptr);
    n += fwrite (plasmamain[m].recomb_simple, sizeof (double), nphot_total, fptr);
    n += fwrite (plasmamain[m].recomb_simple_upweight, sizeof (double), nphot_total, fptr);
    n += fwrite (plasmamain[m].kbf_use, sizeof (double), nphot_total, fptr);
  }

/* Now write out the macro atom info */

  if (geo.nmacro)
  {
    n += fwrite (macromain, sizeof (macro_dummy), NPLASMA, fptr);
    for (m = 0; m < NPLASMA; m++)
    {
      n += fwrite (macromain[m].jbar, sizeof (double), size_Jbar_est, fptr);
      n += fwrite (macromain[m].jbar_old, sizeof (double), size_Jbar_est, fptr);
      n += fwrite (macromain[m].gamma, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].gamma_old, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].gamma_e, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].gamma_e_old, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].alpha_st, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].alpha_st_old, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].alpha_st_e, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].alpha_st_e_old, sizeof (double), size_gamma_est, fptr);
      n += fwrite (macromain[m].recomb_sp, sizeof (double), size_alpha_est, fptr);
      n += fwrite (macromain[m].recomb_sp_e, sizeof (double), size_alpha_est, fptr);
      n += fwrite (macromain[m].matom_emiss, sizeof (double), nlevels_macro, fptr);
      n += fwrite (macromain[m].matom_abs, sizeof (double), nlevels_macro, fptr);
    }
  }

  fclose (fptr);

  Log_silent
    ("wind_write sizes: NPLASMA %d size_Jbar_est %d size_gamma_est %d size_alpha_est %d nlevels_macro %d\n",
     NPLASMA, size_Jbar_est, size_gamma_est, size_alpha_est, nlevels_macro);

  return (n);

}

/*

   wind_read (filename)

   History
	11dec	ksl	Updated so returns -1 if it cannot open the windsave file.  This
			was done to enable one to handle missing files differently in
			different cases
	14jul	nsh	Code added to read in variable length arrays in plasma structure
	15aug	ksl	Updated to read domain structure
	15oct	ksl	Updated to read disk and qdisk stuctures
*/


/**********************************************************/
/** 
 * @brief      Read back the windsavefile 
 *
 * @param [in] char  filename[]   The full name of the windsave file
 * @return     The number of successful reads, or -1 if the file cannot 
 * be opened
 *
 * @details
 * 
 * The routine reads in both the windsave file and the
 * associated atomic data files for a model. It also reads the
 * disk and qdisk structures.
 *
 *
 * ### Notes ###
 *
 * ### Programming Comment ### 
 * This routine calls wind_complete. This looks superfluous, since 
 * wind_complete and its subsidiary routines but it
 * also appears harmless.  ksl 
 *
 **********************************************************/

int
wind_read (filename)
     char filename[];
{
  FILE *fptr;
  int ndom;
  int n, m;
  char header[LINELENGTH];
  char version[LINELENGTH];
  struct stat file_stat;        // Used to check the atomic data exists

  if ((fptr = fopen (filename, "r")) == NULL)
  {
    return (-1);
  }

  n = fread (header, sizeof (header), 1, fptr);
  sscanf (header, "%*s %s", version);
  Log ("Reading Windfile %s created with sirocco version %s with sirocco version %s\n", filename, version, VERSION);

  /* Now read in the geo structure */

  n += fread (&geo, sizeof (geo), 1, fptr);

  /* Read the atomic data file.  This is necessary to do here in order to establish the 
   * values for the dimensionality of some of the variable length structures, associated 
   * with macro atoms, especially but likely to be a good idea ovrall
   */

  if (stat (geo.atomic_filename, &file_stat))
  {
    if (system ("Setup_Py_Dir"))
    {
      Error ("Unable to open %s or create link for atomic data\n", geo.atomic_filename);
      Exit (1);
    }
  }

  get_atomic_data (geo.atomic_filename);

/* Now allocate space for the wind array */

  NDIM2 = geo.ndim2;
  NPLASMA = geo.nplasma;

  n += fread (zdom, sizeof (domain_dummy), geo.ndomain, fptr);
  for (ndom = 0; ndom < geo.ndomain; ++ndom)
  {
    allocate_domain_wind_coords (ndom);
    n += fread (zdom[ndom].wind_x, sizeof (double), zdom[ndom].ndim, fptr);
    n += fread (zdom[ndom].wind_z, sizeof (double), zdom[ndom].mdim, fptr);
    n += fread (zdom[ndom].wind_midx, sizeof (double), zdom[ndom].ndim, fptr);
    n += fread (zdom[ndom].wind_midz, sizeof (double), zdom[ndom].mdim, fptr);
    if (zdom[ndom].coord_type == CYLVAR)
    {
      cylvar_allocate_domain (ndom);
      n += fread (zdom[ndom].wind_z_var, sizeof (double), zdom[ndom].ndim * zdom[ndom].mdim, fptr);
      n += fread (zdom[ndom].wind_midz_var, sizeof (double), zdom[ndom].ndim * zdom[ndom].mdim, fptr);
    }
  }

  calloc_wind (NDIM2);
  n += fread (wmain, sizeof (wind_dummy), NDIM2, fptr);

  /* Read the disk and qdisk structures */

  n += fread (&disk, sizeof (disk), 1, fptr);
  n += fread (&qdisk, sizeof (disk), 1, fptr);

  calloc_plasma (NPLASMA);

  n += fread (plasmamain, sizeof (plasma_dummy), NPLASMA, fptr);

  /*Allocate space for the dynamically allocated plasma arrays */

  calloc_dyn_plasma (NPLASMA);

  /* Read in the dynamically allocated plasma arrays */

  for (m = 0; m < NPLASMA; m++)
  {

    n += fread (plasmamain[m].density, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].partition, sizeof (double), nions, fptr);

    n += fread (plasmamain[m].ioniz, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].recomb, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].inner_recomb, sizeof (double), nions, fptr);

    n += fread (plasmamain[m].scatters, sizeof (int), nions, fptr);
    n += fread (plasmamain[m].xscatters, sizeof (double), nions, fptr);

    n += fread (plasmamain[m].heat_ion, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].cool_rr_ion, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].cool_dr_ion, sizeof (double), nions, fptr);
    n += fread (plasmamain[m].lum_rr_ion, sizeof (double), nions, fptr);

    n += fread (plasmamain[m].levden, sizeof (double), nlte_levels, fptr);
    n += fread (plasmamain[m].recomb_simple, sizeof (double), nphot_total, fptr);
    n += fread (plasmamain[m].recomb_simple_upweight, sizeof (double), nphot_total, fptr);
    n += fread (plasmamain[m].kbf_use, sizeof (double), nphot_total, fptr);
  }


  /*Allocate space for macro-atoms and read in the data */

  if (geo.nmacro > 0)
  {
    calloc_macro (NPLASMA);
    n += fread (macromain, sizeof (macro_dummy), NPLASMA, fptr);
    calloc_estimators (NPLASMA);
    calloc_matom_matrix (NPLASMA);

    for (m = 0; m < NPLASMA; m++)
    {
      n += fread (macromain[m].jbar, sizeof (double), size_Jbar_est, fptr);
      n += fread (macromain[m].jbar_old, sizeof (double), size_Jbar_est, fptr);
      n += fread (macromain[m].gamma, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].gamma_old, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].gamma_e, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].gamma_e_old, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].alpha_st, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].alpha_st_old, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].alpha_st_e, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].alpha_st_e_old, sizeof (double), size_gamma_est, fptr);
      n += fread (macromain[m].recomb_sp, sizeof (double), size_alpha_est, fptr);
      n += fread (macromain[m].recomb_sp_e, sizeof (double), size_alpha_est, fptr);
      n += fread (macromain[m].matom_emiss, sizeof (double), nlevels_macro, fptr);
      n += fread (macromain[m].matom_abs, sizeof (double), nlevels_macro, fptr);

      /* Force recalculation of kpkt_rates and matrix rates */

      macromain[m].kpkt_rates_known = FALSE;
      macromain[m].matrix_rates_known = FALSE;
    }

  }

  fclose (fptr);

  wind_complete ();

  Log ("Read geometry and wind structures from windsavefile %s\n", filename);

  return (n);

}




/**********************************************************/
/** 
 * @brief      A driver routine that calls coordinate-system specific routines
 * that complete the description of the wind
 *
 * @return     Always returns 0
 *
 * @details
 *
 * For the most point, the various routines that are called
 * just recalculate some of the various arrays used for 
 * finding the boundaries of an individual cell.
 *
 * These basically are just 1-d versions of the coordinate
 * grids at the edge and mid-points of each grid cell
 *
 * ### Notes ###
 *
 * The need to call this routine whenever a windsave file is read
 * in could be obviated if all of the variable length arrays
 * in domains were actually written out.  But since these
 * variables are easy to calculate again, it is done here.
 * 
 *
 **********************************************************/

void
wind_complete ()
{
  int ndom;

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    if (zdom[ndom].coord_type == SPHERICAL)
    {
      spherical_wind_complete (ndom, wmain);
    }
    else if (zdom[ndom].coord_type == CYLIND)
    {
      cylind_wind_complete (ndom, wmain);
    }
    else if (zdom[ndom].coord_type == RTHETA)
    {
      rtheta_wind_complete (ndom, wmain);
    }
    else if (zdom[ndom].coord_type == CYLVAR)
    {
      cylvar_wind_complete (ndom, wmain);
    }
    else
    {
      Error ("wind_complete: Don't know how to complete coord_type %d\n", zdom[ndom].coord_type);
      Exit (0);
    }
  }
}


/**********************************************************/
/** 
 * @brief      Save all the spectra to a binary file
 *
 * @param [in] char  filename[]   The file to write to 
 * @return     The number of successful writes
 *
 * @details
 *
 * ### Notes ###
 * 
 * The program exits if one cannot write the file
 *
 **********************************************************/

int
spec_save (filename)
     char filename[];
{

  FILE *fptr;
  char header[LINELENGTH];
  int count;
  int i;

  if ((fptr = fopen (filename, "w")) == NULL)
  {
    Error ("spec_save: Unable to open %s\n", filename);
    Exit (EXIT_FAILURE);
  }

  sprintf (header, "Version %s  nspectra %d NWAVE_IONIZ %d NWAVE_EXTRACT %d NWAVE_MAX %d\n", VERSION, nspectra, NWAVE_IONIZ,
           NWAVE_EXTRACT, NWAVE_MAX);

  count = (int) fwrite (header, sizeof (header), 1, fptr);
  count += (int) fwrite (xxspec, sizeof (spectrum_dummy), nspectra, fptr);

  for (i = 0; i < nspectra; ++i)
  {
    count += (int) fwrite (xxspec[i].f, sizeof (*xxspec[i].f), NWAVE_MAX, fptr);
    count += (int) fwrite (xxspec[i].lf, sizeof (*xxspec[i].lf), NWAVE_MAX, fptr);
    count += (int) fwrite (xxspec[i].f_wind, sizeof (*xxspec[i].f_wind), NWAVE_MAX, fptr);
    count += (int) fwrite (xxspec[i].lf_wind, sizeof (*xxspec[i].lf_wind), NWAVE_MAX, fptr);
  }

  fclose (fptr);

  return (count);
}



/**********************************************************/
/** 
 * @brief      Read the binary file containing all the spectra
 *
 * @param [in] char  filename[]   The name of the file
 * @return     The number of spectra that were read in
 *
 * @details
 *
 * The routine allocates space for the spectra and reads
 * then in.
 *
 * ### Notes ###
 *
 * The first line of the file contains the Python version
 * and the number of spectra to be read in
 * 
 * The program exits if the file does not exist
 *
 **********************************************************/

int
spec_read (filename)
     char filename[];
{
  FILE *fptr;
  int nhead, nwave_ioniz_check;
  int count;
  int i;
  char header[LINELENGTH];
  char version[LINELENGTH];

  if ((fptr = fopen (filename, "r")) == NULL)
  {
    Error ("spec_read: Unable to open %s\n", filename);
    Exit (1);
  }

  count = (int) fread (header, sizeof (header), 1, fptr);

  nhead = sscanf (header, "%*s %s %*s %d %*s %d %*s %d %*s %d", version, &nspectra, &nwave_ioniz_check, &NWAVE_EXTRACT, &NWAVE_MAX);
  if (nhead != 5)
  {
    Error ("Incorrect header format in %s\n", files.specsave);
    Exit (EXIT_FAILURE);
  }
  if (nwave_ioniz_check != (int) NWAVE_IONIZ)
  {
    Error ("The current NWAVE_IONIZ (%d) value is incompatible with the spec_save file which has NWAVE_IONIZ = %d\n", NWAVE_IONIZ,
           nwave_ioniz_check);
    Exit (EXIT_FAILURE);
  }

  Log ("Reading specfile %s with %d spectra and %d wavelength bins, created with sirocco version %s and currently using sirocco version %s\n",
       filename, nspectra, NWAVE_EXTRACT, version, VERSION);

  /* First allocate space */

  xxspec = calloc (nspectra, sizeof (spectrum_dummy));
  if (xxspec == NULL)
  {
    Error ("spectrum_init: Could not allocate memory for %d spectra\n", nspectra);
    Exit (EXIT_FAILURE);
  }
  count += (int) fread (xxspec, sizeof (spectrum_dummy), nspectra, fptr);

  /* Now read the rest of the file */

  spectrum_allocate (nspectra);
  for (i = 0; i < nspectra; ++i)
  {
    count += (int) fread (xxspec[i].f, sizeof (*xxspec[i].f), NWAVE_MAX, fptr);
    count += (int) fread (xxspec[i].lf, sizeof (*xxspec[i].lf), NWAVE_MAX, fptr);
    count += (int) fread (xxspec[i].f_wind, sizeof (*xxspec[i].f_wind), NWAVE_MAX, fptr);
    count += (int) fread (xxspec[i].lf_wind, sizeof (*xxspec[i].lf_wind), NWAVE_MAX, fptr);
  }

  fclose (fptr);

  Log ("Read spec structures from specfile %s\n", filename);

  return (count);

}
