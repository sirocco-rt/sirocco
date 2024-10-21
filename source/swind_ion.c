/***********************************************************/
/** @file  swind_ion.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief   This file contains various subroutines of swind.  
 * It is not part of sirocco!
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
 * @brief      Print informaition about a specific ion
 *
 * @param [in] WindPtr W   The entire wind          
 * @param [in] int element An elemeent
 * @param [in] int istate A particular ion          
 * @param [in] int iswitch  A intenger indicating what specifically
 * to display about the ion
 * @param [in] char rootname []   The rootname of the windsave file
 * @param [in] int ochoice A parameter indcating whether results should
 * be written to a file          
 * @return   Always returns 0  
 *
 * @details
 *	calculates and 
 *	displays information for a specific element and ionization state.  The information
 *	displayed depends upon the value of iswitch  
 *
 * The choices for iswitch are:
 *  * 0 ion fractions
 *  * 1 ion denisities
 *  * 2 number of scatters involving this ion 
 *  * 3 the scatterd flux for this ion
 *
 *
 * ### Notes ###
 *
 **********************************************************/
int
ion_summary (w, element, istate, iswitch, rootname, ochoice)
     WindPtr w;
     int element, istate;
     int iswitch;
     char rootname[];
     int ochoice;
{
  int nion, nelem;
  int n;
  char choice[LINELENGTH], iname[LINELENGTH];
  char name[LINELENGTH];
  char filename[LINELENGTH];
  double x;
  int nplasma;
  double nh;

  x = 0;

  nion = 0;
  while (nion < nions && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;
  if (nion == nions)
  {
    Log ("Error--element %d ion %d not found in define_wind\n", element, istate);
    return (-1);
  }
  nelem = 0;
  while (nelem < nelements && ele[nelem].z != element)
    nelem++;

  strcpy (name, "");


  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    nplasma = w[n].nplasma;
    if (w[n].inwind >= 0)
    {
      if (iswitch == 0)
      {
        sprintf (name, "Element %d (%s) ion %d fractions\n", element, ele[nelem].name, istate);
        aaa[n] = plasmamain[nplasma].density[nion];
        nh = rho2nh * plasmamain[nplasma].rho;
        aaa[n] /= (nh * ele[nelem].abun);
      }
      else if (iswitch == 1)
      {
        sprintf (name, "Element %d (%s) ion %d density\n", element, ele[nelem].name, istate);
        aaa[n] = plasmamain[nplasma].density[nion];
      }
      else if (iswitch == 2)
      {
        sprintf (name, "Element %d (%s) ion %d  #scatters\n", element, ele[nelem].name, istate);
        aaa[n] = plasmamain[nplasma].scatters[nion];
      }
      else if (iswitch == 3)
      {
        sprintf (name, "Element %d (%s) ion %d scattered flux\n", element, ele[nelem].name, istate);
        aaa[n] = plasmamain[nplasma].xscatters[nion];
      }
      else
      {
        Error ("ion_summary : Unknown switch %d \n", iswitch);
        exit (0);
      }
    }
  }

  display (name);




  /* Store the appropriate values in a place where it does not matter */
  if (ochoice)
  {
    for (n = 0; n < NDIM2; n++)
    {
      nplasma = w[n].nplasma;
      if (w[n].inwind >= 0)
      {
        if (iswitch == 0)
          x /= ((plasmamain[nplasma].density[0] + plasmamain[nplasma].density[1]) * ele[nelem].abun);
        else if (iswitch == 1)
        {
          x = plasmamain[nplasma].density[nion];
          x = log10 (x);
        }
        else if (iswitch == 2)
        {
          x = plasmamain[nplasma].scatters[nion];
        }
        else if (iswitch == 3)
        {
          x = plasmamain[nplasma].xscatters[nion];
        }
        else
        {
          Error ("ion_summary : Unknown switch %d \n", iswitch);
          exit (0);
        }


      }
      else
        x = 0.0;
      w[n].x[1] = x;

    }

    strcpy (filename, rootname);
    if (iswitch == 0)
      strcpy (choice, ".ion");
    else if (iswitch == 1)
      strcpy (choice, ".ionc");
    else if (iswitch == 2)
      strcpy (choice, ".ions");
    else if (iswitch == 3)
      strcpy (choice, ".iona");
    else
    {
      Error ("ion_summary : Unknown switch %d \n", iswitch);
      exit (0);
    }

    strcat (choice, ele[nelem].name);
    sprintf (iname, "%d", istate);
    strcat (choice, iname);

    strcat (filename, choice);
    write_array (filename, ochoice);
  }

  return (0);
}


/**********************************************************/
/** 
 * @brief      Calculates tau_ave for a particular transition 
 * of an element and ion assuming one is in resonance
 *
 * @param [in] WindPtr W   The entire wind          
 * @param [in] int element An elemeent
 * @param [in] double freq  Frequency of the line of interest
 * @param [in] int istate A particular ion          
 * @param [in] char rootname []   The rootname of the windsave file
 * @param [in] int ochoice A parameter indcating whether results should
 * be written to a file          
 * @return   Always returns 0  
 *
 * @details
 *
 *
 * ### Notes ###
 *
 * IMPORTANT - This routine does not search for the atomic data associated
 * with line, It assumes a nominal value.
 *
 **********************************************************/
int
tau_ave_summary (w, element, istate, freq, rootname, ochoice)
     WindPtr w;
     int element, istate;
     double freq;
     char rootname[];
     int ochoice;
{
  int nion, nelem;
  int n;
  char choice[LINELENGTH], iname[LINELENGTH];
  char name[LINELENGTH];
  char filename[LINELENGTH];
  int nplasma;



/* Find the CIV ion */
  nion = 0;
  while (nion < nions && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;
  if (nion == nions)
  {
    Log ("Error--element %d ion %d not found in define_wind\n", element, istate);
    return (-1);
  }
  nelem = 0;
  while (nelem < nelements && ele[nelem].z != element)
    nelem++;

  strcpy (name, "");
  sprintf (name, "Element %d (%s) ion %d fractions\n", element, ele[nelem].name, istate);

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = PI_E2_OVER_M * plasmamain[nplasma].density[nion] / freq / w[n].dvds_ave;
    }
  }

  display (name);

  /* Store the appropriate values in a place where it does not matter */
  if (ochoice)
  {
    for (n = 0; n < NDIM2; n++)
    {
      if (w[n].inwind >= 0)
      {
        nplasma = w[n].nplasma;
        w[n].x[1] = plasmamain[nplasma].density[nion] / (0.87 * plasmamain[nplasma].ne * ele[nelem].abun);
      }
      else
        w[n].x[1] = 0;
    }

    strcpy (filename, rootname);

    strcpy (choice, ".ion");
    strcat (choice, ele[nelem].name);
    sprintf (iname, "%d", istate);
    strcat (choice, iname);

    strcat (filename, choice);

    write_array (filename, ochoice);
  }

  return (0);
}



/**********************************************************/
/** 
 * @brief      Calculate the luminoisty for a specific line
 *
 * @param [in] WindPtr W   The entire wind          
 * @param [in] char rootname []   The rootname of the windsave file
 * @param [in] int ochoice A parameter indcating whether results should
 * be written to a file          
 * @return   Always returns 0  
 *
 * @details
 *  For a given line from a preset selection, outputs luminosity.
 *  Can cope with matom/non-matom lines.
 *
 *
 * File rootname.line[ELEM][ION].dat
 *
 *
 * ### Notes ###
 *
 *
 * History:
 *  31062016  SWM modified to add H-alpha, matoms
 *
 **********************************************************/
int
line_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int nion, nelem;
  int element, istate, iline, levu, levl, i_matom_search;
  int n;
  double x, lambda;
  char choice[LINELENGTH], iname[LINELENGTH];
  char name[LINELENGTH];
  char filename[LINELENGTH];
  int nline;
  double freq_search;

  double d1, d2, z, q, energy, rb, tot, omega;
  int nplasma;

  z = q = 0;
  iline = 0;
  lambda = 0;
  i_matom_search = 0;
  rdint ("line (0=C-IV, 1=Ha, 2=Hb, 3=Matom", &iline);
  switch (iline)
  {
  case 0:                      //Carbon-IV
    element = 6;
    istate = 4;
    lambda = 1548.1949e-8;
    break;
  case 1:                      //Hydrogen Alpha
    element = 1;
    istate = 1;
    lambda = 6562.7097e-8;
    levu = 3;
    levl = 2;
    break;
  case 2:                      //Hydrogen Beta
    element = 1;
    istate = 1;
    lambda = 4861.363e-8;
    levu = 4;
    levl = 2;
    break;
  case 3:                      //Generic matom
    i_matom_search = 1;
    element = 1;
    istate = 1;
    levu = 2;
    levl = 1;
    rdint ("Element", &element);
    rdint ("Ion", &istate);
    rdint ("Upper level", &levu);
    rdint ("Lower level", &levl);
    break;
  default:
    Error ("line_summary: Not a valid line.");
    exit (0);

  }

/* Convert wavelength to energy and frequency */
  freq_search = VLIGHT / lambda;
  energy = HC / lambda;

/* Find the ion */
  if (i_matom_search)
  {
    printf ("Searching for matom line...\n");
    nline = 0;
    while (nline < nlines && !(lin_ptr[nline]->z == element && lin_ptr[nline]->istate == istate
                               && lin_ptr[nline]->levu == levu && lin_ptr[nline]->levl == levl))
    {
      nline++;
    }
    if (nline == nlines)
    {
      Error ("line_summary: Could not find line in linelist\n");
      exit (0);
    }
    nelem = 0;
    while (nelem < nelements && ele[nelem].z != element)
      nelem++;
    if (nelem == nelements)
    {
      Log ("line_summary: Could not find element %d", element);
      return (-1);
    }
    nion = 0;
    while (nion < nions && !(ion[nion].z == element && ion[nion].istate == istate))
      nion++;
    if (nion == nions)
    {
      Log ("Error--element %d ion %d not found in define_wind\n", element, istate);
      return (-1);
    }
  }
  else
  {
    nion = 0;
    while (nion < nions && !(ion[nion].z == element && ion[nion].istate == istate))
      nion++;
    if (nion == nions)
    {
      Log ("Error--element %d ion %d not found in define_wind\n", element, istate);
      return (-1);
    }
    nelem = 0;
    while (nelem < nelements && ele[nelem].z != element)
      nelem++;
    if (nelem == nelements)
    {
      Log ("line_summary: Could not find element %d", element);
      return (-1);
    }
    nline = 0;
    freq_search = VLIGHT / lambda;

    while (fabs (1. - lin_ptr[nline]->freq / freq_search) > 0.0001 && nline < nlines)
      nline++;
    if (nline == nlines)
    {
      Error ("line_summary: Could not find line in linelist\n");
      exit (0);
    }
  }


  rdint ("line_transfer(0=pure.abs,1=pure.scat,2=sing.scat,3=escape.prob, 4=off, diagnostic)", &geo.line_mode);
  if (geo.line_mode == LINE_MODE_ABSORB)
    Log ("Pure_abs in line heating/cooling\n");
  else if (geo.line_mode == LINE_MODE_SCAT)
    Log ("Pure_scat_in line heating/cooling\n");
  else if (geo.line_mode == LINE_MODE_SINGLE_SCAT)
    Log ("Single scat for line heating/cooling\n");
  else if (geo.line_mode == LINE_MODE_ESC_PROB)
    Log ("Escape probabilities for line heating/cooling\n");
  else if (geo.line_mode == 4)
    Log ("No line transfer; diagnostic mode only\n");
  else
  {
    Log ("Unknown line mode\n");
    return (0);
  }

  strcpy (name, "");
  if (lin_ptr[nline]->macro_info == TRUE)
  {
    sprintf (name, "%d Luminosity %d (%s) ion %d fractions\n", nline, element, ele[nelem].name, istate);
  }
  else
  {
    sprintf (name, "%d Luminosity %d (%s) ion %d matom %d-%d fractions\n", nline, element, ele[nelem].name, istate,
             lin_ptr[nline]->levu, lin_ptr[nline]->levl);
  }

  tot = 0.0;
  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0.0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;

      if (lin_ptr[nline]->macro_info == TRUE)
      {
        /* If this is a matom line */
        d2 = den_config (&plasmamain[nplasma], lin_ptr[nline]->nconfigu);
        d1 = den_config (&plasmamain[nplasma], lin_ptr[nline]->nconfigl);
      }
      else
      {
        /* If this is not a matom line */
        two_level_atom (lin_ptr[nline], &plasmamain[nplasma], &d1, &d2);
      }


      if (geo.line_mode == 4)
      {
        /* diagnostic mode that just returns n2 * A21 * vol * hnu */
        x = (d2) * a21 (lin_ptr[nline]) * PLANCK * lin_ptr[nline]->freq * w[n].vol;
      }
      else
      {
        /* the below code is essentially duplicated from lum_lines() in lines.c, see #643 */
        x = lin_ptr[nline]->gu / lin_ptr[nline]->gl * d1 - d2;
        z = exp (-H_OVER_K * lin_ptr[nline]->freq / plasmamain[nplasma].t_e);
        q = 1. - scattering_fraction (lin_ptr[nline], &plasmamain[nplasma]);
        x *= q * a21 (lin_ptr[nline]) * z / (1. - z);
        x *= PLANCK * lin_ptr[nline]->freq * plasmamain[nplasma].vol;

        /* Include effects of line trapping */
        if (geo.line_mode == LINE_MODE_ESC_PROB)
          x *= p_escape (lin_ptr[nline], &plasmamain[nplasma]);
      }

      tot += x;
      aaa[n] = x;
    }
  }

  display (name);

  tot = 2. * tot;               // Why is there a factor of 2 here??? ksl

  Log ("The total %s ion %d luminosity (flux) is %8.2g (%8.2g)\n", ele[nelem].name, istate, tot, tot / (4 * PI * 1e4 * PC * PC));

  /* Store the appropriate values in a place where it does not matter */
  if (ochoice)
  {
    for (n = 0; n < NDIM2; n++)
    {
      // Here is the calculation of the effective collisions strength
      if (w[n].inwind >= 0)
      {
        nplasma = w[n].nplasma;
        omega = 5.13 * pow (plasmamain[nplasma].t_e / 1.e5, 0.18);
        rb = 8.629e-6 * exp (-energy / (BOLTZMANN * plasmamain[nplasma].t_e)) / sqrt (plasmamain[nplasma].t_e) * omega;
        w[n].x[1] = plasmamain[nplasma].density[nion] * plasmamain[nplasma].ne * rb * energy * w[n].vol;
      }
      else
        w[n].x[1] = 0;
    }


    strcpy (filename, rootname);
    strcpy (choice, ".line");
    strcat (choice, ele[nelem].name);
    sprintf (iname, "%d", istate);
    strcat (choice, iname);
    if (lin_ptr[nline]->macro_info == TRUE)
    {
      sprintf (iname, ".%d-%d", lin_ptr[nline]->levu, lin_ptr[nline]->levl);
      strcat (choice, iname);
    }
    strcat (filename, choice);
    write_array (filename, ochoice);
  }

  return (0);
}



/**********************************************************/
/** 
 * @brief      Get the total emission of the cells
 *
 * @param [in] WindPtrYW   The entire wind          
 * @param [in] char rootname []   The rootname of the windsave file
 * @param [in] int ochoice A parameter indcating whether results should
 * be written to a file          
 * @return   Always returns 0  
 *
 * @details
 *
 * Calls total_emission
 *
 *
 * ### Notes ###
 *
 **********************************************************/
int
total_emission_summary (rootname, ochoice)
     char rootname[];
     int ochoice;
{
  double tot;
  int n;
  char filename[LINELENGTH];


  tot = 0.0;
  for (n = 0; n < NPLASMA; n++)
  {
    aaa[n] = 0;
    tot += aaa[n] = total_emission (&plasmamain[n], 0.0, VERY_BIG);
  }

  display ("Calculated thermal luminosities of cells");

  Log ("The total luminosity due to thermal process is %8.2g\n", tot);



  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".abs");
    write_array (filename, ochoice);
  }


  return (0);
}


/**********************************************************/
/**    
 * @brief      Find the electron temperature where heating matches
 * the cooling
 *
 * @param [in] WindPtr   The entire wind          
 * @param [in] char rootname []   The rootname of the windsave file
 * @param [in] int ochoice A parameter indcating whether results should
 * be written to a file          
 * @return   Always returns 0  
 *
 * @details
 *
 *
 * ### Notes ###
 *
 **********************************************************/
int
modify_te (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  double x;
  double t_e, calc_te ();
  char filename[LINELENGTH];
  int nplasma;


  for (n = 0; n < NDIM2; n++)
  {
    nplasma = w[n].nplasma;
    aaa[n] = 0;
    if (w[n].inwind >= 0 && (x = plasmamain[nplasma].heat_tot) > 1.0)
    {
      aaa[n] = t_e = calc_te (&plasmamain[nplasma], MIN_TEMP, 1.2 * plasmamain[nplasma].t_r);
    }
  }

  display ("Calculate t_e");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".te_mod");
    write_array (filename, ochoice);
  }

  return (0);
}



/**********************************************************/
/** 
 * @brief      Write out the emission measure for a specific
 * ionization state of an element
 *
 * @param [in] WindPtr W   The entire wind          
 * @param [in] int element The z for an element     
 * @param [in] int istate  The ion state of an element
 * @param [in] char rootname []   The rootname of the windsave file
 * @param [in] int ochoice A parameter indcating whether results should
 * be written to a file          
 * @return   Always returns 0  
 *
 * @details
 *
 *
 * ### Notes ###
 *
 **********************************************************/

int
partial_measure_summary (w, element, istate, rootname, ochoice)
     WindPtr w;
     int element, istate;
     char rootname[];
     int ochoice;
{
  int nion, nelem;
  int n;
  double total;
  char choice[LINELENGTH], iname[LINELENGTH];
  char name[LINELENGTH];
  char filename[LINELENGTH];
  int nplasma;



  total = 0;
  nion = 0;
  while (nion < nions && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;
  if (nion == nions)
  {
    Log ("Error--element %d ion %d not found in define_wind\n", element, istate);
    return (-1);
  }
  nelem = 0;
  while (nelem < nelements && ele[nelem].z != element)
    nelem++;

  strcpy (name, "");
  sprintf (name, "Element %d (%s) ion %d fractions", element, ele[nelem].name, istate);

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    nplasma = w[n].nplasma;
    if (plasmamain[nplasma].ne > 1.0 && w[n].inwind >= 0)
    {
      total += aaa[n] = plasmamain[nplasma].density[nion] * plasmamain[nplasma].ne * w[n].vol;
    }
  }

  display ("Calculate t_e");

  Log ("The partial emission measure ne*nion*V is %8.2g\n", total);

  /* Store the appropriate values in a place where it does not matter */
  if (ochoice)
  {
    for (n = 0; n < NDIM2; n++)
    {
      nplasma = w[n].nplasma;
      if (plasmamain[nplasma].ne > 1.0 && w[n].inwind >= 0)
        w[n].x[1] = plasmamain[nplasma].density[nion] / (0.87 * plasmamain[nplasma].ne * ele[nelem].abun);
      else
        w[n].x[1] = 0;
    }

    strcpy (filename, rootname);
    strcpy (choice, ".part");
    strcat (choice, ele[nelem].name);
    sprintf (iname, "%d", istate);
    strcat (choice, iname);

    strcat (filename, choice);
    write_array (filename, ochoice);
  }

  return (0);
}


/**********************************************************/
/** 
 * @brief      Write out the collision strengths (and related information)
 * for all of the lines 
 *
 * @param [in] WindPtr W   The entire wind          
 * @param [in] char rootname []   The rootname of the windsave file
 * @param [in] int ochoice A parameter indcating whether results should
 * be written to a file          
 * @return   Always returns 0  
 *
 * @details
 *
 *
 * ### Notes ###
 *
 **********************************************************/

int
collision_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int nline, int_te;
  double t_e, qup, qdown, A, wavelength;
  char filename[LINELENGTH], suffix[LINELENGTH];
  FILE *fopen (), *fptr = NULL;

  t_e = 10000.0;

  /* Input from user to request temperature */
  rddoub ("electron temperature for calculation:", &t_e);
  int_te = (int) t_e;           // for file label.

  if (ochoice)
  {
    /* open filename root.coll.dat */
    strcpy (filename, rootname);
    sprintf (suffix, ".t%d.coll.dat", int_te);
    strcat (filename, suffix);
    fptr = fopen (filename, "w");

    Log ("\nWriting collision strengths to file %s...\n\n", filename);

    fprintf (fptr, "# Collision strengths at electron temperature %.1fK\n", t_e);
    fprintf (fptr, "# For atomic data file %s\n", geo.atomic_filename);
    fprintf (fptr, "line wavelength z istate levu levl q12 q21 a21 macro_info\n");
  }
  else
  {
    Log ("Collision strengths at electron temperature %.1fK\n", t_e);
    Log ("line wavelength z istate levu levl q12 q21 a21 macro_info\n");
  }

  nline = 0;

  while (nline < nlines)
  {
    wavelength = VLIGHT / lin_ptr[nline]->freq / ANGSTROM;

    qup = q12 (lin_ptr[nline], t_e);
    qdown = q21 (lin_ptr[nline], t_e);
    A = a21 (lin_ptr[nline]);

    if (ochoice)
    {
      fprintf (fptr, "%d %8.4e %d %d %d %d %8.4e %8.4e %8.4e %d\n",
               nline, wavelength, lin_ptr[nline]->z,
               lin_ptr[nline]->istate, lin_ptr[nline]->levu, lin_ptr[nline]->levl, qup, qdown, A, lin_ptr[nline]->macro_info);
    }
    else
    {
      Log ("%d %8.4e %d %d %d %d %8.4e %8.4e %8.4e %d\n",
           nline, wavelength, lin_ptr[nline]->z,
           lin_ptr[nline]->istate, lin_ptr[nline]->levu, lin_ptr[nline]->levl, qup, qdown, A, lin_ptr[nline]->macro_info);
    }
    nline++;
  }

  if (ochoice)
    fclose (fptr);

  return (0);
}
