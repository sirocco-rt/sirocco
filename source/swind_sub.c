
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	swind is a program which can be used to display various parameters of a wind
	as calculated by sirocco.  This file contains many of the subroutiens
	used

Arguments:


Returns:

Description:
	
		
Notes:

History:
 	97jun	ksl	Coding on swind began.
	05jul	ksl	56d -- Moved all of the subroutines to swind_sub so that it was
			easier to make templates.

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**************************************************************************


  Synopsis:  

	This routine controls what is displayed on the screen.  
	

  Description:	

	There are basically two options, determined by the variable "determine"
	If determine==1, then zoom sets up the variables swind_min, swind_max and
	swind_delta so that the the entire wind is displayed, but it is
	subsampled in the x direction.
 
	If determine!=1, then a section of the wind is displayed at full resolution

  Arguments:  

  Returns:

  Notes:

  History:
 	07jul	ksl	Made modifications to allow for the possibility that
 			the wind has fewer than 10 elemetns in the x direction

 ************************************************************************/

int
zoom (direction)
     int direction;
{
  int center;
  int ndim;


  ndim = zdom[current_domain].ndim;

  if (direction == 1)
  {                             /* then unzoom */
    Log ("Showing selected positions throughout wind\n");
    swind_min = 0;
    swind_max = ndim;
    swind_delta = ndim / 10;
    /*
     * Allow for the possibility that the wind has an xdimension
     * < 10
     */
    if (swind_delta < 1)
      swind_delta = 1;
  }
  else
  {
    Log ("Select part of wind to display\n");
    center = 5;
    rdint ("Center_x", &center);
    swind_min = center - 5;
    if (swind_min < 0)
    {
      Log ("zoom: this choice of center is needlessly close to origin, adjusting swind_min to 0\n");
      swind_min = 0;
    }
    swind_max = swind_min + 10;
    if (swind_max > ndim)
    {
      swind_max = ndim;
      Log ("zoom: this choice of swind_max is lager than NDIM, adusting swind_max to NDIM");
    }
    swind_delta = 1;
  }
  return (0);
}






/**************************************************************************


  Synopsis:  

	overview

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/
int
overview (w, rootname)
     WindPtr w;
     char rootname[];
{
  int n;
  double heating, lines, ff, photo;

  heating = lines = ff = photo = 0.0;

  for (n = 0; n < NPLASMA; n++)
  {
    heating += plasmamain[n].heat_tot;
    lines += plasmamain[n].heat_lines;
    photo += plasmamain[n].heat_photo;
    ff += plasmamain[n].heat_ff;
  }
  Log (" Total cooling  %8.2e heating %8.2e\n", geo.cool_tot_ioniz, heating);
  Log (" Total emission %8.2e heating %8.2e\n", geo.lum_tot_ioniz, heating);
  Log ("    ff emission %8.2e heating %8.2e\n", geo.lum_ff_ioniz, ff);
  Log ("    fb emission %8.2e heating %8.2e\n", geo.cool_rr_ioniz, photo);
  Log ("  line emission %8.2e heating %8.2e\n", geo.lum_lines_ioniz, lines);
  return (0);
}





/**************************************************************************


  Synopsis:  
	Summary of everything at a given position 

  Description:	

  Arguments:  

  Returns:

  Notes:

  This has been modified to work with domains, but it is not obvious that
  this is what one wants, because the position is fixed

  History:

 ************************************************************************/
int
position_summary (w)
     WindPtr w;
{
  double x[3], v[3];
  struct photon p;
  int n;
  int nplasma;
  int inwind, ndom;

  x[0] = geo.rmax / 5;
  x[1] = 0.0;
  x[2] = x[0];

a:Log ("Input x=0,y=0,z=0 to return to main routine\n");
  rddoub ("x", &x[0]);
  rddoub ("y", &x[1]);
  rddoub ("z", &x[2]);

  if (length (x) == 0.0)
    return (0);

  inwind = where_in_wind (x, &ndom);
  if (inwind != W_ALL_INWIND)
  {
    Log ("Position %8.2e  %8.2e %8.2e is not in an active region of grid %d %d\n", x[0], x[1], x[2], inwind, ndom);
    ndom = 0;
  }

  n = where_in_grid (ndom, x);
  nplasma = wmain[n].nplasma;

  Log ("Position %8.2e  %8.2e %8.2e  Cell %5d\n", x[0], x[1], x[2], n);

  Log ("Vertex position %8.2e  %8.2e %8.2e  rtheta %8.2e %8.2e \n", w[n].x[0], w[n].x[1], w[n].x[2], w[n].r, w[n].theta);
  Log ("Center position %8.2e  %8.2e %8.2e  rtheta %8.2e %8.2e \n", w[n].xcen[0], w[n].xcen[1], w[n].xcen[2], w[n].rcen, w[n].thetacen);

  Log ("Electron density: %8.2e  rho %8.2e\n", plasmamain[nplasma].ne, plasmamain[nplasma].rho);
  Log ("Vel cell: %8.2e  %8.2e %8.2e\n", w[n].v[0], w[n].v[1], w[n].v[2]);

  p.x[0] = x[0];
  p.x[1] = x[1];
  p.x[2] = x[2];

  vwind_xyz (ndom, &p, v);
  Log ("Velocity: %8.2e  %8.2e %8.2e\n", v[0], v[1], v[2]);


  goto a;

  return (0);
}




/**************************************************************************


  Synopsis:  
	A summary of the energy absorbed in a cell 

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

	10nov	ksl	What was here previously was bizarre as the filename was
			being added to all the time and this caused ultimately
			a segmentation problem


 ************************************************************************/
int
abs_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  double x, xtot;
  char c;
  char name[LINELENGTH];
  char filename[LINELENGTH];
  int nplasma;

  c = 't';
  xtot = 0.0;

  printf ("Absorption tot=t, lines=l,f=ff,b=fb,h=hydrogen,i=he1,j=he2,z=heavy elements\n");
  rdchar ("Choice", &c);


  strcpy (filename, rootname);

  switch (c)
  {
  case 't':                    /* Total  Absorption */
    strcpy (name, "Total Absorbtion");
    strcat (filename, "heat_tot");
    break;
  case 'f':                    /* ff */
    strcpy (name, "Free Free Absorbtion");
    strcat (filename, "heat_ff");
    break;
  case 'b':                    /* Photoionization */
    strcpy (name, "Total photoionization Heating");
    strcat (filename, "heat_photo");
    break;
  case 'l':                    /* Line  heating */
    strcpy (name, "Resonance Line Heating");
    strcat (filename, "heat_lines");
    break;
  case 'h':                    /* H photoionization */
    strcpy (name, "H  Photoionization Heating");
    strcat (filename, "heat_h");
    break;
  case 'i':                    /* He1 photo */
    strcpy (name, "He I Photoionization Heating");
    strcat (filename, "heat_he1");
    break;
  case 'j':                    /* He2 photo */
    strcpy (name, "He 2 Photoionization Heating");
    strcat (filename, "heat_he2");
    break;
  case 'z':                    /* Metal photo */
    strcpy (name, "Metal Photoionization Heating");
    strcat (filename, "heat_z");
    break;
  default:
    printf ("Not a valid choice\n");
    return (0);
  }



  for (n = 0; n < NDIM2; n++)
  {
    x = 0.0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      switch (c)
      {
      case 't':
        {                       /* Total heating */
          x = plasmamain[nplasma].heat_tot;
          break;
      case 'f':                /* ff heating */
          x = plasmamain[nplasma].heat_ff;
          break;
      case 'b':                /* photoionization heating */
          x = plasmamain[nplasma].heat_photo;
          break;
      case 'l':                /* Line heating */
          x = plasmamain[nplasma].heat_lines;
          break;
      case 'h':                /* H heating */
          x = plasmamain[nplasma].heat_ion[0];
          break;
      case 'i':                /* He1 heating */
          x = plasmamain[nplasma].heat_ion[2];
          break;
      case 'j':                /* He2 heating */
          x = plasmamain[nplasma].heat_ion[3];
          break;
      case 'z':                /* Line heating of high z elements */
          x = plasmamain[nplasma].heat_z;
          break;
      default:
          printf ("Not a valid choice\n");
          return (0);
        }
        xtot += x;

      }
      aaa[n] = x;
    }
  }

  display (name);
  /* Log ("Component heating %8.3e\n", xtot * swind_delta);
     JM130624 v76: Removed swind_delta term as gives too high heating */
  Log ("Component heating %8.3e\n", xtot);

  if (ochoice)
  {
    write_array (filename, ochoice);
  }
  return (0);

}

/**************************************************************************


  Synopsis:  
  A summary of shock heating

  Description:  

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/

int
shock_heating_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  double tot;
  double shock_heating ();
  char filename[LINELENGTH];



  tot = 0.0;
  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      tot += aaa[n] = shock_heating (&w[n]);
    }
  }

  if (geo.nonthermal == 1)
  {
    display ("Shock heating");
    Log ("The total shock heating is %8.2g\n", tot);
  }
  else
  {
    display ("This is only potential shock heating - it was switched off in the model");
    Log ("The total shock heating is %8.2g\n", tot);
  }

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".shock_heating");
    write_array (filename, ochoice);
  }
  return (0);

}





/**************************************************************************


  Synopsis:  
	A summary of adiabatic cooling

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/

int
adiabatic_cooling_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  double tot;
  double adiabatic_cooling ();
  char filename[LINELENGTH];
  double t_e;



  tot = 0.0;
  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      t_e = plasmamain[w[n].nplasma].t_e;
      // ksl - I could not determine what the next line was supposed to do
      // num_recomb (&plasmamain[w[n].nplasma], t_e);
      tot += aaa[n] = adiabatic_cooling (&w[n], t_e);
    }
  }

  if (geo.adiabatic == 1)
  {
    display ("Adiabatic cooling");
    Log ("The total adiabatic cooling is %8.2g\n", tot);
  }
  else
  {
    display ("This is only potential adiabatic cooling - it was switched off in the model");
    Log ("The total adiabatic cooling is %8.2g\n", tot);
  }

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".adiabatic");
    write_array (filename, ochoice);
  }
  return (0);

}





/**************************************************************************


  Synopsis:  
	summary of the lum of a cell.  The routine
	simply reads variables that are already contained
	in the plasma structure


  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/
int
lum_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  double x, xtot;
  char c;
  char name[LINELENGTH];
  char filename[LINELENGTH];
  int nplasma;

  c = 't';
  strcpy (filename, rootname);

  printf ("lum tot=t,r=rad, lines=l,f=ff,b=fb,h=hydrogen,i=he1,j=he2,z=heavy elements\n");
  rdchar ("Choice", &c);
  switch (c)
  {
  case 't':                    /* Total luminosity */
    strcpy (name, "Total Luminosity");
    strcat (filename, ".tot_lum");
    break;
  case 'r':                    /* Radiative energo loss total */
    strcpy (name, "Total Radiative Luminosity");
    strcat (filename, ".rad_lum");
    break;
  case 'f':                    /* Radiative energo loss total */
    strcpy (name, "Free Free Luminosity");
    strcat (filename, ".ff_lum");
    break;
  case 'b':                    /* Radiative energo loss total */
    strcpy (name, "Free Bound (Total Recombination) Luminosity");
    strcat (filename, ".fb_lum");
    break;
  case 'l':                    /* Line lumosity */
    strcpy (name, "Line Luminosity");
    strcat (filename, ".line_lum");
    break;
  case 'h':                    /* H lumosity */
    strcpy (name, "H  Recombination Luminosity");
    strcat (filename, ".h1_recomb_lum");
    break;
  case 'i':                    /* Line lumosity */
    strcpy (name, "He I Recombination Luminosity");
    strcat (filename, ".he1_recom_lum");
    break;
  case 'j':                    /* Line lumosity */
    strcpy (name, "He 2 Recombination Luminosoity");
    strcat (filename, ".he2_recomb_lum");
    break;
  case 'z':                    /* Line lumosity */
    strcpy (name, "Metal Recombination Luminosity");
    strcat (filename, ".z_recomb_lum");
    break;
  default:
    printf ("Not a valid choice\n");
    return (0);
  }



  xtot = 0.0;


  for (n = 0; n < NDIM2; n++)
  {
    x = 0.0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      switch (c)
      {
      case 't':                /* Total luminosity */
        x = plasmamain[nplasma].lum_tot_ioniz;
        break;
      case 'r':                /* Radiative energo loss total */
        x = plasmamain[nplasma].lum_tot_ioniz;
        break;
      case 'f':                /* Radiative energo loss total */
        x = plasmamain[nplasma].lum_ff_ioniz;
        break;
      case 'b':                /* Radiative energo loss total */
        x = plasmamain[nplasma].cool_rr_ioniz;
        break;
      case 'l':                /* Line luminosity */
        x = plasmamain[nplasma].lum_lines_ioniz;
        break;
      case 'h':                /* H luminosity */
        x = plasmamain[nplasma].cool_rr_ion[0];
        break;
      case 'i':                /* Line luminosity */
        x = plasmamain[nplasma].cool_rr_ion[2];
        break;
      case 'j':                /* Line luminosity */
        x = plasmamain[nplasma].cool_rr_ion[3];
        break;
      case 'z':                /* Line luminosity */
        x = plasmamain[nplasma].cool_rr_metals_ioniz;
        break;
      default:
        printf ("Not a valid choice\n");
        return (0);
      }
      xtot += x;

    }
    aaa[n] = x;
  }

  display (name);

  /* Log ("Luminosity total %8.3e\n", xtot * swind_delta);
     JM130624 v76: Removed swind_delta term as gives too high heating */
  Log ("Luminosity total %8.3e\n", xtot);

  if (ochoice)
  {
    write_array (filename, ochoice);
  }
  return (0);

}





/**************************************************************************


  Synopsis:  
	A summary of the number of photoionizations in  a cell

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/
int
photo_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n, ion;
  char filename[LINELENGTH];
  int nplasma;

  ion = 0;

  rdint ("ion(<5)", &ion);
  if (!(0 <= ion && ion < NIONIZ))
  {
    Log ("ion out of range\n");
    return (0);
  }
  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].ioniz[ion] * plasmamain[nplasma].density[ion];
    }
  }
  display ("No of ionizations per second in cell");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".photo_sum");
  }
  return (0);

}





/**************************************************************************


  Synopsis:  
	A summary of the number of recombinations in  a cell

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/
int
recomb_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int ion;
  char filename[LINELENGTH];
  int nplasma;

  ion = 0;
  rdint ("ion(<5)", &ion);
  if (!(0 <= ion && ion < NIONIZ))
  {
    Log ("recomb_summary: ioniz out of range\n");
    return (0);
  }
  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      num_recomb (&plasmamain[nplasma], plasmamain[nplasma].t_e, 1);
      aaa[n] = plasmamain[nplasma].recomb[ion];
    }
  }

  display (" Number of recombinations per second in cell");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".recomb");
    write_array (filename, ochoice);
  }
  return (0);

}






/**************************************************************************


  Synopsis:  
	A summary of electron densities

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/

int
electron_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;


  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0.0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].ne;
    }
  }
  display ("Electron densities");
  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".ne");
    write_array (filename, ochoice);
  }
  return (0);

}






/**************************************************************************


  Synopsis:  
A summary of rho 

  Description:	

  Arguments:  

  Returns:

  Notes:
111002	ksl	Added to try to diagnose what was going 
		on with the torus

  History:

 ************************************************************************/
int
rho_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;


  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0.0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].rho;
    }
  }
  display ("Rho (gm/cm**3)");
  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".rho");
    write_array (filename, ochoice);
  }
  return (0);

}






/**************************************************************************


  Synopsis:  
	A summary of rho 

  Description:	

  Arguments:  

  Returns:

  Notes:
	Note that because of limitations in the way that display 
	works cell numbers greater than 99 are not displayed as
	integers unfortunately


  History:
	111002	ksl	Added to try to diagnose what was going 
			on with the torus

 ************************************************************************/
int
plasma_cell (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;


  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0.0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = nplasma;
    }
  }
  display ("Plasma cell number");
  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".pnum");
    write_array (filename, ochoice);
  }
  return (0);

}








/**************************************************************************


  Synopsis:  
	A summary of the average frequency

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/
int
freq_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].ave_freq;
    }
  }
  display ("Average freqency");
  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".ave_freq");
    write_array (filename, ochoice);

  }

  return (0);

}






/**************************************************************************


  Synopsis:  
	A summary of the number of photons which passed through a cell.

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:
	111002	ksl	Modified to be able to display photons of
 			various types

 ************************************************************************/
int
nphot_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;
  int ichoice;
  char string[LINELENGTH];



  ichoice = 0;
  while (rdint ("nphot(all=0,star=1,bl=2,disk=3,wind=4,agn=5,other=return)", &ichoice) != EOF)
  {
    if (ichoice < 0 || ichoice > 5)
    {
      return (0);
    }

    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {
        nplasma = w[n].nplasma;
        if (ichoice == 0)
        {
          aaa[n] = plasmamain[nplasma].ntot;
          strcpy (string, "Nphot tot per cell");
        }
        else if (ichoice == 1)
        {
          aaa[n] = plasmamain[nplasma].ntot_star;
          strcpy (string, "Nphot star per cell");
        }
        else if (ichoice == 2)
        {
          aaa[n] = plasmamain[nplasma].ntot_bl;
          strcpy (string, "Nphot bl per cell");
        }
        else if (ichoice == 3)
        {
          aaa[n] = plasmamain[nplasma].ntot_disk;
          strcpy (string, "Nphot disk per cell");
        }
        else if (ichoice == 4)
        {
          aaa[n] = plasmamain[nplasma].ntot_wind;
          strcpy (string, "Nphot wind per cell");
        }
        else if (ichoice == 5)
        {
          aaa[n] = plasmamain[nplasma].ntot_agn;
          strcpy (string, "Nphot agn per cell");
        }
        else
        {
          Error ("Unknown choice, try again\n");
        }
      }
    }
    display (string);

    if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".nphot");
      write_array (filename, ochoice);

    }
  }
  return (0);




}





/**************************************************************************


  Synopsis:  
	A summary of the electron temperatures

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/
int
temp_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].t_e;
    }
  }
  display ("T_e");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".te");
    write_array (filename, ochoice);
  }
  return (0);

}



/**************************************************************************


  Synopsis:  
	A summary of the radiation temperatures

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/
int
temp_rad (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;


  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].t_r;
    }
  }
  display ("T_rad");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".tr");
    write_array (filename, ochoice);

  }
  return (0);

}

/**************************************************************************


  Synopsis:  
	A summary of the radiative weights

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/

int
weight_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;
  char filename[LINELENGTH];

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].w;
    }
  }
  display ("Radiative weights");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".w");
    write_array (filename, ochoice);

  }

  return (0);

}





/**************************************************************************


  Synopsis:  
	A summary of the velocities

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/
int
velocity_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  double x;
  char name[LINELENGTH];
  char filename[LINELENGTH];
  int ichoice;

  ichoice = 0;

  /* Note that EOF actually causes an exit via rdpar */
  while (rdint ("|v|=0;v_x=1,v_y=2,v_z=3,return=other", &ichoice) != EOF)
  {
    if (ichoice < 0 || ichoice > 3)
      return (0);

    if (ichoice == 0)
      strcpy (name, "Velocity");
    else if (ichoice == 1)
      strcpy (name, "V_x");
    else if (ichoice == 2)
      strcpy (name, "V_y");
    else
      strcpy (name, "V_z");

    for (n = 0; n < NDIM2; n++)
    {
      x = 0;
      if (w[n].inwind >= 0)
      {
        if (ichoice == 0)
          x = sqrt (w[n].v[0] * w[n].v[0] + w[n].v[1] * w[n].v[1] + w[n].v[2] * w[n].v[2]);
        else if (ichoice == 1)
          x = w[n].v[0];
        else if (ichoice == 2)
          x = w[n].v[1];
        else
          x = w[n].v[2];
      }
      aaa[n] = x;
    }
    display (name);


    if (ochoice)
    {
      strcpy (filename, rootname);

      if (ichoice == 0)
      {
        strcat (filename, ".vel");
        write_array (filename, ochoice);
      }
      if (ichoice == 1)
      {
        strcat (filename, ".vx");
        write_array (filename, ochoice);
      }
      if (ichoice == 2)
      {
        strcat (filename, ".vy");
        write_array (filename, ochoice);
      }
      if (ichoice == 3)
      {
        strcat (filename, ".vz");
        write_array (filename, ochoice);
      }
    }
  }

  return (0);

}






/**************************************************************************


  Synopsis:  
	A summary of the velocities

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/
int
mo_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int ichoice;
  char name[LINELENGTH];
  char filename[LINELENGTH];
  double x;
  PlasmaPtr xplasma;


  ichoice = 0;

  /* Note that EOF actually causes an exit via rdpar */
  while (rdint ("|F_rad|=0;Frac_x=1,Frad_y=2,Frad_z=3,return=other", &ichoice) != EOF)
  {
    if (ichoice < 0 || ichoice > 4)
      return (0);

    if (ichoice == 0)
      strcpy (name, "|F_rad|");
    else if (ichoice == 1)
      strcpy (name, "F_rad_x");
    else if (ichoice == 2)
      strcpy (name, "F_rad_y");
    else
      strcpy (name, "F_rad_z");


    for (n = 0; n < NDIM2; n++)
    {
      xplasma = &plasmamain[w[n].nplasma];
      x = 0;
      if (w[n].inwind >= 0)
      {
        if (ichoice == 0)
          x =
            sqrt (xplasma->dmo_dt[0] * xplasma->dmo_dt[0] +
                  xplasma->dmo_dt[1] * xplasma->dmo_dt[1] + xplasma->dmo_dt[2] * xplasma->dmo_dt[2]);
        else if (ichoice == 1)
          x = xplasma->rad_force_es[0];
        else if (ichoice == 2)
          x = xplasma->rad_force_es[1];
        else
          x = xplasma->rad_force_es[2];
      }
      aaa[n] = x;
    }
    display (name);

    if (ochoice)
    {
      if (ichoice == 0)
      {
        strcpy (filename, rootname);
        strcat (filename, ".f_rad_mod");
        write_array (filename, ochoice);
      }
      else if (ichoice == 1)
      {
        strcpy (filename, rootname);
        strcat (filename, ".f_rad_x");
        write_array (filename, ochoice);
      }
      else if (ichoice == 2)
      {
        strcpy (filename, rootname);
        strcat (filename, ".f_rad_y");
        write_array (filename, ochoice);
      }
      else
      {
        strcpy (filename, rootname);
        strcat (filename, ".f_rad_z");
        write_array (filename, ochoice);
      }
    }

  }
  return (0);

}



/**************************************************************************


  Synopsis:  
	A summary of the volumes of each cell 


  Description:	

  Arguments:  

  Returns:

  Notes:

  History:
	080811	ksl	Add lines from Stuart's version of this
			routine to bring this version into 
			compliance with it.



 ************************************************************************/
int
vol_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;

  char filename[LINELENGTH];

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      aaa[n] = w[n].vol;
    }
  }
  display ("Volumes");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".vol");
    write_array (filename, ochoice);
  }

  return (0);

}


/***********************************************************
                                       Space Telescope Science Institute

Synopsis:
	wind_element provides a detailed look at a single cell

Arguments:


Returns:

Description:

	Note that this routine does not calculate anything.  It
	simple reports on variables in the wind array.
	
		
Notes:

History:
	080804	ksl	60b -- Added reporting of partition function
	101106	ksl	69 -- Changed variable name in rdint to make
			it more obvious what was being discused here.

**************************************************************/

int
wind_element (w)
     WindPtr w;
{
  PlasmaPtr xplasma;
  int m, n, i, j, nn, mm;
  int first, last;
  int ndom;

  ndom = w->ndom;



  n = 50;
a:printf ("There are %i wind elements in this model\n", NDIM2);
  rdint ("Wind.array.element", &n);

  if (n < 0)
    goto b;
  else if (n > NDIM2)
  {
    printf ("No, there are %i wind elements, not %i\n", NDIM2, n);
    goto a;
  }

  wind_n_to_ij (ndom, n, &i, &j);
  xplasma = &plasmamain[w[n].nplasma];

  Log
    ("Element %d (%d,%d)  inwind %d plasma cell %d ntot %d nioniz %d nrad %d\n",
     n, i, j, w[n].inwind, xplasma->nplasma, xplasma->ntot, xplasma->nioniz, xplasma->nrad);
  Log ("xyz %8.2e %8.2e %8.2e vel %8.2e %8.2e %8.2e\n", w[n].x[0], w[n].x[1], w[n].x[2], w[n].v[0], w[n].v[1], w[n].v[2]);
  Log ("r theta %12.6e %12.6e \n", w[n].rcen, w[n].thetacen / RADIAN);

  Log ("rho %8.2e nh %8.2e ne %8.2e t_r %8.2e t_e %8.2e w %8.2e vol %8.2e\n",
       xplasma->rho, xplasma->rho * rho2nh, xplasma->ne, xplasma->t_r, xplasma->t_e, xplasma->w, w[n].vol);

  if (w[n].inwind < 0)
    Log ("\n# Cell is not inwind, expect all zeros to follow\n\n");
  /*70d - added compton - ksl */
  /*70g compton removed from luminosity reporting, it is now a cooling mechanism but does not produce photons
     DR cooling also added in to report */
  Log
    ("t_e  %8.2e lum_tot %8.2e lum_lines  %8.2e lum_ff  %8.2e lum_rr     %8.2e \n",
     xplasma->t_e, xplasma->lum_tot_ioniz, xplasma->lum_lines_ioniz, xplasma->lum_ff_ioniz, xplasma->lum_rr_ioniz);
  Log
    ("t_e %8.2e cool_tot %8.2e lum_lines  %8.2e lum_ff  %8.2e cool_rr     %8.2e cool_comp %8.2e cool_adiab %8.2e cool_DR %8.2e cool_DI %8.2e\n",
     xplasma->t_e, xplasma->cool_tot_ioniz + xplasma->cool_comp_ioniz + xplasma->cool_adiabatic_ioniz + xplasma->cool_dr_ioniz,
     xplasma->lum_lines_ioniz, xplasma->lum_ff_ioniz, xplasma->cool_rr_ioniz, xplasma->cool_comp_ioniz, xplasma->cool_adiabatic_ioniz,
     xplasma->cool_dr_ioniz, xplasma->cool_di_ioniz);
  Log ("t_r %8.2e heat_tot %8.2e heat_lines %8.2e heat_ff %8.2e heat_photo %8.2e heat_auger %8.2e heat_comp %8.2e heat_icomp %8.2e\n",
       xplasma->t_r, xplasma->heat_tot, xplasma->heat_lines, xplasma->heat_ff, xplasma->heat_photo, xplasma->heat_auger, xplasma->heat_comp,
       xplasma->heat_ind_comp);



  Log ("Recombination cooling   HII>HI %8.2e HeII>HeI %8.2e HeIII>HeII %8.2e Metals %8.2e\n", xplasma->cool_rr_ion[0],
       xplasma->cool_rr_ion[2], xplasma->cool_rr_ion[3], xplasma->cool_rr_metals);
  Log ("Photoionization heating HI>HII %8.2e HeI>HeII %8.2e HeII>HeIII %8.2e Metals %8.2e\n", xplasma->heat_ion[0], xplasma->heat_ion[2],
       xplasma->heat_ion[3], xplasma->heat_z);



  Log ("The ratio of rad (total) cooling to heating is %8.2f (%8.2f) \n",
       xplasma->lum_tot_ioniz / xplasma->heat_tot,
       (xplasma->lum_tot_ioniz + xplasma->cool_adiabatic_ioniz + xplasma->cool_comp_ioniz + xplasma->cool_dr_ioniz) / xplasma->heat_tot);
  Log ("Adiabatic cooling %8.2e is %8.2g of total cooling\n",
       xplasma->cool_adiabatic_ioniz,
       xplasma->cool_adiabatic_ioniz / (xplasma->lum_tot + xplasma->cool_adiabatic + xplasma->cool_comp_ioniz + xplasma->cool_dr_ioniz));
  /*70g NSH compton and DR cooling are now reported seperately. */
  Log ("Compton cooling   %8.2e is %8.2g of total cooling\n",
       xplasma->cool_comp_ioniz,
       xplasma->cool_comp_ioniz / (xplasma->lum_tot + xplasma->cool_adiabatic + xplasma->cool_comp_ioniz + xplasma->cool_dr_ioniz));
  Log ("DR cooling        %8.2e is %8.2g of total cooling\n", xplasma->cool_dr_ioniz,
       xplasma->cool_dr_ioniz / (xplasma->lum_tot + xplasma->cool_adiabatic + xplasma->cool_comp_ioniz + xplasma->cool_dr_ioniz));
  Log ("Number of ionizing photons in cell nioniz %d\n", xplasma->nioniz);
  Log ("Log Ionization parameter in this cell U %4.2f xi %4.2f\n", log10 (xplasma->ip), log10 (xplasma->xi));   //70h NSH computed ionizaion parameter
  Log ("ioniz %8.2e %8.2e %8.2e %8.2e %8.2e\n",
       xplasma->ioniz[0], xplasma->ioniz[1], xplasma->ioniz[2], xplasma->ioniz[3], xplasma->ioniz[4]);
  Log
    ("Convergence status: whole %d converging %d t_r %8.2e t_e %8.2e hc %8.2e \n",
     xplasma->converge_whole, xplasma->converging, xplasma->converge_t_r, xplasma->converge_t_e, xplasma->converge_hc);

  Log ("Densities:\n");
  for (nn = 0; nn < 5; nn++)
  {
    first = ele[nn].firstion;
    last = first + ele[nn].nions;
    Log ("%-5s ", ele[nn].name);
    for (m = first; m < last; m++)
      Log (" %8.2e", xplasma->density[m]);
    Log ("\n");
  }


  Log ("Partition function:\n");
  for (nn = 0; nn < 5; nn++)
  {
    first = ele[nn].firstion;
    last = first + ele[nn].nions;
    Log ("%-5s ", ele[nn].name);
    for (m = first; m < last; m++)
      Log (" %8.2e", xplasma->partition[m]);
    Log ("\n");
  }


  Log ("Non LTE Level densities:\n");


  mm = 0;
  for (nn = 0; nn < 10; nn++)
  {
    while (ion[mm].nlte <= 0 && mm < nions)
      mm++;
    if (mm == nions)
      break;

    first = ion[mm].first_levden;
    last = first + ion[mm].nlte;
    Log ("ion %3d %3d", ion[mm].z, ion[mm].istate);
    for (m = first; m < last; m++)
      Log (" %8.2e", xplasma->levden[m]);
    Log ("\n");
    mm++;
  }

  Log ("Spectral model details:\n");
  for (nn = 0; nn < geo.nxfreq; nn++)
  {
    Log ("numin= %9.2e (%9.2e) numax= %9.2e (%9.2e) Model= %2d PL_log_w= %9.2e PL_alpha= %9.2e Exp_w= %9.2e EXP_temp= %9.2e\n",
         xplasma->fmin_mod[nn], geo.xfreq[nn], xplasma->fmax_mod[nn], geo.xfreq[nn + 1], xplasma->spec_mod_type[nn], xplasma->pl_log_w[nn],
         xplasma->pl_alpha[nn], xplasma->exp_w[nn], xplasma->exp_temp[nn]);
  }

  Log ("Flux:\n");
  Log ("F_vis_w = %9.2e  F_vis_phi = %9.2e  F_vis_z = %9.2e \n", xplasma->F_vis[0], xplasma->F_vis[1], xplasma->F_vis[2]);
  Log ("F_UV_w  = %9.2e  F_UV_phi  = %9.2e  F_UV_z  = %9.2e \n", xplasma->F_UV[0], xplasma->F_UV[1], xplasma->F_UV[2]);
  Log ("F_Xray_w= %9.2e  F_Xray_phi= %9.2e  F_Xray_z= %9.2e \n", xplasma->F_Xray[0], xplasma->F_Xray[1], xplasma->F_Xray[2]);




  goto a;

b:return (0);

}





/**************************************************************************


  Synopsis:  
	tau_h_summary (w, rootname, ochoice)

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/
int
tau_h_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = 6.e-18 * plasmamain[nplasma].density[0] * pow (w[n].vol, 0.333);

    }
  }
  display ("tau_h summary");

  return (0);

}





/**************************************************************************


  Synopsis:  
	coolheat_summary (w, rootname, ochoice)

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/
int
coolheat_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].lum_tot_ioniz / plasmamain[nplasma].heat_tot;

    }
  }
  display ("Cooling over heating");

  return (0);

}




/**************************************************************************


  Synopsis:  

  complete_file summary produces a standardised set of output files 
  from the wind_save_file

  Description:	

  At present(1409), we we print temperatures, electron densities,
  convergence information and ion denisities and fractions for
  C III, IV, V, 
  N IV, V, VI
  O V, VI, VII
  Si III, IV

  Arguments:  

  Returns:

  Notes:

  History:

  1409	ksl	Eliminated several files that were of limited use
  		to users.  The deleted files are printed out if
		  DEBUG is set.
  1411  JM debug is now deprecated so replaced with FULL_ION_SUMMARTY

 ************************************************************************/

#define FULL_ION_SUMMARY 0


int
complete_file_summary (w, root, ochoice)
     WindPtr w;
     char root[];
     int ochoice;
{
  temp_summary (w, root, ochoice);
  temp_rad (w, root, ochoice);
  electron_summary (w, root, ochoice);
  convergence_summary (w, root, ochoice);
  ion_summary (w, 6, 3, 0, root, ochoice);
  ion_summary (w, 6, 4, 0, root, ochoice);
  ion_summary (w, 6, 5, 0, root, ochoice);
  ion_summary (w, 7, 4, 0, root, ochoice);
  ion_summary (w, 7, 5, 0, root, ochoice);
  ion_summary (w, 7, 6, 0, root, ochoice);
  ion_summary (w, 8, 4, 0, root, ochoice);
  ion_summary (w, 8, 5, 0, root, ochoice);
  ion_summary (w, 8, 6, 0, root, ochoice);
  ion_summary (w, 8, 7, 0, root, ochoice);
  ion_summary (w, 14, 3, 0, root, ochoice);
  ion_summary (w, 14, 4, 0, root, ochoice);
  ion_summary (w, 14, 5, 0, root, ochoice);

  ion_summary (w, 6, 3, 1, root, ochoice);
  ion_summary (w, 6, 4, 1, root, ochoice);
  ion_summary (w, 6, 5, 1, root, ochoice);
  ion_summary (w, 7, 4, 1, root, ochoice);
  ion_summary (w, 7, 5, 1, root, ochoice);
  ion_summary (w, 7, 6, 1, root, ochoice);
  ion_summary (w, 8, 4, 1, root, ochoice);
  ion_summary (w, 8, 5, 1, root, ochoice);
  ion_summary (w, 8, 6, 1, root, ochoice);
  ion_summary (w, 8, 7, 1, root, ochoice);
  ion_summary (w, 14, 3, 1, root, ochoice);
  ion_summary (w, 14, 4, 1, root, ochoice);
  ion_summary (w, 14, 5, 1, root, ochoice);

#if FULL_ION_SUMMARY
  ion_summary (w, 6, 3, 2, root, ochoice);
  ion_summary (w, 6, 4, 2, root, ochoice);
  ion_summary (w, 6, 5, 2, root, ochoice);
  ion_summary (w, 7, 4, 2, root, ochoice);
  ion_summary (w, 7, 5, 2, root, ochoice);
  ion_summary (w, 7, 6, 2, root, ochoice);
  ion_summary (w, 8, 4, 2, root, ochoice);
  ion_summary (w, 8, 5, 2, root, ochoice);
  ion_summary (w, 8, 6, 2, root, ochoice);
  ion_summary (w, 8, 7, 2, root, ochoice);
  ion_summary (w, 14, 3, 2, root, ochoice);
  ion_summary (w, 14, 4, 2, root, ochoice);
  ion_summary (w, 14, 5, 2, root, ochoice);


  ion_summary (w, 6, 3, 3, root, ochoice);
  ion_summary (w, 6, 4, 3, root, ochoice);
  ion_summary (w, 6, 5, 3, root, ochoice);
  ion_summary (w, 7, 4, 3, root, ochoice);
  ion_summary (w, 7, 5, 3, root, ochoice);
  ion_summary (w, 7, 6, 3, root, ochoice);
  ion_summary (w, 8, 4, 3, root, ochoice);
  ion_summary (w, 8, 5, 3, root, ochoice);
  ion_summary (w, 8, 6, 3, root, ochoice);
  ion_summary (w, 8, 7, 3, root, ochoice);
  ion_summary (w, 14, 3, 3, root, ochoice);
  ion_summary (w, 14, 4, 3, root, ochoice);
  ion_summary (w, 14, 5, 3, root, ochoice);
#endif

  return (0);
}

/* A summary of the regions in the wind */

int
wind_reg_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;


  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = -9;
    if (w[n].inwind >= 0)
    {
      aaa[n] = w[n].inwind;

    }
  }
  display ("Regions of the wind");

  return (0);

}


/* A summary of the dvds_ave */

int
dvds_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  char filename[LINELENGTH], suffix[LINELENGTH];
  int n, ichoice;
  struct photon p;
  //double v1[3]



  rdint ("dvds_ave (0) dvdx (1) dvdy (2) dvdz (3) component LoS (4):", &ichoice);

  if (ichoice == 0)
  {
    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {
        aaa[n] = w[n].dvds_ave;

      }
    }
    display ("Average dvds");
    if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".dvds");
      write_array (filename, ochoice);
    }
  }

  else if (ichoice > 0 && ichoice < 4)
  {


    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {
        p.lmn[0] = 0.0;
        p.lmn[1] = 0.0;
        p.lmn[2] = 0.0;
        stuff_v (w[n].xcen, p.x);

        p.lmn[ichoice - 1] = 1.0;
        aaa[n] = dvwind_ds_cmf (&p);

      }
    }
    display ("Average dvds");
    if (ochoice)
    {
      strcpy (filename, rootname);
      sprintf (suffix, ".dvds_%i", ichoice - 1);
      strcat (filename, suffix);
      write_array (filename, ochoice);
    }
  }

  else
  {
    get_los_dvds (w, rootname, ochoice);
  }


  return (0);
}

/* A summary of inner shell ionization */
/* NSH - this code removed May 18 
int
inner_shell_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;


  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0.0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].gamma_inshl[0];
    }
  }
  display ("Inner_shell");
  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".is");
    write_array (filename, ochoice);
  }
  return (0);

}
*/

/* A summary of the Ionization parameter - might not always be present */

int
IP_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = ((plasmamain[nplasma].ip));
    }
  }
  display ("Ionization parameter");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".IP");
    write_array (filename, ochoice);

  }


  /* JM added printout for xi too */
  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = ((plasmamain[nplasma].xi));
    }
  }
  display ("Xi Ionization parameter");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".xi");
    write_array (filename, ochoice);

  }



  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = ((plasmamain[nplasma].ip_direct));
    }
  }
  display ("Log Ionization parameter (direct)");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".IP_direct");
    write_array (filename, ochoice);

  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = ((plasmamain[nplasma].ip_scatt));
    }
  }
  display ("Log Ionization parameter (scattered)");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".IP_scatt");
    write_array (filename, ochoice);

  }




  return (0);

}



/* A summary of the Sim alpha parameter - might not always be present.

   1108	ksl	Adapted for new version of banded alpha
   1208 nsh	Changed names - reference to sim removed to make way for exponential estimators
   1401 nsh	Added more information - this now writes out all the spectral model parameters
 */

int
alpha_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n, m;
  char filename[LINELENGTH];
  int nplasma;
  char word[LINELENGTH];



  for (m = 0; m < geo.nxfreq; m++)
  {
    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {
        nplasma = w[n].nplasma;
        aaa[n] = plasmamain[nplasma].pl_alpha[m];
      }
    }

    strcpy (word, "");
    sprintf (word, ".pl_alpha%03d", m);
    display (word);

    if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, word);
      write_array (filename, ochoice);
    }
  }

  for (m = 0; m < geo.nxfreq; m++)
  {
    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {
        nplasma = w[n].nplasma;
        aaa[n] = plasmamain[nplasma].pl_log_w[m];
      }
    }

    strcpy (word, "");
    sprintf (word, ".pl_log_w%03d", m);
    display (word);

    if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, word);
      write_array (filename, ochoice);
    }
  }


  for (m = 0; m < geo.nxfreq; m++)
  {
    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {
        nplasma = w[n].nplasma;
        aaa[n] = plasmamain[nplasma].exp_temp[m];
      }
    }

    strcpy (word, "");
    sprintf (word, ".exp_temp%03d", m);
    display (word);

    if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, word);
      write_array (filename, ochoice);
    }
  }

  for (m = 0; m < geo.nxfreq; m++)
  {
    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {
        nplasma = w[n].nplasma;
        aaa[n] = plasmamain[nplasma].exp_w[m];
      }
    }

    strcpy (word, "");
    sprintf (word, ".exp_w%03d", m);
    display (word);

    if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, word);
      write_array (filename, ochoice);
    }
  }

  for (m = 0; m < geo.nxfreq; m++)
  {
    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {
        nplasma = w[n].nplasma;
        aaa[n] = plasmamain[nplasma].spec_mod_type[m];
      }
    }

    strcpy (word, "");
    sprintf (word, ".spec_mod_type%03d", m);
    display (word);

    if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, word);
      write_array (filename, ochoice);
    }
  }



  for (m = 0; m < geo.nxfreq; m++)
  {
    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {
        nplasma = w[n].nplasma;
        aaa[n] = plasmamain[nplasma].fmin_mod[m];
      }
    }

    strcpy (word, "");
    sprintf (word, ".fmin_mod%03d", m);
    display (word);

    if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, word);
      write_array (filename, ochoice);
    }
  }

  for (m = 0; m < geo.nxfreq; m++)
  {
    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {
        nplasma = w[n].nplasma;
        aaa[n] = plasmamain[nplasma].fmax_mod[m];
      }
    }

    strcpy (word, "");
    sprintf (word, ".fmax_mod%03d", m);
    display (word);

    if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, word);
      write_array (filename, ochoice);
    }
  }

  return (0);

}



/* A summary of frequency banded radiation density - a crude spectrum for a cell.

 1301 nsh	Written
 */

int
J_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int i, n;
  char filename[LINELENGTH];
  char number[2], line_number[10];
  int nplasma;
  int uplvl, llvl, njump, lu, ll;
  struct lines *line_ptr;

  njump = 0;
  line_ptr = NULL;

  i = 1;

  rdint ("Band number for J or macro atom J (0), or backup (-1)", &i);


  while (i >= 0)
  {
    if (i == 0)
    {
      printf ("H alpha is 3->2 in this notation!\n");
      rdint ("Upper level macro atom", &uplvl);
      rdint ("Lower level macro atom", &llvl);

      /* Convert 'user levels' into actually levels. i.e. 1 is 0! */
      uplvl = uplvl - 1;
      llvl = llvl - 1;



      /* now we need to find the jbar estimator and line 
         pointer corresponding to this transition */
      njump = 0;
      printf ("Level %i has %i upwards jumps %i downwards jumps\n", llvl + 1, xconfig[llvl].n_bbu_jump, xconfig[llvl].n_bbd_jump);

      while (njump < xconfig[llvl].n_bbu_jump)
      {
        line_ptr = &line[xconfig[llvl].bbu_jump[njump]];
        lu = line_ptr->nconfigu;
        ll = line_ptr->nconfigl;

        if (ll == llvl && lu == uplvl)
          break;
        njump++;
      }

      if (njump >= xconfig[llvl].n_bbu_jump)
      {
        Error ("Couldn't find this transition, try something else!\n");
        return (0);
      }
    }

    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {
        nplasma = w[n].nplasma;
        if (i == 0)
          aaa[n] = macromain[nplasma].jbar_old[xconfig[llvl].bbu_indx_first + njump];
        else
          aaa[n] = (plasmamain[nplasma].xj[i]);
      }
    }

    printf ("Line wavelength is %.2f\n", (VLIGHT / line_ptr->freq) / ANGSTROM);
    printf ("Line freq is %8.4e\n", line_ptr->freq);
    printf ("njump %i llvl %i uplvl %i nres %i", njump, llvl, uplvl, xconfig[llvl].bbu_jump[njump]);
    display ("J in cell");
    //printf ("i=%i", i);
    sprintf (number, "%i", i);

    if (ochoice)
    {
      strcpy (filename, rootname);
      if (i == 0)
      {
        sprintf (line_number, "%ito%i", uplvl, llvl);
        strcat (filename, ".Jbar_");
        strcat (filename, line_number);
      }
      else
      {
        strcat (filename, ".J_band");
        strcat (filename, number);
      }
      write_array (filename, ochoice);

    }

    rdint ("Band number for J or macro atom J (0), or backup (-1)", &i);
  }
  return (0);

}


int
J_scat_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;



  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = (plasmamain[nplasma].j);
    }
  }
  display ("J in cell");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".J_tot");
    write_array (filename, ochoice);
  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = (plasmamain[nplasma].j_direct);
    }
  }
  display ("J in cell from direct photons");
  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".J_direct");
    write_array (filename, ochoice);
  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = (plasmamain[nplasma].j_scatt);
    }
  }
  display ("J in cell from scattered photons");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".J_scatt");
    write_array (filename, ochoice);
  }


  return (0);

}






//Split of photons from different sources in the cell.

int
phot_split (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = (plasmamain[nplasma].ntot_wind);
    }
  }
  display ("Wind photons in cell");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".ntot_wind");
    write_array (filename, ochoice);

  }
  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = (plasmamain[nplasma].ntot_agn);
    }
  }
  display ("AGN photons in cell");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".ntot_agn");
    write_array (filename, ochoice);

  }
  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = (plasmamain[nplasma].ntot_disk);
    }
  }
  display ("Disk photons in cell");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".ntot_disk");
    write_array (filename, ochoice);

  }
  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = (plasmamain[nplasma].ntot_star);
    }
  }
  display ("Stellar photons in cell");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".ntot_star");
    write_array (filename, ochoice);

  }
  return (0);

}

int
thompson (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;
  double ne;
  char filename[LINELENGTH];

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      ne = plasmamain[nplasma].ne;
      aaa[n] = (THOMPSON * ne) * pow (w[n].vol, 0.333);
    }
  }
  display ("Thompson optical depths");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".thomp");
    write_array (filename, ochoice);
  }

  return (0);

}



int
nscat_split (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;
  char filename[LINELENGTH];

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].nscat_es;
    }
  }
  display ("Thompson scatters");

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].nscat_res;
    }
  }
  display ("Resonant scatters");






  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".thomp");
    write_array (filename, ochoice);
  }

  return (0);

}

int
convergence_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;
  char filename[LINELENGTH];

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].converge_whole;
    }
  }
  display ("Convergence (0=converged.  Higher numbers indicate one or more convergence tests failed)");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".conv");
    write_array (filename, ochoice);
  }

  return (0);

}


/* 

   1112	ksl	Write out arrays of use in evaluating what is going on in convergence of
   		various models
*/

int
convergence_all (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;
  char filename[LINELENGTH];




  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].converge_t_e;
    }
  }
  display ("t_e Convergence (0=converged)");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".conv_te");
    write_array (filename, ochoice);
  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].converge_hc;
    }
  }
  display ("hc Convergence (0=converged)");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".conv_hc");
    write_array (filename, ochoice);
  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].converge_t_r;
    }
  }
  display ("t_r Convergence (0=converged)");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".conv_tr");
    write_array (filename, ochoice);
  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].converge_whole;
    }
  }
  display ("Convergence (0=converged)");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".conv_whole");
    write_array (filename, ochoice);
  }


  return (0);
}

/* 

   1306	nsh	Write out information pertaining to the models used in each cell/frequency band
*/

int
model_bands (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n, m;
  int nplasma;
  char filename[LINELENGTH];
  char word[LINELENGTH];



  for (m = 0; m < geo.nxfreq; m++)
  {
    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {
        nplasma = w[n].nplasma;
        aaa[n] = plasmamain[nplasma].nxtot[m];
      }
    }

    strcpy (word, "");
    sprintf (word, ".nxtot%03d", m);
    display (word);

    if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, word);
      write_array (filename, ochoice);
    }
  }

  if (geo.nxfreq == 0)
  {
    printf ("This model does not use bands calculating the ionization\n");
    return (0);
  }


  /* Now write out the value of xj in each band */

  for (m = 0; m < geo.nxfreq; m++)
  {
    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {
        nplasma = w[n].nplasma;
        aaa[n] = plasmamain[nplasma].xj[m];
      }
    }

    strcpy (word, "");
    sprintf (word, ".xj%03d", m);
    display (word);

    if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, word);
      write_array (filename, ochoice);
    }
  }


  /* Now write out the value of xave_freq in each band */

  for (m = 0; m < geo.nxfreq; m++)
  {
    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {
        nplasma = w[n].nplasma;
        aaa[n] = plasmamain[nplasma].xave_freq[m];
      }
    }

    strcpy (word, "");
    sprintf (word, ".xave_freq%03d", m);
    display (word);

    if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, word);
      write_array (filename, ochoice);
    }
  }


  /* Now write out the value of nxtot in each band */

  for (m = 0; m < geo.nxfreq; m++)
  {
    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {
        nplasma = w[n].nplasma;
        aaa[n] = plasmamain[nplasma].nxtot[m];
      }
    }

    strcpy (word, "");
    sprintf (word, ".nxtot%03d", m);
    display (word);

    if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, word);
      write_array (filename, ochoice);
    }
  }

  return (0);
}



/* A summary of adiabatic cooling */
int
heatcool_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;
  char filename[LINELENGTH];
//OLD  float x;

//OLD  x = wind_luminosity (0.0, VERY_BIG, MODE_CMF_TIME);

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].heat_tot;
      if (w[n].div_v < 0.0)     // add in if it is negative and hence a heating term
      {
        aaa[n] += -1.0 * (plasmamain[nplasma].cool_adiabatic_ioniz);
      }
    }
  }
  display ("Total heating");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".heat_tot");
    write_array (filename, ochoice);
  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].heat_lines;
    }
  }
  display ("Line heating");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".heat_lines");
    write_array (filename, ochoice);
  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].heat_ff;
    }
  }
  display ("FF heating");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".heat_ff");
    write_array (filename, ochoice);
  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].heat_comp;
    }
  }
  display ("Compton heating");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".heat_comp");
    write_array (filename, ochoice);
  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].heat_ind_comp;
    }
  }
  display ("Induced Compton heating");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".heat_ind_comp");
    write_array (filename, ochoice);
  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].heat_photo;
    }
  }
  display ("Photo heating");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".heat_photo");
    write_array (filename, ochoice);
  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].lum_lines_ioniz;
    }
  }
  display ("Line Luminosity");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".lum_lines");
    write_array (filename, ochoice);
  }



  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].cool_adiabatic_ioniz;
    }
  }
  display ("Adiabatic Luminosity");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".cool_adiabatic");
    write_array (filename, ochoice);
  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].lum_ff_ioniz;
    }
  }
  display ("Free Free Luminosity");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".lum_ff");
    write_array (filename, ochoice);
  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].cool_comp_ioniz;
    }
  }
  display ("Compton Luminosity");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".cool_comp");
    write_array (filename, ochoice);
  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].cool_dr_ioniz;
    }
  }
  display ("DR Luminosity");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".cool_dr");
    write_array (filename, ochoice);
  }

  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].cool_rr_ioniz;
    }
  }
  display ("FB Luminosity");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".cool_rr");
    write_array (filename, ochoice);
  }





  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] =
        plasmamain[nplasma].cool_rr_ioniz + plasmamain[nplasma].cool_dr_ioniz +
        plasmamain[nplasma].cool_comp_ioniz + plasmamain[nplasma].lum_ff_ioniz + plasmamain[nplasma].lum_lines_ioniz;
      if (w[n].div_v >= 0.0)    //only add in if it is treated as a cooling term
      {
        aaa[n] += plasmamain[nplasma].cool_adiabatic_ioniz;
      }
    }
  }
  display ("Total Cooling");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".lum_total");
    write_array (filename, ochoice);
  }


  for (n = 0; n < NDIM2; n++)
  {
    aaa[n] = 0;
    if (w[n].inwind >= 0)
    {
      nplasma = w[n].nplasma;
      aaa[n] = plasmamain[nplasma].lum_tot_ioniz;
    }
  }
  display ("Total Radiating Luminosity");

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".lum_tot");
    write_array (filename, ochoice);
  }

  return (0);
}








/**************************************************************************


  Synopsis:  
  A summary of the important quantites in a given cell.
  Primarily used for easily making plots

  Description:  

  Arguments:
    w WindPtr
    rootname filename of root pf 
    ochoice whether to save to file  

  Returns:

  Notes:

  History:
  1411 JM coded

************************************************************************/

int
complete_physical_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n, np;
  char filename[LINELENGTH];
  double he1den, he2den, he3den;
  double h1den, h2den, c3den, c4den, c5den;
  double n5den, o6den, si4den;
  int frac_choice;
  int ii, jj;
  double vtot;
  FILE *fptr, *fopen ();
  PlasmaPtr xplasma;
  int ndom;

  rdint ("Save ions as densities (0) or fractions? (1)", &frac_choice);

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".complete");
    fptr = fopen (filename, "w");
  }
  else
    printf ("This mode is recommended purely for file output\n");


  /* JM 1411 -- First we have to write out some headers so that 
     astropy can read the output */

  printf ("n\tnplasma\tinwind\ti\tj\tx\tz\tv\tvx\tvy\tvz\tdvds_ave\tvol\t \
rho\tne\tte\ttr\tnphot\tw\tave_freq\tIP\tconv\tconv_tr\tconv_te\tconv_hc\t \
cool_tot\tlum_tot\tlum_rr\tcool_rr\tlum_ff\tlum_lines\tcool_adiabatic\tcool_comp\tcool_dr\tcool_DI \
heat_tot\theat_photo\theat_auger\theat_lines\theat_ff\theat_comp\theat_ind_comp\theat_shock\t \
ionH1\tionH2\tionHe1\tionHe2\tionHe3\tionC3\tionC4\tionC5\tionN5\tionO6\tionSi4\n");

  if (ochoice)
    fprintf (fptr, "n\tnplasma\tinwind\ti\tj\tx\tz\tr\ttheta\tv\tvx\tvy\tvz\tdvds_ave\tvol\t \
rho\tne\tte\ttr\tnphot\tw\tave_freq\tIP\tXi\tconv\tconv_tr\tconv_te\tconv_hc\t \
cool_tot\tlum_tot\tlum_rr\tcool_rr\tlum_ff\tlum_lines\tcool_adiabatic\tcool_comp\tcool_dr\tcool_DI \t\
heat_tot\theat_photo\theat_auger\theat_lines\theat_ff\theat_comp\theat_ind_comp\theat_shock\t \
ionH1\tionH2\tionHe1\tionHe2\tionHe3\tionC3\tionC4\tionC5\tionN5\tionO6\tionSi4\n");


  Log ("swind_sub does not work yet\n");
  ndom = 0;
  np = 0;
  for (n = 0; n < NDIM2; n++)
  {
    wind_n_to_ij (ndom, n, &ii, &jj);

    if (w[n].inwind >= 0)
    {
      np = w[n].nplasma;

      vtot = sqrt (w[n].v[0] * w[n].v[0] + w[n].v[1] * w[n].v[1] + w[n].v[2] * w[n].v[2]);

      xplasma = &plasmamain[np];

      /* find the density of the main ions (or fractions if frac_choice == 1) */
      h1den = get_density_or_frac (xplasma, 1, 1, frac_choice);
      h2den = get_density_or_frac (xplasma, 1, 2, frac_choice);
      he1den = get_density_or_frac (xplasma, 2, 1, frac_choice);
      he2den = get_density_or_frac (xplasma, 2, 2, frac_choice);
      he3den = get_density_or_frac (xplasma, 2, 3, frac_choice);
      c3den = get_density_or_frac (xplasma, 6, 3, frac_choice);
      c4den = get_density_or_frac (xplasma, 6, 4, frac_choice);
      c5den = get_density_or_frac (xplasma, 6, 5, frac_choice);
      n5den = get_density_or_frac (xplasma, 7, 5, frac_choice);
      o6den = get_density_or_frac (xplasma, 8, 6, frac_choice);
      si4den = get_density_or_frac (xplasma, 14, 4, frac_choice);

      /* printf("%i %i %i %i %i %8.4e %8.4e %8.4e %8.4e %8.4e  %8.4e %8.4e %8.4e %8.4e \
         %8.4e %8.4e %8.4e %i %8.4e %8.4e %8.4e %i %8.4e %8.4e %8.4e %8.4e \
         %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e \
         %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e \
         %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e\n",
         n, np, w[n].inwind, ii, jj, w[n].x[0], w[n].x[2], vtot, w[n].v[0], w[n].v[1], w[n].v[2], w[n].dvds_ave, w[n].vol, 
         plasmamain[np].rho, plasmamain[np].ne, plasmamain[np].t_e, plasmamain[np].t_r, plasmamain[np].ntot,
         plasmamain[np].w, plasmamain[np].ave_freq, plasmamain[np].ip, plasmamain[np].converge_whole, 
         plasmamain[np].converge_t_r, plasmamain[np].converge_t_e, plasmamain[np].converge_hc, 
         plasmamain[np].cool_tot_ioniz, plasmamain[np].lum_tot, plasmamain[np].cool_rr, 
         plasmamain[np].lum_ff, plasmamain[np].lum_lines, plasmamain[np].cool_adiabatic, 
         plasmamain[np].cool_comp, plasmamain[np].cool_dr, plasmamain[np].heat_tot, plasmamain[np].heat_photo, 
         plasmamain[np].heat_lines , plasmamain[np].heat_ff , plasmamain[np].heat_comp, plasmamain[np].heat_ind_comp,
         h1den, h2den, he1den, he2den, he3den, c3den, c4den, c5den, n5den, o6den, si4den);
       */

      if (ochoice)
        fprintf (fptr, "%i %i %i %i %i %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e \
            %8.4e %8.4e %8.4e %i %8.4e %8.4e %8.4e %8.4e %i %8.4e %8.4e %8.4e %8.4e \
            %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e\
            %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e \
            %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e\n", n, np, w[n].inwind, ii, jj, w[n].x[0], w[n].x[2], w[n].rcen, w[n].thetacen / RADIAN, vtot, w[n].v[0], w[n].v[1], w[n].v[2], w[n].dvds_ave, w[n].vol, plasmamain[np].rho, plasmamain[np].ne, plasmamain[np].t_e, plasmamain[np].t_r, plasmamain[np].ntot, plasmamain[np].w, plasmamain[np].ave_freq, plasmamain[np].ip, plasmamain[np].xi, plasmamain[np].converge_whole, plasmamain[np].converge_t_r, plasmamain[np].converge_t_e, plasmamain[np].converge_hc, plasmamain[np].cool_tot_ioniz + plasmamain[np].cool_comp_ioniz + plasmamain[np].cool_adiabatic_ioniz + plasmamain[np].cool_dr_ioniz, plasmamain[np].lum_tot_ioniz, plasmamain[np].lum_rr_ioniz, plasmamain[np].cool_rr_ioniz, plasmamain[np].lum_ff_ioniz, plasmamain[np].lum_lines_ioniz, plasmamain[np].cool_adiabatic_ioniz, plasmamain[np].cool_comp_ioniz, plasmamain[np].cool_dr_ioniz, plasmamain[np].cool_di_ioniz, plasmamain[np].heat_tot, plasmamain[np].heat_photo, plasmamain[np].heat_auger, plasmamain[np].heat_lines, plasmamain[np].heat_ff, plasmamain[np].heat_comp, plasmamain[np].heat_ind_comp, plasmamain[np].heat_shock, h1den, h2den, he1den, he2den, he3den, c3den, c4den, c5den, n5den, o6den, si4den);
    }
    else
    {
      /* if we aren't inwind then print out a load of zeroes */

      /* printf("%i %i %i %i %i %8.4e %8.4e 0.0 0.0 0.0 0.0 0.0 0.0 0.0 \
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 \
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 \
         0.0 0.0 0.0 0.0 0.0 0.0 \
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n",
         n, np, w[n].inwind, ii, jj, w[n].x[0], w[n].x[2]);
       */

      if (ochoice)
        fprintf (fptr, "%i %i %i %i %i %8.4e %8.4e 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 \
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 \
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 \
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 \
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n", n, np, w[n].inwind, ii, jj, w[n].x[0], w[n].x[2]);
    }
  }


  if (ochoice)
  {
    fclose (fptr);
    printf ("\nSaved summary of physical quantites in %s, use py_read_output.py to read\n", filename);
  }

  return (0);


}

/**************************************************************************


  Synopsis:  
  complete_ion_summary outputs a file with all of the ion fractions for a given cell


  History:
  1602 NSH coded

************************************************************************/



int
complete_ion_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;

{
  char cell[5];
  PlasmaPtr xplasma;
  FILE *fptr = NULL, *fopen ();
  char filename[LINELENGTH];


  int n, mm;
  n = 50;
a:printf ("There are %i wind elements in this model\n", NDIM2);
  rdint ("Wind.array.element", &n);

  if (n < 0)
    goto b;
  else if (n > NDIM2)
  {
    printf ("No, there are %i wind elements, not %i\n", NDIM2, n);
    goto a;
  }
  printf ("OK, using cell %i\n", n);
  xplasma = &plasmamain[w[n].nplasma];

  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, "_cell_");
    sprintf (cell, "%04d", n);
    strcat (filename, cell);
    printf ("opening file %s\n", filename);

    strcat (filename, ".ions");
    printf ("opening file %s\n", filename);
    fptr = fopen (filename, "w");
  }

  printf ("ion z n(ion) n(ion)/n(H)\n");


  if (ochoice)
    fprintf (fptr, "ion z n(ion) n(ion)/n(H)\n");


  for (mm = 0; mm < nions; mm++)
  {

    printf ("%i %i %e %e\n", mm, ion[mm].z, xplasma->density[mm], xplasma->density[mm] / (xplasma->rho * rho2nh));
    if (ochoice)
    {
      fprintf (fptr, "%i %i %e %e\n", mm, ion[mm].z, xplasma->density[mm], xplasma->density[mm] / (xplasma->rho * rho2nh));
    }

  }
  goto a;

b:return (0);

}





/**************************************************************************


  Synopsis:  
  get_density_or_frac works out the density of an ion element, istate,
  in a plasma cell xplasma.

  if frac_choice is 1, return an ion fraction

  History:
  1411 JM coded

************************************************************************/


double
get_density_or_frac (xplasma, element, istate, frac_choice)
     PlasmaPtr xplasma;
     int element;
     int istate;
     int frac_choice;
{
  int nion, nelem;
  double nh, density;

  /* find the ion and element in the list */
  nion = find_ion (element, istate);

  nelem = find_element (element);

  /* get density of ion */
  density = xplasma->density[nion];

  /* we want an ion fraction, not a density, so divide by nh */
  if (frac_choice)
  {
    nh = xplasma->density[0] + xplasma->density[1];
    density /= ele[nelem].abun * nh;
  }

  return (density);
}


/**************************************************************************


  Synopsis:  
  find_ion is a little routine which finds which number ion in the list corresponds
  to element with istate. e.g. for CIV, element = 6 and istate = 4.

  History:
  1411 JM coded

************************************************************************/


int
find_ion (element, istate)
     int element;
     int istate;
{
  int nion;

  nion = 0;

  while (nion < nions && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;

  return nion;
}

/**************************************************************************


  Synopsis:  
  find_ion is a little routine which finds which number element in the list corresponds
  to element with z == element. e.g. for CIV, element = 6.

  History:
  1411 JM coded

************************************************************************/


int
find_element (element)
     int element;
{
  int n;

  n = 0;

  while (n < nelements && ele[n].z != element)
    n++;

  return n;
}



/**************************************************************************


  Synopsis:  
  get_los_dvds finds the gradient along the LoS 
  of the velocity projected along the same line of sight.

  History:
  1501 JM coded

************************************************************************/

int
get_los_dvds (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  struct photon p;
  struct photon ptest;
  int n;
  double ds, dvds, v1[3], v2[3], xtest[3];
  double lmn[3], diff[3], phase;
  int vchoice, sight_choice;
  double obs_angle, rzero, r;
  char filename[LINELENGTH], suffix[LINELENGTH];

  int ndom;

  ndom = w->ndom;

  vchoice = 0;
  phase = 0;
  obs_angle = 80.0;
  sight_choice = 0;

  rdint ("use component along LoS (0), or magnitude (1):", &sight_choice);
  rdint ("real (0) poloidal (1) or rotational (2) or back(-1):", &vchoice);

  while (vchoice >= 0)
  {
    rddoub ("angle in deg:", &obs_angle);
    rddoub ("phase (0 to 1):", &phase);

    lmn[0] = sin (obs_angle / RADIAN) * cos (-phase * 360. / RADIAN);
    lmn[1] = sin (obs_angle / RADIAN) * sin (-phase * 360. / RADIAN);
    lmn[2] = cos (obs_angle / RADIAN);

    // if (vchoice == 0) 
    //   strcpy (vstring, "");
    // else if (vchoice == 1) 
    //   strcpy (vstring, "rot");
    // else if (vchoice == 2) 
    //   strcpy (vstring, "pol");


    for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].inwind >= 0)
      {

        stuff_v (w[n].xcen, p.x);
        stuff_v (lmn, p.lmn);
        stuff_phot (&p, &ptest);

        vsub (p.x, w[n].x, diff);
        ds = 0.001 * length (diff);
        move_phot (&ptest, ds);

        if (vchoice == 0)
          dvds = dvwind_ds_cmf (&p);

        /* next choice is for turning off rotational velocity */
        else if (vchoice == 1)
        {
          model_velocity (ndom, p.x, v1);
          model_velocity (ndom, ptest.x, v2);
          v1[1] = 0.0;
          v2[1] = 0.0;

          /* calculate the relevant gradient */
          if (sight_choice == 0)
            dvds = fabs (dot (v1, p.lmn) - dot (v2, p.lmn)) / ds;
          else
            dvds = fabs (length (v1) - length (v2)) / ds;
        }

        /* next choice is for turning rotational velocity only */
        else
        {
          r = sqrt (p.x[0] * p.x[0] + p.x[1] * p.x[1]);
          rzero = sv_find_wind_rzero (ndom, p.x);
          v1[0] = v1[2] = 0.0;
          v1[1] = sqrt (GRAV * geo.mstar * rzero) / r;


          r = sqrt (ptest.x[0] * ptest.x[0] + ptest.x[1] * ptest.x[1]);
          rzero = sv_find_wind_rzero (ndom, ptest.x);
          v2[0] = v2[2] = 0.0;
          v2[1] = sqrt (GRAV * geo.mstar * rzero) / r;

          if (p.x[1] != 0.0)
          {
            project_from_cyl_xyz (p.x, v1, xtest);
            stuff_v (xtest, v1);
          }
          if (ptest.x[1] != 0.0)
          {
            project_from_cyl_xyz (ptest.x, v2, xtest);
            stuff_v (xtest, v2);
          }

          /* calculate the relevant gradient */
          if (sight_choice == 0)
            dvds = fabs (dot (v1, p.lmn) - dot (v2, p.lmn)) / ds;
          else
            dvds = fabs (length (v1) - length (v2)) / ds;
        }

        aaa[n] = dvds;

      }
    }

    printf ("vchoice %i y coord %8.4e direction cosines %.2f %.2f %.2f\n", vchoice, p.x[1], p.lmn[0], p.lmn[1], p.lmn[2]);
    display ("dvds along a LoS");

    if (ochoice)
    {
      strcpy (filename, rootname);
      sprintf (suffix, ".dv%i_ds_A%.1f_P%.2f", vchoice, obs_angle, phase);
      strcat (filename, suffix);
      write_array (filename, ochoice);
    }

    rdint ("real (0) poloidal (1) or rotational (2) or back(-1):", &vchoice);

  }

  return (0);
}


/**********************************************************/
/** 
 * @brief Prints grid boundaries to file
 *
 * @param [in] w              Pointer to wind array
 * @param [in] rootname       Root name of simulation
 * @param [in] ochoice        Whether or not to write out 
 * @return     0
 *
 * Outputs the boundaries for the grids for each domain
 *
 * ###Notes###
 * 6/15 - Written by SWM
***********************************************************/
int
grid_summary (WindPtr w, char rootname[], int ochoice)
{
  char filename[LINELENGTH], suffix[LINELENGTH];
  FILE *fopen (), *fptr;
  int i, j;

  printf ("Outputting grid boundaries to file.\n");

  for (j = 0; j < geo.ndomain; j++)
  {
    strcpy (filename, rootname);
    sprintf (suffix, ".dom%d.grid_x.txt", j);
    strcat (filename, suffix);
    fptr = fopen (filename, "w");
    for (i = 0; i <= zdom[j].ndim; i++)
    {
      fprintf (fptr, "%10.5g\n", zdom[j].wind_x[i]);
    }
    fclose (fptr);

    strcpy (filename, rootname);
    sprintf (suffix, ".dom%d.grid_z.txt", j);
    strcat (filename, suffix);
    fptr = fopen (filename, "w");
    for (i = 0; i <= zdom[j].mdim; i++)
    {
      fprintf (fptr, "%10.5g\n", zdom[j].wind_z[i]);
    }
    fclose (fptr);

  }
  return (0);
}






int
flux_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n, np;
  char filename[LINELENGTH];
  int ii, jj;
  FILE *fptr, *fopen ();
//OLD  PlasmaPtr xplasma;
  int ndom, m;

  np = 0;


  if (ochoice)
  {
    strcpy (filename, rootname);
    strcat (filename, ".flux_summary");
    fptr = fopen (filename, "w");
  }
  else
    printf ("This mode is recommended purely for file output\n");


  /* JM 1411 -- First we have to write out some headers so that 
     astropy can read the output */




  if (ochoice)
  {
    fprintf (fptr, "n\tnplasma\tinwind\ti\tj\tx\tz\tr\ttheta ");

    for (m = 0; m < geo.nxfreq; m++)
    {
      fprintf (fptr, "\tF_w%i\tF_p%i\tF_z%i ", m, m, m);

    }
    fprintf (fptr, "\n");
  }

  Log ("swind_sub does not work yet\n");
  ndom = 0;
  for (n = 0; n < NDIM2; n++)
  {
    wind_n_to_ij (ndom, n, &ii, &jj);

    if (w[n].inwind >= 0)
    {
      np = w[n].nplasma;
//OLD      xplasma = &plasmamain[np];
      if (ochoice)
      {
        fprintf (fptr, "%i %i %i %i %i %8.4e %8.4e %8.4e %8.4e ", n, np, w[n].inwind, ii, jj, w[n].x[0], w[n].x[2], w[n].rcen,
                 w[n].thetacen / RADIAN);
        fprintf (fptr, "%8.4e %8.4e %8.4e ", plasmamain[np].F_vis[0], plasmamain[np].F_vis[1], plasmamain[np].F_vis[2]);
        fprintf (fptr, "%8.4e %8.4e %8.4e ", plasmamain[np].F_UV[0], plasmamain[np].F_UV[1], plasmamain[np].F_UV[2]);
        fprintf (fptr, "%8.4e %8.4e %8.4e ", plasmamain[np].F_Xray[0], plasmamain[np].F_Xray[1], plasmamain[np].F_Xray[2]);


        fprintf (fptr, "\n");
      }

    }
    else
    {
      /* if we aren't inwind then print out a load of zeroes */

      /* printf("%i %i %i %i %i %8.4e %8.4e 0.0 0.0 0.0 0.0 0.0 0.0 0.0 \
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 \
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 \
         0.0 0.0 0.0 0.0 0.0 0.0 \
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n",
         n, np, w[n].inwind, ii, jj, w[n].x[0], w[n].x[2]);
       */

      if (ochoice)
      {
        fprintf (fptr, "%i %i %i %i %i %8.4e %8.4e 0.0 0.0 ", n, np, -2, ii, jj, w[n].x[0], w[n].x[2]);
        for (m = 0; m < geo.nxfreq; m++)
        {
          fprintf (fptr, "0.0 0.0 0.0 ");

        }
        fprintf (fptr, "\n");
      }
    }
  }


  if (ochoice)
  {
    fclose (fptr);
    printf ("\nSaved flux details\n");
  }

  return (0);


}
