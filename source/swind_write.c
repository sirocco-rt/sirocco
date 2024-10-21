
/***********************************************************/
/** @file  swind_write.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  
 * These are subroutines used by swind to write ascii
 * files of various variables for use with smongo or normal display
 * software.
 *
 ***********************************************************/

#include    <math.h>
#include    <stdio.h>
#include    <strings.h>
#include	<string.h>
#include	<stdlib.h>

#include	"atomic.h"
#include	"sirocco.h"
#include    "templates.h"

#define ODIM 256
int init_write_array = 1;
float aout[ODIM][ODIM];


/**********************************************************/
/** 
 * @brief      writes out the arrays in a way which can be read easily
 *    into a routine which can do contouring 
 *
 * @param [in] char  filename[]   The name of the output file
 * @param [in] int  choice   0-> don't write array, 1 -> write without intepolation, 2-> write with interpolation
 * @return     Always returns 0
 *
 * @details
 * This is a general purpose routine used by swind to write 
 * ascii data files.  If a file is written, it can be written
 * at the positions of of the windsave grid, or interpolated
 * onto a linear coordinate grid.  The latter was intended
 * to make it easier to make an "image" of the output file.  
 *
 * ### Notes ###
 * The output files are written in a format that should be
 * readable with the astropy ascii tables routines.
 *
 **********************************************************/

int
write_array (filename, choice)
     char filename[];
     int choice;
{
  //Dynamical allocation is allowed, although I generally avoid it -- 05apr ksl
  float r, z;
  float rmin, rmax, zmin, zmax;
  int ii, jj;
  FILE *fopen (), *fptr;
  char outfile[LINELENGTH];
  char extra[LINELENGTH];

  double length ();
  double xx[3];
  int i;
  int nn, nnn[4], nelem;
  double frac[4];
  int ndom, ndim, mdim, nstart, ndomain;
  ndom = current_domain;
  ndim = zdom[ndom].ndim;
  mdim = zdom[ndom].mdim;
  nstart = zdom[ndom].nstart;
//OLD  float rin[ndim], zin[mdim];



  if (choice == 0)              // Then we don't really want output files
  {
    return (0);
  }

// Open the appropriate file

  strcpy (outfile, filename);
  if (geo.ndomain > 1)
  {
    sprintf (extra, ".%d", current_domain);
    strcat (outfile, extra);
  }
  strcat (outfile, ".dat");     // Add standard extension to filenames, i.e. the one used by tecplot
  fptr = fopen (outfile, "w");


  /* Write out the header information for the file */
  fprintf (fptr, "# TITLE= \"%s\"\n", outfile);
  if (zdom[ndom].coord_type == SPHERICAL)
  {
    fprintf (fptr, "# Coord_Sys SPHERICAL\n");
  }
  else if (zdom[ndom].coord_type == CYLIND)
  {
    fprintf (fptr, "# Coord_Sys CYLIND\n");
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {
    fprintf (fptr, "# Coord_Sys RTHETA\n");
  }
  else if (zdom[ndom].coord_type == CYLVAR)
  {
    fprintf (fptr, "# Coord_Sys CYLVAR\n");
  }
  else
  {
    Error ("write_array: Unknown coordinaate system type: %d\n", zdom[ndom].coord_type);
  }


// Put the r and z coord. grid  into easier to understand arrays
//OLD  for (i = 0; i < ndim; i++)
//OLD  {
//OLD    zin[i] = wmain[nstart + i].x[2];
//OLD    rin[i] = wmain[nstart + i * mdim].x[0];
//OLD  }


/* In swind the filenames are set to start with z if the output coordinates
are linear, and x otherwise.  This is not particularly transparent ?? ksl */

  if (choice == 1)
  {                             // print out the original array elements

    if (zdom[ndom].coord_type == SPHERICAL)
    {

      /* JM 1411 -- added header for reading with astropy ascii module */
      fprintf (fptr, "r var inwind i\n");

      for (i = 0; i < ndim; i++)
      {
        fprintf (fptr, "%8.2e %8.2e %3d %3d \n", wmain[nstart + i].r, aaa[nstart + i], wmain[nstart + i].inwind, i);
      }
    }
    else
    {

      /* JM 1411 -- added header for reading with astropy ascii module */
      fprintf (fptr, "x z var inwind i j\n");

      for (i = 0; i < ndim * mdim; i++)
      {
        wind_n_to_ij (ndom, nstart + i, &ii, &jj);
        fprintf (fptr, "%8.4e %8.4e %8.5e %3d %3d %3d\n",
                 wmain[nstart + i].xcen[0], wmain[nstart + i].xcen[2], aaa[nstart + i], wmain[nstart + i].inwind, ii, jj);

      }
    }

  }
  else if (choice == 2)
  {                             // Then regrid to a linear grid
/* Wind array is organized so that it increments up along z and then over */

/* Want linearly spaced output array */
    rmin = zmin = 0;
    rmax = zmax = geo.rmax;

    xx[1] = 0.0;
    for (ii = 0; ii < ODIM; ii++)
    {
      xx[0] = r = rmin + ii * (rmax - rmin) / (ODIM - 1);
      for (jj = 0; jj < ODIM; jj++)
      {
        xx[2] = z = zmin + jj * (zmax - zmin) / (ODIM - 1);
        if (where_in_wind (xx, &ndomain) == W_ALL_INWIND)
        {                       // Then the position is in the wind region
          coord_fraction (ndom, 0, xx, nnn, frac, &nelem);
          for (nn = 0; nn < nelem; nn++)
          {
            aout[ii][jj] += aaa[nnn[nn]] * frac[nn];
          }

        }
      }
    }

    /* Now print out the array */
    // fprintf (fptr, "VARIABLES= \"X\" \"Z\" \"Var1\" \n");
    // fprintf (fptr, "ZONE I=%d J=%d DATAPACKING=POINT\n", ODIM, ODIM);
    fprintf (fptr, "# Resampled outputs\n");
    for (jj = 0; jj < ODIM; jj++)
    {
      /* JM 1411 -- added header for reading with astropy ascii module */
      fprintf (fptr, "r z var\n");

      z = zmin + (zmax - zmin) * jj / (ODIM - 1);
      for (ii = 0; ii < ODIM; ii++)
      {
        r = rmin + (rmax - rmin) * ii / (ODIM - 1);
        fprintf (fptr, "%8.2e %8.2e %8.2e\n", r, z, aout[ii][jj]);
      }
    }
    /*Finished printing out the array */
  }
  else
  {
    Error ("write_array: Unknown choice %d\n", choice);
  }

  fclose (fptr);
  return (0);
}





/**********************************************************/
/** 
 * @brief      Print results to the screen
 *
 * @param [in] char  name[]   The name associated with the array to be printed
 * @return     Always returns 0
 *
 * @details
 * This is a generalized display routine that outputs the results of a query
 * to the screen for inspection by the user.
 *
 * ### Notes ###
 * For cylindrical coordinates, one goes up in z, and over in x. To cylindrical coordinates
 * 	we want to display so that each colum represents a constant line on the axis
 * 	w[mdim][ndim].  So for a system with mdim=20 and ndim=25.  There are 25 elements
 * 	in the z axis and 20 elements in the xaxis
 * 
 * 	Therefore a row, constant z, is displayed by by incrementing by 20 or mdim
 * 
 * 	For rtheta coordinates, one goes around in theta, up in r.  The fasted moving
 * 	coordinate is r (when thinking of 1-d versions of the array.  Unless we
 * 	are going to write a separate routine, the simplest thing to do
 * 	is to make each column represent a constant r
 * 
 * 	w[mdim][ndim]  So, in spherical polar coorcinates, a system with mdim 20 has
 * 	20 angles, and 25 radii.  
 * 
 * 	So for spherical polar, constant r is displayed by incrementing by 1.  As a result
 * 	it is unclear that one can easily use the same routine for the two situations
 * 	since you seem to be incrementing the opposite axes, since in the one case one 
 * 	wants MDIM rows and other case one wants NDIM rows.  
 * 
 * 	It would be possible if you plot theta lines in each row, but what this means
 * 	is that the first row is closest to the z axis
 *
 **********************************************************/

int
display (name)
     char name[];
{
  int i, j, n;
//OLD  int ndom, ndim, mdim, nstart;
  int ndom, mdim, nstart;

  ndom = current_domain;
//OLD  ndim = zdom[ndom].ndim;
  mdim = zdom[ndom].mdim;
  nstart = zdom[ndom].nstart;

  Log ("Check me now %d %d\n", swind_min, swind_max);


  Log ("\n %s \n", name);
  Log ("z/theta \\x/r");
  for (i = swind_min; i < swind_max; i += swind_delta)
    if (zdom[ndom].coord_type == 1 || swind_project == 1)
      Log ("%8.2e ", wmain[nstart + i * mdim].x[0]);
    else
      Log ("%8.2e ", wmain[nstart + i * mdim].rcen);

  Log ("\n");

  for (j = 0; j < mdim; j++)
  {
    if (zdom[ndom].coord_type == 1 || swind_project == 1)
      Log ("%8.2e ", wmain[nstart + j].x[2]);
    else
      Log ("%8.2e ", wmain[nstart + j].thetacen);

    for (i = swind_min; i < swind_max; i += swind_delta)
    {
      n = nstart + i * mdim + j;
      Log ("%8.2g ", aaa[n]);
    }
    Log ("\n");
  }

  return (0);
}
