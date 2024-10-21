/***********************************************************/
/** @file   paths.c
 * @author SWM
 * @date   July, 2015
 * @brief  Photon pathing functions.
 *
 * File containing photon path tracking functions.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/** @var *double reverb_path_bin
 * @brief	Array of path bin boundaries
 *
 * Set in reverb_init(), contains the boundaries for each 
 * path bin.
 *
 * ###Notes###
 * 10/15	-	Written by SWM
***********************************************************/
double *reverb_path_bin;

/**********************************************************/
/** 
 * @brief	Allocates the arrays for a path histogram
 *
 * @param [in,out] wind		Pointer to parent wind cell
 * @return 					Pointer to onstructed histogram
 *
 * Allocates path bins for a passed wind cell and returns a
 * pointer to the allocated space.
 *
 * ###Notes###
 * 9/3/15	-	Written by SWM
***********************************************************/
Wind_Paths_Ptr
wind_paths_constructor (WindPtr wind)
{
  Wind_Paths_Ptr paths = (Wind_Paths_Ptr) calloc (sizeof (wind_paths_dummy), 1);

  if (paths == NULL)
  {
    Error ("wind_paths_constructor: Could not allocate memory for cell %d\n", wind->nwind);
    Exit (0);
  }

  paths->ad_path_flux = (double *) calloc (sizeof (double), geo.reverb_path_bins);
  paths->ai_path_num = (int *) calloc (sizeof (int), geo.reverb_path_bins);
  paths->ad_path_flux_cent = (double *) calloc (sizeof (double), geo.reverb_path_bins);
  paths->ai_path_num_cent = (int *) calloc (sizeof (int), geo.reverb_path_bins);
  paths->ad_path_flux_disk = (double *) calloc (sizeof (double), geo.reverb_path_bins);
  paths->ai_path_num_disk = (int *) calloc (sizeof (int), geo.reverb_path_bins);
  paths->ad_path_flux_wind = (double *) calloc (sizeof (double), geo.reverb_path_bins);
  paths->ai_path_num_wind = (int *) calloc (sizeof (int), geo.reverb_path_bins);

  if (paths->ad_path_flux == NULL || paths->ai_path_num_wind == NULL)
  {
    Error ("wind_paths_constructor: Could not allocate memory for cell %d bins\n", wind->nwind);
    Exit (0);
  }
  return (paths);
}

/**********************************************************/
/** 
 * @brief	Initialises everything for reverb mapping
 *
 * @param [in,out] wind		Top-level wind cell array
 * @param [in] nangles		Number of observers
 * @return 					0
 *
 * Reports back to user the reverb mode running, and then
 * calls the appropriate setup routines.
 *
 * @see wind_paths_init()
 *
 * ###Notes###
 * 3/15	-	Written by SWM
***********************************************************/
int
reverb_init (WindPtr wind)
{
  char linelist[LINELENGTH];
  char xlinelist[LINELENGTH];
  int i, n;
  double x[3];
  double r_rad_min = 99e99, r_rad_max = 0.0, r_rad_min_log, r_rad_max_log, r_delta;

  if (geo.reverb == REV_WIND || geo.reverb == REV_MATOM)
  {                             //Initialise the arrays that handle photon path data and wind paths
    //Set the convergence fraction required before we start tracking photon paths for matoms
    geo.fraction_converged = 0.0;
    geo.reverb_fraction_converged = 0.90;

    if (geo.reverb_vis == REV_VIS_DUMP || geo.reverb_vis == REV_VIS_BOTH)
    {                           //If we're dumping out arrays of paths per cell
      for (i = 0; i < geo.reverb_dump_cells; i++)
      {                         //For each x/z position, find the n and put it in the array
        x[0] = geo.reverb_dump_cell_x[i];
        x[1] = 0.0;
        x[2] = geo.reverb_dump_cell_z[i];
        wind_x_to_n (x, &n);
        if (n < 0 || n > geo.ndim2)
        {
          Error ("reverb_init: Trying to dump info from a cell outside the grid\n");
        }
        geo.reverb_dump_cell[i] = n;
      }
    }

    for (i = 0; i < geo.ndomain; i++)
    {                           //Find the smallest and largest domain scales in the problem
      if (zdom[i].rmin < r_rad_min)
        r_rad_min = zdom[i].rmin;
      if (zdom[i].rmax > r_rad_max)
        r_rad_max = zdom[i].rmax;
    }

    //Allocate the array of bins
    reverb_path_bin = (double *) calloc (sizeof (double), geo.reverb_path_bins + 1);
    r_rad_min_log = log (r_rad_min);
    r_rad_max_log = log (r_rad_max * 10.0);
    r_delta = (r_rad_max_log - r_rad_min_log) / (double) geo.reverb_path_bins;
    for (i = 0; i <= geo.reverb_path_bins; i++)
    {                           //Set up the bounds for these bins
      reverb_path_bin[i] = exp (r_rad_min_log + i * r_delta);
    }

    wind_paths_init (wind);

    if (geo.reverb == REV_MATOM)
    {                           //If this is matom mode, detail the line numbers being tracked 
      sprintf (linelist, "reverb_init: Macro-atom line path tracking is enabled for lines %d", geo.reverb_line[0]);
      for (i = 1; i < geo.reverb_lines; i++)
      {
        sprintf (xlinelist, "%.100s, %d", linelist, geo.reverb_line[i]);
        sprintf (linelist, "%.200s", xlinelist);
      }
      Log ("%s\n", linelist);

      for (i = 0; i < nlines_macro; i++)
      {                         //Record the h-alpha line index for comparison purposes
        if (line[i].z == 1 && line[i].istate == 1 && line[i].levu == 3 && line[i].levl == 2)
        {
          geo.nres_halpha = line[i].where_in_list;
          break;
        }
      }
    }
    else if (geo.reverb == REV_WIND)
      Log ("reverb_init: Wind cell-based path tracking is enabled\n");
  }
  else if (geo.reverb == REV_PHOTON)
    Log ("reverb_init: Photon-based path tracking is enabled.\n");


  return (0);
}

/**********************************************************/
/** 
 * @brief	Initialises wind path structures
 *
 * @param [in,out] wind		Top-level wind cell array
 * @return 					0
 *
 * Iterates over each wind cell, and declares a single 
 * generic 'wind path' array for binning the paths of
 * incident photons in. Also declares a 'wind path' array
 * for each specific line of interest in matom mode.
 *
 * ###Notes###
 * 10/2/15	-	Written by SWM
***********************************************************/
int
wind_paths_init (WindPtr wind)
{
  int i, j;

  for (i = 0; i < geo.ndim2; i++)
  {                             //For each entry in the wind array
    wind[i].paths = (Wind_Paths_Ptr) wind_paths_constructor (&wind[i]);
    wind[i].line_paths = (Wind_Paths_Ptr *) calloc (sizeof (Wind_Paths_Ptr), geo.reverb_lines);
    for (j = 0; j < geo.reverb_lines; j++)
    {                           //For each line tracked on each cell
      wind[i].line_paths[j] = (Wind_Paths_Ptr) wind_paths_constructor (&wind[i]);
    }
  }
  return (0);
}

/****************************************************************/
/** 
 * @brief		Following a line emission, increments cell paths
 * 
 * @param [in, out] wind 		Wind cell to register photon in.
 * @param [in]		pp 			Photon to register in wind cell
 * @param [in]		nres		
 * @return 						0
 *  
 * When given a wind cell and photon info in matom mode, iterates 
 * over the list of tracked lines to see if the transition is in 
 * it. If so adds the weight to the path array for that line.
 *
 * ###Notes###
 * 27/2/15	-	 Written by SWM
*****************************************************************/
int
line_paths_add_phot (WindPtr wind, PhotPtr pp, int *nres)
{
  int i, j;

  if (geo.reverb_disk == REV_DISK_IGNORE && pp->origin_orig == PTYPE_DISK)
    return (0);
  if (*nres > nlines || *nres < 0)
    return (0);                 //This is a continuum photon

  for (i = 0; i < geo.reverb_lines; i++)
  {                             //Iterate over each tracked line
    if (lin_ptr[*nres]->where_in_list == geo.reverb_line[i])
    {                           //If the passed line exists within the tracked line array
      for (j = 0; j < geo.reverb_path_bins; j++)
      {                         //Iterate over the path bins
        if (pp->path >= reverb_path_bin[j] && pp->path <= reverb_path_bin[j + 1])
        {                       //If the photon's path lies in this bin's bounds, record it
          //printf("DEBUG: Added to line %d in cell %d - path %g weight %g\n",geo.reverb_line[i], wind->nwind, pp->path, pp->w);
          wind->line_paths[i]->ad_path_flux[j] += pp->w;
          wind->line_paths[i]->ai_path_num[j]++;
          switch (pp->origin)
          {
          case PTYPE_STAR:
          case PTYPE_AGN:
          case PTYPE_BL:
            wind->line_paths[i]->ad_path_flux_cent[j] += pp->w;
            wind->line_paths[i]->ai_path_num_cent[j]++;
            break;
          case PTYPE_DISK:
            wind->line_paths[i]->ad_path_flux_disk[j] += pp->w;
            wind->line_paths[i]->ai_path_num_disk[j]++;
            break;
          default:
            wind->line_paths[i]->ad_path_flux_wind[j] += pp->w;
            wind->line_paths[i]->ai_path_num_wind[j]++;
            break;
          }
          //Exit out of this loop
          return (0);
        }
      }
    }
  }
  return (0);
}

/****************************************************************/
/** 
 * @brief		Registers a photon passage with a wind cell.
 * 
 * @param [in, out] wind 	Wind cell to register photon in.
 * @param [in] pp			Photon to register in wind cell.
 * @return 					0
 *  
 * When given a wind cell and photon, adds the photon's weight to
 * the appropriate delay bin.
 *
 * ###Notes### 
 * 27/2/15 	-	Written by SWM 2/15.
*****************************************************************/
int
wind_paths_add_phot (WindPtr wind, PhotPtr pp)
{
  int i;
  if (geo.reverb_disk == REV_DISK_IGNORE && pp->origin_orig == PTYPE_DISK)
    return (0);

  for (i = 0; i < geo.reverb_path_bins; i++)
  {                             //For each bin
    if (pp->path >= reverb_path_bin[i] && pp->path <= reverb_path_bin[i + 1])
    {                           //If the path falls within its bounds, add photon weight
      wind->paths->ad_path_flux[i] += pp->w;
      wind->paths->ai_path_num[i]++;

      switch (pp->origin)
      {
      case PTYPE_STAR:
      case PTYPE_AGN:
      case PTYPE_BL:
        wind->paths->ad_path_flux_cent[i] += pp->w;
        wind->paths->ai_path_num_cent[i]++;
        break;
      case PTYPE_DISK:
        wind->paths->ad_path_flux_disk[i] += pp->w;
        wind->paths->ai_path_num_disk[i]++;
        break;
      default:
        wind->paths->ad_path_flux_wind[i] += pp->w;
        wind->paths->ai_path_num_wind[i]++;
        break;
      }
      return (0);
    }
  }
  return (0);
}

/**********************************************************/
/** 
 * @brief	Generates path for a 'wind' photon in photon mode
 *
 * @param [in,out] pp 	Photon to set path of
 * @return 				0
 *
 * Finds the straight-line distance between photon and the 
 * outer star radius, sets minimum path to that. Used in 
 * photon mode, and in all modes for non-wind starting paths.
 *
 * ###Notes###
 * 20/8/15	-	Written by SWM
***********************************************************/
int
simple_paths_gen_phot (PhotPtr pp)
{
  pp->path = length (pp->x);
  return (0);
}

/**********************************************************/
/** 
 * @brief	Draws a random path from a path histogram
 *
 * @param [in] PathPtr	Path histogram pointer
 * @return 				Path
 *
 * Picks a random path bin, weighted by the flux in each in
 * this cell, then assigns a path from within that bin 
 * (from a uniform random distribution)
 *
 * ###Notes###
 * 26/2/15	-	Written by SWM
 * 24/7/15	-	Removed frequency
***********************************************************/
double
r_draw_from_path_histogram (Wind_Paths_Ptr PathPtr)
{
  double r_rand, r_total, r_bin_min, r_bin_rand, r_path, r_bin_max;
  int i_path = -1;

  r_total = 0.0;
//  r_rand = PathPtr->d_flux * rand () / MAXRAND; DONE
  r_rand = PathPtr->d_flux * random_number (0.0, 1.0);
  i_path = -1;

  //printf("DEBUG: r_rand %g out of total %g\n",r_rand, PathPtr->d_flux);
  while (r_rand > r_total)
  {
    r_total += PathPtr->ad_path_flux[++i_path];
  }

  //Assign photon path to a random position within the bin.
  r_bin_min = reverb_path_bin[i_path - 1];
  r_bin_max = reverb_path_bin[i_path];
//  r_bin_rand = (rand () / MAXRAND) * (r_bin_max - r_bin_min); DONE
  r_bin_rand = random_number (0.0, 1.0) * (r_bin_max - r_bin_min);
  r_path = r_bin_min + r_bin_rand;
  return (r_path);
}


/**********************************************************/
/** 
 * @brief	Generates path for a wind photon
 *
 * @param [in] wind		Wind cell to spawn in
 * @param [in,out] pp 	Photon to set path of
 * @return 				0
 *
 * Picks a random path bin, weighted by the flux in each in
 * this cell, then assigns a path from within that bin 
 * (from a uniform random distribution)
 *
 * @see r_draw_from_path_histogram()
 * @see simple_paths_gen_phot()
 *
 * ###Notes###
 * 26/2/15	-	Written by SWM
 * 24/7/15	-	Removed frequency
***********************************************************/
int
wind_paths_gen_phot (WindPtr wind, PhotPtr pp)
{
  if (geo.ioniz_or_extract == CYCLE_IONIZ)
  {
    simple_paths_gen_phot (pp);
  }
  else if (wind->paths->i_num == 0)
  {                             //If there's no path data registered in this cell, default to simple
    Error ("wind_paths_gen_phot: No path data in cell %d at r=%g, z=%g\n",
           wind->nwind, sqrt (wind->x[0] * wind->x[0] + wind->x[1] * wind->x[1]), wind->x[2]);
    simple_paths_gen_phot (pp);
  }
  else
  {                             //Otherwise, draw a path for the photon from the cell's histogram
    pp->path = r_draw_from_path_histogram (wind->paths);
  }
  return (0);
}

/**********************************************************/
/** 
 * @brief	Generates path for a macro-atom line photon
 *
 * @param [in] wind			Wind cell to spawn in
 * @param [in,out] pp 		Photon to set path of
 * @patam [in] nres			Matom line to generate for
 * @return 					0
 *
 * If the line is being tracked, pick a random path from its 
 * path histogram. If not, default to the regular path 
 * histogram. 
 *
 * @see r_draw_from_path_histogram()
 * @see simple_paths_gen_phot()
 * @see wind_paths_gen_phot()
 *
 * ###Notes###
***********************************************************/
int
line_paths_gen_phot (WindPtr wind, PhotPtr pp, int nres)
{
  int i;
  if (geo.ioniz_or_extract == CYCLE_IONIZ)
  {
    simple_paths_gen_phot (pp);
  }
  else if (wind->paths->i_num == 0)
  {                             //If there's no path data registered in this cell, default to simple
    //Error ("line_paths_gen_phot: No path data in cell %d at r=%g, z=%g\n",
    //         wind->nwind, sqrt(wind->x[0]*wind->x[0] + wind->x[1]*wind->x[1]), wind->x[2]);
    simple_paths_gen_phot (pp);
  }
  else if (nres < 0 || nres >= nlines || lin_ptr[nres]->macro_info == FALSE)
  {                             //If this line is invalid, continuum or non-matom then default to wind
    pp->path = r_draw_from_path_histogram (wind->paths);
  }
  else
  {                             //Iterate over array to see if this line is tracked. If so, use that 
    //array index to find the path array for the line (i.e. the n^th line
    //being tracked is the n^th line path histogram).
    for (i = 0; i < geo.reverb_lines; i++)
    {
      if (lin_ptr[nres]->where_in_list == geo.reverb_line[i])
      {                         //Line identified using its position in nres as unique ID
        if (wind->line_paths[i]->i_num > 0)
        {                       //If there photons recorded in this histogram
          pp->path = r_draw_from_path_histogram (wind->line_paths[i]);
        }
        else
        {                       //If there are no photons in this histogram, log and default
          //to using the wind path histogram.
          //Error("line_paths_gen_phot: No path data for line %d in cell %d at r=%g, z=%g\n",
          // wind->nwind, nres, sqrt(wind->x[0]*wind->x[0] + wind->x[1]*wind->x[1]), wind->x[2]);
          pp->path = r_draw_from_path_histogram (wind->paths);
        }
        return (0);
      }
    }
    //If the line isn't being tracked, default to wind
    pp->path = r_draw_from_path_histogram (wind->paths);
  }
  return (0);
}


/****************************************************************/
/** 
 * @brief	Evaluates individual wind cell paths
 *
 * @param [in,out] wind		Wind cell to evaluate
 *
 * Records the total flux in the cell, as well as making a 
 * simple 'average path' calculation.
 *
 * @see wind_paths_evaluate()
 *
 * ###Notes###
 * 26/2/15	-	Written by SWM
 * 24/7/15	-	Removed frequency
*****************************************************************/
int
wind_paths_evaluate_single (Wind_Paths_Ptr paths)
{
  int i;
  paths->d_flux = 0.0;
  paths->d_path = 0.0;
  paths->i_num = 0;

  for (i = 0; i < geo.reverb_path_bins; i++)
  {                             //For each path bin, add its contribution to total flux & avg path
    paths->d_flux += paths->ad_path_flux[i];
    paths->i_num += paths->ai_path_num[i];
    paths->d_path += paths->ad_path_flux[i] * (reverb_path_bin[i] + reverb_path_bin[i + 1]) / 2.0;
  }

  //If there was any data in this cell, calculate avg. path
  if (paths->i_num > 0)
    paths->d_path /= paths->d_flux;

  return (0);
}


/****************************************************************/
/** 
 * @brief	Evaluates wind path details for a cycle
 *
 * @param [in,out] wind	Wind to evaluate
 * @return 				0
 *
 * Iterates over each cell in the wind.
 * 
 * @see wind_paths_evaluate_single()
 *
 * ###Notes###
 * 26/2/15	-	Written by SWM
 * 24/7/15	-	Removed frequency
*****************************************************************/
int
wind_paths_evaluate (WindPtr wind, int i_rank)
{
  int i, j;
  for (i = 0; i < geo.ndim2; i++)
  {                             //For each cell in the wind
    if (wind[i].inwind >= 0)
    {                           //If this is a wind cel;, evaluate each of the path histograms
      wind_paths_evaluate_single (wind[i].paths);
      for (j = 0; j < geo.reverb_lines; j++)
      {
        wind_paths_evaluate_single (wind[i].line_paths[j]);
      }
    }
  }

  if (geo.reverb_vis == REV_VIS_VTK || geo.reverb_vis == REV_VIS_BOTH)
  {                             //If we're visualising the delays in 3d
    Log ("wind_paths_evaluate: Outputting %d .vtk visualisations\n", geo.ndomain);
    for (i = 0; i < geo.ndomain; i++)
    {                           //For each domain, dump a vtk file of the grid & info
      wind_paths_output_vtk (wind, i);
    }
  }
  if (geo.reverb_vis == REV_VIS_DUMP || geo.reverb_vis == REV_VIS_BOTH)
  {                             //Dump the path delay information for certain tracked cells to file
    wind_paths_output_dump (wind, i_rank);
  }

  Log (" Completed evaluating wind path arrays.");

  return (0);
}


/****************************************************************/
/** 
 * @brief	Dumps wind path arrays for a wind cell
 *
 * @param [in] wind		Wind cell to dump
 * @return 				0
 *
 * Outputs a file "root.wind_paths_[nwind].csv" containing the path
 * length bins and total flux in each for each line histogram (and
 * the regular wind histogram) in the given cell.
 *
 * ###Notes###
 * 10/15	-	Written by SWM
*****************************************************************/
int
wind_paths_dump (WindPtr wind, int rank_global)
{
  FILE *fopen (), *fptr;
  char c_file[LINELENGTH];
  int j, k;

  //Setup file name and open the file
  sprintf (c_file, "%.100s.wind_paths_%d.%d.csv", files.root, wind->nwind, rank_global);
  fptr = fopen (c_file, "w");

  //Print out metadata header specifying the domain and position
  //As nwind is unique but not descriptive, and X:Y:Domain would be too long
  fprintf (fptr, "'Wind domain', %d,'X',%g,'Y',%g\n", wind->ndom, wind->xcen[0], wind->xcen[2]);

  //Now print out the header for the actual path data for this wind cell 
  fprintf (fptr, "'Path Bin', 'Heating', 'Heating Central', 'Heating Disk', 'Heating Wind'");
  for (j = 0; j < geo.reverb_lines; j++)
  {
    fprintf (fptr, ", 'Line %d', 'Line %d Central', 'Line %d Disk', 'Line %d Wind'",
             geo.reverb_line[j], geo.reverb_line[j], geo.reverb_line[j], geo.reverb_line[j]);
  }
  fprintf (fptr, "\n");

  for (k = 0; k < geo.reverb_path_bins; k++)
  {                             //For each path bin, print the 'wind' weight 
    fprintf (fptr, "%g, %g, %g, %g, %g", reverb_path_bin[k],
             wind->paths->ad_path_flux[k],
             wind->paths->ad_path_flux_cent[k], wind->paths->ad_path_flux_disk[k], wind->paths->ad_path_flux_wind[k]);

    for (j = 0; j < geo.reverb_lines; j++)
    {                           //For each tracked line, print the weight in this bin
      fprintf (fptr, ", %g, %g, %g, %g",
               wind->line_paths[j]->ad_path_flux[k],
               wind->line_paths[j]->ad_path_flux_cent[k],
               wind->line_paths[j]->ad_path_flux_disk[k], wind->line_paths[j]->ad_path_flux_wind[k]);
    }
    fprintf (fptr, "\n");
  }
  fclose (fptr);
  return (0);
}


/****************************************************************/
/** 
 * @brief	Iterates through the wind, dumping cells of interest
 *
 * @param [in] wind		Wind array to dump from
 * @return 				0
 *
 * Outputs a file "root.ngrid.csv" containing the path length bins
 * and total flux in each for each line histogram (and the regular
 * wind histogram) in the given cell.
 *
 * ###Notes###
 * 10/15	-	Written by SWM
*****************************************************************/
int
wind_paths_output_dump (WindPtr wind, int i_rank)
{
  int i, d, n;
  double x[3];
  for (i = 0; i < geo.reverb_dump_cells; i++)
  {                             //For each location we want to dump the details for
    x[0] = geo.reverb_dump_cell_x[i];
    x[1] = 0.0;
    x[2] = geo.reverb_dump_cell_z[i];

    for (d = 0; d < geo.ndomain; d++)
    {                           //For each domain, check if this position is within it
      n = where_in_grid (d, x);
      if (n >= 0)
      {                         //If it is, then dump the delay information 
        wind_paths_dump (&wind[n], i_rank);
      }
    }
  }
  return (0);
}


/****************************************************************/
/** 
 * @brief		Given trz index in wind, returns vtk data index.
 * 
 * @param [in] i 	Theta index of cell.
 * @param [in] j 	Radius index of cell.
 * @param [in] k 	Height index of cell.
 * @param [in] i_top Whether the point is above or below the disk
 * @return			Index for vtk poly.
 * 
 * When given the position of a cell in the wind, returns the
 * corresponding index in the flat VTK data arrays. RZ uses the
 * uses the standard wind location, theta is set as defined in
 * geo. Used in wind_paths_output() only.
 *
 * ###Notes### Written by SWM 4/15.
 *****************************************************************/
int
wind_paths_point_index (int i, int j, int k, int i_top, DomainPtr dom)
{
  int n;
  n = i * 2 * (geo.reverb_angle_bins + 1) * dom->ndim + j * 2 * (geo.reverb_angle_bins + 1) + k * 2 + i_top;
  return (n);
}

/****************************************************************/
/** 
 * @brief		Given rtt index in wind, returns vtk data index.
 * 
 * @param [in] i 	Theta index of cell.
 * @param [in] j 	Radius index of cell.
 * @param [in] k 	Height index of cell.
 * @return			Index for vtk poly.
 * 
 * When given the position of a cell in the wind, returns the
 * corresponding index in the flat VTK data arrays. RZ uses the
 * uses the standard wind location, theta is set as defined in
 * geo. Used in wind_paths_output() only.
 *
 * ###Notes### Written by SWM 4/15.
 *****************************************************************/
int
wind_paths_sphere_point_index (int i, int j, int k)
{
  int n;
  n = i * (geo.reverb_angle_bins + 1) * (geo.reverb_angle_bins + 1) + j * (geo.reverb_angle_bins + 1) + k;
  return (n);
}


/****************************************************************/
/** 
 * @brief		Outputs wind path information to vtk.
 * 
 * @param [in] wind 		Pointer to wind array.
 * @param [in] c_file_in	Name of input file.
 * @return 					0
 *  
 * When given a wind containing position and delay map information
 * generated using REVERB_WIND, outputs a 3d model of the wind to
 * file in ASCII .vtk format.
 *
 * ###Notes### Written by SWM 4/15.
*****************************************************************/
int
wind_paths_output_vtk (WindPtr wind, int ndom)
{
  FILE *fopen (), *fptr;
  char c_file[LINELENGTH];
  int i, j, k, n, i_obs, i_cells, i_points;
  double r_azi, r_inc, r_x, r_y, r_z, r_err;
  PhotPtr p_test = calloc (sizeof (p_dummy), 1);
  DomainPtr dom;

  //Get output filename
  sprintf (c_file, "%.100s.%d.wind_paths.vtk", files.root, ndom);

  if ((fptr = fopen (c_file, "w")) == NULL)
  {                             //If this file can't be opened, error out
    Error ("wind_paths_output_vtk: Unable to open %s for writing\n", c_file);
    Exit (0);
  }
  Log ("Outputting wind path information to file '%s'.\n", c_file);

  dom = &zdom[ndom];

  //Setup the number of vertexes and cells within this mesh
  //Notably, we skip the outside cell as it's empty
  //For the spherical case, we generate theta and thi bins using the same resolution
  switch (dom->coord_type)
  {
  case SPHERICAL:
    i_cells = 1 * (dom->ndim - 1) * geo.reverb_angle_bins * geo.reverb_angle_bins;
    i_points = 1 * dom->ndim * (geo.reverb_angle_bins + 1) * (geo.reverb_angle_bins + 1);
    break;
  case CYLIND:
    i_cells = 2 * (dom->ndim - 1) * (dom->mdim - 1) * geo.reverb_angle_bins;
    i_points = 2 * dom->ndim * dom->mdim * (geo.reverb_angle_bins + 1);
    break;
  default:
    free (p_test);
    fclose (fptr);
    Error ("wind_paths_output_vtk: Mesh format not yet supported (CYLVAR/RTHETA)");
    return (0);
  }

  //Write out header
  fprintf (fptr, "# vtk DataFile Version 2.0\n");
  fprintf (fptr, "Wind file data\nASCII\n");

  //Write out positions of corners of each wind cell as vertexes
  fprintf (fptr, "DATASET UNSTRUCTURED_GRID\n");
  fprintf (fptr, "POINTS %d float\n", i_points);
  if (dom->coord_type == SPHERICAL)
  {
    for (i = 0; i < dom->ndim; i++)
    {
      n = dom->nstart + i;
      for (j = 0; j <= geo.reverb_angle_bins; j++)
      {
        r_inc = j * (PI / (double) geo.reverb_angle_bins) - PI / 2.0;

        for (k = 0; k <= geo.reverb_angle_bins; k++)
        {
          r_azi = k * (PI / (double) geo.reverb_angle_bins);
          r_x = wind[n].x[0] * sin (r_inc) * cos (r_azi);
          r_y = wind[n].x[0] * sin (r_inc) * sin (r_azi);
          r_z = wind[n].x[0] * cos (r_inc);
          fprintf (fptr, "%10.5g %10.5g %10.5g\n", r_x, r_y, r_z);
        }
      }
    }
  }
  else
  {
    for (i = 0; i < dom->ndim; i++)
    {
      for (j = 0; j < dom->mdim; j++)
      {
        wind_ij_to_n (ndom, i, j, &n);
        for (k = 0; k <= geo.reverb_angle_bins; k++)
        {
          r_azi = k * (PI / (double) geo.reverb_angle_bins);
          r_x = wind[n].x[0] * cos (r_azi);
          r_y = wind[n].x[0] * sin (r_azi);
          fprintf (fptr, "%10.5g %10.5g %10.5g\n", r_x, r_y, wind[n].x[2]);
          fprintf (fptr, "%10.5g %10.5g %10.5g\n", r_x, r_y, -wind[n].x[2]);
        }
      }
    }
  }
  fprintf (fptr, "\n");

  //Write out the vertexes comprising each cell
  fprintf (fptr, "CELLS %d %d\n", i_cells, 9 * i_cells);
  if (dom->coord_type == SPHERICAL)
  {
    for (i = 0; i < dom->ndim - 1; i++)
    {
      n = dom->nstart + i;
      for (j = 0; j < geo.reverb_angle_bins; j++)
      {
        r_inc = j * (PI / (double) geo.reverb_angle_bins) - PI / 2.0;

        for (k = 0; k < geo.reverb_angle_bins; k++)
        {
          fprintf (fptr, "8 %d %d %d %d %d %d %d %d\n",
                   wind_paths_sphere_point_index (i, j, k),
                   wind_paths_sphere_point_index (i, j, k + 1),
                   wind_paths_sphere_point_index (i, j + 1, k + 1),
                   wind_paths_sphere_point_index (i, j + 1, k),
                   wind_paths_sphere_point_index (i + 1, j, k),
                   wind_paths_sphere_point_index (i + 1, j, k + 1),
                   wind_paths_sphere_point_index (i + 1, j + 1, k + 1), wind_paths_sphere_point_index (i + 1, j + 1, k));
        }
      }
    }
  }
  else
  {
    for (i = 0; i < dom->ndim - 1; i++)
    {
      for (j = 0; j < dom->mdim - 1; j++)
      {
        for (k = 0; k < geo.reverb_angle_bins; k++)
        {
          fprintf (fptr, "8 %d %d %d %d %d %d %d %d\n",
                   wind_paths_point_index (i, j, k, 1, dom),
                   wind_paths_point_index (i, j, k + 1, 1, dom),
                   wind_paths_point_index (i, j + 1, k + 1, 1, dom),
                   wind_paths_point_index (i, j + 1, k, 1, dom),
                   wind_paths_point_index (i + 1, j, k, 1, dom),
                   wind_paths_point_index (i + 1, j, k + 1, 1, dom),
                   wind_paths_point_index (i + 1, j + 1, k + 1, 1, dom), wind_paths_point_index (i + 1, j + 1, k, 1, dom));
          fprintf (fptr, "8 %d %d %d %d %d %d %d %d\n",
                   wind_paths_point_index (i, j, k, 0, dom),
                   wind_paths_point_index (i, j, k + 1, 0, dom),
                   wind_paths_point_index (i, j + 1, k + 1, 0, dom),
                   wind_paths_point_index (i, j + 1, k, 0, dom),
                   wind_paths_point_index (i + 1, j, k, 0, dom),
                   wind_paths_point_index (i + 1, j, k + 1, 0, dom),
                   wind_paths_point_index (i + 1, j + 1, k + 1, 0, dom), wind_paths_point_index (i + 1, j + 1, k, 0, dom));
        }
      }
    }
  }

  //Write the type for each cell (would be unnecessary if STRUCTURED_GRID used)
  //But this allows for more flexible expansion
  fprintf (fptr, "CELL_TYPES %d\n", i_cells);
  for (i = 0; i < i_cells; i++)
    fprintf (fptr, "12\n");
  fprintf (fptr, "\n");

  //Write out the arrays containing the various properties in the appropriate order
  fprintf (fptr, "CELL_DATA %d\n", i_cells);

  fprintf (fptr, "SCALARS phot_count float 1\n");
  fprintf (fptr, "LOOKUP_TABLE default\n");
  if (dom->coord_type == SPHERICAL)
  {
    for (i = 0; i < dom->ndim - 1; i++)
    {
      n = dom->nstart + i;
      for (j = 0; j < geo.reverb_angle_bins; j++)
      {
        for (k = 0; k < geo.reverb_angle_bins; k++)
        {
          fprintf (fptr, "%d\n", wind[n].paths->i_num);
        }
      }
    }
  }
  else
  {
    for (i = 0; i < dom->ndim - 1; i++)
    {
      for (j = 0; j < dom->mdim - 1; j++)
      {
        wind_ij_to_n (ndom, i, j, &n);
        for (k = 0; k < geo.reverb_angle_bins; k++)
        {
          fprintf (fptr, "%d\n", wind[n].paths->i_num);
          fprintf (fptr, "%d\n", wind[n].paths->i_num);
        }
      }
    }
  }

  //Use statistical error
  fprintf (fptr, "SCALARS path_errors float 1\n");
  fprintf (fptr, "LOOKUP_TABLE default\n");
  if (dom->coord_type == SPHERICAL)
  {
    for (i = 0; i < dom->ndim - 1; i++)
    {
      n = dom->nstart + i;
      for (j = 0; j < geo.reverb_angle_bins; j++)
      {
        for (k = 0; k < geo.reverb_angle_bins; k++)
        {
          if (wind[n].paths->i_num > 0)
          {
            r_err = sqrt ((double) wind[n].paths->i_num) / (double) wind[n].paths->i_num;
            fprintf (fptr, "%g\n", r_err);
          }
          else
          {
            fprintf (fptr, "-1\n");
          }
        }
      }
    }
  }
  else
  {
    for (i = 0; i < dom->ndim - 1; i++)
    {
      for (j = 0; j < dom->mdim - 1; j++)
      {
        wind_ij_to_n (ndom, i, j, &n);
        for (k = 0; k < geo.reverb_angle_bins; k++)
        {
          if (wind[n].paths->i_num > 0)
          {
            r_err = sqrt ((double) wind[n].paths->i_num) / (double) wind[n].paths->i_num;
            fprintf (fptr, "%g\n", r_err);
            fprintf (fptr, "%g\n", r_err);
          }
          else
          {
            fprintf (fptr, "-1\n");
            fprintf (fptr, "-1\n");
          }

        }
      }
    }
  }

  fprintf (fptr, "SCALARS path_rel_diff_from_direct float 1\n");
  fprintf (fptr, "LOOKUP_TABLE default\n");
  if (dom->coord_type == SPHERICAL)
  {
    for (i = 0; i < dom->ndim - 1; i++)
    {
      n = dom->nstart + i;
      for (j = 0; j < geo.reverb_angle_bins; j++)
      {
        for (k = 0; k < geo.reverb_angle_bins; k++)
        {
          if (wind[n].paths->i_num > 0)
          {
            double f_diff;
            f_diff = wind[n].paths->d_path;
            f_diff -= (sqrt (wind[n].xcen[0] * wind[n].xcen[0] +
                             wind[n].xcen[1] * wind[n].xcen[1] + wind[n].xcen[2] * wind[n].xcen[2]) - geo.rstar);
            f_diff = fabs (f_diff);
            f_diff /= wind[n].paths->d_path;

            fprintf (fptr, "%g\n", f_diff);
          }
          else
          {
            fprintf (fptr, "-1\n");
          }
        }
      }
    }
  }
  else
  {
    for (i = 0; i < dom->ndim - 1; i++)
    {
      for (j = 0; j < dom->mdim - 1; j++)
      {
        wind_ij_to_n (ndom, i, j, &n);
        for (k = 0; k < geo.reverb_angle_bins; k++)
        {
          if (wind[n].paths->i_num > 0)
          {
            double f_diff;
            f_diff = wind[n].paths->d_path;
            f_diff -= (sqrt (wind[n].xcen[0] * wind[n].xcen[0] +
                             wind[n].xcen[1] * wind[n].xcen[1] + wind[n].xcen[2] * wind[n].xcen[2]) - geo.rstar);
            f_diff = fabs (f_diff);
            f_diff /= wind[n].paths->d_path;

            fprintf (fptr, "%g\n", f_diff);
            fprintf (fptr, "%g\n", f_diff);
          }
          else
          {
            fprintf (fptr, "-1\n");
            fprintf (fptr, "-1\n");
          }
        }
      }
    }
  }

  fprintf (fptr, "SCALARS path_average float 1\n");
  fprintf (fptr, "LOOKUP_TABLE default\n");
  if (dom->coord_type == SPHERICAL)
  {
    for (i = 0; i < dom->ndim - 1; i++)
    {
      n = dom->nstart + i;
      for (j = 0; j < geo.reverb_angle_bins; j++)
      {
        for (k = 0; k < geo.reverb_angle_bins; k++)
        {
          if (wind[n].paths->i_num > 0)
          {
            fprintf (fptr, "%g\n", wind[n].paths->d_path);
          }
          else
          {
            fprintf (fptr, "-1\n");
          }
        }
      }
    }
  }
  else
  {
    for (i = 0; i < dom->ndim - 1; i++)
    {
      for (j = 0; j < dom->mdim - 1; j++)
      {
        wind_ij_to_n (ndom, i, j, &n);
        for (k = 0; k < geo.reverb_angle_bins; k++)
        {
          if (wind[n].paths->i_num > 0)
          {
            fprintf (fptr, "%g\n", wind[n].paths->d_path);
            fprintf (fptr, "%g\n", wind[n].paths->d_path);
          }
          else
          {
            fprintf (fptr, "-1\n");
            fprintf (fptr, "-1\n");
          }
        }
      }
    }
  }


  for (i_obs = 0; i_obs < geo.nangles; i_obs++)
  {
    fprintf (fptr, "SCALARS path_%s float 1\n", xxspec[MSPEC + i_obs].name);
    fprintf (fptr, "LOOKUP_TABLE default\n");
    stuff_v (xxspec[MSPEC + i_obs].lmn, p_test->lmn);

    if (dom->coord_type == SPHERICAL)
    {
      for (i = 0; i < dom->ndim - 1; i++)
      {
        n = dom->nstart + i;
        for (j = 0; j < geo.reverb_angle_bins; j++)
        {
          r_inc = ((double) j + 0.5) * (PI / (double) geo.reverb_angle_bins) - PI / 2.0;

          for (k = 0; k < geo.reverb_angle_bins; k++)
          {
            if (wind[n].paths->i_num > 0)
            {
              r_azi = ((double) k + 0.5) * (PI / (double) geo.reverb_angle_bins);
              p_test->x[0] = wind[n].xcen[0] * sin (r_inc) * cos (r_azi);
              p_test->x[1] = wind[n].xcen[0] * sin (r_inc) * sin (r_azi);
              p_test->x[2] = wind[n].xcen[0] * cos (r_inc);
              fprintf (fptr, "%g\n", wind[n].paths->d_path + delay_to_observer (p_test));
            }
            else
            {
              fprintf (fptr, "-1\n");
            }
          }
        }
      }
    }
    else
    {
      for (i = 0; i < dom->ndim - 1; i++)
      {
        for (j = 0; j < dom->mdim - 1; j++)
        {
          wind_ij_to_n (ndom, i, j, &n);
          for (k = 0; k < geo.reverb_angle_bins; k++)
          {
            if (wind[n].paths->i_num > 0)
            {
              r_azi = ((double) k + 0.5) * (PI / (double) geo.reverb_angle_bins);
              p_test->x[0] = wind[n].xcen[0] * cos (r_azi);
              p_test->x[1] = wind[n].xcen[0] * sin (r_azi);
              p_test->x[2] = wind[n].xcen[2];
              fprintf (fptr, "%g\n", wind[n].paths->d_path + delay_to_observer (p_test));
              p_test->x[2] = -wind[n].xcen[2];
              fprintf (fptr, "%g\n", wind[n].paths->d_path + delay_to_observer (p_test));
            }
            else
            {
              fprintf (fptr, "-1\n");
              fprintf (fptr, "-1\n");
            }
          }
        }
      }
    }
  }

  free (p_test);
  fclose (fptr);
  return (0);
}
