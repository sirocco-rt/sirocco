/***********************************************************/
/** @file  setup_reverb.c
 * @author swm
 * @date   January, 2018
 *
 * @brief  Routines for reading in the settings for reverb
 *
 * Contains the function for reading in reverberation mapping
 * and visualisation settings.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief      Gets parameters from input file for reverb
 * and visualisation settings.
 *
 * @return    0
 *
 * @details
 * Reads in reverberation mapping and basic visualisation for
 * reverberation mapping settings from file. These settings
 * are fully-documented in the Sphinx docs, and briefly here:
 *
 * ### reverb.type ###
 * Sets whether or not to do reverb mapping, and if so how to
 * assign photon starting paths for non-CO photons.
 *   'none': Off
 *   'photon': starting paths set by distance to the CO.
 *   'wind': photons generated in the wind assigned starting
 *     paths from the distribution of paths heating the wind
 *     cell they were spawned in.
 *   'matom': as wind but photons generated in a line assigned
 *     starting paths from the distribution of paths of photons
 *     that de-excited into that line. For lines in matom_lines.
 *
 * ### reverb.matom_lines, reverb.matom_line ###
 * The number of macro-atom lines to track (above). Internal
 * line number! TODO: Convert to elem:ion:lvlu:lvll format.
 *
 * ### reverb.disk_type ###
 * How the starting paths of photons from the disk are assigned.
 *   'correlated': Set by distance to the CO.
 *   'uncorrelated': Set to 0. Not recommended.
 *   'ignored': Disk photons do not contribute to 'wind' and
 *     'matom' distributions.
 *
 * ### reverb.path_bins ###
 * How many bins to store wind & matom path distributions in.
 * Typically 1000.
 *
 * ### reverb.visualisation ###
 *   'none': No visualisation
 *   'vtk': Output a .vtk file showing the mean paths in each cell,
 *     for visualising in something like VisIt or Paraview.
 *   'dump': Output a series of files dumping the path distributions
 *     for cells set by reverb.dump_cells.
 *   'both': Output both.
 *
 * ### reverb.dump_cells, reverb.dump_cell_x, reverb.dump_cell_z ###
 * The number of cells to dump the path distributions for, and the
 * x/z coordinate pairs for each cell.
 *
 * ### reverb.filter_lines, reverb.filter_line ###
 * Whether or not to filter the output list of photons. If -1,
 * excludes all continuum photons. If >0, provide that many
 * line numbers.
 *
 * ### Notes ###
 * 6/5/18 - Documented by SWM
 * 17/2/19 - Updated to use rdchoice
 **********************************************************/
int
get_meta_params (void)
{

  int i, j, k, z, istate, levl, levu;
  char trackline[LINELENGTH];
  char values[LINELENGTH], answer[LINELENGTH];

  rdpar_comment ("Parameters for Reverberation Modeling (if needed)");

  // ========== DEAL WITH BASIC REVERB TYPE ==========
  strcpy (answer, "none");
  sprintf (values, "%d,%d,%d,%d", REV_NONE, REV_PHOTON, REV_WIND, REV_MATOM);
  geo.reverb = rdchoice ("Reverb.type(none,photon,wind,matom)", values, answer);

  // ========== DEAL WITH DISK SETTINGS ==========
  if (geo.disk_type != DISK_NONE && geo.reverb != REV_NONE)
  {
    strcpy (answer, "correlated");
    sprintf (values, "%d,%d,%d", REV_DISK_CORRELATED, REV_DISK_UNCORRELATED, REV_DISK_IGNORE);
    geo.reverb_disk = rdchoice ("Reverb.disk_type(correlated,uncorrelated,ignore)", values, answer);
  }

  // ========== DEAL WITH VISUALISATION SETTINGS ==========
  if (geo.reverb == REV_WIND || geo.reverb == REV_MATOM)
  {
    //If this requires further parameters, set defaults
    geo.reverb_lines = 0;
    geo.reverb_path_bins = 1000;
    geo.reverb_angle_bins = 100;
    geo.reverb_dump_cells = 0;
    geo.reverb_vis = REV_VIS_NONE;

    rdint ("Reverb.path_bins", &geo.reverb_path_bins);

    strcpy (answer, "none");
    sprintf (values, "%d,%d,%d,%d", REV_VIS_NONE, REV_VIS_VTK, REV_VIS_DUMP, REV_VIS_BOTH);
    geo.reverb_vis = rdchoice ("Reverb.visualisation(none,vtk,dump,both)", values, answer);

    if (geo.reverb_vis == REV_VIS_VTK || geo.reverb_vis == REV_VIS_BOTH)
      //If we're producing a 3d visualisation, select bins. This is just for aesthetics
      rdint ("Reverb.angle_bins(for_vtk)", &geo.reverb_angle_bins);
    if (geo.reverb_vis == REV_VIS_DUMP || geo.reverb_vis == REV_VIS_BOTH)
    {                           //If we;re dumping path arrays, read in the number of cells to dump them for
      rdint ("Reverb.dump_cells(number)", &geo.reverb_dump_cells);
      geo.reverb_dump_cell_x = (double *) calloc (geo.reverb_dump_cells, sizeof (double));
      geo.reverb_dump_cell_z = (double *) calloc (geo.reverb_dump_cells, sizeof (double));
      geo.reverb_dump_cell = (int *) calloc (geo.reverb_dump_cells, sizeof (int));
      for (k = 0; k < geo.reverb_dump_cells; k++)
      {                         //For each we expect, read a paired cell coord as "[i]:[j]". May need to use swind to find indexes.
        rdline ("Reverb.dump_cell(x:z_position)", trackline);
        if (sscanf (trackline, "%lf:%lf", &geo.reverb_dump_cell_x[k], &geo.reverb_dump_cell_z[k]) == EOF)
        {                       //If this line is malformed, warn the user
          Error ("Reverb.dump_cell: Invalid position line '%s'\n \
            Expected format '[x]:[z]'\n", trackline);
          Exit (0);
        }
      }
    }
  }

  // ========== DEAL WITH MATOM LINES ==========
  if (geo.reverb == REV_MATOM)
  {                             //If this is macro-atom mode
    if (geo.rt_mode != RT_MODE_MACRO)
    {                           //But we're not actually working in matom mode...
      Error ("reverb.type: Invalid reverb mode.\n \
      Macro-atom mode selected but macro-atom scattering not on.\n");
      Exit (0);
    }

    //Read in the number of lines to be tracked and allocate space for them
    rdint ("Reverb.matom_lines(number)", &geo.reverb_lines);
    geo.reverb_line = (int *) calloc (geo.reverb_lines, sizeof (int));
    if (geo.reverb_lines < 1)
    {                           //If this is <1, then warn the user and quit
      Error ("Reverb.matom_lines: \
      Must specify 1 or more lines to watch in macro-atom mode.\n");
      Exit (0);
    }

    for (i = 0; i < geo.reverb_lines; i++)
    {                           //Finally, for each line we expect, read it in
      rdline ("Reverb.matom_line(line_index)", trackline);
      if (sscanf (trackline, "%d:%d:%d:%d", &z, &istate, &levu, &levl) == EOF)
      {                         //If this line is malformed, warn the user
        Error ("Reverb.matom_line: Malformed line '%s'\n \
          Expected format '[z]:[istate]:[upper level]:[lower level]'\n", trackline);
        Exit (0);
      }
      else
      {                         //Otherwise, sift through the line list to find what this transition corresponds to
        for (j = 0; j < nlines_macro; j++)
        {                       //And record the line position in geo for comparison purposes
          if (line[j].z == z && line[j].istate == istate && line[j].levu == levu && line[j].levl == levl)
          {                     //We're matching z, ionisation state, and upper and lower level transitions
            geo.reverb_line[i] = line[j].where_in_list;
          }
        }
      }
    }
  }
  else if (geo.reverb == REV_WIND)
  {                             //For wind mode...
    if (geo.wind_radiation == 0)
    {                           //Warn if this data is being gathered but not used (can be useful for debug)
      Error ("reverb.type: Wind radiation is off but wind-based path tracking is enabled!\n");
    }
  }

  // ========== DEAL WITH LINE CULLING ==========
  if (geo.reverb != REV_NONE)
  {
    //Should we filter any lines out?
    //If -1, blacklist continuum, if >0 specify lines as above and whitelist
    //-2 dumps everything to disk, can be big!
    //Automatically include matom_lines
    rdint ("Reverb.filter_lines(0=off,-1=continuum,-2=dumpall,>0=count)", &geo.reverb_filter_lines);
    if (geo.reverb_filter_lines > 0)
    {                           //If we're given a whitelist, allocate temp storage (up to 256 lines!)
      int temp[256], bFound;
      for (i = 0; i < geo.reverb_filter_lines; i++)
      {                         //For each provided line, read in
        rdint ("Reverb.filter_line(line_index)", &temp[i]);
      }
      if (geo.reverb == REV_MATOM)
      {                         //If we're in matom mode, check if those lines have already been included
        for (i = 0; i < geo.reverb_lines; i++)
        {                       //For each matom line
          bFound = 0;
          for (j = 0; j < geo.reverb_filter_lines; j++)
          {                     //Check if it's in the filter list
            if (geo.reverb_line[i] == temp[j])
              bFound = 1;
          }
          if (!bFound)
          {                     //If it's not, add it to the filter list and increment the total lines
            temp[geo.reverb_filter_lines++] = geo.reverb_line[i];
          }
        }
      }
      //Allocate enough space for the filter list
      geo.reverb_filter_line = calloc (geo.reverb_filter_lines, sizeof (int));
      for (i = 0; i < geo.reverb_filter_lines; i++)
      {                         //Populate the filter list from the temp list
        geo.reverb_filter_line[i] = temp[i];
      }
    }
  }
  return (0);
}
