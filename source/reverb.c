/***********************************************************/
/** @file   reverb.c
 * @author SWM
 * @date   July, 2015
 * @brief  Reverberation mapping functions.
 *
 * File containing reverberation mapping functions.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "atomic.h"
#include "sirocco.h"

char delay_dump_file[LINELENGTH];
int delay_dump_bank_size = 65535, delay_dump_bank_curr = 0;
int *delay_dump_spec;
PhotPtr delay_dump_bank;


/**********************************************************/
/** 
 * @brief	Calculates the delay to the observer plane
 *
 * @param [in] pp			Pointer to test photon
 * @return 					Distance from photon to plane
 *
 * Sets up the observer plane, then calculates the distance
 * from the current photon to it. Uses ds_to_plane() for
 * reliability.
 *
 * ###Notes###
 * 9/14	-	Written by SWM
***********************************************************/
double
delay_to_observer (PhotPtr pp)
{
  plane_dummy observer;
  observer.x[0] = pp->lmn[0] * geo.rmax;
  observer.x[1] = pp->lmn[1] * geo.rmax;
  observer.x[2] = pp->lmn[2] * geo.rmax;
  observer.lmn[0] = pp->lmn[0];
  observer.lmn[1] = pp->lmn[1];
  observer.lmn[2] = pp->lmn[2];
  return (pp->path + ds_to_plane (&observer, pp, FALSE));
}

/**********************************************************/
/** 
 * @brief	Prepares delay dump output file
 *
 * @param [in] filename		File root for run
 * @param [in] restart_stat If this is a restart run
 * @param [in] i_rank		Parallel rank
 * @return 					0
 *
 * Sets up filenames, allocates bank for temporary storage
 * of photons that are to be dumped, and if it's not a resume
 * run prints out a header for the file. Filename is set per
 * thread. The file is then built up in batches using
 * delay_dump() in increments of #delay_dump_bank_size.
 *
 * ###Notes###
 * 9/14	-	Written by SWM
***********************************************************/
int
delay_dump_prep (int restart_stat)
{
  FILE *fptr;
  char s_time[LINELENGTH];
  int i;

  //Get output filename
  if (rank_global > 0)
  {
    sprintf (delay_dump_file, "%.100s.delay_dump%d", files.root, rank_global);
  }
  else
  {
    sprintf (delay_dump_file, "%.100s.delay_dump", files.root);
  }

  //Allocate and zero dump files and set extract status
  delay_dump_bank = (PhotPtr) calloc (sizeof (p_dummy), delay_dump_bank_size);
  delay_dump_spec = (int *) calloc (sizeof (int), delay_dump_bank_size);
  for (i = 0; i < delay_dump_bank_size; i++)
    delay_dump_spec[i] = 0;

  if (restart_stat == TRUE)
  {                             //Check whether the output file already has a header
    Log ("delay_dump_prep: Resume run, skipping writeout\n");
    return (0);
  }


  if ((fptr = fopen (delay_dump_file, "w")) != NULL)
  {                             //If this isn't a continue run, prep the output file
    if (rank_global > 0)
    {
      fprintf (fptr, "# Delay dump file for slave process %d\n", rank_global);
    }
    else
    {                           // Construct and write a header string for the output file
      fprintf (fptr, "# Python Version %s\n", VERSION);
      get_time (s_time);
      fprintf (fptr, "# Date	%s\n#  \n", s_time);
      fprintf (fptr, "#\n# %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n", "Np", "Freq.", "Lambda",
               "Weight", "LastX", "LastY", "LastZ", "Scat.", "RScat.", "Delay", "Spec.", "Orig.", "Res.", "LineRes.");
    }
    fclose (fptr);
    Log ("delay_dump_prep: Thread %d successfully prepared file '%s' for writing\n", rank_global, delay_dump_file);
  }
  else
  {
    Error ("delay_dump_prep: Thread %d failed to open file '%s' due to error %d: %s\n", rank_global, delay_dump_file, errno,
           strerror (errno));
  }
  return (0);
}

/**********************************************************/
/** 
 * @brief	Finishes dumping tracked photons to file
 *
 * @return 					0
 *
 * Dumps the remaining tracked photons to file, frees memory.
 *
 * ###Notes###
 * 6/15	-	Written by SWM
***********************************************************/
int
delay_dump_finish (void)
{
  Log ("delay_dump_finish: Dumping %d photons to file\n", delay_dump_bank_curr - 1);
  if (delay_dump_bank_curr > 0)
  {
    delay_dump (delay_dump_bank, delay_dump_bank_curr - 1);
  }
  free (delay_dump_bank);
  free (delay_dump_spec);
  return (0);
}

/**********************************************************/
/** 
 * @brief	Prepares delay dump output file
 *
 * @param [in] i_ranks		Number of parallel processes
 * @return 					0
 *
 * Collects all the delay dump files together at the end. 
 * Called by the master thread. Uses 'cat' for simplicity.
 *
 * ###Notes###
 * 6/15	-	Written by SWM
***********************************************************/
int
delay_dump_combine (int i_ranks)
{
  FILE *fopen ();               //, *f_base, *f_cat;
  char c_call[LINELENGTH];      //, c_cat[LINELENGTH], c_char;
  //int i;
/*
	f_base = fopen(delay_dump_file, 'a');
	for(i=1;i<iRanks;i++)
	{
		sprintf(c_cat, "%s%d", c_cat, i);
		if((f_cat = fopen(c_cat, 'r')) == NULL)
			Error("delay_dump_combine: Missing file %s%d", c_cat,i);
		else
			while((c_char = fgetc(f_cat)) != EOF) fputc(c_char, f_base);
		fclose(f_cat);
	}
	fclose(f_base);
*/
  //Yes this is done as a system call and won 't work on Windows machines. Lazy solution!
  sprintf (c_call, "cat %.50s[0-9]* >> %.100s", delay_dump_file, delay_dump_file);
  if (system (c_call) < 0)
  {
    Error ("delay_dump_combine: Error calling system command '%s'", c_call);
  }
  else
  {
    sprintf (c_call, "rm %.50s[0-9]*", delay_dump_file);
    if (system (c_call) < 0)
    {
      Error ("delay_dump_combine: Error calling system command '%s'", c_call);
    }
  }
  return (0);
}

/**********************************************************/
/** 
 * @brief	Dumps tracked photons to file
 *
 * @param [in] np			Pointer to photon array tp dump
 * @param [in] np			Number photons in p
 * @param [in] iExtracted	Array of extract flags for p
 * @return 					0
 *
 * Sifts through the photons in the delay_dump file, checking
 * if they've scattered or were generated in the wind and so
 * contribute to the delay map. Uses the same filters as 
 * the spectra_create() function for scatters & angles.
 *
 * ###Notes###
 * 6/15	-	Written by SWM
***********************************************************/
int
delay_dump (PhotPtr p, int np)
{
  FILE *fopen (), *fptr;
  int nphot, mscat, mtopbot, i, subzero;
  double delay;
  subzero = 0;

  Log ("delay_dump: Dumping %d photons\n", np);
  /*
   * Open a file for writing the spectrum
   */
  if ((fptr = fopen (delay_dump_file, "a")) == NULL)
  {
    Error ("delay_dump: Unable to reopen %s for writing\n", delay_dump_file);
    Exit (0);
  }
  for (nphot = 0; nphot < np; nphot++)
  {
    /*
     * Complicated if statement to allow one to choose
     * whether to construct the spectrum from all photons
     * or just from photons which have scattered a
     * specific number of times.  01apr13--ksl-Modified
     * if statement to change behavior on negative
     * numbers to say that a negative number for mscat
     * implies that you accept any photon with |mscat| or
     * more scatters
     */
    i = delay_dump_spec[nphot];
    if (((mscat = xxspec[i].nscat) >= MAXSCAT ||
         p[nphot].nscat == mscat ||
         (mscat < 0 && p[nphot].nscat >= (-mscat))) && ((mtopbot = xxspec[i].top_bot) == 0 || (mtopbot * p[nphot].x[2]) > 0))
    {
      delay = (delay_to_observer (&p[nphot]) - geo.rmax) / VLIGHT;
      if (delay < 0)
        subzero++;

      fprintf (fptr, "%-12d %-12.5g %-12.7g %-12.5g %-12.5g %-12.5g %-12.5g %-12d %-12d %-12.5g %-12d %-12d %-12d %-12d\n",
               p[nphot].np, p[nphot].freq, VLIGHT * 1e8 / p[nphot].freq, p[nphot].w, p[nphot].x[0], p[nphot].x[1], p[nphot].x[2],
               p[nphot].nscat, p[nphot].nrscat, delay, i - MSPEC, p[nphot].origin, p[nphot].nres, p[nphot].line_res);
    }
  }

  if (subzero > 0)
  {
    Error ("delay_dump: %d photons with <0 delay found! Increase path bin resolution to minimise this error.", subzero);
  }
  fclose (fptr);
  return (0);
}

/**********************************************************/
/** 
 * @brief	Preps a single photon to be dumped
 *
 * @param [in] pp			Pointer to extracted photon
 * @param [in] i_spec		Spectrum p extracted to
 * @return 					0
 *
 * Takes a photon and copies it to the staging arrays for 
 * delay dumping, to be output to file later.
 *
 * ###Notes###
 * 6/15	-	Written by SWM
***********************************************************/
int
delay_dump_single (PhotPtr pp, int i_spec)
{
  if (geo.reverb_filter_lines == -1 && pp->nres == -1)
  {
    //If we're filtering out continuum photons and this is a continuum photon, throw it away.
    return (1);
  }
  else if (geo.reverb_filter_lines > 0)
  {
    //If we're filtering to *only* photons of given lines, is this one of them? If not, throw away
    int i, bFound = 0;
    for (i = 0; i < geo.reverb_filter_lines; i++)
      if (pp->nres == geo.reverb_filter_line[i])
        bFound = 1;
    if (!bFound)
      return (1);
  }

  stuff_phot (pp, &delay_dump_bank[delay_dump_bank_curr]);      //Bank single photon in temp array
  delay_dump_spec[delay_dump_bank_curr] = i_spec;       //Record photon spectrum too
  if (delay_dump_bank_curr == delay_dump_bank_size - 1) //If temp array is full
  {
    delay_dump (delay_dump_bank, delay_dump_bank_size);
    delay_dump_bank_curr = 0;   //Dump to file, zero array position
  }
  else
  {
    delay_dump_bank_curr++;
  }
  return (0);
}
