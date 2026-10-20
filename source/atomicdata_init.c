
/***********************************************************/
/** @file  atommicdata_init.c
 * @author ksl
 * @date   July, 2020  
 *
 * @brief  Initialize the various structures that contain atomic
 * data in Sirocco 
 * and other similar programs
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "log.h"
// If routines are added cproto > atomic_proto.h should be run
#include "atomic_proto.h"

#define LINELENGTH 400
#define MAXWORDS    20

#define NXSECTION_INIT 100      /* Initial number of records in phot_top, inner_cross and auger_macro; they double as needed */

/* Initial values of log_sigma, the log of the last cross section, for the two tables */
#define PHOT_TOP_LOG_SIGMA_INIT    (-1.0)
#define INNER_CROSS_LOG_SIGMA_INIT 0.0

static int phot_top_max = 0;    /* The number of records currently allocated in phot_top */
static int inner_cross_max = 0; /* The number of records currently allocated in inner_cross */
static int auger_macro_max = 0; /* The number of records currently allocated in auger_macro */


/**********************************************************/
/**
 * @brief      Grow a table of cross-section records so that it
 * holds record n
 *
 * @param [in] Topbase_phot *  table   The table, which may be NULL
 * @param [in, out] int *  nmax   The number of records allocated in the table
 * @param [in] int  n   The index of the record that must exist
 * @param [in] double  log_sigma   The initial value of log_sigma for new records
 * @param [in] char *  name   The name of the table, for error messages
 * @return     The (possibly moved) table
 *
 * @details
 * The table is doubled in size until it holds record n.  New records
 * are initialized to the values used before the tables were grown
 * dynamically, with no cross-section points allocated.
 *
 * ### Notes ###
 * Because realloc may move the table, nothing may hold a pointer into
 * it while the atomic data are being read.  The frequency-ordered
 * pointers phot_top_ptr and inner_cross_ptr are made after reading.
 *
 **********************************************************/

static Topbase_phot *
grow_xsections (Topbase_phot *table, int *nmax, int n, double log_sigma, char *name)
{
  int nnew, i;

  if (n < *nmax)
    return (table);

  nnew = (*nmax > 0) ? *nmax : NXSECTION_INIT;
  while (nnew <= n)
    nnew *= 2;

  table = (Topbase_phot *) realloc (table, nnew * sizeof (Topbase_phot));
  if (table == NULL)
  {
    Error ("grow_xsections: There is a problem in allocating memory for %d records of %s\n", nnew, name);
    exit (0);
  }

  for (i = *nmax; i < nnew; i++)
  {
    memset (&table[i], 0, sizeof (Topbase_phot));
    table[i].nlev = (-1);
    table[i].uplev = (-1);
    table[i].nion = (-1);       //the ion to which this cross section belongs
    table[i].n_elec_yield = -1; //index into the electron yield array (inner shell only)
    table[i].n = -1;            //n of the shell (inner shell only)
    table[i].l = -1;            //l of the subshell (inner shell only)
    table[i].z = (-1);          //atomic number
    table[i].np = (-1);         //number of points in the fit
    table[i].macro_info = (-1); //Initialise - don't know if using Macro Atoms or not: set to -1 (SS)
    table[i].f = (-1);          //last frequency
    table[i].log_f = (-1);      //log of last frequency
    table[i].sigma = 0.0;       //last cross section
    table[i].log_sigma = log_sigma;     //log of last cross section
  }

  *nmax = nnew;
  return (table);
}


/**********************************************************/
/**
 * @brief      Free a table of cross-section records and the
 * points of each record
 *
 * @param [in, out] Topbase_phot **  table   The table, set to NULL on return
 * @param [in, out] int *  nmax   The number of records allocated, set to 0 on return
 * @return     Always returns 0
 *
 **********************************************************/

static int
free_xsections (Topbase_phot **table, int *nmax)
{
  int n;

  for (n = 0; n < *nmax; n++)
    free ((*table)[n].freq);

  free (*table);
  *table = NULL;
  *nmax = 0;
  return (0);
}


/**********************************************************/
/**
 * @brief      Allocate the cross-section points of one record
 *
 * @param [in, out] TopPhotPtr  xptr   The record
 * @param [in] int  np   The number of points
 * @return     Always returns 0
 *
 * @details
 * freq, log_freq, x and log_x are allocated as one block of 4*np
 * doubles, so freeing freq frees all four.  Any points the record
 * already had are freed first, since a cross section may be replaced
 * by one with a different number of points.
 *
 **********************************************************/

static int
alloc_xsection_points (TopPhotPtr xptr, int np)
{
  double *block;

  if (np < 1)
  {
    Error ("alloc_xsection_points: A cross section cannot have %d points\n", np);
    exit (0);
  }

  free (xptr->freq);
  block = (double *) calloc (4 * np, sizeof (double));
  if (block == NULL)
  {
    Error ("alloc_xsection_points: There is a problem in allocating memory for %d points\n", np);
    exit (0);
  }

  xptr->freq = block;
  xptr->log_freq = block + np;
  xptr->x = block + 2 * np;
  xptr->log_x = block + 3 * np;
  return (0);
}


/**********************************************************/
/**
 * @brief      Make sure phot_top record n exists, with room for np
 * cross-section points
 *
 * @param [in] int  n   The index of the record in phot_top
 * @param [in] int  np   The number of points in the cross section
 * @return     A pointer to the record
 *
 * @details
 * phot_top may move, so the returned pointer is only valid until the
 * next call.  Fields other than the points of an existing record are
 * left as they are.
 *
 **********************************************************/

TopPhotPtr
alloc_phot_top_record (int n, int np)
{
  phot_top = grow_xsections (phot_top, &phot_top_max, n, PHOT_TOP_LOG_SIGMA_INIT, "phot_top");
  alloc_xsection_points (&phot_top[n], np);
  return (&phot_top[n]);
}


/**********************************************************/
/**
 * @brief      Make sure inner_cross record n exists, with room for np
 * cross-section points
 *
 * @param [in] int  n   The index of the record in inner_cross
 * @param [in] int  np   The number of points in the cross section
 * @return     A pointer to the record
 *
 * @details
 * As alloc_phot_top_record, but for the inner shell cross sections.
 *
 **********************************************************/

TopPhotPtr
alloc_inner_cross_record (int n, int np)
{
  inner_cross = grow_xsections (inner_cross, &inner_cross_max, n, INNER_CROSS_LOG_SIGMA_INIT, "inner_cross");
  alloc_xsection_points (&inner_cross[n], np);
  return (&inner_cross[n]);
}


/**********************************************************/
/**
 * @brief      Grow the table of macro-atom Auger records so that it
 * holds record n
 *
 * @param [in] AugerPtr  table   The table, which may be NULL
 * @param [in, out] int *  nmax   The number of records allocated in the table
 * @param [in] int  n   The index of the record that must exist
 * @return     The (possibly moved) table
 *
 * @details
 * As grow_xsections, but for auger_macro.  Records are referred to
 * by index (xconfig[].iauger) while the data are read, so the table
 * may move.
 *
 **********************************************************/

static AugerPtr
grow_auger_macro (AugerPtr table, int *nmax, int n)
{
  int nnew, i, j;

  if (n < *nmax)
    return (table);

  nnew = (*nmax > 0) ? *nmax : NXSECTION_INIT;
  while (nnew <= n)
    nnew *= 2;

  table = (AugerPtr) realloc (table, nnew * sizeof (auger_dummy));
  if (table == NULL)
  {
    Error ("grow_auger_macro: There is a problem in allocating memory for %d Auger records\n", nnew);
    exit (0);
  }

  for (i = *nmax; i < nnew; i++)
  {
    memset (&table[i], 0, sizeof (auger_dummy));
    table[i].z = -1;
    table[i].nion = -1;
    table[i].istate = -1;
    table[i].nconfig = -1;
    table[i].iauger = -1;
    table[i].nauger = -1;
    table[i].Avalue_auger = 0.0;
    for (j = 0; j < NAUGER_ELECTRONS; j++)
    {
      table[i].nconfig_target[j] = -1;
      table[i].branching_ratio[j] = 0.0;
    }
  }

  *nmax = nnew;
  return (table);
}


/**********************************************************/
/**
 * @brief      Make sure auger_macro record n exists
 *
 * @param [in] int  n   The index of the record in auger_macro
 * @return     A pointer to the record
 *
 * @details
 * auger_macro may move, so the returned pointer is only valid until
 * the next call.
 *
 **********************************************************/

AugerPtr
alloc_auger_macro_record (int n)
{
  auger_macro = grow_auger_macro (auger_macro, &auger_macro_max, n);
  return (&auger_macro[n]);
}


/**********************************************************/
/**
 * @brief      Free the photoionization and inner shell cross sections
 *
 * @return     Always returns 0
 *
 * @details
 * Frees phot_top and inner_cross, the points of each record, and the
 * frequency-ordered pointers phot_top_ptr and inner_cross_ptr, and
 * sets them all to NULL.
 *
 **********************************************************/

int
free_xsection_data ()
{
  free_xsections (&phot_top, &phot_top_max);
  free_xsections (&inner_cross, &inner_cross_max);
  free (phot_top_ptr);
  phot_top_ptr = NULL;
  free (inner_cross_ptr);
  inner_cross_ptr = NULL;
  return (0);
}


/**********************************************************/
/**
 * @brief      routine to initialze the data structures 
 *  	defined in "atomic.h"
 *
 *
 * @return     Always returns 0
 *
 * @details
 *
 *
 * ### Notes ###
 *
 * This functionality was moved from get_atomicdata basically
 * to shorten that very long routine in July 2020.  Some
 * counters are still set in get_atomicdata
 *
 *
 **********************************************************/

int
init_atomic_data ()
{

  int n, i, j;
  int n1;


/* Allocate structures for storage of data */

  if (ele != NULL)
  {
    free (ele);
  }
  ele = (ElemPtr) calloc (sizeof (ele_dummy), NELEMENTS);

  if (ele == NULL)
  {
    Error ("There is a problem in allocating memory for the element structure\n");
    exit (0);
  }
  else
  {
    Log_silent
      ("Allocated %10d bytes for each of %6d elements of   elements totaling %10.0f Mb \n",
       sizeof (ele_dummy), NELEMENTS, 1.e-6 * NELEMENTS * sizeof (ele_dummy));
  }


  if (ion != NULL)
  {
    free (ion);
  }
  ion = (IonPtr) calloc (sizeof (ion_dummy), NIONS);

  if (ion == NULL)
  {
    Error ("There is a problem in allocating memory for the ion structure\n");
    exit (0);
  }
  else
  {
    Log_silent
      ("Allocated %10d bytes for each of %6d elements of       ions totaling %10.1f Mb \n",
       sizeof (ion_dummy), NIONS, 1.e-6 * NIONS * sizeof (ion_dummy));
  }


  if (xconfig != NULL)
  {
    free (xconfig);
  }
  xconfig = (ConfigPtr) calloc (sizeof (config_dummy), NLEVELS);

  if (xconfig == NULL)
  {
    Error ("There is a problem in allocating memory for the config structure\n");
    exit (0);
  }
  else
  {
    Log_silent
      ("Allocated %10d bytes for each of %6d elements of     config totaling %10.1f Mb \n",
       sizeof (config_dummy), NLEVELS, 1.e-6 * NLEVELS * sizeof (config_dummy));
  }

  if (line != NULL)
  {
    free (line);
  }
  line = (LinePtr) calloc (sizeof (line_dummy), NLINES);

  if (line == NULL)
  {
    Error ("There is a problem in allocating memory for the line structure\n");
    exit (0);
  }
  else
  {
    Log_silent
      ("Allocated %10d bytes for each of %6d elements of       line totaling %10.1f Mb \n",
       sizeof (line_dummy), NLINES, 1.e-6 * NLINES * sizeof (line_dummy));
  }

  /* auger_macro is grown by alloc_auger_macro_record as the data are read */

  free (auger_macro);
  auger_macro = NULL;
  auger_macro_max = 0;
  auger_macro = grow_auger_macro (auger_macro, &auger_macro_max, 0);


  /* The photoionization and inner shell cross sections start small and are
     grown by alloc_phot_top_record and alloc_inner_cross_record as the data
     are read.  Free anything left from a previous call first. */

  free_xsection_data ();

  phot_top = grow_xsections (phot_top, &phot_top_max, 0, PHOT_TOP_LOG_SIGMA_INIT, "phot_top");
  inner_cross = grow_xsections (inner_cross, &inner_cross_max, 0, INNER_CROSS_LOG_SIGMA_INIT, "inner_cross");

  /* Initialize variables */



  phot_freq_min = VERY_BIG;
  inner_freq_min = VERY_BIG;

  for (n = 0; n < NELEMENTS; n++)
  {
    strcpy (ele[n].name, "none");
    ele[n].z = (-1);
    ele[n].abun = (-1);
    ele[n].firstion = (-1);
    ele[n].nions = (-1);
    ele[n].istate_max = (-1);
  }

  for (n = 0; n < NIONS; n++)
  {
    ion[n].z = (-1);
    ion[n].istate = (-1);
    ion[n].nelem = (-1);
    ion[n].ip = (-1);
    ion[n].g = (-1);
    ion[n].log_g = (-1);        //Log version of multiplicity
    ion[n].nmax = (-1);
    ion[n].firstlevel = (-1);
    ion[n].nlevels = (-1);
    ion[n].first_nlte_level = (-1);
    ion[n].first_levden = (-1);
    ion[n].nlte = (-1);
    ion[n].phot_info = (-1);
    ion[n].macro_info = (-1);   //Initialise - don't know if using Macro Atoms or not: set to -1 (SS)
    ion[n].ntop_first = 0;      // The fact that ntop_first and ntop  are initialized to 0 and not -1 is important
    ion[n].ntop_ground = 0;     //NSH 0312 initialize the new pointer for GS cross sections
    ion[n].ntop = 0;
    ion[n].nxphot = (-1);
    ion[n].lev_type = (-1);     // Initialise to indicate we don't know what types of configurations will be read
    ion[n].drflag = 0;          //Initialise to indicate as far as we know, there are no dielectronic recombination parameters associated with this ion.
    ion[n].total_rrflag = 0;    //Initialise to say this ion has no Badnell total recombination data
    ion[n].nxtotalrr = -1;      //Initialise the pointer into the bad_t_rr structure.
    ion[n].bad_gs_rr_t_flag = 0;        //Initialise to say this ion has no Badnell ground state recombination data
    ion[n].bad_gs_rr_r_flag = 0;        //Initialise to say this ion has no Badnell ground state recombination data
    ion[n].nxbadgsrr = -1;      //Initialise the pointer into the bad_gs_rr structure.
    ion[n].dere_di_flag = 0;    //Initialise to say this ion has no Dere DI rate data
    ion[n].nxderedi = -1;       //Initialise the pointer into the Dere DI rate structure
    ion[n].n_inner = 0;         //Initialise the pointer to say we have no inner shell ionization cross sections
    ion[n].n_ch_ex = -1;        //Initialise the pointer to the charge exchange 

    for (i = 0; i < N_INNER; i++)
      ion[n].nxinner[i] = -1;   //Inintialise the inner shell pointer array
  }


  /* The phot_top and inner_cross records are initialized by grow_xsections
     as they are allocated */

  for (n = 0; n < NIONS * N_INNER; n++) //Initialise atomic arrasy with dimension NIONS*NINNER
  {
    inner_elec_yield[n].nion = inner_fluor_yield[n].nion = (-1);
    inner_elec_yield[n].n = inner_fluor_yield[n].n = (-1);
    inner_elec_yield[n].l = inner_fluor_yield[n].l = (-1);
    inner_elec_yield[n].z = inner_fluor_yield[n].z = (-1);
    inner_elec_yield[n].I = inner_elec_yield[n].Ea = 0.0;
    inner_fluor_yield[n].freq = inner_fluor_yield[n].yield = 0.0;
    for (j = 0; j < 10; j++)
      inner_elec_yield[n].prob[j] = 0.0;
  }


  for (i = 0; i < NLEVELS; i++)
  {
    xconfig[i].n_bbu_jump = 0;  // initialising the number of jumps from each level to 0. (SS)
    xconfig[i].n_bbd_jump = 0;
    xconfig[i].n_bfu_jump = 0;
    xconfig[i].n_bfd_jump = 0;
    xconfig[i].iauger = -1;
    xconfig[i].nauger = 0;
  }

  for (n = 0; n < NLINES; n++)
  {
    line[n].freq = -1;
    line[n].f = -1;
    line[n].nion = -1;
    line[n].gl = line[n].gu = 0;
    line[n].el = line[n].eu = 0.0;
    line[n].macro_info = -1;
    line[n].coll_index = -999;
  }

  /* The auger_macro records are initialized by grow_auger_macro as they are allocated */

/* The following lines initialise the dielectronic recombination structure */
  for (n = 0; n < NIONS; n++)
  {
    drecomb[n].nion = -1;
    drecomb[n].nparam = -1;     //the number of parameters - it varies from ion to ion
    for (n1 = 0; n1 < MAX_DR_PARAMS; n1++)
    {
      drecomb[n].c[n1] = 0.0;
      drecomb[n].e[n1] = 0.0;
    }
  }


/* The following lines initialise the badnell total radiative recombination rate structure */
  for (n = 0; n < NIONS; n++)
  {
    total_rr[n].nion = -1;
    for (n1 = 0; n1 < T_RR_PARAMS; n1++)
    {
      total_rr[n].params[n1] = 0.0;
    }

  }


/* The following lines initialise the badnell total radiative recombination rate structure */
  for (n = 0; n < NIONS; n++)
  {
    bad_gs_rr[n].nion = -1;
    for (n1 = 0; n1 < BAD_GS_RR_PARAMS; n1++)
    {
      bad_gs_rr[n].temps[n1] = 0.0;
      bad_gs_rr[n].rates[n1] = 0.0;
    }

  }


/* The following lines initialise the dere direct ionization rate struture */
  for (n = 0; n < NIONS; n++)
  {
    dere_di_rate[n].nion = -1;
    dere_di_rate[n].xi = 0.0;
    dere_di_rate[n].min_temp = 1e99;
    for (n1 = 0; n1 < DERE_DI_PARAMS; n1++)
    {
      dere_di_rate[n].temps[n1] = 0.0;
      dere_di_rate[n].rates[n1] = 0.0;
    }

  }


/* The following lines initialise the Sutherland gaunt factors */
  gaunt_n_gsqrd = 0;            //The number of sets of scaled temperatures we have data for
  for (n = 0; n < MAX_GAUNT_N_GSQRD; n++)
  {
    gaunt_total[n].log_gsqrd = 0.0;
    gaunt_total[n].gff = 0.0;
    gaunt_total[n].s1 = 0.0;
    gaunt_total[n].s2 = 0.0;
    gaunt_total[n].s3 = 0.0;
  }

/* The following lines initialise the Sutherland gaunt factors */
  n_charge_exchange = 0;        //The number of sets of scaled temperatures we have data for
  for (n = 0; n < MAX_CHARGE_EXCHANGE; n++)
  {
    charge_exchange[n].nion1 = charge_exchange[n].nion2 = -1;
    charge_exchange[n].a = charge_exchange[n].b = charge_exchange[n].c = charge_exchange[n].d = 0.0;
    charge_exchange[n].tmin = -1e99;
    charge_exchange[n].tmax = 1e99;
    charge_exchange[n].energy_defect = 0.0;
    charge_exchange[n].delta_e_ovr_k = 0.0;
  }


/* The following lines initialise the collision strengths */
  for (n = 0; n < NLINES; n++)
  {
    coll_stren[n].n = -1;       //Internal index
    coll_stren[n].lower = -1;   //The lower energy level - this is in Chianti notation and is currently unused
    coll_stren[n].upper = -1;   //The upper energy level - this is in Chianti notation and is currently unused
    coll_stren[n].energy = 0.0; //The energy of the transition
    coll_stren[n].gf = 0.0;
    coll_stren[n].hi_t_lim = 0.0;       //The high temperature limit
    coll_stren[n].n_points = 0; //The number of points in the spline fit
    coll_stren[n].type = -1;    //The type of fit, this defines how one computes the scaled temperature and scaled coll strength
    coll_stren[n].scaling_param = 0.0;  //The scaling parameter C used in the Burgess and Tully calculations
    for (n1 = 0; n1 < N_COLL_STREN_PTS; n1++)
    {
      coll_stren[n].sct[n1] = 0.0;      //The scaled temperature points in the fit
      coll_stren[n].scups[n1] = 0.0;
    }
  }


  return (0);
}
