
/***********************************************************/
/** @file  atommicdata_init.c
 * @author ksl
 * @date   July, 2020  
 *
 * @brief  Initialize the various structures that contain atomic
 * data in Python 
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


  if (config != NULL)
  {
    free (config);
  }
  config = (ConfigPtr) calloc (sizeof (config_dummy), NLEVELS);

  if (config == NULL)
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


  /*This initializes the top_phot array - it is used for all ionization processes so some elements
     are only used in some circumstances
   */

  for (n = 0; n < NLEVELS; n++)
  {
    phot_top[n].nlev = (-1);
    phot_top[n].uplev = (-1);
    phot_top[n].nion = (-1);    //the ion to which this cross section belongs
    phot_top[n].n_elec_yield = -1;      //pointer to the electron yield array (for inner shell)
    phot_top[n].n = -1;         //pointer to shell (inner shell)
    phot_top[n].l = -1;         //pointer to l subshell (inner shell only)
    phot_top[n].z = (-1);       //atomic number
    phot_top[n].np = (-1);      //number of points in the fit
    phot_top[n].macro_info = (-1);      //Initialise - don't know if using Macro Atoms or not: set to -1 (SS)
    for (j = 0; j < NCROSS; j++)        //initialise the crooss sectiond
    {
      phot_top[n].freq[j] = (-1);
      phot_top[n].log_freq[j] = (-1);
      phot_top[n].x[j] = (-1);
      phot_top[n].log_x[j] = (-1);

    }
    phot_top[n].f = (-1);       //last frequency
    phot_top[n].log_f = (-1);   //log of last frequency    
    phot_top[n].sigma = 0.0;    //last cross section
    phot_top[n].log_sigma = -1.0;       //log of last cross section

  }


  for (n = 0; n < NIONS * N_INNER; n++) //Initialise atomic arrasy with dimension NIONS*NINNER
  {
    inner_cross[n].nlev = (-1);
    inner_cross[n].uplev = (-1);
    inner_cross[n].nion = inner_elec_yield[n].nion = inner_fluor_yield[n].nion = (-1);
    inner_cross[n].n_elec_yield = -1;
    inner_cross[n].n = inner_elec_yield[n].n = inner_fluor_yield[n].n = (-1);
    inner_cross[n].l = inner_elec_yield[n].l = inner_fluor_yield[n].l = (-1);
    inner_cross[n].z = inner_elec_yield[n].z = inner_fluor_yield[n].z = (-1);
    inner_elec_yield[n].I = inner_elec_yield[n].Ea = 0.0;
    inner_fluor_yield[n].freq = inner_fluor_yield[n].yield = 0.0;
    for (j = 0; j < 10; j++)
      inner_elec_yield[n].prob[j] = 0.0;
    inner_cross[n].np = (-1);
    inner_cross[n].macro_info = (-1);   //Initialise - don't know if using Macro Atoms or not: set to -1 (SS)
    for (j = 0; j < NCROSS; j++)
    {
      inner_cross[n].freq[j] = (-1);
      inner_cross[n].x[j] = (-1);
      inner_cross[n].log_freq[j] = (-1);
      inner_cross[n].log_x[j] = (-1);
    }
    inner_cross[n].f = (-1);
    inner_cross[n].sigma = 0.0;
    inner_cross[n].log_f = (-1);
    inner_cross[n].log_sigma = 0.0;
  }


  for (i = 0; i < NLEVELS; i++)
  {
    config[i].n_bbu_jump = 0;   // initialising the number of jumps from each level to 0. (SS)
    config[i].n_bbd_jump = 0;
    config[i].n_bfu_jump = 0;
    config[i].n_bfd_jump = 0;
  }

  for (n = 0; n < NLINES; n++)
  {
    line[n].freq = -1;
    line[n].f = 0;
    line[n].nion = -1;
    line[n].gl = line[n].gu = 0;
    line[n].el = line[n].eu = 0.0;
    line[n].macro_info = -1;
    line[n].coll_index = -999;
  }

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
