
/***********************************************************/
/** @file  gridwind.c
 * @author ksl,sss,jm
 * @date   April, 2018
 *
 * @brief  This file contains routines for allocating space for most
 * of the structures needed by Python.  These include 
 * Wind and Plasma structures, and if using macro atoms, the
 * various structures rwquired for this.  The file also contains
 * a routine that maps grid cells in the wind to those in the 
 * Plasma structure, and vice versa.
 *
 *
 * The hierachy of stuctures that describe a model in Python 
 * are domains (which describe a region of the wind), wmain 
 * which contains basic information about the cells in each domain,
 * and plasmamain which contains more detailed information about
 * cells which are actually in the wind.  
 *
 * The cells in wmain are created on a simple grid, often
 * a cylindrical region, but no all of this region may be
 * in the wind, as the case in a binconical flow.  To save
 * memory, one usually only creates elements of plasmamain
 * that correspond to wind cells that are actually in the 
 * wind.
 *
 * The Plasma structure includes space for storing all of
 * the information needed to calculate ionization and equilibrium 
 * in the two-level approximation.  However, additional structures
 * are required when operating in macro atom mode.  Routines
 * for allocating the structures associated with macroatoms
 * are also included here.
 *
 * Note that a coniderable effort has been made to minimize the
 * total size of the various structures, and so PlasmaMain 
 * contains a number of variables that are actually pointers to
 * other structures that are allocated dynamically depeding on
 * the atomic data.  This is also true in the case of the main
 * macro atom structure MacroMain.  
 *
 * Note that some structures are allocated in other portions of the
 * code, e.g get_atomic_data, but all of the structures that 
 * are needed to set up the wind are here. 
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/** 
 * @brief      Create a map between wind and plasma cells
 *
 * @return     Always returns 0
 *
 * @details
 * The routine fills variables, nplasma in wmain, and nwind in plasmamain
 * that map elements of wmain to those in plasmamain anc vice versus
 *
 *
 * ### Notes ###
 * Normally, there are fewer plasma elements than wind elements, since
 * plasma elements are created only for those wind elements which have 
 * finite volume in the wind.  For each domain, a coordinate grid is
 * created that covers a region of space, e.g a portion of a cylindrical
 * grid.  However some cells in this region are empty or matter, e.g
 * in a bi-conical flow.  
 * 
 * wmain and plasmamain will already have been allocated memory by
 * the time this routine is called.
 *
 *
 *
 **********************************************************/

int
create_wind_and_plasma_cell_maps ()
{
  int i, j;
  j = 0;

  for (i = 0; i < NDIM2; i++)
  {
    wmain[i].nwind = i;
    if (wmain[i].vol > 0)
    {
      wmain[i].nplasma = j;
      plasmamain[j].nplasma = j;
      plasmamain[j].nwind = i;
      j++;
      if (wmain[i].inwind < 0)
      {
        Error
          ("create_wind_and_plasma_cell_maps: wind cell %d (nplasma %d) has volume but not flagged as in wind! Critical error, could cause undefined behaviour. Exiting.\n",
           i, j);
        Exit (0);
      }
    }
    else
    {
      wmain[i].nplasma = NPLASMA;
      if (wmain[i].inwind >= 0)
      {
        Error
          ("create_wind_and_plasma_cell_maps: wind cell %d has zero volume but flagged inwind! Critical error, could cause undefined behaviour. Exiting.\n",
           i);
        Exit (0);
      }
    }
  }
  if (j != NPLASMA)
  {
    Error ("create_wind_and_plasma_cell_maps: Problems with matching cells -- Expected %d Got %d\n", NPLASMA, j);
    Exit (0);
  }

  plasmamain[NPLASMA].nplasma = NPLASMA;
  plasmamain[NPLASMA].nwind = -1;
  return (0);
}





/**********************************************************/
/** 
 * @brief      Allocate space for the wind domain
 *
 * @param [in] nelem   The number of elements in the wind domain to allocate
 * @return     Returns 0 unless there is insufficient space, in which
 * case the routine will call the program to exits
 *
 * @details
 *
 * ### Notes ###
 * The wind is conisists of grids of cells that paper the active region
 * of the calculation.  Each domain has a number of these grid cells.
 * nelem is the total number of such cells.
 *
 **********************************************************/

int
calloc_wind (nelem)
     int nelem;
{

  if (wmain != NULL)
  {
    free (wmain);
  }

  wmain = (WindPtr) calloc (nelem + 1, sizeof (wind_dummy));

  if (wmain == NULL)
  {
    Error ("There is a problem in allocating memory for the wind structure\n");
    Exit (0);
  }
  else
  {
    Log
      ("Allocated %10d bytes for each of %5d elements of             totaling %10.1f Mb\n",
       sizeof (wind_dummy), nelem + 1, 1.e-6 * (nelem + 1) * sizeof (wind_dummy));
  }

  return (0);
}






/**********************************************************/
/** 
 * @brief      Allocate memory for plasmamain, which contains temperature, densities, etc.
 * for cells which are in the wind
 *
 * @param [in, out] int  nelem   The number of elements of plasmamain to allocate
 * @return     Returns 0, unless the program is unable to allocate the requested memory 
 * in which case the program exits
 *
 * @details
 *
 * ### Notes ###
 *
 * This only allocates elements.  It does not populate them
 * with any information.
 *
 * The routine also allocates space for storing photons associated with
 * each plasma cell.  This is used for creating wind photons in a cell
 * in cases where one wants to created multiple photons of a certain 
 * type, and generating the cdf for this is time consuming, but generating
 * photons once you have the cdf is not.  It is currently used for fb processes
 * only.
 *
 **********************************************************/

int
calloc_plasma (nelem)
     int nelem;
{

  if (plasmamain != NULL)
  {
    free (plasmamain);
  }

  /*Allocate one extra element to store data where there is no volume */

  plasmamain = (PlasmaPtr) calloc (sizeof (plasma_dummy), (nelem + 1));
  geo.nplasma = nelem;

  if (plasmamain == NULL)
  {
    Error ("There is a problem in allocating memory for the plasma structure\n");
    Exit (0);
  }
  else
  {
    Log
      ("Allocated %10d bytes for each of %5d elements of      plasma totaling %10.1f Mb \n",
       sizeof (plasma_dummy), (nelem + 1), 1.e-6 * (nelem + 1) * sizeof (plasma_dummy));
  }

  /* Now allocate space for storing photon frequencies -- 57h */
  if (photstoremain != NULL)
  {
    free (photstoremain);
  }
  photstoremain = (PhotStorePtr) calloc (sizeof (photon_store_dummy), (nelem + 1));

  if (photstoremain == NULL)
  {
    Error ("There is a problem in allocating memory for the photonstore structure\n");
    Exit (0);
  }
  else
  {
    Log
      ("Allocated %10d bytes for each of %5d elements of photonstore totaling %10.1f Mb \n",
       sizeof (photon_store_dummy), (nelem + 1), 1.e-6 * (nelem + 1) * sizeof (photon_store_dummy));
  }

  /* Repeat above for matom storage photon frequencies -- 82h */
  if (matomphotstoremain != NULL)
  {
    free (matomphotstoremain);
  }
  matomphotstoremain = (MatomPhotStorePtr) calloc (sizeof (matom_photon_store_dummy), (nelem + 1));

  if (matomphotstoremain == NULL)
  {
    Error ("There is a problem in allocating memory for the matomphotonstore structure\n");
    Exit (0);
  }
  else
  {
    Log
      ("Allocated %10d bytes for each of %5d elements of matomphotonstore totaling %10.1f Mb \n",
       sizeof (matom_photon_store_dummy), (nelem + 1), 1.e-6 * (nelem + 1) * sizeof (matom_photon_store_dummy));
  }

  return (0);
}



/**********************************************************/
/** 
 * @brief      Check that one is not trying to get information from
 * the dummy Plasma element that is associated with wind elements that
 * are not really in the wind
 *
 * @param [in] PlasmaPtr  xplasma   A pointer to the plasma element
 * @param [in] char  message[]   A message which is printed out if 
 * one has tried to access the dummy Plsma element
 * @return     0 (FALSE) if one is not in the empty plasma cell, 1 
 * (TRUE) if you have gotten there.
 *
 * @details
 * Wind cells that have no associated volume are assigned to
 * the same plasma element.  One should never be asked to 
 * generate a photon or indeeed have a scattering event in
 * such wind cells.  This little routine is just a diagnostic
 * routine that throws an error message if one has gotten
 * into such a situation.
 *
 * ### Notes ###
 *
 **********************************************************/

int
check_plasma (xplasma, message)
     PlasmaPtr xplasma;
     char message[];
{
  if (xplasma->nplasma == NPLASMA)
  {
    Error ("check_plasma: In Dummy Plasma Cell when probably don't want to be:  %s \n", message);
    return (TRUE);
  }
  else
    return (FALSE);
}




/**********************************************************/
/** 
 * @brief      Allocate memory to store information for macro atoms for
 * each Plasma element
 *
 * @param [in] int  nelem   The number of cells that are actually included
 * in the wind. 
 * @return     Always returns 0, unless memory cannot allocate the
 * requested memory in which case the program exits
 *
 * @details
 * This routine allocates memory for the structure macromain, in the
 * situation where the program is being run in macro-atom mode.
 *
 * ### Notes ###
 * The size of calloc macro depends on the number cells in the wind,
 * but many of the elements of macro main are simply stucture pointers
 * that are allocated by calloc_estimators.
 * 
 *
 **********************************************************/

int
calloc_macro (nelem)
     int nelem;
{

  /* JM 1502 -- commented out this if loop because we want 
     the macro structure to be allocated regardless in geo.rt_mode = RT_MODE_MACRO. see #138 */

  if (macromain != NULL)
  {
    free (macromain);
  }

  //Allocate one extra element to store data where there is no volume

  macromain = (MacroPtr) calloc (sizeof (macro_dummy), (nelem + 1));
  geo.nmacro = nelem;

  if (macromain == NULL)
  {
    Error ("calloc_macro: There is a problem in allocating memory for the macro structure\n");
    Exit (0);
  }
  else if (nlevels_macro > 0 || geo.nmacro > 0)
  {
    Log
      ("Allocated %10d bytes for each of %5d elements of macro totaling %10.1f Mb \n",
       sizeof (macro_dummy), (nelem + 1), 1.e-6 * (nelem + 1) * sizeof (macro_dummy));
  }
  else
  {
    Log ("calloc_macro: Allocated no space for macro since nlevels_macro==0\n");
  }

  return (0);
}




/**********************************************************/
/** 
 * @brief      Dynamically allocate various arrays in macromain
 *
 * @param [in] int  nelem   The number of elements in macromain
 * @return     Returns 0, unless the desired memory cannot be 
 * allocated in which the routine exits after an error message
 *
 * @details
 * To minimize the total amount of memory required by Python 
 * in macromode, allocate memory for various arrays in macromain
 * which depend upon the number of ions which are treated as
 * macro atoms
 *
 * ### Notes ###
 *
 **********************************************************/

int
calloc_estimators (nelem)
     int nelem;
{
  int n;

  if (nlevels_macro == 0 && geo.nmacro == 0)
  {
    geo.nmacro = 0;
    Log_silent ("Allocated no space for MA estimators since nlevels_macro==0 and geo.nmacro==0\n");
    return (0);
  }
  //Allocate one extra element to store data where there is no volume


  size_Jbar_est = 0;
  size_gamma_est = 0;
  size_alpha_est = 0;

  for (n = 0; n < nlevels_macro; n++)
  {
    Log_silent
      ("calloc_estimators: level %d has n_bbu_jump %d  n_bbd_jump %d n_bfu_jump %d n_bfd_jump %d\n",
       n, xconfig[n].n_bbu_jump, xconfig[n].n_bbd_jump, xconfig[n].n_bfu_jump, xconfig[n].n_bfd_jump);
    xconfig[n].bbu_indx_first = size_Jbar_est;
    size_Jbar_est += xconfig[n].n_bbu_jump;
    xconfig[n].bfu_indx_first = size_gamma_est;
    size_gamma_est += xconfig[n].n_bfu_jump;
    xconfig[n].bfd_indx_first = size_alpha_est;
    size_alpha_est += xconfig[n].n_bfd_jump;
  }




  Log ("calloc_estimators: size_Jbar_est %d size_gamma_est %d size_alpha_est %d\n", size_Jbar_est, size_gamma_est, size_alpha_est);


  for (n = 0; n < nelem; n++)
  {
    if ((macromain[n].jbar = calloc (sizeof (double), size_Jbar_est)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].jbar_old = calloc (sizeof (double), size_Jbar_est)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].gamma = calloc (sizeof (double), size_gamma_est)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].gamma_old = calloc (sizeof (double), size_gamma_est)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].gamma_e = calloc (sizeof (double), size_gamma_est)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].gamma_e_old = calloc (sizeof (double), size_gamma_est)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].alpha_st = calloc (sizeof (double), size_gamma_est)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].alpha_st_old = calloc (sizeof (double), size_gamma_est)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].alpha_st_e = calloc (sizeof (double), size_gamma_est)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].alpha_st_e_old = calloc (sizeof (double), size_gamma_est)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].recomb_sp = calloc (sizeof (double), size_alpha_est)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].recomb_sp_e = calloc (sizeof (double), size_alpha_est)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].matom_emiss = calloc (sizeof (double), nlevels_macro)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].matom_abs = calloc (sizeof (double), nlevels_macro)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].cooling_bf = calloc (sizeof (double), nphot_total)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].cooling_bf_col = calloc (sizeof (double), nphot_total)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }

    if ((macromain[n].cooling_bb = calloc (sizeof (double), nlines)) == NULL)
    {
      Error ("calloc_estimators: Error in allocating memory for MA estimators\n");
      Exit (0);
    }
  }



  if (nlevels_macro > 0 || geo.nmacro > 0)
  {
    Log
      ("Allocated %10.1f Mb for MA estimators \n",
       1.e-6 * (nelem + 1) * (2. * nlevels_macro + 2. * size_alpha_est + 8. * size_gamma_est + 2. * size_Jbar_est) * sizeof (double));

  }
  else
  {
    Log ("Allocated no space for macro estimators since nlevels_macro==0\n");
  }

  return (0);
}




/**********************************************************/
/** 
 * @brief      This subroutine allocates space for dynamic arrays in the plasma 
 * 	structure
 *
 * @param [in] int  nelem  the number of plasma cells (actually NPLASMA+1) 
 * to allow for empty cell
 * @return     Returns 0, unless the memory cannot be allocated in which case
 * the proram exits
 *
 * @details
 * This subroutine allocates space for variable length arrays in the plasma structure.
 *
 * ### Notes ###
 * Arrays sized to the number of ions are largest,
 * and dominate the size of nplasma, so these were first to be 
 * dynamically allocated.
 *
 **********************************************************/

int
calloc_dyn_plasma (nelem)
     int nelem;
{
  int n;

/*  Loop over all elements in the plasma array, adding one for an empty cell 
 *  used for extrapolations.
 */

  for (n = 0; n < nelem + 1; n++)
  {
    if ((plasmamain[n].density = calloc (sizeof (double), nions)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for density\n");
      Exit (0);
    }
    if ((plasmamain[n].partition = calloc (sizeof (double), nions)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for partition\n");
      Exit (0);
    }
    if ((plasmamain[n].ioniz = calloc (sizeof (double), nions)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for ioniz\n");
      Exit (0);
    }
    if ((plasmamain[n].recomb = calloc (sizeof (double), nions)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for recomb\n");
      Exit (0);
    }
    if ((plasmamain[n].scatters = calloc (sizeof (int), nions)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for scatters\n");
      Exit (0);
    }
    if ((plasmamain[n].xscatters = calloc (sizeof (double), nions)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for xscatters\n");
      Exit (0);
    }
    if ((plasmamain[n].heat_ion = calloc (sizeof (double), nions)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for heat_ion\n");
      Exit (0);
    }
    if ((plasmamain[n].heat_inner_ion = calloc (sizeof (double), nions)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for heat_ion\n");
      Exit (0);
    }
    if ((plasmamain[n].cool_rr_ion = calloc (sizeof (double), nions)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for cool_rr_ion\n");
      Exit (0);
    }
    if ((plasmamain[n].lum_rr_ion = calloc (sizeof (double), nions)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for lum_rr_ion\n");
      Exit (0);
    }
    if ((plasmamain[n].inner_recomb = calloc (sizeof (double), nions)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for inner_recomb\n");
      Exit (0);
    }
    if ((plasmamain[n].inner_ioniz = calloc (sizeof (double), n_inner_tot)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for inner_ioniz\n");
      Exit (0);
    }
    if ((plasmamain[n].cool_dr_ion = calloc (sizeof (double), nions)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for lum_inner_recomb\n");
      Exit (0);
    }


    if ((plasmamain[n].levden = calloc (sizeof (double), nlte_levels)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for levden\n");
      Exit (0);
    }

    if ((plasmamain[n].recomb_simple = calloc (sizeof (double), nphot_total)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for recomb_simple\n");
      Exit (0);
    }

    if ((plasmamain[n].recomb_simple_upweight = calloc (sizeof (double), nphot_total)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for recomb_simple_upweight\n");
      Exit (0);
    }


    if ((plasmamain[n].kbf_use = calloc (sizeof (double), nphot_total)) == NULL)
    {
      Error ("calloc_dyn_plasma: Error in allocating memory for kbf_use\n");
      Exit (0);
    }
  }

  Log
    ("Allocated %10d bytes for each of %5d elements variable length plasma arrays totaling %10.1f Mb \n",
     sizeof (double) * nions * 14, (nelem + 1), 1.e-6 * (nelem + 1) * sizeof (double) * (nions * 14 + nlte_levels + nphot_total * 2));

  return (0);
}



/**********************************************************/
/** 
 * @brief      Dynamically allocate the macro-atom matrixes
 *
 * @param [in] int  nelem   The number of elements in macromain
 * @return     Returns 0, unless the desired memory cannot be 
 * allocated in which the routine exits after an error message
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

int
calloc_matom_matrix (nelem)
     int nelem;
{
  int nrows = nlevels_macro + 1;
  int n;
  int nmatrices_allocated = 0;
  if (nlevels_macro == 0 && geo.nmacro == 0)
  {
    geo.nmacro = 0;
    Log_silent ("Allocated no space for MA matrix since nlevels_macro==0 and geo.nmacro==0\n");
    return (0);
  }

  for (n = 0; n < nelem; n++)
  {
    if (macromain[n].store_matom_matrix == TRUE)
    {
      allocate_macro_matrix (&macromain[n].matom_matrix, nrows);
      nmatrices_allocated += 1;
    }
  }

  if (nlevels_macro > 0 && nmatrices_allocated > 0)
  {
    Log ("Allocated %10.1f Mb for MA matrix \n", 1.e-6 * (nmatrices_allocated + 1) * (nrows * nrows) * sizeof (double));
  }

  return (0);
}

/**********************************************************/
/**
 * @brief  Allocate memory for a square matom_matrix array
 *
 * @param [in, out]  double ***  matrix_addr  The address to the double pointer
 * @param [in]       int         matrix_size  The size of the square matrix
 *
 * @details
 *
 * This will allocate a square matrix of size matrix_size * matrix_size. To
 * use this function,
 *
 * double ***matrix;
 * int matrix_size = 10;
 * allocate_macro_matrix(&matrix, matrix_size);
 *
 * To free memory allocated by this function,
 *
 * free(matrix[0]);
 * free(matrix);
 *
 **********************************************************/

void
allocate_macro_matrix (double ***matrix_addr, int matrix_size)
{
  /* We're doing a trick here to allocate a contiguous chunk of memory
   * for a 2d array. The first step is we allocate nrow pointers, which
   * will be the rows of the matrix. */
  *matrix_addr = calloc (matrix_size, sizeof (double *));
  if (matrix_addr == NULL)
  {
    Error ("allocate_macro_matrix: unable to allocate rows for macro matrix_addr");
    Exit (EXIT_FAILURE);
  }

  /* Now allocate the memory required for the entire matrix on the first row */
  (*matrix_addr)[0] = calloc (matrix_size * matrix_size, sizeof (double));
  if ((*matrix_addr)[0] == NULL)
  {
    Error ("allocate_macro_matrix: unable to allocate elements for macro matrix_addr\n");
    Exit (EXIT_FAILURE);
  }

  /* The final step is to reshape the big allocation on the first row into
   * smaller chunks by using pointer arithmetic to point the pointer in the
   * first allocation to some offset into second allocation. We're basically
   * moving memory around manually. */
  for (int row = 1; row < matrix_size; ++row)
  {
    (*matrix_addr)[row] = (*matrix_addr)[row - 1] + matrix_size;
  }
}
