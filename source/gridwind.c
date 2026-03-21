
/***********************************************************/
/** @file  gridwind.c
 * @author ksl,sss,jm
 * @date   April, 2018
 *
 * @brief  This file contains routines for allocating space for most
 * of the structures needed by Sirocco.  These include 
 * Wind and Plasma structures, and if using macro atoms, the
 * various structures rwquired for this.  The file also contains
 * a routine that maps grid cells in the wind to those in the 
 * Plasma structure, and vice versa.
 *
 *
 * The hierachy of stuctures that describe a model in Sirocco 
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
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/* Convenience macros to pass MPI_Win addresses or NULL depending on MPI mode */
#ifdef MPI_ON
#define PLASMA_WIN(field) plasma_block_ptrs.field
#define MACRO_WIN(field)  macro_block_ptrs.field
#else
/* Use a dummy variable to satisfy the void* parameter when MPI is off */
static int _dummy_win;
#define PLASMA_WIN(field) _dummy_win
#define MACRO_WIN(field)  _dummy_win
#endif


/**********************************************************/
/**
 * @brief  Allocate a contiguous block, using MPI shared memory when available.
 *
 * @param [in]     count       Number of elements to allocate
 * @param [in]     elem_size   Size of each element in bytes
 * @param [out]    ptr         Pointer set to the allocated block (shared or private)
 * @param [out]    win         MPI_Win handle (set only when MPI shared memory is used)
 * @param [in]     use_shared  If TRUE, use MPI_Win_allocate_shared (state/derived);
 *                             if FALSE, use regular calloc (est arrays)
 * @return         0 on success, exits on failure
 *
 * @details
 * When MPI-3 shared memory is available and use_shared is TRUE, only the
 * node leader (node_rank == 0) allocates memory; other ranks on the same
 * node query the leader's pointer via MPI_Win_shared_query.  This means
 * one physical copy per node for read-only (state) and broadcast (derived)
 * data.
 *
 * For estimator arrays (use_shared == FALSE), each rank gets its own
 * private allocation via regular calloc, since estimators are accumulated
 * independently per rank.
 *
 * In non-MPI builds, always uses calloc.
 **********************************************************/

static int
alloc_block_double (long count, double **ptr, void *win_ptr, int use_shared)
{
#ifdef MPI_ON
  if (use_shared && np_mpi_global > 1)
  {
    MPI_Win *win = (MPI_Win *) win_ptr;
    MPI_Aint block_size = (node_rank == 0) ? count * (MPI_Aint) sizeof (double) : 0;
    MPI_Win_allocate_shared (block_size, sizeof (double), MPI_INFO_NULL, node_comm, ptr, win);

    if (node_rank != 0)
    {
      MPI_Aint sz;
      int disp;
      MPI_Win_shared_query (*win, 0, &sz, &disp, ptr);
    }

    /* Zero-initialize the shared block (only leader needs to, but barrier ensures visibility) */
    if (node_rank == 0)
    {
      memset (*ptr, 0, count * sizeof (double));
    }
    MPI_Barrier (node_comm);

    if (*ptr == NULL)
    {
      Error ("alloc_block_double: MPI_Win_allocate_shared returned NULL\n");
      Exit (0);
    }
    return (0);
  }
#endif

  (void) win_ptr;
  (void) use_shared;
  *ptr = calloc (count, sizeof (double));
  if (*ptr == NULL)
  {
    Error ("alloc_block_double: calloc failed for %ld doubles\n", count);
    Exit (0);
  }
  return (0);
}


/**********************************************************/
/**
 * @brief  Allocate a contiguous int block, using MPI shared memory when available.
 *
 * @details Same as alloc_block_double but for int arrays.
 **********************************************************/

static int
alloc_block_int (long count, int **ptr, void *win_ptr, int use_shared)
{
#ifdef MPI_ON
  if (use_shared && np_mpi_global > 1)
  {
    MPI_Win *win = (MPI_Win *) win_ptr;
    MPI_Aint block_size = (node_rank == 0) ? count * (MPI_Aint) sizeof (int) : 0;
    MPI_Win_allocate_shared (block_size, sizeof (int), MPI_INFO_NULL, node_comm, ptr, win);

    if (node_rank != 0)
    {
      MPI_Aint sz;
      int disp;
      MPI_Win_shared_query (*win, 0, &sz, &disp, ptr);
    }

    if (node_rank == 0)
    {
      memset (*ptr, 0, count * sizeof (int));
    }
    MPI_Barrier (node_comm);

    if (*ptr == NULL)
    {
      Error ("alloc_block_int: MPI_Win_allocate_shared returned NULL\n");
      Exit (0);
    }
    return (0);
  }
#endif

  (void) win_ptr;
  (void) use_shared;
  *ptr = calloc (count, sizeof (int));
  if (*ptr == NULL)
  {
    Error ("alloc_block_int: calloc failed for %ld ints\n", count);
    Exit (0);
  }
  return (0);
}


/**********************************************************/
/**
 * @brief  Free a block, using MPI_Win_free for shared blocks or free() for private.
 *
 * @param [in,out]  ptr          Pointer to set to NULL after freeing
 * @param [in]      win_ptr      MPI_Win handle (or NULL for non-shared)
 * @param [in]      is_shared    TRUE if block was allocated with MPI shared memory
 **********************************************************/

static void
free_block (void **ptr, void *win_ptr, int is_shared)
{
  if (*ptr == NULL)
    return;

#ifdef MPI_ON
  if (is_shared && np_mpi_global > 1)
  {
    MPI_Win *win = (MPI_Win *) win_ptr;
    MPI_Win_free (win);
    *ptr = NULL;
    return;
  }
#endif

  (void) win_ptr;
  (void) is_shared;
  free (*ptr);
  *ptr = NULL;
}


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
calloc_wind (int nelem)
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
calloc_plasma (int nelem)
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
check_plasma (PlasmaPtr xplasma, char message[])
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
calloc_macro (int nelem)
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
 * To minimize the total amount of memory required by Sirocco 
 * in macromode, allocate memory for various arrays in macromain
 * which depend upon the number of ions which are treated as
 * macro atoms
 *
 * ### Notes ###
 *
 **********************************************************/

int
calloc_estimators (int nelem)
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


  /* Allocate contiguous blocks for all macro atom dynamic arrays.
   * State and derived arrays use MPI shared memory (one copy per node).
   * Estimator arrays are always private (each rank accumulates independently). */

  int use_shared = FALSE;
#ifdef MPI_ON
  use_shared = (np_mpi_global > 1) ? TRUE : FALSE;
#endif

  /* state arrays — shared across node-local ranks */
  alloc_block_double ((long) nelem * size_Jbar_est, &macro_block_ptrs.jbar_old_block, &MACRO_WIN (win_jbar_old), use_shared);
  alloc_block_double ((long) nelem * size_gamma_est, &macro_block_ptrs.gamma_old_block, &MACRO_WIN (win_gamma_old), use_shared);
  alloc_block_double ((long) nelem * size_gamma_est, &macro_block_ptrs.gamma_e_old_block, &MACRO_WIN (win_gamma_e_old), use_shared);
  alloc_block_double ((long) nelem * size_gamma_est, &macro_block_ptrs.alpha_st_old_block, &MACRO_WIN (win_alpha_st_old), use_shared);
  alloc_block_double ((long) nelem * size_gamma_est, &macro_block_ptrs.alpha_st_e_old_block, &MACRO_WIN (win_alpha_st_e_old), use_shared);

  /* est arrays — always private per rank */
  alloc_block_double ((long) nelem * size_Jbar_est, &macro_block_ptrs.jbar_block, NULL, FALSE);
  alloc_block_double ((long) nelem * size_gamma_est, &macro_block_ptrs.gamma_block, NULL, FALSE);
  alloc_block_double ((long) nelem * size_gamma_est, &macro_block_ptrs.gamma_e_block, NULL, FALSE);
  alloc_block_double ((long) nelem * size_gamma_est, &macro_block_ptrs.alpha_st_block, NULL, FALSE);
  alloc_block_double ((long) nelem * size_gamma_est, &macro_block_ptrs.alpha_st_e_block, NULL, FALSE);
  alloc_block_double ((long) nelem * size_alpha_est, &macro_block_ptrs.recomb_sp_block, NULL, FALSE);
  alloc_block_double ((long) nelem * size_alpha_est, &macro_block_ptrs.recomb_sp_e_block, NULL, FALSE);
  alloc_block_double ((long) nelem * nlevels_macro, &macro_block_ptrs.matom_abs_block, NULL, FALSE);
  alloc_block_double ((long) nelem * nphot_total, &macro_block_ptrs.cooling_bf_block, NULL, FALSE);
  alloc_block_double ((long) nelem * nphot_total, &macro_block_ptrs.cooling_bf_col_block, NULL, FALSE);
  alloc_block_double ((long) nelem * nlines, &macro_block_ptrs.cooling_bb_block, NULL, FALSE);

  /* derived arrays — shared across node-local ranks */
  alloc_block_double ((long) nelem * nlevels_macro, &macro_block_ptrs.matom_emiss_block, &MACRO_WIN (win_matom_emiss), use_shared);

#ifdef MPI_ON
  macro_block_ptrs.shared_memory_active = use_shared;
#endif

  /* Point each cell's pointers into the contiguous blocks */
  for (n = 0; n < nelem; n++)
  {
    macromain[n].state.jbar_old = macro_block_ptrs.jbar_old_block + n * size_Jbar_est;
    macromain[n].state.gamma_old = macro_block_ptrs.gamma_old_block + n * size_gamma_est;
    macromain[n].state.gamma_e_old = macro_block_ptrs.gamma_e_old_block + n * size_gamma_est;
    macromain[n].state.alpha_st_old = macro_block_ptrs.alpha_st_old_block + n * size_gamma_est;
    macromain[n].state.alpha_st_e_old = macro_block_ptrs.alpha_st_e_old_block + n * size_gamma_est;

    macromain[n].est.jbar = macro_block_ptrs.jbar_block + n * size_Jbar_est;
    macromain[n].est.gamma = macro_block_ptrs.gamma_block + n * size_gamma_est;
    macromain[n].est.gamma_e = macro_block_ptrs.gamma_e_block + n * size_gamma_est;
    macromain[n].est.alpha_st = macro_block_ptrs.alpha_st_block + n * size_gamma_est;
    macromain[n].est.alpha_st_e = macro_block_ptrs.alpha_st_e_block + n * size_gamma_est;
    macromain[n].est.recomb_sp = macro_block_ptrs.recomb_sp_block + n * size_alpha_est;
    macromain[n].est.recomb_sp_e = macro_block_ptrs.recomb_sp_e_block + n * size_alpha_est;
    macromain[n].est.matom_abs = macro_block_ptrs.matom_abs_block + n * nlevels_macro;
    macromain[n].est.cooling_bf = macro_block_ptrs.cooling_bf_block + n * nphot_total;
    macromain[n].est.cooling_bf_col = macro_block_ptrs.cooling_bf_col_block + n * nphot_total;
    macromain[n].est.cooling_bb = macro_block_ptrs.cooling_bb_block + n * nlines;

    macromain[n].derived.matom_emiss = macro_block_ptrs.matom_emiss_block + n * nlevels_macro;
  }



  if (nlevels_macro > 0 || geo.nmacro > 0)
  {
    double macro_shared_bytes = (double) nelem * sizeof (double) * (size_Jbar_est + 4.0 * size_gamma_est + nlevels_macro);
    double macro_private_bytes = (double) nelem * sizeof (double) *
      (size_Jbar_est + 4.0 * size_gamma_est + 2.0 * size_alpha_est + nlevels_macro + 2.0 * nphot_total + nlines);
    Log
      ("Macro-atom memory per rank: dynamic shared %.1f MB (one copy/node), dynamic private %.1f MB (per rank)\n",
       1.e-6 * macro_shared_bytes, 1.e-6 * macro_private_bytes);

  }
  else
  {
    Log ("Allocated no space for macro estimators since nlevels_macro==0\n");
  }

  return (0);
}




/**********************************************************/
/**
 * @brief      Allocate contiguous blocks for all dynamic plasma arrays.
 *
 * @param [in] int  nelem  the number of plasma cells (NPLASMA; one extra
 * element is added internally for the empty/dummy cell)
 * @return     Returns 0, unless the memory cannot be allocated in which case
 * the program exits
 *
 * @details
 * Allocates contiguous blocks for variable-length arrays in the plasma structure,
 * then points each cell's sub-struct pointers into the blocks at the correct offset.
 *
 * In MPI builds with np_mpi_global > 1, blocks are allocated using MPI-3
 * shared memory so that only one physical copy exists per node:
 *
 * - **State arrays** (density, partition, levden, recomb_simple, etc.) —
 *   shared across node-local ranks (read-only during photon transport).
 * - **Estimator arrays** (ioniz, heat_ion, heat_inner_ion, inner_ioniz) —
 *   always private per rank (each rank accumulates independently).
 * - **Derived arrays** (recomb, cool_rr_ion, lum_rr_ion, cool_dr_ion,
 *   inner_recomb) — shared across node-local ranks.
 * - **scatters, xscatters** — always private per rank despite being
 *   derived quantities, because they are incremented during photon
 *   transport and would otherwise create a race condition in shared memory.
 *
 * In non-MPI builds or with a single MPI rank, all blocks use regular calloc.
 *
 * ### Notes ###
 * Arrays sized to the number of ions are largest,
 * and dominate the memory footprint of the plasma grid.
 *
 **********************************************************/

int
calloc_dyn_plasma (int nelem)
{
  int n;
  int nelem_alloc = nelem + 1;  /* One extra element for the empty/dummy cell */
  long nalloc_ions = (long) nelem_alloc * nions;
  long nalloc_nlte = (long) nelem_alloc * nlte_levels;
  long nalloc_phot = (long) nelem_alloc * nphot_total;
  long nalloc_inner = (long) nelem_alloc * n_inner_tot;
  int use_shared = FALSE;
  int was_shared = FALSE;

#ifdef MPI_ON
  use_shared = (np_mpi_global > 1) ? TRUE : FALSE;
  was_shared = plasma_block_ptrs.shared_memory_active;
#endif

  /* Free any previously allocated blocks */
  if (plasma_block_ptrs.density_block != NULL)
  {
    free_block ((void **) &plasma_block_ptrs.density_block, &PLASMA_WIN (win_density), was_shared);
    free_block ((void **) &plasma_block_ptrs.partition_block, &PLASMA_WIN (win_partition), was_shared);
    free_block ((void **) &plasma_block_ptrs.levden_block, &PLASMA_WIN (win_levden), was_shared);
    free_block ((void **) &plasma_block_ptrs.recomb_simple_block, &PLASMA_WIN (win_recomb_simple), was_shared);
    free_block ((void **) &plasma_block_ptrs.recomb_simple_upweight_block, &PLASMA_WIN (win_recomb_simple_upweight), was_shared);
    free_block ((void **) &plasma_block_ptrs.kbf_use_block, &PLASMA_WIN (win_kbf_use), was_shared);
    free_block ((void **) &plasma_block_ptrs.ioniz_block, NULL, FALSE);
    free_block ((void **) &plasma_block_ptrs.heat_ion_block, NULL, FALSE);
    free_block ((void **) &plasma_block_ptrs.heat_inner_ion_block, NULL, FALSE);
    free_block ((void **) &plasma_block_ptrs.inner_ioniz_block, NULL, FALSE);
    free_block ((void **) &plasma_block_ptrs.recomb_block, &PLASMA_WIN (win_recomb), was_shared);
    free_block ((void **) &plasma_block_ptrs.scatters_block, NULL, FALSE);
    free_block ((void **) &plasma_block_ptrs.xscatters_block, NULL, FALSE);
    free_block ((void **) &plasma_block_ptrs.cool_rr_ion_block, &PLASMA_WIN (win_cool_rr_ion), was_shared);
    free_block ((void **) &plasma_block_ptrs.lum_rr_ion_block, &PLASMA_WIN (win_lum_rr_ion), was_shared);
    free_block ((void **) &plasma_block_ptrs.cool_dr_ion_block, &PLASMA_WIN (win_cool_dr_ion), was_shared);
    free_block ((void **) &plasma_block_ptrs.inner_recomb_block, &PLASMA_WIN (win_inner_recomb), was_shared);
    free_block ((void **) &plasma_block_ptrs.state_xbands_dblock, &PLASMA_WIN (win_state_xbands_d), was_shared);
    free_block ((void **) &plasma_block_ptrs.state_spec_mod_type_block, &PLASMA_WIN (win_state_spec_mod_type), was_shared);
    free_block ((void **) &plasma_block_ptrs.derived_persist_force_block, &PLASMA_WIN (win_derived_persist_force), was_shared);
    free_block ((void **) &plasma_block_ptrs.derived_persist_angle_block, &PLASMA_WIN (win_derived_persist_angle), was_shared);
  }

  /* Allocate contiguous blocks for all dynamic plasma arrays.
   * State and derived arrays use MPI shared memory (one copy per node).
   * Estimator arrays are always private (each rank accumulates independently). */

  /* state arrays — shared across node-local ranks */
  alloc_block_double (nalloc_ions, &plasma_block_ptrs.density_block, &PLASMA_WIN (win_density), use_shared);
  alloc_block_double (nalloc_ions, &plasma_block_ptrs.partition_block, &PLASMA_WIN (win_partition), use_shared);
  alloc_block_double (nalloc_nlte, &plasma_block_ptrs.levden_block, &PLASMA_WIN (win_levden), use_shared);
  alloc_block_double (nalloc_phot, &plasma_block_ptrs.recomb_simple_block, &PLASMA_WIN (win_recomb_simple), use_shared);
  alloc_block_double (nalloc_phot, &plasma_block_ptrs.recomb_simple_upweight_block, &PLASMA_WIN (win_recomb_simple_upweight), use_shared);
  /* kbf_use is int* but historically allocated/written as doubles for binary compat */
  alloc_block_double (nalloc_phot, &plasma_block_ptrs.kbf_use_block, &PLASMA_WIN (win_kbf_use), use_shared);

  /* est arrays — always private per rank */
  alloc_block_double (nalloc_ions, &plasma_block_ptrs.ioniz_block, NULL, FALSE);
  alloc_block_double (nalloc_ions, &plasma_block_ptrs.heat_ion_block, NULL, FALSE);
  alloc_block_double (nalloc_ions, &plasma_block_ptrs.heat_inner_ion_block, NULL, FALSE);
  alloc_block_double (nalloc_inner, &plasma_block_ptrs.inner_ioniz_block, NULL, FALSE);

  /* derived arrays — shared across node-local ranks (except scatters/xscatters
   * which are written during photon transport and must be private per rank) */
  alloc_block_double (nalloc_ions, &plasma_block_ptrs.recomb_block, &PLASMA_WIN (win_recomb), use_shared);
  alloc_block_int (nalloc_ions, &plasma_block_ptrs.scatters_block, NULL, FALSE);
  alloc_block_double (nalloc_ions, &plasma_block_ptrs.xscatters_block, NULL, FALSE);
  alloc_block_double (nalloc_ions, &plasma_block_ptrs.cool_rr_ion_block, &PLASMA_WIN (win_cool_rr_ion), use_shared);
  alloc_block_double (nalloc_ions, &plasma_block_ptrs.lum_rr_ion_block, &PLASMA_WIN (win_lum_rr_ion), use_shared);
  alloc_block_double (nalloc_ions, &plasma_block_ptrs.cool_dr_ion_block, &PLASMA_WIN (win_cool_dr_ion), use_shared);
  alloc_block_double (nalloc_ions, &plasma_block_ptrs.inner_recomb_block, &PLASMA_WIN (win_inner_recomb), use_shared);

  /* state fixed-size arrays — combined contiguous blocks (shared) */
  {
    int state_xbands_stride = 6 * NXBANDS + 2 * (NXBANDS + 1);
    alloc_block_double ((long) nelem_alloc * state_xbands_stride, &plasma_block_ptrs.state_xbands_dblock,
                        &PLASMA_WIN (win_state_xbands_d), use_shared);
  }
  alloc_block_int ((long) nelem_alloc * NXBANDS, &plasma_block_ptrs.state_spec_mod_type_block,
                   &PLASMA_WIN (win_state_spec_mod_type), use_shared);

  /* derived fixed-size arrays — combined contiguous blocks (shared) */
  alloc_block_double ((long) nelem_alloc * 6 * NFORCE_DIRECTIONS, &plasma_block_ptrs.derived_persist_force_block,
                      &PLASMA_WIN (win_derived_persist_force), use_shared);
  alloc_block_double ((long) nelem_alloc * 3 * NFLUX_ANGLES, &plasma_block_ptrs.derived_persist_angle_block,
                      &PLASMA_WIN (win_derived_persist_angle), use_shared);

#ifdef MPI_ON
  plasma_block_ptrs.shared_memory_active = use_shared;
#endif

  /* Point each cell's pointers into the contiguous blocks at the right offset */
  for (n = 0; n < nelem_alloc; n++)
  {
    plasmamain[n].state.density = plasma_block_ptrs.density_block + n * nions;
    plasmamain[n].state.partition = plasma_block_ptrs.partition_block + n * nions;
    plasmamain[n].state.levden = plasma_block_ptrs.levden_block + n * nlte_levels;
    plasmamain[n].state.recomb_simple = plasma_block_ptrs.recomb_simple_block + n * nphot_total;
    plasmamain[n].state.recomb_simple_upweight = plasma_block_ptrs.recomb_simple_upweight_block + n * nphot_total;
    plasmamain[n].state.kbf_use = (int *) (plasma_block_ptrs.kbf_use_block + n * nphot_total);

    plasmamain[n].est.ioniz = plasma_block_ptrs.ioniz_block + n * nions;
    plasmamain[n].est.heat_ion = plasma_block_ptrs.heat_ion_block + n * nions;
    plasmamain[n].est.heat_inner_ion = plasma_block_ptrs.heat_inner_ion_block + n * nions;
    plasmamain[n].est.inner_ioniz = plasma_block_ptrs.inner_ioniz_block + n * n_inner_tot;

    plasmamain[n].derived.recomb = plasma_block_ptrs.recomb_block + n * nions;
    plasmamain[n].derived.scatters = plasma_block_ptrs.scatters_block + n * nions;
    plasmamain[n].derived.xscatters = plasma_block_ptrs.xscatters_block + n * nions;
    plasmamain[n].derived.cool_rr_ion = plasma_block_ptrs.cool_rr_ion_block + n * nions;
    plasmamain[n].derived.lum_rr_ion = plasma_block_ptrs.lum_rr_ion_block + n * nions;
    plasmamain[n].derived.cool_dr_ion = plasma_block_ptrs.cool_dr_ion_block + n * nions;
    plasmamain[n].derived.inner_recomb = plasma_block_ptrs.inner_recomb_block + n * nions;

    /* State spectral block: stride = 6*NXBANDS + 2*(NXBANDS+1) = 162 */
    {
      int state_xbands_stride = 6 * NXBANDS + 2 * (NXBANDS + 1);
      double *sbase = plasma_block_ptrs.state_xbands_dblock + n * state_xbands_stride;
      plasmamain[n].state.f1 = sbase;
      plasmamain[n].state.f2 = sbase + (NXBANDS + 1);
      plasmamain[n].state.pl_alpha = sbase + 2 * (NXBANDS + 1);
      plasmamain[n].state.pl_log_w = sbase + 2 * (NXBANDS + 1) + NXBANDS;
      plasmamain[n].state.exp_temp = sbase + 2 * (NXBANDS + 1) + 2 * NXBANDS;
      plasmamain[n].state.exp_w = sbase + 2 * (NXBANDS + 1) + 3 * NXBANDS;
      plasmamain[n].state.fmin_mod = sbase + 2 * (NXBANDS + 1) + 4 * NXBANDS;
      plasmamain[n].state.fmax_mod = sbase + 2 * (NXBANDS + 1) + 5 * NXBANDS;
    }
    plasmamain[n].state.spec_mod_type = (enum spec_mod_type_enum *) (plasma_block_ptrs.state_spec_mod_type_block + n * NXBANDS);

    /* Derived persistent force block: stride = 6 * NFORCE_DIRECTIONS */
    {
      double *dbase = plasma_block_ptrs.derived_persist_force_block + n * 6 * NFORCE_DIRECTIONS;
      plasmamain[n].derived.F_vis_persistent = dbase;
      plasmamain[n].derived.F_UV_persistent = dbase + NFORCE_DIRECTIONS;
      plasmamain[n].derived.F_Xray_persistent = dbase + 2 * NFORCE_DIRECTIONS;
      plasmamain[n].derived.rad_force_es_persist = dbase + 3 * NFORCE_DIRECTIONS;
      plasmamain[n].derived.rad_force_ff_persist = dbase + 4 * NFORCE_DIRECTIONS;
      plasmamain[n].derived.rad_force_bf_persist = dbase + 5 * NFORCE_DIRECTIONS;
    }

    /* Derived persistent angle block: stride = 3 * NFLUX_ANGLES */
    {
      double *dbase = plasma_block_ptrs.derived_persist_angle_block + n * 3 * NFLUX_ANGLES;
      plasmamain[n].derived.F_UV_ang_theta_persist = dbase;
      plasmamain[n].derived.F_UV_ang_phi_persist = dbase + NFLUX_ANGLES;
      plasmamain[n].derived.F_UV_ang_r_persist = dbase + 2 * NFLUX_ANGLES;
    }
  }

  /* Report memory breakdown: shared (one copy per node) vs private (per rank) vs base struct */
  {
    double shared_bytes =
      (double) nelem_alloc * sizeof (double) * (2.0 * nions + nlte_levels + 3.0 * nphot_total + 5.0 * nions + n_inner_tot);
    double private_bytes = (double) nelem_alloc * sizeof (double) * (3.0 * nions + n_inner_tot) + (double) nelem_alloc * sizeof (int) * nions;  /* scatters(int) + xscatters(double) already in shared_bytes above... */

    /* Recalculate properly:
     * Shared state: density(nions) + partition(nions) + levden(nlte) + recomb_simple(nphot) + recomb_simple_upweight(nphot) + kbf_use(nphot) = 2*nions + nlte + 3*nphot
     *   + state_xbands_d(6*NXBANDS + 2*(NXBANDS+1)) + spec_mod_type(NXBANDS ints)
     * Shared derived: recomb(nions) + cool_rr_ion(nions) + lum_rr_ion(nions) + cool_dr_ion(nions) + inner_recomb(nions) = 5*nions
     *   + persist_force(6*NFORCE_DIRECTIONS) + persist_angle(3*NFLUX_ANGLES)
     * Private est: ioniz(nions) + heat_ion(nions) + heat_inner_ion(nions) + inner_ioniz(n_inner) = 3*nions + n_inner
     * Private derived: scatters(nions, int) + xscatters(nions, double) = nions*(4+8) */
    shared_bytes = (double) nelem_alloc *(sizeof (double) * (2.0 * nions + nlte_levels + 3.0 * nphot_total + 5.0 * nions
                                                             + 6.0 * NXBANDS + 2.0 * (NXBANDS + 1)
                                                             + 6.0 * NFORCE_DIRECTIONS + 3.0 * NFLUX_ANGLES) + sizeof (int) * NXBANDS);
    private_bytes = (double) nelem_alloc *(sizeof (double) * (3.0 * nions + n_inner_tot + nions) + sizeof (int) * nions);
    double base_struct_bytes = (double) nelem_alloc * sizeof (plasma_dummy);

    Log
      ("Plasma memory per rank: base struct %.1f MB, dynamic shared %.1f MB (one copy/node), dynamic private %.1f MB (per rank)\n",
       1.e-6 * base_struct_bytes, 1.e-6 * shared_bytes, 1.e-6 * private_bytes);
    Log
      ("  nions=%d, nlte_levels=%d, nphot_total=%d, n_inner_tot=%d, nelem=%d\n",
       nions, nlte_levels, nphot_total, n_inner_tot, nelem_alloc - 1);
  }

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
calloc_matom_matrix (int nelem)
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
    if (macromain[n].state.store_matom_matrix == TRUE)
    {
      allocate_macro_matrix (&macromain[n].derived.matom_matrix, nrows);
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
