/***********************************************************/
/**
 * @file define_wind.c
 * @author EJP
 * @date November 2023
 *
 * @brief  Initialise the wind in parallel
 *
 * @details
 *
 * All the functions, other than define_wind, are static.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

/**********************************************************/
/**
 * @brief Calculate an approximate value for the mass loss of the wind
 *
 * @details
 *
 * ### Notes ###
 *
 * This checks if m_dot is correct. The calculation is very approximate. It only
 * calculates m_dot flowing up through the grid, i.e. any material that leaks
 * out of the side will not be counted. A much better approach is needed, e.g.
 * to calculate the mass flow in a spherical shell at the edge of the disk. In
 * addition, it  uses a very inexact way to determine which grid points are in
 * the wind.
 *
 * This routine is only guaranteed to work with spherical coordinate systems,
 * but will still output some numbers for other coordinate systems. It should be
 * possible to get the correct value for all coordinate systems by integrating
 * using interpolated values, rather than the analytical form used below.
 *
 * It's unclear if this really works with domains, but the following code
 * compiles.
 *
 **********************************************************/

static void
calculate_mdot_wind (void)
{
  int i;
  int n;
  int ndom;
  int ndim;
  int mdim;
  int nstart;
  int nplasma;

  double mdotwind;
  double mdotbase;
  double rr;

  WindPtr w = wmain;

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    if (zdom[ndom].coord_type == CYLIND)
    {
      Log ("Warning: wind2d.c: next lines may not make sense with multiple domains\n");
      ndim = zdom[ndom].ndim;
      mdim = zdom[ndom].mdim;
      nstart = zdom[ndom].nstart;

      mdotwind = 0;
      mdotbase = 0;
      for (i = 0; i < ndim - 1; i++)
      {
        n = nstart + i * mdim + (mdim / 2);
        /* rr is router * router - rinner * rinner for this shell */
        rr = w[n + mdim].x[0] * w[n + mdim].x[0] + w[n + mdim].x[1] * w[n + mdim].x[1] - (w[n].x[0] * w[n].x[0] + w[n].x[1] * w[n].x[1]);
        if (w[nstart + i * mdim].inwind == W_ALL_INWIND)
        {
          nplasma = w[nstart + i * mdim].nplasma;
          mdotbase += plasmamain[nplasma].rho * PI * rr * w[nstart + i * mdim].v[2];
        }
        if (w[n].inwind == W_ALL_INWIND)
        {
          nplasma = w[n].nplasma;
          mdotwind += plasmamain[nplasma].rho * PI * rr * w[n].v[2];
        }
      }

      /* Factor of 2 to allow for wind on both sides of disk */
      Log ("m dot wind: Domain %d Desired %g   Base %g Calculated %g\n", ndom, zdom[ndom].wind_mdot, 2. * mdotbase, 2. * mdotwind);
      mdot_wind (w, 1.e6, zdom[ndom].rmax / 2.);
    }
    else
    {
      Error ("wind2d.c: Not currently able to calculate mdot wind for coord type %d in domain %d\n", zdom[ndom].coord_type, ndom);
    }
  }
}

/**********************************************************/
/**
 * @brief Initialise the basic properties for the macro atom grid
 *
 * @details
 *
 * This should only be called after the wind grid has been created, and when
 * macro atoms are being used. If no macro atoms are present, then there will
 * be a failed allocation.
 *
 **********************************************************/

static void
create_macro_grid (void)
{
  int n_plasma;
  int n_start;
  int n_stop;

  if (geo.rt_mode != RT_MODE_MACRO)     /* Let's be extra safe */
  {
    Error ("Trying to initialise the macro grid without macro atoms\n");
    return;
  }

  calloc_macro (NPLASMA);
  calloc_estimators (NPLASMA);

  /* At this point in time, there is no need to parallelise this step as all
   * we need to do is set the value of two integers for each macro cell. The
   * overhead associated with parallelisation and communication surely far
   * exceeds the time saved doing this in parallel */
  n_start = 0;
  n_stop = NPLASMA;

  for (n_plasma = n_start; n_plasma < n_stop; ++n_plasma)
  {
    macromain[n_plasma].store_matom_matrix = modes.store_matom_matrix;
    macromain[n_plasma].matom_transition_mode = geo.matom_transition_mode;
  }

  calloc_matom_matrix (NPLASMA);
}

/**********************************************************/
/**
 * @brief Initialise the spectral models for each band
 *
 * @details
 *
 **********************************************************/

static void
set_spectral_models (PlasmaPtr cell)
{
  int n_band;

  for (n_band = 0; n_band < NXBANDS; n_band++)
  {
    cell->spec_mod_type[n_band] = SPEC_MOD_FAIL;        /*NSH 120817 - setting this to
                                                           a negative number means that at the outset, we assume we do not have a
                                                           suitable model for the cell */
    cell->exp_temp[n_band] = geo.tmax;  /*NSH 120817 - as an initial guess,
                                           set this number to the hottest part of the model -
                                           this should define where any exponential dropoff becomes important */
    cell->exp_w[n_band] = 0.0;  /* 120817 Who knows what this should be! */
    cell->pl_alpha[n_band] = geo.alpha_agn;     /*Awind2d: For domains an initial guess we assume the whole wind is
                                                   optically thin and so the spectral index for a PL illumination will be the
                                                   same everywhere.  */
    cell->pl_log_w[n_band] = -1e99;     /*131114 - a tiny weight - just to fill the variable */
    cell->fmin_mod[n_band] = 1e99;      /* Set the minium model frequency to the max frequency in the band - means it will never be used which is correct at this time - there is no model */
    cell->fmax_mod[n_band] = 1e-99;     /* Set the maximum model frequency to the min frequency in the band */
  }
}

/**********************************************************/
/**
 * @brief Set the temperature of a plasma cell.
 *
 * @details
 *
 * The temperature is either gotten from a hydro file, from an import file or
 * from the parameter file. The temperature in the parameter file is the
 * radiation temperature. In all winds other than imported models, the electron
 * temperature is either fixed (set to the temperature in the parameter file) or
 * set using the Lucy approximation: T_e = 0.9 T_r.
 *
 **********************************************************/

static void
set_plasma_temperature (PlasmaPtr cell, int ndom, double xcen[3])
{
  DomainPtr dom = &zdom[ndom];

  if (dom->wind_type == HYDRO)
  {
    cell->t_r = hydro_temp (xcen);
  }
  else if (dom->wind_type == IMPORT)
  {
    cell->t_r = import_temperature (ndom, xcen, FALSE);
  }
  else                          /* Taken from parameter file */
  {
    cell->t_r = dom->twind;
  }

  /* Initialize variables having to do with convergence in initial stages */
  cell->gain = 0.5;
  cell->dt_e_old = 0.0;
  cell->dt_e = 0.0;

  /* The next lines set the electron temperature to 0.9 times the radiation temperature (which is a bit
     odd since the input is the wind temperature, but is taken to be the radiation temperature). If we have
     a fixed temperature calculation,then the wind temperature is set to be the wind temperature so the
     user gets what they are expecting */
  if (dom->wind_type == IMPORT)
  {
    cell->t_e = import_temperature (ndom, xcen, TRUE);
  }
  else if (modes.fixed_temp == FALSE && modes.zeus_connect == FALSE)
  {
    cell->t_e = cell->t_e_old = 0.9 * cell->t_r;        // Lucy guess
  }
  else                          //If we want to fix the temperature, we set it to tr which has previously been set to twind
  {
    cell->t_e = cell->t_e_old = cell->t_r;
  }
}

/**********************************************************/
/**
 * @brief Initialise the properties of the plasma grid
 *
 * @details
 *
 * This has to be done after creation of the wind grid, otherwise NPLASMA will
 * not be populated nor will we baby able to get some properties from the
 * wind grid (such as the velocity).
 *
 **********************************************************/

static void
create_plasma_grid (void)
{
  int nion;
  int n_plasma;
  int ndom;
  int nwind;
  int ierr;
  int n_start;
  int n_stop;
  int n_cells_rank;
  double rrstar;
  double xcen[3];
  WindPtr w = wmain;

  if (NPLASMA <= 0)             /* Let's be extra safe */
  {
    Error ("Invalid value for number of plasma cells %d\n", NPLASMA);
    Exit (EXIT_FAILURE);
  }

#ifdef MPI_ON
  n_cells_rank = get_parallel_nrange (rank_global, NPLASMA, np_mpi_global, &n_start, &n_stop);
#else
  n_start = 0;
  n_stop = NPLASMA;
  n_cells_rank = NPLASMA;
#endif

  calloc_plasma (NPLASMA);
  calloc_dyn_plasma (NPLASMA);
  create_wind_and_plasma_cell_maps ();

  for (n_plasma = n_start; n_plasma < n_stop; n_plasma++)
  {
    nwind = plasmamain[n_plasma].nwind;
    ndom = wmain[nwind].ndom;
    stuff_v (w[nwind].xcen, xcen);

    plasmamain[n_plasma].xgamma = wmain[nwind].xgamma_cen;
    plasmamain[n_plasma].rho = model_rho (ndom, xcen) / (zdom[ndom].fill * plasmamain[n_plasma].xgamma);
    plasmamain[n_plasma].vol = w[nwind].vol * zdom[ndom].fill;

    /* This is where we initialise the spectral models for the wind. */
    set_spectral_models (&plasmamain[n_plasma]);

    /* Initialise cell temperatures */
    set_plasma_temperature (&plasmamain[n_plasma], ndom, xcen);

    /* Calculate an initial guess for the weight of the PL spectrum (constant / area of a sphere / 4pi) */
    rrstar = 1.0 - (geo.rstar * geo.rstar) / (xcen[0] * xcen[0] + xcen[1] * xcen[1] + xcen[2] * xcen[2]);
    if (rrstar > 0)
    {
      plasmamain[n_plasma].w = 0.5 * (1 - sqrt (rrstar));
    }
    else
    {                           /* Allow for possibility that grid point is inside star */
      plasmamain[n_plasma].w = 0.5;
    }

    /* Determine the initial ion abundances, using either LTE or fixed concentrations */
    if (geo.ioniz_mode != IONMODE_FIXED)
    {                           /* Carry out an LTE determination of the ionization (using T_r) */
      ierr = ion_abundances (&plasmamain[n_plasma], IONMODE_LTE_TR);
    }
    else
    {                           /* Set the concentrations to specified values */
      ierr = ion_abundances (&plasmamain[n_plasma], IONMODE_FIXED);
    }
    if (ierr != 0)
    {
      Error
        ("wind_define after ion_abundances: cell %d rho %8.2e t_r %8.2 t_e %8.2e w %8.2e\n",
         n_plasma, plasmamain[n_plasma].rho, plasmamain[n_plasma].t_r, plasmamain[n_plasma].t_e, plasmamain[n_plasma].w);
    }

    /* Initialize arrays for tracking number of scatters in the cell */
    for (nion = 0; nion < nions; nion++)
    {
      plasmamain[n_plasma].scatters[nion] = 0;
      plasmamain[n_plasma].xscatters[nion] = 0;
    }

    /* Calculate adiabatic and non-thermal heating contributions */
    if (geo.adiabatic)
    {
      plasmamain[n_plasma].cool_adiabatic = adiabatic_cooling (&w[nwind], plasmamain[n_plasma].t_e);
    }
    else
    {
      plasmamain[n_plasma].cool_adiabatic = 0.0;
    }
    if (geo.nonthermal)
    {
      nwind = plasmamain[n_plasma].nwind;
      plasmamain[n_plasma].heat_shock = shock_heating (&w[nwind]);
    }
    else
    {
      plasmamain[n_plasma].heat_shock = 0.0;
    }
  }

  /* zero the counters which record diagnostics for mean_intensity */
  nerr_Jmodel_wrong_freq = 0;
  nerr_no_Jmodel = 0;

  broadcast_plasma_grid (n_start, n_stop, n_cells_rank);
}

/**********************************************************/
/**
 * @brief  Create the wind (and etc.) coordinate grid for the given
 *         wind or coordinate system specified
 *
 * @details
 *
 * This function will create the coordinate grid depending on the type of the
 * wind or the coordinate type. For the imported, shell and hydro wind types,
 * Python has special functions to make the coordinate grid. For all other wind
 * types, Python uses a more generic approach to create the coordinate grid.
 *
 **********************************************************/

static void
make_coordinate_grid (void)
{
  int ndom;

  /* This is done on domain-by-domain basis first */
  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    /* Checking first for two unique cases of wind type */
    if (zdom[ndom].wind_type == IMPORT)
    {
      import_make_grid (ndom, wmain);
    }
    else if (zdom[ndom].wind_type == SHELL)
    {
      shell_make_grid (ndom, wmain);
    }
    else if (zdom[ndom].wind_type == HYDRO)
    {
      rtheta_make_hydro_grid (ndom, wmain);
    }
    /* Now we can check for coord types to define the coordinates */
    else if (zdom[ndom].coord_type == SPHERICAL)
    {
      spherical_make_grid (ndom, wmain);
    }
    else if (zdom[ndom].coord_type == CYLIND)
    {
      cylind_make_grid (ndom, wmain);
    }
    else if (zdom[ndom].coord_type == RTHETA)
    {
      rtheta_make_grid (ndom, wmain);
    }
    else if (zdom[ndom].coord_type == CYLVAR)
    {
      cylvar_make_grid (ndom, wmain);
    }
    else
    {
      Error ("make_coordinate_grid: unknown wind or coordinate type\n");
      Exit (EXIT_FAILURE);
    }
  }
}

/**********************************************************/
/**
 * @brief Compute the volumes of each cell in the wind grid.
 *
 * @details
 *
 * Compute the volume of each cell given the specific coordinate system for
 * each domain.
 *
 **********************************************************/

static void
calculate_cell_volume (WindPtr cell)
{
  int ndom;

  ndom = cell->ndom;

  if (zdom[ndom].coord_type == SPHERICAL)
  {
    spherical_cell_volume (cell);
  }
  else if (zdom[ndom].coord_type == CYLIND)
  {
    cylind_cell_volume (cell);
  }
  else if (zdom[ndom].coord_type == CYLVAR)
  {
    cylvar_cell_volume (cell);
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {
    if (zdom[ndom].wind_type == HYDRO)
    {
      rtheta_hydro_cell_volume (cell);
    }
    else
    {
      rtheta_cell_volume (cell);
    }
  }
  else
  {
    Error ("calculate_cell_volume: unknown coordinate type %d for cell %d in domain %d\n", zdom[ndom].coord_type, cell->nwind, ndom);
    Exit (EXIT_FAILURE);
  }
}

/**********************************************************/
/**
 * @brief Complete the creation of the wind grid
 *
 * @return
 *
 * @details
 *
 * The main purpose of this function is to find out how many plasma cells there
 * should be (cells with non-zero volume), and to check that the properties of
 * the wind cells make sense.
 *
 **********************************************************/

static void
complete_wind_grid_creation (void)
{
  int i;
  int j;
  int n;
  int ndom;
  int n_vol;
  int n_inwind;
  int n_part;

  WindPtr w = wmain;

  /* Count the number of cells which have volume and catch cases where we have
   * made a mistake and made a grid/model with no cells in the wind or cells
   * without volume */
  NPLASMA = 0;
  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    n_vol = 0;
    n_inwind = 0;
    n_part = 0;
    for (n = zdom[ndom].nstart; n < zdom[ndom].nstop; n++)
    {
      if (w[n].vol > 0.0)
      {
        n_vol++;
      }
      if (w[n].inwind == W_ALL_INWIND)
      {
        n_inwind++;
      }
      if (w[n].inwind == W_PART_INWIND)
      {
        n_part++;
      }
    }
    NPLASMA += n_vol;
    Log
      ("wind2d: For domain %d there are %3d cells of which %d are in inwind, %d partially in_wind, & %d with pos. vol\n",
       ndom, zdom[ndom].ndim2, n_inwind, n_part, n_vol);

    if (n_inwind == 0)
    {
      Error ("wind2d: There are no wind cells in domain %d.  This seems unreasonble\n", ndom);
      if (n_part > 0)
      {
        Log ("wind2d: partially in wind cells have been found, so consider increasing resolution\n");
      }
      Exit (EXIT_FAILURE);
    }

    /* Next we need to check that each cell which has volume has in the wind or
     * not. This can sometimes happen with odd wind parameters. */
    if (zdom[ndom].coord_type != SPHERICAL && zdom[ndom].wind_type != IMPORT)
    {
      for (n = zdom[ndom].nstart; n < zdom[ndom].nstop; n++)
      {
        n_inwind = check_corners_inwind (n);

        if (w[n].vol == 0 && n_inwind > 0)
        {
          wind_n_to_ij (ndom, n, &i, &j);
          Error ("wind2d: Cell %3d (%2d,%2d) in domain %d has %d corners in wind, but zero volume\n", n, i, j, ndom, n_inwind);
          w[n].inwind = W_IGNORE;
        }
        if (w[n].inwind == W_PART_INWIND && n_inwind == 4)
        {
          wind_n_to_ij (ndom, n, &i, &j);
          Error ("wind2d: Cell %3d (%2d,%2d) in domain %d has 4 corners in wind, but is only partially in wind\n", n, i, j, ndom);
        }
      }
    }
    else if (zdom[ndom].coord_type == SPHERICAL)
    {
      Log ("wind2d: Not checking corners_in_wind for SPHERICAL coordinates used in domain %d\n", ndom);
    }
    else if (zdom[ndom].wind_type == IMPORT)
    {
      Log ("wind2d: Not checking corners_in_wind for imported model, domain %d\n", ndom);
    }
  }

  /* This is a final check to make sure everything is reasonable. What this
   * basically does is check for super-luminal cells, as well as NaNs or INFs
   * for the various properties of the wind */
  wind_check (wmain, -1);
}

/**********************************************************/
/**
 * @brief Initialise the properties of the wind grid
 *
 * @details
 *
 * This function deals with all of the steps to create the grid coordinate
 * system and to define the volume and velocity of each cell. Parts of the
 * domain structure are also defined in here, as well as the number of plasma
 * cells (e.g. cells with non-zero volume).
 *
 * ### Notes ###
 *
 * To speed up grid creation, make_coordinate_grid(), calculate_cell_volume()
 * and define_wind_velocities() are all parallelised.
 *
 * TODO: loops over ndomain need to be refactored into loops over NDIM2
 *
 **********************************************************/

static void
create_wind_grid (void)
{
  int n;
  int ndom;
  int n_start;
  int n_stop;
  int n_cells_rank;
  double v_cen[3];
  WindPtr cell;

  /* Set up indices for starting and ending positions of each wind
   * domain in the wind grid. Note that zdom[ndom].ndim and zdom[ndom].mdim
   * should have been already established in `get_grid_params`. The wind grid
   * hasn't been allocated just yet, as we are counting the number of cells in
   * the loop which need to be allocated */
  n = 0;
  for (ndom = 0; ndom < geo.ndomain; ++ndom)
  {
    zdom[ndom].nstart = n;
    n += zdom[ndom].ndim * zdom[ndom].mdim;
    zdom[ndom].nstop = n;
    NDIM2 += zdom[ndom].ndim * zdom[ndom].mdim;
  }
  geo.ndim2 = NDIM2;

  calloc_wind (NDIM2);

  /* Assign the domain for each cell in the wind grid */
  int offset = 0;
  for (ndom = 0; ndom < geo.ndomain; ++ndom)
  {
    for (n = zdom[ndom].nstart; n < zdom[ndom].nstop; ++n)
    {
      wmain[n].ndom = ndom;
      wmain[n].inwind = W_NOT_ASSIGNED;
      wmain[n].dfudge = DFUDGE;
      wmain[n].nwind = n + offset;
      wmain[n].nwind_dom = n;
    }
    offset += zdom[ndom].ndim;
  }

  /* The first thing we need to do is to create the coordinate grid. We'll do
   * this in serial, as it's very difficult to untangle this process into
   * something done in parallel due to data dependencies and how certain/special
   * wind types are set up in differing ways. At the moment, it is not worth the
   * overhead or human cost to do this stage (and the next) in parallel */
  make_coordinate_grid ();

  /* With the grid defined, we need to "complete" it by populating some 1d
   * array of coordinates which are used in the next step to compute some other
   * properties, namely the volume and boundaries of the grid. This bit is
   * also tangled up in other serial parts of the code and is also very short.
   * As above, it is impractical and not worthwhile to parallelise this step */
  wind_complete ();

#ifdef MPI_ON
  n_cells_rank = get_parallel_nrange (rank_global, NDIM2, np_mpi_global, &n_start, &n_stop);
#else
  n_start = 0;
  n_stop = NDIM2;
  n_cells_rank = NDIM2;
#endif

  /* The next stages are done in parallel, as calculating volumes and some
   * velocity gradients is expensive. Within this loop we:
   *  - create the coordinate grid
   *  - compute the cell volumes
   *  - define the velocity, velocity gradients and divergence at cell vertex
   *  - compute the average and max dv/ds in each cell
   *  - initialise the gamma factor in each cell
   * */
  for (n = n_start; n < n_stop; ++n)
  {
    cell = &wmain[n];
    calculate_cell_volume (cell);

    if (zdom[cell->ndom].wind_type != IMPORT)
    {
      model_velocity (cell->ndom, cell->x, cell->v);
    }

    model_vgrad (cell->ndom, cell->x, cell->v_grad);
    wind_div_v (cell->ndom, cell);

    /* Now we do two expensive calculations to figure out the direction of
     * the largest velocity gradient in the cell as well as the angle average
     * velocity gradient of the cell */
    calculate_cell_dvds_ave (cell->ndom, cell); /* Defined at cell center */
    calculate_cell_dvds_max (cell->ndom, cell); /* Defined at cell vertex */

    if (rel_mode == REL_MODE_FULL)
    {
      cell->xgamma = calculate_gamma_factor (cell->v);
      model_velocity (cell->ndom, cell->xcen, v_cen);
      cell->xgamma_cen = calculate_gamma_factor (v_cen);
      cell->vol *= cell->xgamma_cen;
    }
    else
    {
      cell->xgamma = 1.0;
      cell->xgamma_cen = 1.0;
    }
  }

  /* Now communicate the cells between ranks */
  broadcast_wind_grid (n_start, n_stop, n_cells_rank);

  /* This wind *should* be fully initialised at this point, so we'll perform
   * some consistency checks to make sure everything is OK. Note that this
   * function is also used to assign the value of NPLASMA (the number of plasma
   * cells in the wind). This isn't really worth doing in parallel, so we have
   * kept it serial for now. */
  complete_wind_grid_creation ();
}

/**********************************************************/
/**
 * @brief Initialise the structures which characterise the wind, including
 *        wmain and plasmamain.
 *
 * @details
 *
 * This is the controlling routine for initialising the wind and related
 * structures, calculating the salient properties (e.g. the volume that is in
 * the wind) and each grid cell. In addition to initialise the properties, the
 * memory allocation is handled in here too.
 *
 * The general flow is first to create the wind grid, followed by the plasma
 * grid and then the macro grid if appropriate.
 *
 **********************************************************/

void
define_wind (void)
{
  int n;

  /* The first thing we need to do is define the wind grid, as this is the base
   * grid everything else shoots off from. The wind grid defines the position
   * of cells, as well as the velocity of the wind and which cells are in
   * or out of the wind */
  create_wind_grid ();

  /* With the wind grid set up, we can now deal with the grid where all the
   * atomic physics and radiative properties are stored. We have to do this
   * after creating the wind grid because that defines the number of cells with
   * wind volume. These wind cells are the ones which have plasma cell
   * counterparts. */
  create_plasma_grid ();

  /* If macro atoms are in use, we need to allocate and initialise the
   * macro grid as well. */
  if (geo.rt_mode == RT_MODE_MACRO)
  {
    create_macro_grid ();
  }

  /* We can now clean up any memory used to import a wind */
  for (n = 0; n < geo.ndomain; n++)
  {
    if (zdom[n].wind_type == IMPORT)
    {
      free_import (zdom[n].coord_type, n);
    }
  }

  /* For purely diagnostic purposes, a very (very) approximate value of
   * mdot_wind is calculated which can be compared to the input in the
   * parameter file */
  calculate_mdot_wind ();
}
