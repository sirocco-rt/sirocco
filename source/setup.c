/***********************************************************/
/** @file  setup.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  Various intialization routines, including a number
 * that read values from the parameter files
 *
 * ### Notes ###
 *
 * Because the input and setup of Python is relatively complex,
 * we have over time moved much of this out of main into
 * subroutines, and we have generally tried to collect
 * related portion of the initialization into different setup
 * routines, such as setup_disk or setup_files.
 *
 * This file contains a collection of these setup routines,
 * that either stand-alone or have not been moved into their
 * own files.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/**
 * @brief      initializes the geo structure to something that is semi-reasonable
 *
 * @return     Always returns 0
 *
 * @details
 * Initial values for all of the variables that are not part of the individual
 * wind descriptions that are actully read(!) into the program should be
 * created here.  The derived values are not needed.
 *
 * ### Notes ###
 *
 * In general, cgs units are the working units for Python.  Thus all
 * intialization should be converted to these units.
 *
 * @bug Currently init_geo is set up for CVs and Stars and not AGN.  We now
 * read in the system type as the first variable. The intent was to
 * allow one to use the system type to particularize how geo (aqnd
 * other variables) were intialized. But this has yet to be carreid
 * out.
 *
 *
 **********************************************************/

int
init_geo ()
{
  geo.ndomain = 0;              /*ndomain is a convenience variable so we do not always
                                   need to write geo.ndomain but it should nearly always
                                   be set to the same value as geo.ndomain */
  geo.run_type = 0;             /* Indicates this is a run from scratch, which includes
                                   the case where we already have a wind model but want
                                   to change some of the parameters.  init_goe should not
                                   be called at all if we are simply continuing a previous
                                   run */
  geo.hydro_domain_number = -1;
  geo.nplasma = 0;
  geo.nmacro = 0;

  if (geo.system_type == SYSTEM_TYPE_CV || geo.system_type == SYSTEM_TYPE_BH)
  {
    geo.binary = TRUE;
  }

/*  The domains have been created but have not been initialized at all */

  zdom[0].coord_type = CYLIND;
  zdom[0].ndim = 30;
  zdom[0].mdim = 30;
  zdom[0].log_linear = COORD_TYPE_LOG;  /* Set intervals to be logarithmic */

  zdom[1].coord_type = CYLIND;
  zdom[1].ndim = 30;
  zdom[1].mdim = 10;
  zdom[1].log_linear = COORD_TYPE_LOG;  /* Set intervals to be logarithmic */


  geo.disk_z0 = geo.disk_z1 = 0.0;
  geo.adiabatic = 1;            // Default is now set so that adiabatic cooling is included in the wind
  geo.auger_ionization = TRUE;  //Default is on.


  geo.run_type = 0;             // Not a restart of a previous run

  geo.star_ion_spectype = geo.star_spectype
    = geo.disk_ion_spectype = geo.disk_spectype = geo.bl_ion_spectype = geo.bl_spectype = SPECTYPE_BB;
  geo.agn_ion_spectype = SPECTYPE_POW;


  geo.rstar = 7e8;
  geo.rstar_sq = geo.rstar * geo.rstar;
  geo.mstar = 0.8 * MSOL;
  geo.m_sec = 0.4 * MSOL;
  geo.period = 3.2 * 3600;
  geo.tstar = 40000;

  geo.ioniz_mode = IONMODE_ML93;        /* default is on the spot and find the best t */
  geo.line_mode = LINE_MODE_ESC_PROB;   /* default is escape probabilites */
  geo.matom_transition_mode = MATOM_MC_JUMPS;   /* default for macro-atom transition mode is jumping not matrix */

  geo.star_radiation = TRUE;    /* 1 implies star will radiate */
  geo.disk_radiation = TRUE;    /* 1 implies disk will radiate */
  geo.bl_radiation = FALSE;     /*1 implies boundary layer will radiate */
  geo.wind_radiation = TRUE;    /* 1 implies wind will radiate */

  geo.disk_type = DISK_FLAT;    /*1 implies existence of a disk for purposes of absorption */
  geo.disk_rad_max = 2.4e10;
  geo.disk_mdot = 1.e-8 * MSOL / YR;
  geo.colour_correction = FCOL_OFF;

  geo.t_bl = 100000.;

  geo.pl_geometry = PL_GEOMETRY_SPHERE; // default to spherical geometry
  geo.lamp_post_height = 0.0;   // should only be used if geo.pl_geometry is PL_GEOMETRY_LAMP_POST
  geo.bubble_size = 0.0;        // should only be used if geo.pl_geometry is PL_GEOMETRY_BUBBLE_

  strcpy (geo.atomic_filename, "data/standard80.dat");
  strcpy (geo.fixed_con_file, "none");

  geo.wcycles = geo.pcycles = 1;
  geo.wcycle = geo.pcycle = 0;

  // Note that geo.model_list is initialized through get_spectype
  geo.model_count = 0;          //The number of models read in

  /* We should set the frame ASAP in the geo struct, so the grid initialisation
   * functions know what frame we're working with */
  if (rel_mode == REL_MODE_FULL)
  {
    geo.frame = CMF_FRAME;
  }
  else
  {
    geo.frame = OBS_FRAME;
  }

  return (0);
}

/// This is to assure that we read model lists in the same order everytime
char get_spectype_oldname[LINELENGTH] = "data/kurucz91.ls";
//int model_count = 0;


/**********************************************************/
/**
 * @brief      Generalized routine to get the spectrum type for a radiation source
 * and if necessary read a set of precalculated spectra
 *
 * @param [in] int  yesno  An integer used to decide whether to ask for a spectrum type
 * @param [in, out] char *  question  The query for a spectrum type for this component
 * @param [in, out] int *  spectype   The type of spectrum to assign
 * @return    the spectype, a number that says what type of spectrum (bb, power law, etc)
 * to generate for this source
 *
 * @details
 *
 * If yesno is non-zero, the user will be asked for the type of spectrum
 * for a radation source.  If yesno is 0, we presume that this radiation source (a star
 * a disk or the wind itself) is not to radiate in the calculation, and the spectrum
 * type is set to SPECTYPE_NONE.
 *
 * If one wants to geneate spectra from a seriels of models (such as synthetic spectra
 * calculated for stars), this routine calls routines to read the models.
 *
 * ### Notes ###
 *
 * This routine is awkwardly constructed. For reasons that are probably
 * historical the routine converts the internally stored values to different
 * values in order to ask what type of spectrum the user wants, and then
 * translates these back to the internal value.
 *
 *
 **********************************************************/

int
get_spectype (yesno, question, spectype)
     int yesno;
     char *question;
     int *spectype;
{
  char model_list[LINELENGTH];
  char one_choice[LINELENGTH];
  char choices[LINELENGTH];
  int get_models ();            // Note: Needed because get_models cannot be included in templates.h
  int i;
  int init_choices (), get_choices ();


  if (yesno)
  {
    init_choices ();            // Populate the spect array

    /* If this is is RUN_TYPE_PREVIOUS or RUN_TYPE_RESTART, we need to check
       whether a model spectrum was read in, and if so set the spectype to
       SPEC_TYPE_MODEL.  Note that this assume that all other spectrum types
       have negative values
     */

    if (*spectype >= 0 && (geo.run_type == RUN_TYPE_RESTART || geo.run_type == RUN_TYPE_PREVIOUS || geo.ioniz_or_extract == CYCLE_IONIZ))
    {
      *spectype = SPECTYPE_MODEL;
    }
    /* Locate the word that corresponds to the spectype that was entered
     */

    for (i = 0; i < zz_spec.n; i++)
    {
      if (*spectype == zz_spec.vals[i])
      {
        strcpy (one_choice, zz_spec.choices[i]);
        break;
      }
    }

    if (i == zz_spec.n)
    {
      Error ("get_spectype: Programming error.  Unknown spectype %d\n", *spectype);
      exit (0);
    }

    get_choices (question, choices, &zz_spec);
    *spectype = rdchoice (question, choices, one_choice);

    if (*spectype == SPECTYPE_MONO)
    {
      double x;
      x = 1000.;
      rddoub ("monochromatic.wavelength", &x);
      geo.mono_freq = VLIGHT * 1e8 / x;
    }


    if (*spectype == SPECTYPE_MODEL)
    {
      if (geo.run_type == RUN_TYPE_PREVIOUS)
      {                         // Continuing an old model
        strcpy (model_list, geo.model_list[geo.model_count]);
      }
      else
      {                         // Starting a new model
        strcpy (model_list, get_spectype_oldname);
      }

      rdstr ("Input_spectra.model_file", model_list);

      for (i = 0; i < geo.model_count; i++)     //See if we have already read in this model
      {
        if (strcmp (model_list, geo.model_list[i]) == 0)
        {
          *spectype = i;
          return (*spectype);
        }
      }
      get_models (model_list, 2, spectype);
      strcpy (geo.model_list[geo.model_count], model_list);     // Copy it to geo
      strcpy (get_spectype_oldname, model_list);        // Also copy it back to the old name

      geo.model_count++;
    }
  }


  else
  {
    *spectype = SPECTYPE_NONE;  // No radiation
  }

  return (*spectype);
}




/**********************************************************/
/**
 * @brief      simply initialises the set of
 * advanced modes stored in the modes structure to a
 * default value.
 *
 * @return   Allways returns 0
 *
 *
 * @details
 *
 * ### Notes ###
 * modes is a structure declared in sirocco.h
 *
 * Most of the modes are initialized to FALSE, which means that
 * activites that could be undertaken for this mode, will
 * not be.
 *
 * Advanced modes are turned off by default, unless the
 * -d flag is invoked at runtime.  If it is invoked,
 * one will be asked whether to turn-on an advanced mode.
 *
 * In most cases, the advanced modes cause extra diagnostic
 * information to be saved.
 *
 * see #111 and #120 for recent work
 *
 **********************************************************/

int
init_advanced_modes ()
{
  modes.iadvanced = FALSE;      // this is controlled by the -d flag, global mode control.
  modes.extra_diagnostics = FALSE;      //  when set, want to save some extra diagnostic info
  modes.save_cell_stats = FALSE;        // save photons statistics by cell
  modes.keep_ioncycle_windsaves = FALSE;        // save wind file each ionization cycle
  modes.keep_ioncycle_spectra = FALSE;  // to save spectrum file each ionization cycle
  modes.track_resonant_scatters = FALSE;        // track resonant scatters
  modes.save_extract_photons = FALSE;   // save details on extracted photons
  modes.adjust_grid = FALSE;    // the user wants to adjust the grid scale
  modes.diag_on_off = FALSE;    // extra diagnostics
  modes.use_debug = FALSE;
  modes.print_dvds_info = FALSE;        // print information on velocity gradients
  modes.quit_after_inputs = FALSE;      // check inputs and quit
  modes.quit_after_wind_defined = FALSE;        // define wind and quit
  modes.fixed_temp = FALSE;     // do not attempt to change temperature 
  modes.zeus_connect = FALSE;   // connect with zeus

  //note write_atomicdata  is defined in atomic.h, rather than the modes structure
  write_atomicdata = FALSE;     // print out summary of atomic data


  modes.keep_photoabs = TRUE;   // keep photoabsorption in final spectrum

  modes.jumps_for_detailed_spectra = FALSE;     //use old jumps mode for calculating macro atom
  //emissivites
  modes.use_upweighting_of_simple_macro_atoms = FALSE;  //use upweighting mode for handling 
  //bf interactions with simple macro atoms

  modes.store_matom_matrix = TRUE;      /* default is to store the macro-atom matrix */

  modes.run_xtest_diagnostics = FALSE;  /* allow special xtest_diagnostics in the various routines */
  modes.partial_cells = PC_ZERO_DEN;    /* Default is to omit partial cells in calculation */

  modes.no_macro_pops_for_ions = FALSE; /* use the ion densities from macro_pops where applicable */

  return (0);
}



/**********************************************************/
/**
 * @brief      get inputs that describe the detailed spectra that
 * one intends to extract.
 *
 * @return   Always returns 0
 *
 * @details
 *
 * ### Notes ###
 *
 * There are inputs here both for a normal extaction and in
 * advanced mode for selecting spectra from photons with specific numbers
 * of scatters, or from a certain reagion.  It is also possible
 * to extact spectra in the 'live or die" mode, in which one
 * simply tracks photons that emerge in a certain inclination
 * range.
 *
 **********************************************************/

int
init_observers ()
{
  int n;
  int ichoice;
  char answer[LINELENGTH];

  geo.nangles = 4;
  geo.angle[0] = 10;
  geo.angle[1] = 30.;
  geo.angle[2] = 60.;
  geo.angle[3] = 80.;
  for (n = 4; n < NSPEC; n++)
    geo.angle[n] = 45;
  for (n = 0; n < NSPEC; n++)
  {
    geo.phase[n] = 0.5;
    geo.scat_select[n] = MAXSCAT;
    geo.top_bot_select[n] = 0;
  }
  geo.swavemin = 850;
  geo.swavemax = 1850;

  rdpar_comment ("The minimum and maximum wavelengths in the final spectra and the number of wavelength bins");
  NWAVE_EXTRACT = NWAVE_IONIZ;
  rdint ("Spectrum.nwave", &NWAVE_EXTRACT);

  if (NWAVE_EXTRACT < 0)
  {
    Error ("The number of bins in the extracted spectra must be positive\n");
    Exit (0);
  }

  if (NWAVE_EXTRACT < NWAVE_MIN)
  {
    Error ("Minimum number of spectrum bins for detailed spectra is %d\n", NWAVE_MIN);
    NWAVE_EXTRACT = NWAVE_MIN;
  }

  /*
   * Determine the maximum number of wavelength bins we need to allocate. The
   * ionization and spectral cycles spectra use a different amount of wavelength
   * bins, but we allocate the same amount of memory for each cycle type
   */

  if (NWAVE_EXTRACT < NWAVE_IONIZ)
  {
    NWAVE_MAX = NWAVE_IONIZ;
  }
  else
  {
    NWAVE_MAX = NWAVE_EXTRACT;
  }

  /*
   * Now read in the read of the variables for the spectra
   */

  if (geo.mono_freq > 0)
  {
    geo.swavemax = 1.5 * VLIGHT * 1e8 / geo.mono_freq;
    geo.swavemin = 0.667 * VLIGHT * 1e8 / geo.mono_freq;
  }

  rddoub ("Spectrum.wavemin(Angstroms)", &geo.swavemin);
  rddoub ("Spectrum.wavemax(Angstroms)", &geo.swavemax);



  if (geo.swavemin > geo.swavemax)
  {
    geo.swavemax = geo.swavemin;
    geo.swavemin = geo.swavemax;
  }

  /* convert wavelengths to frequencies and store for use
     in computing macro atom and k-packet emissivities. */

  geo.sfmin = VLIGHT / (geo.swavemax * 1.e-8);
  geo.sfmax = VLIGHT / (geo.swavemin * 1.e-8);

  if (geo.mono_freq > 0)
  {
    if (geo.mono_freq < geo.sfmin || geo.mono_freq > geo.sfmax)
    {
      Error ("Monochromatic frequency %e out of range of frequency limits %e %e\n", geo.mono_freq, geo.sfmin, geo.sfmax);
      Exit (1);
    }
  }


  geo.matom_radiation = 0;      //initialise for ionization cycles - don't use pre-computed emissivities for macro-atom levels/ k-packets.

  rdpar_comment ("The observers and their location relative to the system");

  if (geo.system_type == SYSTEM_TYPE_STAR)
  {
    geo.nangles = 1;
    geo.angle[0] = 45;
  }

  rdint ("Spectrum.no_observers", &geo.nangles);

  if (geo.nangles < 1 || geo.nangles > NSPEC)
  {
    Error ("no_observers %d should not be > %d or <0\n", geo.nangles, NSPEC);
    Exit (0);
  }


  for (n = 0; n < geo.nangles; n++)
    rddoub ("Spectrum.angle(0=pole)", &geo.angle[n]);

  /* Phase 0 in this case corresponds to
   * an extraction direction which is in the xz plane
   */

  if (geo.system_type == SYSTEM_TYPE_CV || geo.system_type == SYSTEM_TYPE_BH)
  {

    for (n = 0; n < geo.nangles; n++)
      rddoub ("Spectrum.orbit_phase(0=inferior_conjunction)", &geo.phase[n]);
  }
  else
    Log ("No phase information needed as system type %i is not a binary\n", geo.system_type);


  strcpy (answer, "extract");
  geo.select_extract = rdchoice ("Spectrum.live_or_die(live.or.die,extract)", "0,1", answer);
  if (geo.select_extract == TRUE)
  {
    Log ("OK, detailed spectra will be created using the extract option\n");
  }
  else
    Log ("OK, detailed spectra will be created using the live or die option\n");

  /* In advanced mode, select spectra with certain numbers of scatterings and or other
   * characteristics
   */

  if (modes.iadvanced)
  {
    strcpy (answer, "no");
    ichoice = rdchoice ("@Spectrum.select_specific_no_of_scatters_in_spectra(yes,no)", ",1,0", answer);

    if (ichoice)
    {
      Log ("OK n>=%d->all; 0<=n<%d -> n scatters; n<0 -> >= |n| scatters\n", MAXSCAT, MAXSCAT);
      for (n = 0; n < geo.nangles; n++)
      {
        rdint ("@Spectrum.select_scatters", &geo.scat_select[n]);
      }
    }
    strcpy (answer, "no");
    ichoice = rdchoice ("@Spectrum.select_photons_by_position(yes,no)", "1,0", answer);
    if (ichoice && !geo.select_extract)
    {
      Error ("setup.c: Cannot select photons by position if in Live.or.die mode. Use Extract mode!\n");
      Exit (0);
    }
    else if (ichoice)
    {
      for (n = 0; n < geo.nangles; n++)
      {
        strcpy (answer, "all");
        geo.top_bot_select[n] = rdchoice ("@Spectrum.select_location(all,below_disk,above_disk,spherical_region)", "0,-1,1,2", answer);


        if (geo.top_bot_select[n] == 2)
        {
          Log ("Warning: Make sure that position will be in wind, or no joy will be obtained\n");
          rddoub ("@Spectrum.select_rho(cm)", &geo.rho_select[n]);
          rddoub ("@Spectrum.select_z(cm)", &geo.z_select[n]);
          rddoub ("@Spectrum.select_azimuth(deg)", &geo.az_select[n]);
          rddoub ("@Spectrum.select_r(cm)", &geo.r_select[n]);

        }
      }
    }
  }

  /* Select the units of the output spectra.  This is always needed.
   * There are 3 basics choices flambda, fnu, and the internal units
   * of the program.  The first two imply that output units are scaled
   * to a distance of 100 pc. The internal units are basically a luminosity
   * within a wavelength/frequency interval. */

  strcpy (answer, "flambda");
  geo.select_spectype = rdchoice ("Spectrum.type(flambda,fnu,basic)", "1,2,3", answer);

  if (geo.select_spectype == 1)
  {
    Log ("OK, generating flambda at 100pc\n");
    geo.select_spectype = SPECTYPE_FLAMBDA;
  }
  else if (geo.select_spectype == 2)
  {
    Log ("OK, generating fnu at 100 pc\n");
    geo.select_spectype = SPECTYPE_FNU;
  }
  else
  {
    Log ("OK, basic Monte Carlo spectrum\n");
    geo.select_spectype = SPECTYPE_RAW;
  }

  return (0);
}


/**********************************************************/
/**
 * @brief      gets information about the number of
 * 	cycles and how many photons there should be per cycle.
 *
 * @return     Generally returns 0
 *
 * @details
 * ??? DESCRIPTION ???
 *
 * ### Notes ###
 * The routine also allocates memory for the photon structure.
 * If the routine is unable to allocate this membory, the routine
 * will exit.
 *
 **********************************************************/

PhotPtr
init_photons ()
{
  PhotPtr p;

  /* Although Photons_per_cycle is really an integer,
     read in as a double so it is easier for input
     (in scientific notation) */


  double nphot = 1e5;
  rddoub ("Photons_per_cycle", &nphot); // NPHOT is photons/cycle
  if ((NPHOT = (int) nphot) <= 0)
  {
    Error ("%1.2e is invalid choice for NPHOT; NPHOT > 0 required.", (double) NPHOT);
    Exit (1);
  }


#ifdef MPI_ON
  Log ("Photons per cycle per MPI task will be %d\n", NPHOT / np_mpi_global);
  NPHOT /= np_mpi_global;
#endif

  rdint ("Ionization_cycles", &geo.wcycles);

  /* On restarts, the spectra that are read in have to be renormalized if
   * the number of spectral cycles has been increased before a restart, and
   * so we need to record this number. If this is not a restart, then
   * geo.pcycles_renorm will not be used.
   */

  geo.pcycles_renorm = geo.pcycles;

  rdint ("Spectrum_cycles", &geo.pcycles);


  /* If both geo.wcycles and geo.pcycles are 0, there is nothing to do, but for
   * except that for diagnostic reasons one migh want to inspect the initial
   * setup in the windsave file.  So now we wait to exit until this actually
   * written out.
   */

  if (geo.wcycles == 0 && geo.pcycles == 0)
  {
    Log ("Both ionization and spectral cycles are set to 0; the rest of the inputs will be gathered\n ");
    Log ("After that, the windsave file will be written to disk but then the program will exit\n");
  }

  /* Allocate the memory for the photon structure now that NPHOT is established */

  photmain = p = (PhotPtr) calloc (sizeof (p_dummy), NPHOT);
  /* If the number of photons per cycle is changed, NPHOT can be less, so we define NPHOT_MAX
   * to the maximum number of photons that one can create.  NPHOT is used extensively with
   * Python.  It is the NPHOT in a particular cycle, in a given thread.
   */

  NPHOT_MAX = NPHOT;


  if (p == NULL)
  {
    Error ("init_photons: There is a problem in allocating memory for the photon structure\n");
    Exit (0);
  }
  else
  {
    /* large photon numbers can cause problems / runs to crash. Report to use (see #209) */
    Log
      ("Allocated %10d bytes for each of %5d elements of photon structure totaling %10.1f Mb \n",
       sizeof (p_dummy), NPHOT, 1.e-6 * NPHOT * sizeof (p_dummy));
    if ((NPHOT * sizeof (p_dummy)) > 1e9)
      Error ("Over 1 GIGABYTE of photon structure allocated. Could cause serious problems.\n");
  }

  return (p);
}





/**********************************************************/
/**
 * @brief      Select the ionization and line transfer modes, along
 * with various other parameters about ionization
 *
 * @return   Always returns 0
 *
 * @details
 * The routine queries for a variety of parameters which describe
 * how ionization and line transfer are handled.  This includes
 * whether surfaces reflect
 *
 * ### Notes ###
 *
 * The routine is not particular well named, and it is not
 * immediately obvious that all of the parmenters that are queried
 * for here are related.
 *
 **********************************************************/

int
init_ionization ()
{
  int thermal_opt;
  char answer[LINELENGTH];



  strcpy (answer, "matrix_bb");
  geo.ioniz_mode =
    rdchoice ("Wind.ionization(on.the.spot,ML93,LTE_tr,LTE_te,fixed,matrix_bb,matrix_pow,matrix_est)", "0,3,1,4,2,8,9,10", answer);

  if (geo.ioniz_mode == IONMODE_FIXED)
  {
    char *a_warning;

    a_warning = "\
\n\
Warning: Fixed concentrations is intended primarily as a diagnostic mode  \n\
for studying the details of radiative transfer.   The calculated ne will   \n\
only agree with what is in the Fixed_concentrations file, if the elements/ions \n\
section of the input data agrees with the elements containted in the  \n\
fixed concentration file. \n\
\n\
";
    printf ("%s", a_warning);

    rdstr ("Wind.fixed_concentrations_file", &geo.fixed_con_file[0]);
  }


  /* get_line_transfer_mode reads in the Line_transfer question from the user,
     then alters the variables geo.line_mode, geo.scatter_mode, geo.rt_mode and geo.macro_simple.
     This is fairly involved and so is a separate routine  */

  get_line_transfer_mode ();

  /* we don't allow the user to use matrix_est ionization mode and macro-atom line transfer, see #875 */
  if ((geo.ioniz_mode == IONMODE_MATRIX_ESTIMATORS) && (geo.rt_mode == RT_MODE_MACRO))
  {
    Error ("matrix_est ionization mode cannot be used with macro-atom line transfer. Exiting.\n");
    Exit (EXIT_FAILURE);
  }


  strcpy (answer, "reflect");
  geo.absorb_reflect = rdchoice ("Surface.reflection.or.absorption(reflect,absorb,thermalized.rerad)", "1,0,2", answer);


  /* Setup options associated with non radiative process that can affect the thermal balance.  At present
   * these are adiabatic heating and an extra heating term explicitly implemented for FU Ori stars.  The
   * default is set to 0.  Adiabatic cooling only
   */

  strcpy (answer, "adiabatic");

  thermal_opt = rdchoice ("Wind_heating.extra_processes(none,adiabatic,nonthermal,both)", "1,0,2,3", answer);

  if (thermal_opt == 0)
  {
    geo.adiabatic = 1;
    geo.nonthermal = 0;
  }

  else if (thermal_opt == 1)
  {
    geo.adiabatic = 0;
    geo.nonthermal = 0;
  }
  else if (thermal_opt == 2)
  {
    geo.adiabatic = 0;
    geo.nonthermal = 1;
  }

  else if (thermal_opt == 3)
  {
    geo.adiabatic = 1;
    geo.nonthermal = 1;
  }
  else
  {
    Error ("Unknown thermal balance mode %d\n", thermal_opt);
    exit (0);
  }

  if (geo.nonthermal)
  {
    /* The shock heating is defined initally as a luminosity to be added to wind
     * but is immediately converted to a luminosity per unit volumne
     *
     * Since nearly all systems that we are dealing with have a star we initialize
     * the amount of extra heating as a fraction of the stellar luminosity
     *
     * See cooling.c shock_heating
     */

    Log
      ("Warning: Non-thermal heating has been selected.  This is a very special option put in place for modelling FU Ori stars, and should be used with extreme caution\n");

    geo.shock_factor = 0.001 * 4 * PI * pow (geo.rstar, 2) * STEFAN_BOLTZMANN * pow (geo.tstar, 4.);
    rddoub ("Wind_heating.extra_luminosity", &geo.shock_factor);
    geo.shock_factor /= (4 * PI * pow (geo.rstar, 3));
    Log ("The non_thermal emissivity at the base is %.2e\n", geo.shock_factor);

    if (geo.rt_mode == RT_MODE_MACRO)
    {
      geo.frac_extra_kpkts = 0.1;
      rddoub ("Wind_heating.kpacket_frac", &geo.frac_extra_kpkts);

    }
  }

  return (0);

}





/**********************************************************/
/**
 * @brief      Sets a global "push_through_distance" designed
 * to prevent photons from being trapped at boundaries
 *
 * @return     dfudge 	the push through distance
 *
 * @details
 *
 * dfudge is intended to be a small positive number that prevents
 * photons from being trapped at the edges of cells or wind cones.
 * It is often added to the distance to an "edge" to push the photon
 * thorough whatever boundary is being calculated.
 *
 * ### Notes ###
 *
 * setup_dfudge is used to modify the variable DFUDGE which
 * can be found in sirocco.h.
 *
 * Originally, DFUDGE was the only number used to push through
 * boundariies, but today DFUDGE is used sparingly, if at all,
 * as push through distances within wind cells are defined
 * differently for each cell.
 *
 * There are two competing factors in defining DFUDGE.  It
 * should be short enough so that the push through distane
 * goes only a small way into s cell.  It should be large
 * enough though that round-off errors do not prevent one
 * from actually getting into a cell.
 *
 **********************************************************/

double
setup_dfudge ()
{
  double dfudge;
  double delta;
  double rmin;
  int ndom;

  rmin = VERY_BIG;
  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    if (rmin > zdom[ndom].rmin)
    {
      rmin = zdom[ndom].rmin;
    }
  }

  delta = geo.rmax - rmin;

  if (delta < 1.e8)
  {
    dfudge = (geo.rmax - rmin) / 1000.0;
  }
  else if (delta < 1e15)
  {
    dfudge = 1e5;
  }
  else
  {
    dfudge = geo.rmax / 1.e10;
  }


  return (dfudge);
}



/**********************************************************/
/**
 * @brief Initialise the atomic data
 *
 * @details
 *
 **********************************************************/

void
setup_atomic_data (const char *atomic_filename)
{
  int rc;                       // Return code from running Setup_Py_Dir
  struct stat file_stat;        // Used to check the atomic data exists
  char answer[LINELENGTH];

  /* read a variable which controls whether to save a summary of atomic data
     this is defined in atomic.h, rather than the modes structure */

  if (modes.iadvanced)
  {
    strcpy (answer, "no");
    write_atomicdata = rdchoice ("@Diag.write_atomicdata(yes,no)", "1,0", answer);
    if (write_atomicdata)
      Log ("You have opted to save a summary of the atomic data\n");
  }

  /*
   * Check that geo.atomic_filename exists - i.e. that the directory is readable
   * and in the directory Python is being executed from. If it isn't - then
   * try to run Setup_Py_Dir. If both fail, then warn the user and exit Python
   */

  if (stat (atomic_filename, &file_stat))
  {
    Log ("Unable to open atomic masterfile %s\n", atomic_filename);
    Log ("Running Setup_Py_Dir to try and fix the situation\n");
    if (rank_global == 0)
    {
      rc = system ("Setup_Py_Dir");
      if (rc)
      {
        Error ("Unable to open %s and run Setup_Py_Dir\n", atomic_filename);
        Exit (1);
      }
    }
  }

  get_atomic_data ((char *) atomic_filename);

  /* throw a fatal error if there are macro-atom levels but rt_mode is non macro */
  if (nlevels_macro > 0 && geo.rt_mode != RT_MODE_MACRO)
  {
    Error ("Fatal error: you specified macro-atom data but standard line transfer. Not supported.\n");
    Log ("Try changing either to a simple data set or to macro-atom line transfer\n");
    Exit (1);
  }

  /* Prevent bf calculation of macro_estimators when no macro atoms are present.   */
  if (nlevels_macro == 0)
  {
    geo.macro_simple = TRUE;    // Make everything simple if no macro atoms -- 57h
  }
  /* initialise the choice of handling for macro pops. */
  if (geo.run_type == RUN_TYPE_PREVIOUS)
  {
    geo.macro_ioniz_mode = MACRO_IONIZ_MODE_ESTIMATORS; // Now that macro atom properties are available for restarts
  }
  else
  {
    geo.macro_ioniz_mode = MACRO_IONIZ_MODE_NO_ESTIMATORS;
  }


}
