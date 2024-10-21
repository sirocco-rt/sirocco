#ifdef MPI_ON
#include "mpi.h"
#endif



#define UV_low 7.4e14           /**< The lower frequency bound of the UV band as defined in IOS 21348
                                  */
#define UV_hi 3e16              /**< The lower frequency bound of the UV band as defined in IOS 21348
                                  */

extern int np_mpi_global;      /**< Global variable which holds the number of MPI processes
                                */

extern int rank_global;

extern int verbosity;          /**< verbosity level for printing out information. 0 low, 10 is high
                                 */

/* Integer representations for the logging level in Python's logging functions,
 * as defined in xlog.c */
#define SHOW_PARALLEL	     1
#define SHOW_ERROR	       2
#define SHOW_LOG  	       3
#define SHOW_DEBUG	       4
#define SHOW_LOG_SILENT    5
#define SHOW_ERROR_SILENT  5

#define TRUE  1
#define FALSE 0


#define USE_GRADIENTS        TRUE       /**< IF true use interpolated velocity gradients to calculate dv_ds */


#define REL_MODE_LINEAR 0       /**< Only make v/c corrections when doing frame transfers */
#define REL_MODE_FULL   1       /**< Make full corrections for special relativity including co-moving frame effects */
#define REL_MODE_SR_FREQ 2      /**< Make full corrects for special relativity frequency shifts, but ignore co-moving frame effects */

extern int rel_mode;                   /**< How doppler effects and co-moving frames are  */

extern int run_xtest;                  /**< Variable if TRUE causes a special test mode to be run */

extern int NDIM2;                      /**< The total number of wind cells in wmain
                                         */
extern int NPLASMA;                    /**< /The number of cells with non-zero volume or the size of plasma structure
                                        */

/* These are tunable parameters that control various aspects of sirocco
 * and the assorted programs.  In some cases they affect the "care" with
 * which a calculation is made, and this can in principle affect the
 * speed of the program.  One day if may be desirable to change some
 * of these parameters with time in the program.  At present they are
 * simply collected here
 *
 * */

extern double DFUDGE;
#define XFUDGE   1e-5           /**<  The scale factor used in setting up cell x cell dfudge
                                  */

#define VCHECK	1.e6            /**<  The maximum allowable error in calculation of the velocity in calculate_ds
                                  */


extern double SMAX_FRAC;               /**< In translate_in_wind, a limit is placed on the maximum distance a
                                   photon can travel in one step.  It is a fraction SMAX_FRAC of the
                                   distance of the photon from the origin.  This had been hardwired to
                                   0.1 for up to 57h.  Changing it to 0.5 speeds up the current version
                                   of the code by as much as a factor of 2 for small sized grids.  This
                                   had been introduced as part of the attempt to assure ourselves that
                                   line shapes were calculated as accurately as possible.  The underlhying
                                   rational for having a maximum disstance is associated with the fact that
                                   we use linear interpolation along the line of sight to establish velocities
                                   and most of our grid cells in 2.5d are actually hoop shaped, which means
                                   one can travel a long distance within a hoop if the direction of the photon
                                   is not more or less radial, but if moving along the hoop.
                                 */
extern double DENSITY_PHOT_MIN;   /**< This constant is a minimum density for the purpose of
                                    * calculating photoionization heating and recombination cooling.
                                    * It is important that heating and cooling
                                    * be calculated self-consistently.  Program speed is somewhat sensitive
                                    * to this parameter, at the 10% level if raised from 1e-3 to 1.
                                    * There is a
                                    * trade-off since lower minima may give better results, especially for
                                    * macro atoms.
                                    */

#define LDEN_MIN        1e-3    /**< The minimum density required for a line to be conidered for scattering
                                   or emission in calculate_ds and lum_lines */

#define DILUTION_FACTOR_MINIMUM 1e-10 /**< The next term globally defines a minimum
                                        * value for the dilution faction */

/* End of "care factor" definition */


#define EPSILON  			1.e-6   /**< A general purpose fairly small number */
#define NSTAT 				10      /**<  JM increased this to ten to allow for adiabatic
                                                  */
#define TAU_MAX				20.     /**<  Sets an upper limit in extract on when
                                                  *  a photon can be assumed to be completely absorbed
                                                 */
#define TAU_MIN                         1e-6    /**< minimum value for tau for p_escape_from_tau function-
                                                  *  below this we set to p_escape_ to 1
                                                 */

#define TMAX_FACTOR			1.5     /**< Factor by which t_e can exceed
                                                  *  t_r in order for absorbed to
                                                  *  match emitted flux */

#define TMAX    5e8             /**< This is the maximum temperature permitted. Tthis was
                                  * introduced following problems with adiabatically heated cells increasing
                                  * forever. The value was suggested by DP as a
                                  * sensible compton teperature for the PK05/P05 Zeus models. */
#define TMIN   100              /**<  A minimum temperature below which various emission processes
                                  * are set to 0 */

#define DELTA_V	 1.      /**< This is the accuracy in velocity space (cm/s) that we sample edges
                           * when producing freebound photons */

#define DANG_LIVE_OR_DIE   0.2  /**< If constructing photons from a live or die run of the code, the
                                  *  angle over which photons will be accepted must be defined */

extern double PHOT_RANGE;      /**< When a variable number of photons are called in different ionization
                                *  cycles this is the log of the difference between NPHOT_MAX
                                *  and the value in the first cycle
                                */
extern int NPHOT_MAX;          /**< The maximum number of photon bundles created per cycle */
extern int NPHOT;              /**< The number of photon bundles created, defined in setup.c */

extern int NWAVE_MAX;
extern int NWAVE_EXTRACT;     /**< The number of wavelength bins for spectra during the spectrum cycles
                                */
extern int NWAVE_NOW;         /**< Either NWAVE_IONIZ or NWAVE_EXTRACT depending on whether in
                                * ionizaiton of spectrum cycles
                                */


#define NWAVE_IONIZ 10000     /**< The number of wavelength bins for spectra during the ionization cycles
                                */
#define NWAVE_MIN   100       /**< The minimum number of wavelength bins in during spectral cycles
                               */
#define MAXSCAT    2000

/* Define the structures */
#include "math_struc.h"


/* Geometry is an actual structure which exists at the beginning of the program
   It carries the variables which define the geometry.  Reasonable values of each of
   these should be defined before it is altered with inputs from the terminal.
   The geometry structure is used to transfer all of the information about a wind
 */

/* Definitions of spectral types, which are all negative because when
 * one reads a spectrum from a list of models these are numbered beginning
 * with zero, see the discussion in get_models.c   080518
 */
#define SPECTYPE_BB      (-1)
#define SPECTYPE_UNIFORM (-2)
#define SPECTYPE_NONE	   (-3)
#define SPECTYPE_POW     (-4)
#define SPECTYPE_CL_TAB  (-5)   // This is to emulate cloudy
#define SPECTYPE_BREM    (-6)
#define SPECTYPE_MONO    (-7)
#define SPECTYPE_BB_FCOL (-8)
#define SPECTYPE_MODEL	 (-99)  // This is just used briefly, before a model number is assigned

/* definitions of types of colour correction */
#define FCOL_OFF  0
#define FCOL_DONE 1

/* Number of model_lists that one can have, should be the same as NCOMPS in models.h */
#define NCOMPS 	10
#ifdef LINELENGTH
#undef LINELENGTH
#endif
#define LINELENGTH 	400

/* This structure contains the information needed for each separate region of space, e.g the
 * wind and the disk
 */

enum coord_type_enum
{ SPHERICAL = 0,                //!< Spherical coordinates
  CYLIND = 1,                   //!< Standard cylindirical coordinates
  RTHETA = 2,                   //!< Polar coordinates
  CYLVAR = 3                    //!< Cylindrical coordinates but the z axis varies with rho
};


/* List of possible wind_types */

#define SV   			0
#define	STAR    		1
#define	HYDRO 			3
#define	CORONA 			4
#define KNIGGE			5
#define	HOMOLOGOUS 		6
#define	SHELL 			9
#define IMPORT                  10
#define	DISK_ATMOS 		11



/* Next define structures that pertain to possible region geometries
   as well as some used for vector operations
*/

/**
  * A structure defining a plane in 3 dimensions
  */

typedef struct plane
{
  double x[3];                  /**< A position included in the plane (usually the "center" */
  double lmn[3];                /**< A unit vector perpendicular to the plane (usually in the "positive" direction */
} plane_dummy, *PlanePtr;


/* These planes define the ends of a cylinder/pillbox which encapsulate the
   secondary with respect to the origin which also the location of the central source.
   The values are defined in the routin binary_basics*/

extern plane_dummy plane_m2_near, plane_m2_far;

/**
  * A structure defining a cone
  *
  * As currently defined cones are only defined azimuthally arounds the z axis and
  * not in 3d.  Thus they can  be specified in terms of an intercept along the
  * z axis and an opening anble or slope
  */
typedef struct cone
{
  double z;                     /**< The place where the cone intersects the z axis */
  double dzdr;                  /**< the slope */
}
cone_dummy, *ConePtr;

extern double velocity_electron[3];     // velocity of the electron when thermal effects are included

/* End of structures which are used to define boundaries to the emission regions */
/*******************DOMAIN structure***********************************************/

#define MAX_DOM			10
#define NDIM_MAX 3000           /**< maximum size of the grid in each dimension */
#define NDIM_MAX2D NDIM_MAX * NDIM_MAX  // Maximum dimensions for 2D importing

/**
  * The structure which defines a wind model
  *
  * In sirocco multiple domains can be created for different
  * portion of the wind.
  *
  */

#define COORD_TYPE_LOG 0
#define COORD_TYPE_LINEAR 1

typedef struct domain
{
  char name[LINELENGTH];
  int wind_type;
  int ndim, mdim, ndim2;        /**< ndim is the size in the x direction, while mdim is the
                                  size in z or theta direction */
  int nstart, nstop;            /**< the beginning and end (-1) location in wmain of this component */
  enum coord_type_enum coord_type;  /**< The type of coordinate system used for this domain */
  int log_linear;               /**< 0 -> the grid spacing will be logarithmic in x and z, 1-> linear */
  double xlog_scale, zlog_scale;        /**< Scale factors for setting up a logarithmic grid, the [1,1] cell
                                           will be located at xlog_scale,zlog_scale */

  /* The next few structures define the boundaries of an emission region */
  struct cone windcone[2];      /**< The cones that define the boundary of winds like SV or kwd */
  struct plane windplane[2];    /**< Planes which define the top and bottom of the wind */
//  double rho_min, rho_max;      /**< The values defining inner and outer cylinders that bound the wind */

  double *wind_x, *wind_z;
  double *wind_midx, *wind_midz;

  ConePtr cones_rtheta;         /**< A ptr to the cones that define boundaries of cells in the theta direction
                                   when rtheta coords  are being used */
/* Next two lines are for cyl_var coordinates.  They are used primarily for locating where a position is
 * in a grid with cyl_var coordinates. See cylvar_where in grid
 */

  double **wind_z_var, **wind_midz_var;

//  double wind_z_var[NDIM_MAX][NDIM_MAX];
//  double wind_midz_var[NDIM_MAX][NDIM_MAX];


  /* Generic parameters for the wind */
  double wind_mdot, stellar_wind_mdot;  /**< Mass loss rate in disk and stellar wind */
  double rmin, rmax;            /**< Spherical extent of the wind */
  double zmin, zmax;            /**<  Vertical extent of the wind, often the same as rmax */
  double wind_rhomin_at_disk, wind_rhomax_at_disk;      /**< Min/Max rho for wind in disk plane */
  double wind_thetamin, wind_thetamax;  /**< Angles defining inner and outer cones of wind, measured from disk plane */
  double mdot_norm;             /**< A normalization factor used in SV wind, and Knigge wind */

  double twind;                 /**< Initial temperature for a domain */

  /* Parameters defining Shlossman & Vitello Wind */
  double sv_lambda;             /**< power law exponent describing from  what portion of disk wind is radiated */
  double sv_rmin, sv_rmax, sv_thetamin, sv_thetamax, sv_gamma;  /**< parameters defining the goemetry of the wind */
  double sv_v_zero;             /**< velocity at base of wind */
  int sv_v_zero_mode;           /**< use fixed initial velocity or multiple of sound speed */
#define FIXED 0
#define SOUND_SPEED 1
  double sv_r_scale, sv_alpha;  /**< the scale length and power law exponent for the velocity law */
  double sv_v_infinity;         /**< the factor by which the velocity at infinity exceeds the excape velocity */


  /* Parameters defining Knigge Wind */
  double kn_dratio;             /**< parameter describing collimation of wind */
  double kn_lambda;             /**< power law exponent describing from  what portion of disk wind is radiated */
  double kn_r_scale, kn_alpha;  /**< the scale length and power law exponent for the velocity law */
  double kn_v_infinity;         /**< the factor by which the velocity at infinity exceeds the excape velocity */
  double kn_v_zero;             /**< NSH 19/04/11 - Added in as the multiple of the sound speed to use as the initial velocity */

  /* Parameters describing Castor and Larmors spherical wind */
  double cl_v_zero, cl_v_infinity, cl_beta;     /**< Power law exponent */
  double cl_rmin;

  /*Parameters defining a corona in a ring above a disk */
  double corona_rmin, corona_rmax;      /**< the minimum and maximu radius of the corona */
  double corona_zmax;                   /**< The maximum vertical extent of the corona */
  double corona_base_density, corona_scale_height;      /**< the density at the base of the corona and the scale height */
  double corona_vel_frac;       /**< the radial velocity of the corona in units of the keplerian velocity */

  double fill; /**< The filling factor for the wind or corona */
}
domain_dummy, *DomainPtr;       // One structure for each domain

extern DomainPtr zdom;          //This is the array pointer that contains the domains
extern int current_domain;      // This integer is used by swind only


/*******************GEOMETRY structure*********************************************/
#define SYSTEM_TYPE_STAR   0
#define SYSTEM_TYPE_CV     1
#define SYSTEM_TYPE_BH     4
#define SYSTEM_TYPE_AGN    2
#define	SYSTEM_TYPE_PREVIOUS   	   3

/* RUN_TYPE differs from SYSTEM_TYPE in that
  it has implications on how the program is run
  wherease SYSTEM_TYPE refers (mainly) to the type
  of system, with the exception
*/

#define RUN_TYPE_NEW       0
#define RUN_TYPE_RESTART   1
#define RUN_TYPE_PREVIOUS  3


/**
 * the geometry structure contains information that applies to all domains, including
 * the basic system geometry, descriptions of the radition sources, and truly
 * global information including how ionization calculations are caried out.
 *
 * Information that is domain specific should be placed directly in the domain
 * structure.
 */

struct geometry
{

#define OBS_FRAME 0
#define CMF_FRAME 1

  int frame;                    /**< Records frame parmeters like density and volumes are stroed */
  int system_type;              /**< See allowed types above. system_type should only be used for setup */
  int binary;                   /**< Indicates whether or not the system is a binary. TRUE or FALSE */

  int ndomain;                  /**< The number of domains in a model */
  int ndim2;                    /**< The total number of windcells in all domains */
  int nplasma, nmacro;          /**< The total number of cells in the plasma and macro structures in all domains */

  /* variables which store the domain numbers of the wind, disk atmosphere.
     Other components should be added here.  Right now we need a wind_domain
     number because the inputs for the disk and a putativel disk atmosphere are
     interrsed.  The first step will be to put this information into alocal variale
     in sirocco.c. We should not have to carry this forward */

  int wind_domain_number;
  /* Ultimately the next variable should not be needed but to avoid a bad
   * interaction with Nick's effort meld zeus calculations with Python
   * I have added a new variable.  It is likely that the domain number for
   * this will always be 0, but one could imagine cases where that might
   * not be the case  ksl -160927
   */

  int hydro_domain_number;      /**< Created for the special case of runs with Zeus
                                  */


  /* This section allows for restarting the program, and adds parameters used
   * in the calculation */

  int wcycle, pcycle;           /**< The number of completed ionization and spectrum cycles */
  int wcycles, pcycles, pcycles_renorm; /**< The number of ionization and spectrum cycles desired, pcycles_renorm
                                         * is only used on restarts.  See spectrum_restart_renormalize
                                         */
#define CYCLE_IONIZ    0
#define CYCLE_EXTRACT  1
  int ioniz_or_extract;         /**<  Set to CYCLE_IONIZ during ionization cycles, set to CYCLE_EXTRACT during calculation of
                                   detailed spectrum.
                                 */


  /* This section stores information which specifies the spectra to be extracted.  Some of the parameters
   * are used only in advanced modes.
   */

#define NSPEC   20
  int nangles;   /**< The number of angles to create spectra for */
  double angle[NSPEC], phase[NSPEC];  /**< The angle and associated binary phase (if relevant) for the extracted spectra */
  int scat_select[NSPEC], top_bot_select[NSPEC];  /**< Variables to constrain the spectra by number of scatters
                                                    * and whether the photons "originate" from above or relow the disk
                                                    * plane
                                                    */
  double rho_select[NSPEC], z_select[NSPEC], az_select[NSPEC], r_select[NSPEC];  /**< Variables which can be used to
                                                                                   * constrain the spectra by position
                                                                                   */
  double swavemin, swavemax, sfmin, sfmax;      // The minimum and maximum wavelengths/freqs for detailed spectra
  int select_extract, select_spectype;  /**< select_extract is TRUE if extract mode, FALSE if Live or Die
                                          * select_spectype indicates what type of spectrum, e.g FLAMBA, will
                                          * be created.
                                          */

/* Begin description of the actual geometry */

  double rmax, rmax_sq;         /**< The maximum distance (and distance**2) to which a photon should be followed */

/* Basic paremeters of the system, as opposed to elements of the wind or winds */

  double mstar, rstar, rstar_sq, tstar, gstar;  /**< Basic parameters for the star (often a WD) in the system */
  double tstar_init;            /**< The temperature of the star, before backscattering is taken into account */
  double lum_star_init, lum_star_back;  /**< The luminosity of the star as determined by tstar_init */

  double tmax;                  /**< the maximum temperature of any element of the model
                                   - used to help estimate things for an exponential representation of the spectrum in a cell */

#define DISK_MISSED 0
#define DISK_HIT_TOP 1
#define DISK_HIT_BOT 2
#define DISK_HIT_EDGE  3

#define DISK_NONE   0
#define DISK_FLAT   1
#define DISK_VERTICALLY_EXTENDED   2
#define DISK_WITH_HOLE      3

  int disk_type;

#define BACK_RAD_ABSORB_AND_DESTROY  0  /**< Disk simply absorbs the radiation and it is lost */
#define BACK_RAD_SCATTER            1   /**< Disk reradiates the radiation immediately via electron scattering */
#define BACK_RAD_ABSORB_AND_HEAT     2  /**< Correct disk temperature for illumination by photons
                                           which hit the dsik.  Disk radiation is absorbed and changes
                                           the temperature of the disk for future ionization cycles
                                         */

  int absorb_reflect;           /**< Controls what happens when a photon hits the disk or star
                                 */

#define DISK_TPROFILE_STANDARD          0       /**< This is a standard Shakura-Sunyaev disk. The profile depends on mstar and mdot_disk
                                                  */
#define DISK_TPROFILE_READIN            1       /**< Here the temperature profile for the disk is simply read in as a function of radius
                                                  */

  int disk_tprofile;            /**<  Variable used to specify a standard accretion disk (0) or
                                   one that has been read in and stored. */
  double disk_mdot;             /**<  mdot of  DISK */
  double disk_rad_min, disk_rad_max;
  double disk_z0, disk_z1;      /**<  For vertically extended disk, z=disk_z0*(r/diskrad)**disk_z1 *diskrad */
  double lum_disk_init, lum_disk_back;  /**<  The intrinsic luminosity of the disk, the back scattered luminosity */
  int run_type;                 /**<  Variable that describes whether this is a continuation of a previous run
                                   Added in order to separate the question of whether we are continuing an old run fro
                                   the type of wind model.  Bascially if run_type is 0, this is a run from scratch,
                                   if SYSTEM_TYPE_PREVIOUS it is an old run     */
  int star_radiation, disk_radiation;   /**<  1 means consider radiation from star, disk,  bl, and/or wind */
  int bl_radiation, wind_radiation, agn_radiation;

  int matom_radiation;          /**<  for use in macro atom computations of detailed spectra
                                   - 1 means use emissivities for BOTH macro atom levels and kpkts. 0 means don't
                                   (which is correct for the ionization cycles. */
  int ioniz_mode;               /**<  describes the type of ionization calculation which will
                                   be carried out.  The various ioniz_modes are defined by #defines IONMODE_MATRIX_BB
                                   etc.  See the documentation in this file for what these mean. */
#define MACRO_IONIZ_MODE_NO_ESTIMATORS  0
#define MACRO_IONIZ_MODE_ESTIMATORS     1
  int macro_ioniz_mode;         /**<  Controls the use of macro atom populations and
                                   ionization fractions. If it is set to MACRO_IONIZ_MODE_ESTIMATOR then macro atom populations
                                   computed from estimators are used. If set to MACRO_IONIZ_MODE_NO_ESTIMATORS then the macro atom
                                   populations are computed as for simple ions. By default it is set to
                                   MACRO_IONIZ_MODE_NO_ESTIMATORS initially and then set to  MACRO_IONIZ_MODE_ESTIMATORS the first time that
                                   Monte Carlo estimators are normalised. */
  int macro_simple;             /**<  As default this is set to FALSE, in which case a full
                                   Macro Atom calculation is performed. If it is set to TRUE it means
                                   that although Macro Atom data has been read in, all lines/continua are treated
                                   using the simplified two-level approximation. Such a calculation should reproduce
                                   similar results to the simple atom case.                 */

#define LINE_MODE_ABSORB      0
#define LINE_MODE_SCAT        1
#define LINE_MODE_SINGLE_SCAT 2
#define LINE_MODE_ESC_PROB    3



  int line_mode;                /**< 0, the atomosphere is a completely absorbing and no photons
                                   will be scattered.  In this mode, assuming the wind is a source
                                   of emission, the emissivity will be the Einstein A coefficient
                                   1, the atmosphere is completely scattering, there will be no
                                   interchange of energy between the photons and the electrons
                                   as a result of radiation transfer
                                   2, then a simple single scattering approximation is applied in which
                                   case the scattered flux is just  A21/(C21+A21).
                                   3, then radiation trapping aka escape probabilities are included as well.
                                   6 - 9,  If set to 6-9 initially, this switches on the macro atom case
                                   and then behaves like LINE_MODE_ESC_PROB.
                                 */
/* Note that the scatter_mode is actually a subsidiary variable of the line_mode.  Choosing a line_mode
 * results in the selection of a scatter_mode */

#define SCATTER_MODE_ISOTROPIC    0
#define SCATTER_MODE_THERMAL      2

  int scatter_mode;             /**< The way in which scattering for resonance lines is treated
                                  ** 0  isotropic
                                  ** 2  thermally broadened anisotropic
                                  */

#define RT_MODE_2LEVEL  1
#define RT_MODE_MACRO   2

  int rt_mode;                  /**<  radiative transfer mode.
                                  ** 2 for Macro Atom method,
                                  ** 1 for non-Macro Atom methods
                                  */

  /* define a global transition mode for macro-atoms */
  /* the global storage mode is set in the modes structure */
  int matom_transition_mode;

  /* Define the choices for calculating the FB, see, e.g. integ_fb */

#define FB_FULL         0       /**<  Calculate fb emissivity including energy associated with the threshold */
#define FB_REDUCED      1       /**<  Calculqate the fb emissivity without the threshold energy */
#define FB_RATE         2       /**<  Calulate the fb recombinarion rate  */


  /* Define the modes for free bound integrals */
#define OUTER_SHELL  1
#define INNER_SHELL  2

  /* The next set of variables defineds frequency bands used a boundaries for accumulating coarse (and fine) spectral information
     about the spectrum of photons in a cell. There are the coarse bands currently used for creating spectral models and there
     are a finer set of frequency intervals used for the fine spectra.  The spectra themselves can be found in the
     Plasma structure for each cell */

#define  NXBANDS 20             /**<  the maximum number of bands (frequency intervals that can be defined for
                                   storing coarse spectra for each plasma cell */

  int nxfreq;                   /**<  the number of frequency intervals actually used */
  double xfreq[NXBANDS + 1];    /**<  the frequency boundaries for the coarse spectra  */

#define NBINS_IN_CELL_SPEC   1000       /**< The number of bins in the cell spectra  */

  double cell_log_freq_min, cell_log_freq_max, cell_delta_lfreq;        /**< Parameters defining freqency intervals for cell spectra.
                                                                           These are defined as logarithmic frequency intervals */


  /* The next set pf variables assign a SPECTYPE (see above) for
     each possible source of radiation in a model.  The value assigned can be different for
     the ionization and detailed spectrum generation part of the code */

  int star_ion_spectype, star_spectype; /**<  The type of spectrum used to create the continuum
                                           for the star in the ionization and final spectrum calculation */
  int disk_ion_spectype, disk_spectype; /**<  Same as above but for the disk */
  int bl_ion_spectype, bl_spectype;     /**<  Same as above but for the boundary layer */
  int agn_ion_spectype, agn_spectype;   /**<  Same as above but for the AGN */

  int colour_correction;    /** colour correction mode */

  /* Searchlight mode is very experimental.  See notes in diag.c */
  double searchlight_x[3], searchlight_lmn[3];/**< location and direction of all photons in the spectral
                                                * cycles when searchlight mode (an advanced option) is
                                                * invoked
                                                */
  double mono_freq;                     /**< The frequency of all photons in the so-called monochromatic mode, which
                                          *  can be use for certain diagnositc experiments
                                          */

  char model_list[NCOMPS][LINELENGTH];  /**<  The file which contains the model names and the associated values for the model */
  int model_count;              /**< The number of distinct models that have been read in */

  double mdot_norm;             /**< A normalization factor used in SV wind, and Knigge wind */
  int adiabatic;                /**< 0-> Do not include adiabatic heating in calculating the cooling of the wind
                                   1-> Use adiabatic heating in calculating the cooling of the wind
                                 */
  int nonthermal;               /**<  0 --> No extra heating due to shocks
                                   1 --> Extra heating due to shocks (etc)  (Added for FU Ori)
                                 */

  double shock_factor;          /**<  A scaling factor used for including an extra heating term (for FU Ori stars
                                 */
  double frac_extra_kpkts;      /**<  in the case that we have extra heating and macro-atoms, the fraction of
                                   photons to reserve for those generated directly by k-packets */

  int auger_ionization;         /**< 0 -> Do not include innershell photoionization /Auger effects; 1-> include them */

/* Initial values for defining wind structure for a planar geometry.  These are currently only used by balance and this
   may not be the best approach generally and depending on where this ends up. Some consolidation is desirable */
  double pl_vol, pl_vmax;
  double pl_t_r, pl_t_e, pl_w;
  double pl_nh;

  /* Variables having to do with heating and cooling */

  double lum_tot, lum_star, lum_disk, lum_bl, lum_wind; /**<  The total luminosities of the disk, star, bl,
                                                          * & wind are actually not used in a fundamental
                                                          * way in the program */
  double lum_agn;               /**< The total luminosity of the AGN or point source at the center */
  double lum_ff, lum_rr, lum_lines;     /**<  The luminosity of the wind as a result of ff, fb, and line radiation */
  double cool_rr;               /**< the cooling rate due to radiative recombination - not the same as the luminosity */
  double cool_comp;             /**< The luminosity of the wind as a result of compton cooling */
  double cool_di;               /**< The direct ionization luminosity */
  double cool_dr;               /**< The luminosity of the wind due to dielectronic recombination */
  double cool_adiabatic;        /**< The cooling of the wind due to adiabatic expansion */
  double heat_adiabatic;        /**< The heating of the wind due to adiabatic heating - split out
                                  * from cool_adiabatic to get an accurate idea of whether it is important */
  double heat_shock;            /**< The amount of extra heating going into the wind due to shock heating. Added for FU Ori project */

  double f1, f2;                /**<  The freguency minimum and maximum for which the band limted luminosities have been calculated */
  double f_tot, f_star, f_disk, f_bl, f_agn, f_wind;    /**<  The integrated specific L between a freq min and max which are
                                                           used to establish the band limited luminosity  of photons of various types.
                                                           For detailed spectra cycles, this will the band limed luminosity between
                                                           the minimum and maximum wavelength of the detailed spectrum */

/* These variables are copies of the lum variables above, and are only calculated during ionization cycles
   This is a bugfix for JM130621, windsave bug */
  double lum_ff_ioniz, lum_rr_ioniz, lum_lines_ioniz;
  double cool_rr_ioniz;
  double cool_comp_ioniz;
  double cool_di_ioniz;         /**<  The direct ionization luminosity */
  double cool_dr_ioniz;
  double cool_adiabatic_ioniz;
  double lum_wind_ioniz, lum_star_ioniz, lum_disk_ioniz, lum_bl_ioniz, lum_tot_ioniz;

  double f_matom, f_kpkt;       /**< Used in computations of detailed spectra - the
                                   energy emitted in the band via k-packets and macro atoms respectively. */

  double n_ioniz, cool_tot_ioniz;

/* The next set of parameters relate to the secondary
 */

  double m_sec, q;              /**<  Mass of the secondary, mass ratio of system */
  double period;                /**<  Period of the systems in seconds */
  double a, l1, l2, phi;        /**<  Separation of primary and secondary, distance of l1 from primary,phi at l1 */
  double l1_from_m2, r2_far;    /**<  Distance to l1 from m2, distance to far side of secondary from primary */
  double r2_width;              /**<  Maximum half width of Roche surface of secondary in the plane of orbit */

  double t_bl;                  /**< temperature of the boundary layer */
  double weight;                /**< weight factor for photons/defined in define_phot */

/* The next set of parameters relate to the central source of an AGN
 */

  double brem_temp;             /**< The temperature of a bremsstrahlung source */
  double brem_alpha;            /**< The exponent of the nu term for a bremstrahlung source */

  double pl_low_cutoff;         /**<  accessible only in advanced mode- see #34. default to zero */

  double alpha_agn;             /**< The power law index of a BH at the center of an AGN.  Note that the luminosity
                                   of the agn is elsewhere in the structure */
  double const_agn;             /**< The constant for the Power law, there are lots of ways of defining the PL which is best? */
  double d_agn;                 /**<  the distance to the agn - only used in balance to calculate the ionization fraction */


  int pl_geometry;              /**<  geometry of X-ray point source */
#define PL_GEOMETRY_SPHERE 0
#define PL_GEOMETRY_LAMP_POST 1
#define PL_GEOMETRY_BUBBLE 2
  double lamp_post_height;      /**<  height of X-ray point source if lamp post */
  double bubble_size;           /**<  size of a bubble if any */

/* The next four variables added by nsh Apr 2012 to allow broken power law to match the cloudy table command */
  double agn_cltab_low;         /**< break at which the low frequency power law ends */
  double agn_cltab_hi;          /**< break at which the high frequency power law cuts in */
  double agn_cltab_low_alpha;   /**< photon index for the low frequency end */
  double agn_cltab_hi_alpha;    /**< photon index for the high frequency end   */

// The next set of parameters describe the input datafiles that are read
  char atomic_filename[132];    /**<  The masterfile for the atomic data */
  char fixed_con_file[132];     /**<  For fixed concentrations, the file specifying concentrations */

  //Added by SWM for tracking C-IV/H-A hotspots
  int nres_halpha;

  /* Variables used for reverberation mapping */

  double fraction_converged, reverb_fraction_converged;
  int reverb_filter_lines, *reverb_filter_line;
  enum reverb_disk_enum
  { REV_DISK_CORRELATED = 0, REV_DISK_UNCORRELATED = 1, REV_DISK_IGNORE = 3 } reverb_disk;
  enum reverb_enum
  { REV_NONE = 0, REV_PHOTON = 1, REV_WIND = 2, REV_MATOM = 3 } reverb;
  enum reverb_vis_enum
  { REV_VIS_NONE = 0, REV_VIS_VTK = 1, REV_VIS_DUMP = 2, REV_VIS_BOTH = 3 } reverb_vis;
  int reverb_wind_cycles;
  int reverb_path_bins, reverb_angle_bins;      //** Number of bins for path arrays, vtk output angular bins */
  int reverb_dump_cells;        /**< Number of cells to dump */
  double *reverb_dump_cell_x, *reverb_dump_cell_z;
  int *reverb_dump_cell;
  int reverb_lines, *reverb_line;       /**< Number of lines to track, and array of line 'nres' values */

  int spec_mod;                 /**< A flag saying we do have spectral models.
                                  * Used in Compton scattering and matrix_ion
                                  */
};


extern struct geometry geo;

/******************************END GEOMETRY STRUCTURE, BEGIN XDISK*********************************/
enum band_definition_enum
{
  T_STAR_BAND = 0,
  MIN_MAX_FREQ_BAND = 1,
  CV_BAND = 2,
  YSO_BAND = 3,
  USER_DEF_BAND = 4,
  CLOUDY_TEST_BAND = 5,
  WIDE_BAND = 6,
  AGN_BAND = 7,
  LOG_USER_DEF_BAND = 8,
  TDE_BB_BAND = 9
};

#define NRINGS	3001            /**< The actual number of rings completely defined
                                   is NRINGS-1 ... or from 0 to NRINGS-2.  This is
                                   because you need an outer radius...but the rest
                                   of this element is not filled in. */

/**
  xdisk is a structure that is used to store information about the disk in a system
*/
struct xdisk
{
  double r[NRINGS];             /**< The inner radius of an annulus */
  double t[NRINGS];             /**< The temperature at the middle of the annulus */
  double g[NRINGS];             /**< The gravity at the middle of the annulus */
  double v[NRINGS];             /**< The velocity at the middle of the annulus */
  double emit[NRINGS];          /**< The radiative energy of the photons emitted from the disk */
  double heat[NRINGS];          /**< The total energy flux of photons hitting each annulus */
  double ave_freq[NRINGS];      /**< The flux weighted average of frequency of photons hitting each annulus */
  double w[NRINGS];             /**< The weight relative to a BB of the photons that hit the disk */
  double t_hit[NRINGS];         /**< The effective T of photons hitting the disk, based on the average frequency */
  int nphot[NRINGS];        /**< The number of disk photons created in each annulus */
  int nhit[NRINGS];         /**< The number of photons which hit each annulus */
};
extern struct xdisk disk, qdisk;   /**<  disk defines zones in the disk which in a specified frequency band emit equal amounts
                                   of radiation. disk gets reinitialized whenever the frequency interval of interest
                                   is changed.  qdisk stores the amount of heating of the disk as a result of
                                   illumination by the star or wind. It's boundaries are fixed throughout a cycle */

#define NBLMODEL 5000

/** the blmodel structure is intended to store a non standard temperature
   and (optionally gravity) profile for the disk

   n_params should be 1 or 2, depending on whether t, or t and g
   are read in
   */

struct blmodel
{
  int n_params;
  int n_blpts;
  double r[NBLMODEL];
  double t[NBLMODEL];
  double g[NBLMODEL];
};
extern struct blmodel blmod;


/*************************WIND_PATHS for Reveberation Mapping *****************************************/
/**    Used for reverperation mapping, Wind paths is defined per cell and contains a binned array holding the spectrum of paths. Layers are
    For each frequency:
      For each path bin:
        What's the total fluxback of all these photons entering the cell?
*/
typedef struct wind_paths
{
  double *ad_path_flux;         /**< Array[by frequency, then path] of total flux of photons with the given v&p
                                  */
  double *ad_path_flux_disk;
  double *ad_path_flux_wind;
  double *ad_path_flux_cent;    /**<  As above, by source */
  int *ai_path_num;             /**< Array [by frequency, then path] of the number of photons in this bin
                                  */
  int *ai_path_num_disk;
  int *ai_path_num_wind;
  int *ai_path_num_cent;        /**<  As above, by source */
  double d_flux, d_path;        /**< Total flux, average path */
  int i_num;                    /**< Number of photons hitting this cell */
} wind_paths_dummy, *Wind_Paths_Ptr;

/******************************WIND STRUCTURE*******************************/
#define NIONIZ	5               /*The number of ions (normally H and He) for which one separately tracks ionization
                                   and recombinations */

/**
  * The structure that defines the wind for an individual domain after the
  * wind model has been intepreted

  * Note that most of the macro atom information is in a separate structure.  This was
  * done to make it easier to control the size of the entire structure   06jul-ksl
  */


typedef struct wind
{
  int ndom;                     /**< The domain associated with this element of the wind */
  int nwind;                    /**< A self-reference to this cell in the wind structure */
  int nwind_dom;                /**< The element number of the wind cell in its wind domain */
  int nplasma;                  /**< A cross refrence to the corresponding cell in the plasma structure */
  double x[3];                  /**< position of inner vertex of cell */
  double xcen[3];               /**< position of the "center" of a cell */
  double r, rcen;               /**< radial location of cell (Used for spherical, spherical polar
                                   coordinates. */
  double theta, thetacen;       /**< Angle of coordinate from z axis in degrees  */
  double dtheta, dr;            /**<  widths of bins, used in hydro import mode */
  struct cone wcone;            /**<  cone structure that defines the bottom edge of the cell in
                                   CYLVAR coordinates */
  double v[3];                  /**< velocity at inner vertex of cell in the observer frame.  For 2d coordinate systems this
                                   is defined in the xz plane */
  double v_grad[3][3];          /**< velocity gradient tensor  at the inner vertex of the cell in the co-moving frame */
  double div_v;                 /**< Divergence of v at center of cell in the co-moving frame */
  double dvds_ave;              /**<  Average value of dvds */
  double dvds_max;              /**< The maximum value of dvds */
  double vol;                   /**<  valid volume of this cell (that is the volume of the cell that is considered
                                   to be in the wind.  This differs from the volume in the Plasma structure
                                   where the volume is the volume that is actually filled with material.
                                   The vol that is stored here after the progam has initialized itself is the
                                   co-moving frame volume.
                                 */
  double xgamma, xgamma_cen;    /**<  1./sqrt(1-beta**2) at x at edge and center of cell */
  double dfudge;                /**<  A number which defines a push through distance for this cell, which replaces the
                                   global variable DFUDGE in many instances */
  enum inwind_enum
  { W_IN_DISK = -5, W_IN_STAR = -4, W_IGNORE = -2, W_NOT_INWIND = -1,
    W_ALL_INWIND = 0, W_PART_INWIND = 1, W_NOT_ASSIGNED = -999
  } inwind;                     /**< Basic information on the nature of a particular cell. */
  Wind_Paths_Ptr paths, *line_paths;    /**<  Path data struct for each cell */
}
wind_dummy, *WindPtr;

extern WindPtr wmain;

/*****************************PLASMA STRUCTURE**************************/
/** Plasma is a structure that contains information about the properties of the
 * plasma in regions of the geometry that are actually included in the wind
 * Note that a number of the arrays are dynamically allocated.
 */

typedef struct plasma
{
  int nwind;                    /**<  A cross reference to the corresponding cell in the  wind structure */
  int nplasma;                  /**<  A self reference to this  in the plasma structure */
  double ne;                    /**<  Electron density in the shell (CMF) */
  double rho;                   /**<  Density at the center of the cell. (CMF) For clumped models, this is rho of the clump */
  double vol;                   /**<  Volume of this cell in CMF frame (more specifically the volume  that is filled with material
                                   which can differs from the valid volume of the cell due to clumping.) */
  double xgamma;                /**<  1./sqrt(1-beta**2) at center of cell */
  double *density;              /**<  The number density of a specific ion in the CMF.  The order of the ions is
                                   the same as read in by the atomic data routines. */
  double *partition;            /**<  The partition function for each  ion.  */
  double *levden;               /* The number density (occupation number?) of a specific level */

  double kappa_ff_factor;       /**<  Multiplicative factor for calculating the FF heating for a photon. */


  double *recomb_simple;        /**<  "alpha_e - alpha" (in Leon's notation) for b-f processes in simple atoms. */
  double *recomb_simple_upweight;       /* multiplicative factor to account for ratio of total to "cooling" energy for b-f processes in simple atoms. */

/* Beginning of macro information */
  double kpkt_emiss;            /**< This is the luminosity produced due to the conversion k-packet -> r-packet in the cell
                                   in the frequency range that is required for the final spectral synthesis. (SS) */

  double kpkt_abs;              /**<  k-packet equivalent of matom_abs. (SS) */

  /* kbf_use and kbf_nuse are set by the routine kbf_need, and they provide indices into the photoinization processes
   * that are "significant" in a plasma cell, based on the density of a particular ion in a cell and the x-section
   * at the photoinization edge.  This process was introduced as a means to speed the program up by ignoring those
   * bf processes that would contribute negligibly to the bf opacity
   */

  int *kbf_use;                 /**<  List of the indices of the photoionization processes to be used for kappa_bf.  */
  int kbf_nuse;                 /**<  Total number of photoionization processes to be used for kappa_bf. (SS) */

/* End of macro information */


  double t_r, t_r_old;          /**< radiation temperature of cell */
  double t_e, t_e_old;          /**< electron temperature of cell */
  double dt_e, dt_e_old;        /**< How much t_e changed in the previous iteration */
  double heat_tot, heat_tot_old;        /**<  heating from all sources */
  double abs_tot;
  double heat_lines, heat_ff;
  double heat_comp;             /**<  The compton heating for the cell */
  double heat_ind_comp;         /**<  The induced compton heatingfor the cell */
  double heat_lines_macro, heat_photo_macro;    /**<  bb and bf heating due to macro atoms. Subset of heat_lines
                                                   and heat_photo. SS June 04. */
  double heat_photo, heat_z;    /**< photoionization heating total and of metals */
  double heat_auger;            /**<  photoionization heating due to inner shell ionizations */
  double heat_ch_ex;
  double abs_photo, abs_auger;  /**<  this is the energy absorbed from the photon due to these processes - different from
                                   the heating rate because of the binding energy */
  double w;                     /**< The dilution factor of the wind */

  int ntot;                     /**< Total number of photon passages */

  /*  counters of the number of photon passages by origin */

  int ntot_star;
  int ntot_bl;
  int ntot_disk;
  int ntot_wind;
  int ntot_agn;


  int nscat_es;                 /**<  The number of electrons scatters in the cell */
  int nscat_res;                /**<  The number of resonant line scatters in the cell */

  double mean_ds;               /**<  Mean photon path length in a cell. */
  int n_ds;                     /**<  Number of times a path lengyh was added; needed to compute mean_ds */
  int nrad;                     /**<  Total number of photons created within the cell */
  int nioniz;                   /**<  Total number of photon passages by photons capable of ionizing H */
  double *ioniz, *recomb;       /**<  Number of ionizations and recombinations for each ion.
                                   The sense is ionization from ion[n], and recombinations
                                   to each ion[n].  */
  double *inner_ioniz, *inner_recomb;
  int *scatters;                /**<  The number of scatters in this cell for each ion. */
  double *xscatters;            /**<  Diagnostic measure of energy scattered out of beam on extract. */
  double *heat_ion;             /**<  The amount of energy being transferred to the electron pool
                                   by this ion via photoionization. */
  double *heat_inner_ion;       /**<  The amount of energy being transferred to the electron pool
                                   by this ion via photoionization. */
  double *cool_rr_ion;          /**<  The amount of energy being released from the electron pool
                                   by this ion via recombination. */
  double *lum_rr_ion;           /**<  The recombination luminosity
                                   by this ion via recombination. */


#define MEAN_INTENSITY_BB_MODEL  1
#define MEAN_INTENSITY_ESTIMATOR_MODEL 2

  double *cool_dr_ion;
  double j, ave_freq;           /**<  Mean (angle-averaged) total intensity, intensity-averaged frequency */

  /* Information related to spectral bands used for modelling */
  double xj[NXBANDS], xave_freq[NXBANDS];       /**<  Frequency limited versions of j and ave_freq */
  double fmin[NXBANDS], fmax[NXBANDS];         /**<  Minimum (Maximum) frequency photon observed in a band -
                                                 * this is incremented during photon flight */
  double fmin_mod[NXBANDS], fmax_mod[NXBANDS];  /**<  Minimum (Maximum) frequency of the band-limited model
                                                  *  after allowing possibility that the observed limit,
                                                  *  is primarily due to photon statistics. See epectral_estimators.c */
  double xsd_freq[NXBANDS];     /**<  The standard deviation of the frequency in the band */
  int nxtot[NXBANDS];           /**<  The total number of photon passages in frequency bands */

  enum spec_mod_type_enum
  {
    SPEC_MOD_PL = 1,
    SPEC_MOD_EXP = 2,
    SPEC_MOD_FAIL = -1
  } spec_mod_type[NXBANDS];     /**<  A switch to say which type of representation we are using for this band in this cell.
                                   Negative means we have no useful representation, 0 means power law, 1 means exponential */

  double pl_alpha[NXBANDS];     /**< Computed spectral index for a power law spectrum representing this cell */
  double pl_log_w[NXBANDS];     /**< This is the log version of the power law weight. It is in an attempt to allow very large
                                   values of alpha to work with the PL spectral model to avoide NAN problems.
                                   The pl_w version can be deleted once testing is complete */


  double exp_temp[NXBANDS];     /**<  The effective temperature of an exponential representation of the radiation field in a cell */
  double exp_w[NXBANDS];        /**<  The prefactor of an exponential representation of the radiation field in a cell */

  double cell_spec_flux[NBINS_IN_CELL_SPEC];    /**< The array where the cell spectra are accumulated. */

#define NFLUX_ANGLES 36 /**< The number of bins into which the directional flux is calculated */


  /*Binned fluxes */
  double F_UV_ang_theta[NFLUX_ANGLES];
  double F_UV_ang_phi[NFLUX_ANGLES];
  double F_UV_ang_r[NFLUX_ANGLES];


  /*A version of the binned flux that is averaged over cycles */
  double F_UV_ang_theta_persist[NFLUX_ANGLES];
  double F_UV_ang_phi_persist[NFLUX_ANGLES];
  double F_UV_ang_r_persist[NFLUX_ANGLES];

  /* The term direct here means from photons which have not been scattered. These are photons which have been
     created by the central object, or the disk, or in the simple case the wind, but which have not undergone
     any kind of interaction which would change their direction
   */
  double j_direct, j_scatt;     /**<  Mean intensity due to direct photons and scattered photons.
                                 Direct photons include photons created in the wind in simple mode. */
  double ip_direct, ip_scatt;   /**<  Ionization parameter due  to direct photons and scattered photons. See ip */
  double max_freq;              /**<   The maximum frequency photon seen in this cell */
  double cool_tot;              /**< The total cooling in a cell */
  /* The total luminosity of all processes in the cell, basically the emissivity of the cell times it volume. Not the same
     as what escapes the cell, since photons can interact within the cell and lose weight or even be destroyed */
  double lum_lines, lum_ff, cool_adiabatic;
  double lum_rr, lum_rr_metals; /**<  the radiative recombination luminosity - not the same as the cooling rate */
  double cool_comp;             /**<  The compton luminosity of the cell */
  double cool_di;               /**<  The direct ionization luminosity */
  double cool_dr;               /**<  The dielectronic recombination luminosity of the cell */
  double cool_rr, cool_rr_metals;       /**< fb luminosity & fb of metals metals */
  double lum_tot, lum_tot_old;  /**<  The specific radiative luminosity in frequencies defined by freqmin
                                   and freqmax.  This will depend on the last call to total_emission */
  double cool_tot_ioniz;
  double lum_lines_ioniz, lum_ff_ioniz, cool_adiabatic_ioniz;
  double lum_rr_ioniz;
  double cool_comp_ioniz;       /**<  The compton luminosity of the cell */
  double cool_di_ioniz;         /**<  The direct ionization luminosity */
  double cool_dr_ioniz;         /**<  The dielectronic recombination luminosity of the cell */
  double cool_rr_ioniz, cool_rr_metals_ioniz;   /**< fb luminosity & fb of metals metals */
  double lum_tot_ioniz;         /**<  The specfic radiative luminosity in frequencies defined by freqmin
                                   and freqmax.  This will depend on the last call to total_emission */
  double heat_shock;            /**<  An extra heating term added to allow for shock heating of the plasma (Implementef for FU Ori Project */

  double bf_simple_ionpool_in, bf_simple_ionpool_out; /**<Varibles to track net flow of energy
                                                        into and from ionization pool
                                                        in BF_SIMPLE_EMISSIVITY_APPROACH
                                                        */
#define N_PHOT_PROC 500
  int n_bf_in[N_PHOT_PROC], n_bf_out[N_PHOT_PROC];
                                                 /**<Counters to track bf excitations and de-exitations.
                                                   */

  double comp_nujnu;            /**<  The integral of alpha(nu)nuj(nu) used to
                                   compute compton cooling-  only needs computing once per cycle
                                 */

#define N_DMO_DT_DIRECTIONS 3
#define NFORCE_DIRECTIONS 4
  /* directional fluxes (in observer frame) in 3 wavebands. - last element contains the  magnitude of flux) */
  double F_vis[NFORCE_DIRECTIONS];
  double F_UV[NFORCE_DIRECTIONS];
  double F_Xray[NFORCE_DIRECTIONS];

  double F_vis_persistent[NFORCE_DIRECTIONS];
  double F_UV_persistent[NFORCE_DIRECTIONS];
  double F_Xray_persistent[NFORCE_DIRECTIONS];

  double dmo_dt[N_DMO_DT_DIRECTIONS];             /**< Radiative force of wind */
  double rad_force_es[NFORCE_DIRECTIONS];       /**< Radiative force of wind - 4th element is sum of magnitudes */
  double rad_force_ff[NFORCE_DIRECTIONS];       /**< Radiative force of wind - 4th element is sum of magnitudes */
  double rad_force_bf[NFORCE_DIRECTIONS];       /**< Radiative force of wind - 4th element is sum of magnitudes */

  double rad_force_es_persist[NFORCE_DIRECTIONS];       /**< Radiative force of wind - 4th element is sum of magnitudes */
  double rad_force_ff_persist[NFORCE_DIRECTIONS];       /**< Radiative force of wind - 4th element is sum of magnitudes */
  double rad_force_bf_persist[NFORCE_DIRECTIONS];       /**< Radiative force of wind - 4th element is sum of magnitudes */

  double gain;                  /**<  The gain being used in iterations of the structure */
  double converge_t_r, converge_t_e, converge_hc;
  /**< Three measures of whether the program believes the grid is converged.  The first two
     are the fractional changes in t_r, t_e between this and the last cycle.
     The third number is the fraction between heating and cooling divided by the sum of the 2
   */

  int trcheck, techeck, hccheck;
  /**< The individual convergence checks used to calculate converge_whole.
     Each of these values is 0 if the fractional change or in the case of the last
     check error is less than a value, currently set to 0.05.  This number is now
     also used to say if the cell is over temperature - it is set to 2 in this case   */

  int converge_whole, converging;
  /**< converge_whole is the sum of the individual convergence checks.  It is 0 if all
     of the convergence checks indicated convergence. converging is an indicator of whether
     the program thought the cell is on the way to convergence 0 implies converging */

#define CELL_CONVERGING 0       /*  converging - temperature is oscillating and decreasing */
#define CELL_NOT_CONVERGING 1   /*  not converging (temperature is shooting off in one direction) */
#define CONVERGENCE_CHECK_PASS 0        /* Cell has passed a convergence check */
#define CONVERGENCE_CHECK_FAIL 1        /* Cell has failed a convergence check */
#define CONVERGENCE_CHECK_OVER_TEMP 2   /* Cell has electron temperature is more than TMAX */

  double ip;                    /**<  Ionization parameter calculated as number of photons over the lyman limit entering a cell, 
                                  divided by the number density of hydrogen for the cell.  This is the definnition used in Cloudy */
  double xi;                    /**<  Ionization parameter as defined by Tarter, Tucker, and Salpeter  1969 (ApJ 156, 943).  
                                  It is the ionizing flux over the number of hydrogen atoms */
} plasma_dummy, *PlasmaPtr;

extern PlasmaPtr plasmamain;

/*******************************PHOTON_STORE*********************************************/
#define NSTORE 10

/** A storage area for photons.  The idea is that it is sometimes time-consuming to create the
cumulative distribution function for a process, but trivial to create more than one photon
of a particular type once one has the cdf,  This appears to be case for f fb photons.  But
the same procedure could be used for line photons */

typedef struct photon_store
{
  int n;                        /**<  This is the photon number that was last used */
  double t, f1, f2, freq[NSTORE];

} photon_store_dummy, *PhotStorePtr;

extern PhotStorePtr photstoremain;

/** A second photon store: this is very similar to photon_store above but for use
   in generating macro atom bf photons from cfds*/
typedef struct matom_photon_store
{
  int n;                        /**<  This is the photon number that was last used */
  double t, nconf, freq[NSTORE];

} matom_photon_store_dummy, *MatomPhotStorePtr;

extern MatomPhotStorePtr matomphotstoremain;
//OLD #define MATOM_BF_PDF 1000       /**< number of points to use in a macro atom bf PDF
//OLD                                   */


/*******************************MACRO STRUCTURE*****************************/
/**
  The stucture used for storing infomration for macro atoms

  The various arrays created here are organized sequentially by macro level
   and so the number of elements in each is the number of macro levels.
*/
typedef struct macro
{
  double *jbar; /**<  This will store the Sobolev mean intensity in transitions which is needed
     for Macro Atom jumping probabilities. The indexing is by configuration (the
     NLTE_LEVELS) and then by the upward bound-bound jumps from that level
     (the NBBJUMPS) (SS) */

  double *jbar_old;

  double *gamma; /**< This is similar to the jbar but for bound-free transitions. It records the
     appropriate photoionisation rate co-efficient. (SS) */

  double *gamma_old;

  double *gamma_e; /**< This is Leon's gamma_e: very similar to gamma but energy weighted. Needed
     for division of photoionisation energy into excitation and k-packets. (SS) */

  double *gamma_e_old;

  double *alpha_st; /**< Same as gamma but for stimulated recombination rather than photoionisation. (SS) */

  double *alpha_st_old;

  double *alpha_st_e; /**< Same as gamma_e but for stimulated recombination rather than photoionisation. (SS) */

  double *alpha_st_e_old;

  double *recomb_sp; /**< Spontaneous recombination. (SS) */

  double *recomb_sp_e; /**< "e" version of the spontaneous recombination coefficient. (SS) */

  double *matom_emiss; /**< This is the luminosity due to the de-activation of macro atoms in the cell
     in the frequency range that is required for the final spectral synthesis. (SS) */

  double *matom_abs; /**< This is the energy absorbed by the macro atom levels - recorded during the ionization
     cycles and used to get matom_emiss (SS) */

  /* This portion of the macro structure  is not written out by windsave */
  int kpkt_rates_known;

  double *cooling_bf;
  double *cooling_bf_col;
  double *cooling_bb;

  /* set of cooling rate stores, which are calculated for each macro atom each cycle,
     and used to select destruction rates for kpkts */
  double cooling_normalisation;
  double cooling_bbtot, cooling_bftot, cooling_bf_coltot;
  double cooling_bb_simple_tot;
  double cooling_ff, cooling_ff_lofreq;
  double cooling_adiabatic;     // this is just cool_adiabatic / vol / ne

#define MATOM_MC_JUMPS 0
#define MATOM_MATRIX   1
  int matom_transition_mode;    /**<  what mode to use for the macro-atom transition probabilities */
  int store_matom_matrix;
  int matrix_rates_known;
  double **matom_matrix;        /**<  array to store transitions probabilities */
} macro_dummy, *MacroPtr;

extern MacroPtr macromain;

//extern int xxxpdfwind;          // When 1, line luminosity calculates pdf

extern int size_Jbar_est, size_gamma_est, size_alpha_est;


//These constants are used in the various routines which compute ionization state
#define SAHA 4.82907e15         /**<  2* (2.*PI*MELEC*k)**1.5 / h**3  (Calculated in constants) */
#define MAXITERATIONS	200     /**< /The number of loops to do to try to converge in ne */
#define FRACTIONAL_ERROR 0.03   /**< The change in n_e which causes a break out of the loop for ne */
#define THETAMAX	 1e4    /**< Used in initial calculation of n_e */
#define MIN_TEMP	100.    /**<  ??? this is another minimum temperature, which is
                                  * used in saha.c and variable_temperature.c )
                                  */

// these definitions are for various ionization modes
#define IONMODE_ML93_FIXTE 0    /**< Lucy Mazzali using existing t_e (no HC balance) */
#define IONMODE_LTE_TR 1        /**< LTE using t_r */
#define IONMODE_LTE_TE 4        /**<  LTE using t_e */
#define IONMODE_FIXED 2         /**<  Hardwired concentrations */
#define IONMODE_ML93 3          /**<  Lucy Mazzali */
#define IONMODE_MATRIX_BB 8     /**<  matrix solver BB model */
#define IONMODE_MATRIX_SPECTRALMODEL 9  /**< matrix solver spectral model based on power laws */
#define IONMODE_MATRIX_ESTIMATORS 10    /**<  matrix solver spectral model based on power laws */

// and the corresponding modes in nebular_concentrations
#define NEBULARMODE_TR 0        /**< LTE using t_r */
#define NEBULARMODE_TE 1        /**< LTE using t_e */
#define NEBULARMODE_ML93 2      /**< ML93 using correction */
#define NEBULARMODE_NLTE_SIM 3  /**< Non_LTE with SS modification (Probably could be removed) */
#define NEBULARMODE_LTE_GROUND 4        /**< A test mode which forces all levels to the GS (Probably could be removed) */
#define NEBULARMODE_PAIRWISE_ML93 6     /**< pairwise ML93 (diffuse BB) */
#define NEBULARMODE_PAIRWISE_SPECTRALMODEL 7    /**< pairwise spectral models (power law or expoentials) */
#define NEBULARMODE_MATRIX_BB 8 /**<  matrix solver BB model */
#define NEBULARMODE_MATRIX_SPECTRALMODEL 9      /**< matrix solver spectral model */
#define NEBULARMODE_MATRIX_ESTIMATORS 10        /**<  matrix solver spectral model */

#define NEBULARMODE_MATRIX_MULTISHOT   11    /**<  matrix solver spectral model based on power laws which
                                          * updates T_e multiple times before arriving at a final
                                          * solution 
                                          */

// modes for the wind_luminosity routine
#define MODE_OBSERVER_FRAME_TIME 0
#define MODE_CMF_TIME 1

/***************************************PHOTON STRUCTURE*********************************/
/**
  * The structure that contains information about indidual photons as they pass through
  * the wind
  */
typedef struct photon
{
  double x[3];                  /**<  The position of packet */
  double lmn[3];                /**<  Direction cosines of the packet */
  double freq, freq_orig;       /**<  current, original frequency (redshifted) of this packet */
  double w, w_orig;             /**<  current and original weight of this packet */
  double tau;                   /**<  optical depth of the photon since its creation or last interaction */

#define N_ISTAT 13              /**<  number of entries in the istat_enum */
  enum istat_enum
  {
    P_INWIND = 0,               /**< in wind, */
    P_SCAT = 1,                 /**< in process of scattering, */
    P_ESCAPE = 2,               /**< Escaped to reach the universe, */
    P_HIT_STAR = 3,             /**< absorbed by photosphere of star, */
    P_TOO_MANY_SCATTERS = 4,    /**< in wind after MAXSCAT scatters */
    P_ERROR = 5,                /**< Trying to scatter a photon in a location where it should not scatter */
    P_ABSORB = 6,               /**< Photoabsorbed within wind */
    P_HIT_DISK = 7,             /**< Banged into disk */
    P_SEC = 8,                  /**< Photon hit secondary */
    P_ADIABATIC = 9,            /**< records that a photon created a kpkt which was destroyed by adiabatic cooling */
    P_ERROR_MATOM = 10,         /**< Some kind of error in processing of a photon which excited a macroattom */
    P_LOFREQ_FF = 11,           /**< records a photon that had too low a frequency */
    P_REPOSITION_ERROR = 12     /**< A photon passed through the disk due to dfudge pushing it through incorrectly */
  } istat;                      /**< status of photon. */

  enum frame
  {
    F_LOCAL = 0,     /**< The photon is in the local frame */
    F_OBSERVER = 1   /**< The photon is in the observer frame */
  } frame;

  int nscat;         /**< Number of scatters for this photon */
  int nrscat;        /**<  number of resonance scatterings */
  int nmacro;        /**<  number of macro atom interactions */

  int nres;          /**< For line scattering, indicates the actual transition;
                                   for continuum scattering, meaning
                                   depends on matom vs non-matom. See headers of emission.c
                                   or matom.c for details. */
  int line_res;      /**<  The line which a photon belongs to. A photon can tagged as a line, then continuum
                                   scatter. line_res will still be tagged as the same line until it scatters off another line. */
  int nnscat;        /**<  Used for the thermal trapping model of
                       * anisotropic scattering to carry the number of
                       * scattering to "extract" when needed for wind
                       * generated photons SS05. */
  int grid;          /**< grid position of the photon in the wind, if
                       * the photon is in the wind.  If the photon is not
                       * in the wind, then -1 implies inside the wind cone and
                       * -2 implies outside the wind */

  enum origin_enum
  { PTYPE_STAR = 0,
    PTYPE_BL = 1,
    PTYPE_DISK = 2,
    PTYPE_WIND = 3,
    PTYPE_AGN = 4,
    PTYPE_STAR_MATOM = 10,
    PTYPE_BL_MATOM = 11,
    PTYPE_DISK_MATOM = 12,
    PTYPE_WIND_MATOM = 13,
    PTYPE_AGN_MATOM = 14,
    PTYPE_DUMMY = -1
  } origin, origin_orig;        /* Where this photon originated.  If the photon has
                                   scattered its "origin" may be changed to "wind". */
  /* note that we add 10 to origin when processed by a macro-atom
     which means we need these values in the enum list.  In making spectra in spectrum_create
     10 is subtracted from the types.  If ever this logic is changed one must the be careful
     that it is fixed in create_spectra as well.

     Comment - ksl - 180712 - The logic for all of this is obscure to me, since we keep track of the
     photons origin separately.  At some point one might want to revisit the necessity for this
   */
  int np;                       /* The photon number, which eases tracking a photon for diagnostic
                                   purposes */
  double path;                  /* The total path length of a photon (used for reverberation calcuations) */
  double ds;                    /* the distance a photon has moved since its creation or last interaction */
}
p_dummy, *PhotPtr;

#define NRES_ES (-1)
#define NRES_FF (-2)
#define NRES_NOT_SET (-3)
#define NRES_BF NLINES

extern PhotPtr photmain;               /**< A pointer to all of the photons that have been created in a subcycle. Added to ease
                                        breaking the main routine of sirocco into separate rooutines for inputs and
                                        running the program */

/**************************************SPECTRUM STRUCTURE ***********************/
    /* The next section defines the spectrum arrays.  The spectrum structure contains
       the selection criteria for each spectrum as well as the array in which to store the
       spectrum.  The first MSPEC spectra are reserved for the total generated spectrum,
       total emitted spectrum, the total spectrum of photons which are scattered and
       the total spectrum of photons which are absorbed.  The remainder of the spectra pertain
       to the spectrum as seen from various directions. Note that the space for the spectra
       are allocated using a calloc statement in spectrum_init.   
       */



#define SPECTYPE_RAW        0   /**< As written this produces L_nu and so to get a Luminosity
                                  * one needs to integrate
                                  */
#define SPECTYPE_FLAMBDA    1
#define SPECTYPE_FNU        2
#define D_SOURCE 100.0          /**<  distance to the source in parsecs for genearating spectra
                                  */

extern int nspectra;           /**< After create_spectrum, the number of elements allocated for s, or
                                 * alternatively the number of spectra one has to work with.  Note that
                                 * general s[0],s[1] and s[2] are the escaping, scattered and absorbed photons,
                                 * while elements higher than this will contain spectra
                                 as seen by different observers */


#define MSPEC               8   /**< The number of standard spectra - i.e. not user defined angles */
#define SPEC_CREATED        0   /**< The spectrum of from external sources with  weights before transmission through the wind */
#define SPEC_CWIND          1   /**< The spectrum created in the wind with rheir original weights */
#define SPEC_EMITTED        2   /**< The emitted spectrum - i.e. photons with their weights changed by transmission through the wind */
#define SPEC_CENSRC         3   /**< The emitted spectrum from photons emitted from the central source (if there is one) */
#define SPEC_DISK           4   /**< The emitted spectrum from photons emitted from the disk (if there is one) */
#define SPEC_WIND           5   /**< The emitted spectrum from photons emitted from the wind itself */
#define SPEC_HITSURF        6   /**< The spectrum for photons which hit the a surface and were absorbed
                                  * - should be zero for when reflection
                                  * is turned on */
#define SPEC_SCATTERED      7   /* The spectrum of photons which were scattered at least
                                 * once in the wind - the weight used is the final
                                 * weight after transmission through the wind */

extern int nscat[MAXSCAT + 1], nres[MAXSCAT + 1], nstat[NSTAT];

typedef struct spectrum
{
  char name[40];
  double freqmin, freqmax, dfreq;
  double lfreqmin, lfreqmax, ldfreq;    /**<  NSH 1302 - values for logarithmic spectra */
  double lmn[3];
  double mmax, mmin;            /**<  Used only in live or die situations, mmax=cos(angle-DANG_LIVE_OR_DIE)
                                  * and mmim=cos(angle+DANG_LIVE_OR_DIE).   In actually defining this
                                  * one has to worry about signs and exactly whether you are above
                                  * or below the plane */
  double renorm;                /**<  Factor used in Live or die option where it is 4*PI/domega;
                                  * 1.0 otherwise */
  int nphot[NSTAT];             /**<  Used to help define what happened to photons when trying
                                  * to construct this particular spectrum */
  int nscat;                    /**<  nscat>MAXSCAT -> select everything
                                  *  0 < nscat < MAXScat select only those photons which have
                                  *  scattered nscat times, number of scattering photons in
                                  * spectrum, if nscat is negative
                                  *  then all photons with more than |nscat| are included.  */
  int top_bot;                  /**<  0 ->select photons regardless of location
                                   >0     -> select only photons whose "last" position is above the disk
                                   <0    -> select only photons whose last position is below the disk */
  double x[3], r;               /**<  The position and radius of a special region from which to extract spectra.
                                   x is taken to be the center of the region and r is taken to be the radius of
                                   the region.   */
  double *f;                    /**<  The spectrum in linear (wavelength or frequency) units */
  double *lf;                   /**<  The specturm in log (wavelength or frequency)  units  */

  double *f_wind;               /**<  The spectrum of photons created in the wind or scattered in the wind. Created for
                                   reflection studies but possibly useful for other reasons as well. */
  double *lf_wind;              /**<  The logarithmic version of this */
}
spectrum_dummy, *SpecPtr;


extern SpecPtr xxspec;

/* Parameters used only by swind
 * swind_projecti	0 -> simply print the various parameters without
 * 			atempting toproject onto a yz plane
 * 			1 -> project onto a yz plane.
 */

extern int swind_min, swind_max, swind_delta, swind_project;
extern double *aaa;             // A pointer to an array used by swind

/***************************CDF STRUCTURE*****************************/
/* This is the structure for storing cumulative distribution functions. The CDFs are
generated from a function which is usually only proportional to the probability density
function or from an array.  It is sometimes useful, e.g. in calculating the reweighting function to
have access to the proper normalization.
*/


extern struct Cdf cdf_ff;
extern struct Cdf cdf_fb;
extern struct Cdf cdf_bb;
extern struct Cdf cdf_brem;



/* Variable used to allow something to be printed out the first few times
   an event occurs */
extern int itest, jtest;

extern char hubeny_list[132];    /**< /Location of listing of files representing hubeny atmospheres
                                   */


/* ***********************XBAND STRUCTURE *********************/
#define NBANDS 20

/** The xbands structure is associated with frequency limits used for photon
  * generation and for calculating heating and cooling
  *
  * xbands is used to control the number of photons generated in each
  * band
  */

struct xbands
{
  double f1[NBANDS], f2[NBANDS];  /**< The miniumum and maximum frequency of each  band */
  double alpha[NBANDS];
  double pl_const[NBANDS];
  double min_fraction[NBANDS];
  double nat_fraction[NBANDS];  /**< The fraction of the accepted luminosity in this band */
  double used_fraction[NBANDS];
  double flux[NBANDS];          /**< The "luminosity" within a band */
  double weight[NBANDS];        /**< The weight/energy for photons created with each  and */
  int nphot[NBANDS];            /**< The number of photons created in each band */
  int nbands;                   /**< Actual number of bands in use */
};

extern struct xbands xband;


/***************************FBSTRUC ***********************************
* The next section contains the freebound structures that can be used for both the
*
* ksl - 211030 most of the FBSTRUC was moved into recomb.c, since that is the only place
* these structures were used
* However kap_bf remains becuase it is used in the macro atom routines
*/

//#define NTEMPS        60      // The number of temperatures which are stored in each fbstruct
                                /* NSH this was increased from 30 to 60 to take account of 3 extra OOM
                                   intemperature we wanted to have in fb */
//#define NFB   20              // The maximum number of frequency intervals for which the fb emission is calculated

//struct fbstruc
//{
//  double f1, f2;
//  double cool[NIONS][NTEMPS];   //cooling rate due to radiative recombination
//  double lum[NIONS][NTEMPS];    //emissivity due to radiative recombinaion
//  double cool_inner[NIONS][NTEMPS];     //cooling rate due to recombinations to inner shells
//}
//freebound[NFB];

//extern double xnrecomb[NIONS][NTEMPS]; // There is only one set of recombination coefficients
//extern double xninnerrecomb[NIONS][NTEMPS];    // There is only one set of recombination coefficients

// extern double fb_t[NTEMPS];
// extern int nfb;                        // Actual number of freqency intervals calculated

/* kap_bf stores opacities for a single cell and as calculated by the routine kappa_bf.
 * It was made an external array to avoid having to pass it between various calling routines
 * but this means that one has to be careful that data is not stale.  It is required for
 * macro-atoms where bf is a scattering process, but not for the simple case.
 */

extern double kap_bf[NLEVELS];



// 12jun nsh - some commands to enable photon logging in given cells. There is also a pointer in the geo

extern FILE *pstatptr;        /**< pointer to a diagnostic file that will contain photon data for given cells
                                */
extern int cell_phot_stats;   /**< 1=do  it, 0=dont do it
                                */
#define  NCSTAT 10            /**< The maximum number of cells we are going to log */
extern int ncstat;            /**<  the actual number we are going to log */
extern int ncell_stats[NCSTAT];   /**< the numbers of the cells we are going to log */


/* Added variables which count number of times two situations occur (See #91) */
extern int nerr_no_Jmodel;
extern int nerr_Jmodel_wrong_freq;



enum partial_cells_enum
{ PC_INCLUDE = 0,   /**< Include partial cells, as has been done historically */
  PC_ZERO_DEN = 1,  /**<  Exclude partial cells, by setting their density to 0 */
  PC_EXTEND = 2     /**< Exclude partial cells, by extending the density of a full cell,
                      * into a partial cell
                      */
};


/***********************ADVANCED_MODES STRUCTURE **********************/
/**
  * Structure that contains the various switches that contoll which
  * which of the various modes for running the code are exercized
  *
  * Most (though perhaps not all) of these modes are accessed
  * via special inputs that are only read in when sirocco
  * is started up with athe -d command line option.
  *
  * All of the various variables take TRUE or FALSE, effectively
  */

struct advanced_modes
{
  int iadvanced;                /**< Set to TRUE when sirocco is invoked with -d flag.
                                 * By itself, simply causes Python ot acces the
                                 * various advanced modes commands in the .pf file
                                 */
  int extra_diagnostics;        /**< when set to TRUE, requests various additional
                                  information about what the user wishes to
                                  * write out
                                  */
  int save_cell_stats;          /**< when TRUE,save photons statistics by cell */
  int keep_ioncycle_windsaves;  /**< when TRUE, saves wind files for each ionization cycle */
  int keep_ioncycle_spectra;    /**< when TRUE, saves the total spectra for each ionization cycle */
  int make_tables;              /**< when TRUE, create tables showing various parameters for each cycle */
  int track_resonant_scatters;  /**< when TRUE, tracks resonant scatters */
  int save_photons;             /**< when TRUE, tracks photons (in photon2d) */
  int save_extract_photons;     /**< when TRUE, saves details on extracted photons */
  int adjust_grid;              /**< when TRUE, allows  wants to adjust the grid scale */
  int diag_on_off;              /**< when TRUE, prints extra diagnostics */
  int use_debug;                /**< when TRUE, prints out debug statements */
  int print_dvds_info;          /**< when TRUE, print out information on the velocity gradients */
  int keep_photoabs;            /**< when TRUE, keep photoabsorption in final spectrum */
  int quit_after_inputs;        /**< when TRUE, quit after inputs are read in.  Note that
                                  ** this flag is set from the command line with the -i option
                                  */
  int quit_after_wind_defined;  /**< when TRUE, quit after the wind grid has been defined and saved */
  int fixed_temp;               /**< do not alter temperature from that set in the parameter file */
  int zeus_connect;             /**< We are connecting to zeus, do not seek new temp and output
                                  * a heating and cooling file
                                  */
  int rand_seed_usetime;        /**< When TRUE use a seed set by the time.  The default case here
                                  * is to seed the random number generated with a fixed seed,
                                  * not based on time.  The decision of which seed to use
                                  * is set by the command line option --rseed*/
  int photon_speedup;           /**< The default is for each ionization cycle to have the same
                                  * number of photons bundles.  However, one can by using
                                  * the command lines swich -p cause the number of photons
                                  * to increase logarithmically in each ionization cyCle
                                  * The variable captures this option.
                                  */
  int save_rng;                 /**< IF TRUE, save the GSL RNG stage.  This flag is set by the
                                 * command line option --rng.  It causes the state of
                                 * random number generator to be written to one or more files
                                 * in a hidden directory .ring_ROOT.  A file is created
                                 * for each thread. The purpose of this is to handle
                                 * restarts of very long runs. */
  int load_rng;                 /**< IF TRUE, load the GSL RNG state.  This flag is also set
                                  * by the command line option --rng.  It causes
                                  * the state of the random numbe generator to be
                                  * read from a file.*/
  int store_matom_matrix;       /**< If TRUE, write the macro-atom matrix ot a file*/
  int jumps_for_detailed_spectra;   /**< If true, use the older deprecated jump method for
                                      * calculating emissivities
                                      * in detailed spectra.  Note that this is not
                                      * the same as using matrices or jumps for ionization
                                      * cyles, which is currently controlled directly
                                      * from the normal .pf file
                                      */
  int use_upweighting_of_simple_macro_atoms; /**< If TURE, use the deprecated method for simple atoms
                                                    *in macro scheme
                                                    */
  int run_xtest_diagnostics;     /**< If TRUE, then xtest is being run, which is
                                   * a special routine to run xtest  specific diagnostic
                                   * tests.
                                   */
  int partial_cells;             /**< Switch to decribe treatment of partial cells. */
  int searchlight;               /**< Switch to invoke search light option. This is
                                  a very experimental diagnostic mode in which photons
                                  originate at a specific place in the grid in the
                                  detailed spectrum stage.  It probably should be
                                  used in conjunction with analysing individual photons
                                  as they pass through the grid.  There are lots of
                                  issues with how the detailed spectra are constucted
                                  that make it less useful than it might seem. */
  int no_macro_pops_for_ions;     /* if true, then use the ion densities from the ionization mode
                                     for macro-atoms, rather than from macro_pops */
};

extern struct advanced_modes modes;


extern FILE *optr;               /**< pointer to a diagnostic file that will contain dvds information */



/***********************FILENAMES STRUCTURE ****************************/
/**
 *  Structure containing all of the file and directory names used
 * in the process of running a model
 */
struct filenames
{
  char root[LINELENGTH];        /**< main rootname */
  char windsave[LINELENGTH];    /**< wind save filename */
  char old_windsave[LINELENGTH];        /**< old windsave name */
  char specsave[LINELENGTH];    /**< spec save filename */
  char diag[LINELENGTH];        /**< diag file */
  char diagfolder[LINELENGTH];  /**< diag folder */
  char input[LINELENGTH];       /**< input name if creating new pf file */
  char new_pf[LINELENGTH];      /**< name of generated pf file */
  char lwspec[LINELENGTH];      /**< .log_spec_tot file (spectra from last ionization cycle
                                  * in log wavelength scale)
                                  */
  char lwspec_wind[LINELENGTH]; /**< .log_spec_tot_wind (same as above but in log units)  */
  char spec[LINELENGTH];        /**< .spec file (extracted spectra on linear scale) */
  char lspec[LINELENGTH];       /**<. spec file (extracted spectra on a log scale) */
  char spec_wind[LINELENGTH];   /**< .spec file (extracted spectra limited to wind photons on a linear scale) */
  char lspec_wind[LINELENGTH];  /**< .spec file (extracted spectra limited to wind photons on a log scale) */
  char disk[LINELENGTH];        /**< disk diag file name */
  char tprofile[LINELENGTH];    /**< non standard tprofile fname */
  char phot[LINELENGTH];        /**< photfile e.g. sirocco.phot */
  char windrad[LINELENGTH];     /**< wind rad file */
  char extra[LINELENGTH];       /**< extra diagnositcs file opened by init_extra_diagnostics */
};

extern struct filenames files;



#define NMAX_OPTIONS 20

/* these two variables are used by xdefine_phot() in photon_gen.c
   to set the mode for get_matom_f()in matom.c and tell it
   whether it has already calculated the matom emissivities or not. */
#define CALCULATE_MATOM_EMISSIVITIES 0
#define USE_STORED_MATOM_EMISSIVITIES 1

/* Used in macro_gov elsewhere to descibe choices between being or going
   to a kpkt or macro atom state */
#define KPKT 2
#define MATOM 1
/* modes for kpkt calculations */
#define KPKT_MODE_CONTINUUM  0  /* only account for k->r processes */
#define KPKT_MODE_ALL        1  /* account for all cooling processes */
#define KPKT_MODE_CONT_PLUS_ADIABATIC 2 /* account for k->r and adiabatic destruction */


/* whether or not to use the implicit/accelerated macro-atom scheme, in which
   a matrix inversion is used in the emissivity calcualtion rather than
   a MC sampling of the transition probabilities */
#define MATOM_MATRIX_EMISSIVITIES  TRUE
#define STORE_B_MATRIX TRUE
#define MATOM_TRANSITION_MODE MATOM_MATRIX

/* Variable introducted to cut off macroatom / estimator integrals
   when exponential function reaches extreme values. Effectivevly a max limit
   imposed on x = hnu/kT terms */

#define ALPHA_MATOM_NUMAX_LIMIT 30      /**< maximum value for h nu / k T to be considered in integrals */
#define ALPHA_FF 100.           /**< maximum h nu / kT to create the free free CDF
                                  */


/* non-radiative heat flow mode */
#define KPKT_NET_HEAT_MODE 0


/* DIAGNOSTIC for understanding problems with imported models
 */


#define BOUND_NONE   0
#define BOUND_INNER_CONE  1
#define BOUND_OUTER_CONE  2
#define BOUND_RMAX 3
#define BOUND_RMIN 4
#define BOUND_ZMIN 5
#define BOUND_ZMAX 6
#define BOUND_INNER_RHO 7
#define BOUND_OUTER_RHO 8

extern int xxxbound;


/** Structure associated with rdchoice.  This
 * structure is required only in cases where one
 * wants to use rdchoice multiple times with different
 * options.  But it can, or course be used for
 * any such call.  It allows one to associate
 * an input word with an output value, using the routine
 * get_choices.  There needs to be one structure
 * for each input variable.  At present, this is only
 * used for the selection of spec_types
 */
#define MAX_RDPAR_CHOICES 20 

typedef struct rdpar_choices
{
  char choices[MAX_RDPAR_CHOICES][LINELENGTH];
  int vals[MAX_RDPAR_CHOICES];
  int n;
} dummy_choices, *ChoicePtr;

extern struct rdpar_choices zz_spec;


/* ************************************************************************** */
/**
 * The constants and structures which require global scope throughout the wind
 * import functions are defined here. This includes the structure which will
 * temporarily hold the imported model data.
 *
 * ************************************************************************** */


/*
 * The following definitions are used to try and improve the readability for
 * some of the code used to read in 1D and 2D grids from an input file. They
 * refer to how many columns are in the data file and what will be read in.
 */

#define READ_NO_TEMP_1D          4
#define READ_ELECTRON_TEMP_1D    (READ_NO_TEMP_1D + 1)
#define READ_BOTH_TEMP_1D        (READ_NO_TEMP_1D + 2)

#define READ_NO_TEMP_2D          9
#define READ_ELECTRON_TEMP_2D    (READ_NO_TEMP_2D + 1)
#define READ_BOTH_TEMP_2D        (READ_NO_TEMP_2D + 2)

/**
 * The Import structure will contain all of the required information for
 * creating a wind grid using an imported model.
 */

struct Import
{
  int init_temperature;         /**<  initialise to t.wind.init if TRUE */
  int ncell;                    /**<  the total number of cells read in */
  int ndim, mdim;               /**<  the number of coordinates in the n and m dimensions */
  int *i, *j;                   /**<  the i (row) and j (column) elements */
  int *inwind;                  /**<  flag for the cell being inwind or not inwind */
  double *x, *z, *r, *theta;    /**<  the x/r or z/theta coordinates of the grid in cgs units */
  double *v_x, *v_y, *v_z;      /**<  the velocity in Cartesian coordinates in cgs units */
  double *v_r;                  /**<  the radial velocity in cgs units */
  double *mass_rho;             /**<  the mass density in cgs units */
  double *t_e, *t_r;            /**<  the electron and radiation temperature in Kelvin */
  double *wind_x, *wind_z;      /**<  the wind grid coordinates */
  double *wind_midx, *wind_midz;        /**<  the wind grid mid points */
};

extern struct Import *imported_model;   // MAX_DOM is defined in sirocco.h and as such import.h has to be included after

/* the functions contained in log., rdpar.c and lineio.c are
   declare separately from templates. This is because some functions
   only use log.h and don't use sirocco.h due to repeated definitions */
#include "log.h"
#include "version.h"
#include "templates.h"

/* We're going to keep the matrix GPU functions seperate from the other templates */
#ifdef CUDA_ON
#include "matrix_gpu.h"
#endif
