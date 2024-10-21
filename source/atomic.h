/* The header that descibes atomic data */

/* The constants have been separated from atomic.h to enable one to create routines that do not 
   require the atomic data, but do need various physical constants */

#include "constants.h"

/* The next term attempts to globally define a minimum density to prevent zero divides in some routines */
#define DENSITY_MIN		1.e-20

/*
   Note that the structure ele array may not be completely filled.  In this structure, the dimension is simply
   the order in which elements are read in, and one may skip elements (which are not of interest).

   So that one can easily find, say carbon, however, *elz contains pointers to ele such that elz[6] is the
   pointer to carbon.  Elements for which there is no data, point to the dummy element structure foo_element

   If you want to cycle through the elements one will act on ele directly.  If you want to know something
   about a specific element *elz may be better, at least that is the idea

 */


#define NELEMENTS		50      /**< Maximum number of elements to consider */
extern int nelements;                  /**< The actual number of ions read from the data file */
#define NIONS		500     /**< Maximum number of ions to consider */
extern int nions;                      /**< The actual number of ions read from the datafile */
#define NLEVELS 	12000   /**<  Maximum number of levels for all elements and ions */
extern int nlevels;                    /**< These are the actual number of levels which were read in */
#define NLTE_LEVELS	12000   /**<  Maximum number of levels to treat explicitly */
extern int nlte_levels;                /**<  Actual number of levels to treat explicityly */


/* AUGER NOTE: recommended to increase NLEVELS_MACRO to at least 500 for Auger macro-atoms */
#define NLEVELS_MACRO   600     /**<  Maximum number of macro atom levels. (SS, June 04) */
extern int nlevels_macro;              /**<  Actual number of macro atom levels. (SS, June 04) */
#define NLINES 		200000  /**<  Maximum number of lines to be read */
extern int nlines;                     /**<  Actual number of lines that were read in */
extern int nlines_macro;               /**<  Actual number of Macro Atom lines that were read in.  New version of get_atomic
                                   data assumes that macro lines are read in before non-macro lines */
#define N_INNER     10          /**< Maximum number of inner shell ionization cross sections per ion */
extern int n_inner_tot;                /**< The actual number of inner shell ionization cross sections in total */


#define LARGE_MACRO_ATOMS FALSE 
/* sometimes we have to compile with larger hardwired values here for large macro-atoms e.g. Fe17to27 dataset */
#if LARGE_MACRO_ATOMS
  #define NBBJUMPS         200    /**<  Maximum number of Macro Atom bound-bound jumps from any one configuration (SS) */
  #define NBFJUMPS         200    /**<  Maximum number of Macro Atom Bound-free jumps from any one configuration (SS) */
  #define NAUGER_MACRO     2000        /* number of Auger processes */
#else
  /* AUGER NOTE: recommended to increase these values to at least 200 for Auger macro-atoms */
  #define NBBJUMPS         100    /**<  Maximum number of Macro Atom bound-bound jumps from any one configuration (SS) */
  #define NBFJUMPS         100    /**<  Maximum number of Macro Atom Bound-free jumps from any one configuration (SS) */
  #define NAUGER_MACRO     200        /* number of Auger processes */
#endif

extern int nauger_macro;        /* the number of auger processes read in associated with macro atoms */
#define NAUGER_ELECTRONS 4

#define MAXJUMPS          1000000       /**<  The maximum number of Macro Atom jumps before emission (if this is exceeded
                                           it gives up (SS) */
#define NAUGER 2                /**< Maximum number of "auger" processes */
extern int nauger;                     /**< Actual number of innershell edges for which autoionization is to be computed */



/** Element contains physical parameters that apply to element as a whole and 
   * provides an index to the various ions associated with that element
   */

typedef struct elements
{
  char name[20];                /**<  Element name */
  int z;                        /**<  Atomic number */
  int firstion;                 /**<  Index into struct ions  ion[firstion] is the lowest ionization state of this ion */
  int lastion;                  /**<  Index into struct ions ion[lastion] is the higest ionization state of this ion */
  int nions;                    /**<  The number of ions actually read in from the data file for this element */
  double abun;                  /**<  Abundance */
  double atomic_weight;         /**<  Atomic weight */
  int istate_max;               /**<  highest ionization stage of element */
}
ele_dummy, *ElemPtr;

extern ElemPtr ele;

extern double rho2nh;                  /**<  Conversion constant from rho to nh the number of H atoms per unit volume */

/** ions is the basic structure describing individual ions.  
  *
  * It is filled from 0 up to nions. 
  *
  * It provides some basic information about each ion, and contints indices
  * to atomic data, such as photioninization x-sections, for each ion.
 */


typedef struct ions
{
  int z;                        /**<  defines the element for this ion */
  int istate;                   /**<  1=neutral, 2 = once ionized, etc. */
  int nelem;                    /**<  index to elements structure */
  double ip;                    /**<  ionization potential of this ion (converted from eV to ergs by get_atomic) */
  double g, log_g;              /**<  multiplicity of ground state, note that this is not totally consistent
                                   with energy levels and this needs to be reconciled  
	                               there is also a log version for use in freebound calculations*/
  int nmax;                     /**<  The maximum number of allowed lte configurations for this ion. 
                                   Used only by get_atomic */
  int n_lte_max;                /**<  The maximum number of allowed nlte configurations for this ion. 
                                   Used only by get_atomic */
  int firstlevel;               /**<  Index into config struc; should also refer to ground state for this ion */
  int nlevels;                  /**<  Actual number of "lte" levels for this ion. "lte" here refers to any level 
                                   in which densities are to be calculated on the fly, instead of the beginning 
                                   or end of an ionization cycle */
  int first_nlte_level;         /**<  Index into config structure for the nlte_levels. There will be nlte 
                                   levels associated with this ion */
  int first_levden;             /**<  Index into the levden array in wind structure, not necessairly 
                                   the same as first_nlte_level  (Name changed 080810 -- 62 */
  int nlte;                     /**<  Actual number of nlte levels for this ion */
  int phot_info;                /**< 
                                  ** -1 means there are no photoionization cross sections for this ion, 
                                  ** 0  means the photoionization is given on an ion basis (e.g. for the 
                                  ** ground state using VFKY and the xphot structure
                                  ** 1 means photoionization is given on a level basis, using topbase value
                                  ** Topbase photoinization trumps VFKY photoionization
				  */
  int macro_info;               /**< Identifies whether ion is to be treated using a Macro Atom approach.
                                   ** set to -1 initially (means not known) 
                                   ** set to 1 if macro atom method is used
                                   ** set to 0 if this ion is a simple ion.            
                                   Note: program will exit if -1 before leaving get_atomicdata
                                 */
  int ntop_first;               /**< Internal index into topbase photionization structure */
  int ntop_ground;              /**< Index to the ground state topbase photoionization state */
  int ntop;                     /**< Number of topbase photoionization cross sections for this ion */
  int nxphot;                   /**< Internal index into VFKY photionionization structure.  There
                                   is only one of these per ion */
  int lev_type;                 /**< The type of configurations used to describe this ion 
                                   ** set to -1 initially (means not known)
                                   ** set to  0 if Kurucz (Note ultimately this may be treated similarly
                                   to topbase.  
                                   ** set to  1 if configurations are described as macro atoms, even if
                                   macro_methods are not used
                                   ** set to  2 if topbase inputs, even if their ar no so-called non-lte
                                   levels, i.e. levels in the levden array
                                 */
  int drflag;                   /**<  The number of dielectronic recombination parameter types read in. 0
                                   probably means it has no data. 2 is good, 3 or 1 means an error has taken
                                   place - this is trapped in get_atomicdata.c */
  int nxdrecomb;                /**<  link into the drecomb structure to give the location of the dielectronic
                                   recombination coefficients for this ion */
  int total_rrflag;             /**< Flag to say wether we have badnell style total radiative rate 
                                   coefficients for this ion */
  int nxtotalrr;                /**< index into the bad_t_rr structure to give the location of the
                                   Badnell fit coefficient for the total radiative recombination rate for 
                                   this ion if it exists */
  int bad_gs_rr_t_flag;         /**< Flag to say whether we have badnell style resolved ground state radiative temperature
                                   data for this ion */
  int bad_gs_rr_r_flag;         /**< Flag to say wether we have badnell style resolved ground state radiative rate 
                                   coefficients for this ion */
  int nxbadgsrr;                /**< index into the bad_gs_rr structure to give the location of the
                                   Badnell fit coefficient for the resolved ground state recombination rate for 
                                   this ion if it exists */
  int dere_di_flag;             /**< Flag to say if we have DERE direct ionization data for this ion */
  int nxderedi;                 /**< index into the dere direct ionization structure to give the location of the data for this ion */
  int nxinner[N_INNER];         /**< index to each of the inner shell cross sections associtated with this ion */
  int n_inner;                  /**< The number of inner shell cross section associated with this ion */
  int n_ch_ex;                  /**< The number of the charge exchange rate that applies to this ion */

}
ion_dummy, *IonPtr;

extern IonPtr ion;


/** The structrued descibing the  energy levels.  

 The number of configurations is given by nlevels 
 */

typedef struct configurations
{
  int z;                        /**< The element associated with this configuration */
  int istate;                   /**< The ion associated with the configuration */
  int nion;                     /**< Internal index to ion structure */
  int nden;                     /**< Internal index to levden array in plasma structure. -1 if no such level */
  int isp;                      /**< Topbase (ISLP) description of angular momentum (2s+1)*100+l*10+p */
  int ilv;                      /**<  Unique Level number (within an ion. Topbase ilv when Topbase record */
  int macro_info;               /**<  Unambigous way to determine this is part of a macro-level 
                                  **  set to -1 initially (means not known)
                                  **  set to 1 if macro atom method is used
                                  ** set to 0 if macro atom method is not used for this configuration 
                                  *
                                  * Note: program will exit if -1 before leaving get_atomicdata
                                  */
  int n_bbu_jump, n_bbd_jump;   /**<  No. of bound-bound upward/downward jumps available to this configuration (SS) */
  int n_bfu_jump, n_bfd_jump;   /**<  No. of bound-free upward/downward jumps available to this configuration (SS) */
  double g, log_g;              /* multiplicity for level and log version for use in freebound integrations */
  double q_num;                 /* principal quantum number.  In Topbase this has non-integer values */
  double ex;                    /*excitation energy of level */
  double rad_rate;              /* Total spontaneous radiative de-excitation rate for level */
  int bbu_jump[NBBJUMPS];       /* list of upwards bb jumps from this configuration */
  int bbd_jump[NBBJUMPS];       /* list of downwards bb jumps from this configuration */
  int bfu_jump[NBFJUMPS];       /* list of upwards bf jumps from this configuration */
  int bfd_jump[NBFJUMPS];       /* list of downwards bf jumps from this configuration (SS) */
  int bbu_indx_first;           /* index to first MC estimator for bb jumps from this configuration */
  int bfu_indx_first;           /* index to first MC estimator for bf jumps from this configuration */
  int bfd_indx_first;           /* index to first rate for downward bf jumps from this configuration */
  int iauger;                   /* if this configuration is a vacancy state in a macro-atom, then this 
                                   indexes where the entry is in the Auger structure. Otherwise, -1. */
  int nauger;                   /* number of possible Auger channels (nominal maximum of 4) */
}
config_dummy, *ConfigPtr;

extern ConfigPtr xconfig;

typedef struct auger
{
  int nion;                     /* The ion no (in sirocco) of the transition */
  int z, istate;                /* element and ion associated with the line */
  int nconfig;                  /* the entry in the config structure where the vacancy state lies */
  int iauger;                   /* the index to this auger entry */
  int nauger;                   /* number of Auger jumps */
  double Avalue_auger;          /* the total A value associated with creating Auger electrons */
  int nconfig_target[NAUGER_ELECTRONS]; /* the configuration of the upper ion we are jumping to */
  double branching_ratio[NAUGER_ELECTRONS];     /* probability of creating 0-3 Auger electrons */
} auger_dummy, *AugerPtr;

extern AugerPtr auger_macro;



/** Structure to describe the lines */

typedef struct lines
{
  int nion;                     /**< The ion no (in sirocco) of the transition */
  int z, istate;                /**< element and ion associated with the line */
  double gl, gu;                /**< multiplicity of lower and upper level respectively */
  int nconfigl, nconfigu;       /**< The configuration no (in sirocco) of the transition */
  int levl, levu;               /**< level no of transition..parallel/redundant with el,eu hopefully */
  int macro_info;               /**<  Identifies whether line is to be treated using a Macro Atom approach.
                                  **  set to -1 (not known initially) 
                                  **  set to 0 if not a macro atom line  
                                  **  set to 1 if a macro atom line  
                                  *
                                  * Note: program will exit if -1 before leaving get_atomicdata
                                  */
  double freq;                  /**<  The frequency of the resonance line */
  double f;                     /**< oscillator strength.  Note: it might be better to keep PI_E2_OVER_MEC flambda times this.
                                   Could do that by initializing */
  double el, eu;                /**<  The energy of the lower and upper levels for the transition */
  double pow;                   /**< The power in the lines as last calculated in total_line_emission */
  int where_in_list;            /**<  Position of line in the line list: i.e. lin_ptr[line[n].where_in_list] points
                                   to the line. Added by SS for use in macro atom method. */
  int down_index;               /**<  This is to map from the line to knowing which macro atom jump it is (and therefore find
                                   the estimator with which it is associated. The estimator is identified by the
                                   upper configuration (nconfigu) and then down_index (for deexcitation) or the lower
                                   configuration (nconfigl) and then up_index. (SS) */
  int up_index;
  int coll_index;               /**<  A link into the collision strength data, if its -999 it means there is no data and van reg should be used */

}
line_dummy, *LinePtr;


extern LinePtr line, lin_ptr[NLINES];  /**<  line[] is the actual structure array that contains all the data, *lin_ptr
                                   is an array which contains a frequency ordered set of ptrs to line */
                                /**<  fast_line (added by SS August 05) is going to be a hypothetical
                                   rapid transition used in the macro atoms to stabilise level populations */
extern struct lines fast_line;

extern int nline_min, nline_max, nline_delt;   /**<  Used to select a range of lines in a frequency band from the lin_ptr array 
                                           in situations where the frequency range of interest is limited, including for defining which
                                           lines come into play for resonant scattering along a line of sight, and in
                                           calculating band_limit luminosities.  The limits are established by the
                                           routine limit_lines.
                                         */


/* coll_stren is the collision strength interpolation data extracted from Chianti */


#define N_COLL_STREN_PTS	20      //The maximum number of parameters in the interpolations
extern int n_coll_stren;

typedef struct coll_stren
{
  int n;                        /**< Internal index */
  int lower;                    /**< The lower energy level - this is in Chianti notation and 
                                  * is currently unused
                                  */
  int upper;                    /**< The upper energy level - this is in Chianti notation and is currently unused
                                  */
  double energy;                /**< The energy of the transition */
  double gf;                    /**< The effective oscillator strength 
                                  *
                                  * NSH thinks this is oscillator strtength x lower level multiplicity
                                  */
  double hi_t_lim;              /**< The high temperature limit */
  int n_points;                 /**< The number of points in the spline fit */
  int type;                     /**< The type of fit, this defines how one computes 
                                  * the scaled temperature and scaled coll strength
                                  */
  float scaling_param;          /**< The scaling parameter C used in the Burgess and Tully calculations
                                  */
  double sct[N_COLL_STREN_PTS]; /**< The scaled temperature points in the fit
                                  */
  double scups[N_COLL_STREN_PTS];       /**< The scaled coll strengths in ythe fit.
                                          */
} Coll_stren, *Coll_strenptr;

extern Coll_stren coll_stren[NLINES];



extern int nxphot;                     /**< The actual number of ions for which there are VFKY photoionization x-sections */
extern double phot_freq_min;           /**< The lowest frequency for which photoionization can occur */
extern double inner_freq_min;          /**< The lowest frequency for which inner shell ionization can take place */

#define NCROSS 3000            /**<  Maximum number of x-sections for a single photionization process */

#if LARGE_MACRO_ATOMS
  #define NTOP_PHOT 600           /**<  Maximum number of photoionisation processes.  */
#else
  /* AUGER NOTE: recommended to increase NTOP_PHOT to 1000 for Auger macro-atoms */
  #define NTOP_PHOT 1000           /**<  Maximum number of photoionisation processes.  */
#endif
extern int ntop_phot;                  /**<  The actual number of TopBase photoionzation x-sections */
extern int nphot_total;                /**<  total number of photoionzation x-sections = nxphot + ntop_phot */

/** 
  * Photionization x-sections in topbase format
  *
  * Note
  *  If the old topbase treatment is to be replaced by Macro Atoms perhaps this
                                   can be dropped - need for backward compatibility? (SS) 
  */
typedef struct topbase_phot
{
  int nlev;                     /**< Internal index to the config structure for the lower state */
  int uplev;                    /**< Internal index to the config structure for the upper state (SS) */
  int nion;                     /**< Internal index to the    ion structure for this x-section */
  int z, istate;
  int np;                       /**< the number of points in the corr section fit */
  int n, l;                     /**< Shell and subshell, used for inner shell */
  int nlast;                    /**<  nlast is an index into the arrays freq and x.  It allows 
                                   a quick check to see whether one needs to do a search
                                   for the correct array elements, or whether this process
                                   can be short circuited */
  int n_elec_yield;             /**< Index to the electron yield array - only used for inner shell ionizations */
//  int n_fluor_yield;            /**< Index to the fluorescent photon yield array - only used for inner shell ionizations */
  int macro_info;               /**< Identifies whether line is to be treated using a Macro Atom approach.
                                  ** set to -1 initially
                                  ** set to 0 if not a macro atom line  
                                  ** set to 1 if a macro atom line  (ksl 04 apr)
                                  *
                                  * Note: Program will exit before leaving get_atomicdata if not initallized to 0 or 1
                                 */
  int down_index;               /**< This is to map from the line to knowing which macro atom jump it is (and therefore find
                                   the estimator with which it is associated. The estimator is identified by the
                                   upper configuration (uplev) and then down_index (for deexcitation) or the lower
                                   configuration (nlev) and then up_index. (SS) */
  int up_index;
  int use;                      /**< It we are to use this cross section. This allows unused VFKY cross sections to sit in the array. */
  double freq[NCROSS], log_freq[NCROSS], x[NCROSS], log_x[NCROSS];
                                                                 /**< frequency and cross sections, plus log versions for FB integrals*/
  double f, log_f, sigma, log_sigma;            /**< last freq, last x-section and log versions*/
} Topbase_phot, *TopPhotPtr;

extern Topbase_phot phot_top[NLEVELS];
extern TopPhotPtr phot_top_ptr[NLEVELS];       /**<  Pointers to phot_top in threshold frequency order - this */

extern Topbase_phot inner_cross[N_INNER * NIONS];  /**< Pointer to inner shell cross sections which use the same structure type */
extern TopPhotPtr inner_cross_ptr[N_INNER * NIONS];  /**< Pointer to inner shell cross sections in frequency order */



/** The structure for electron yield data for inner shell ionization from Kaastra and Mewe */
typedef struct inner_elec_yield
{
  int nion;                     /**< Index to the ion which was the parent of the inner shell ionization */
  int z, istate;
  int n, l;                     /**< Quantum numbers of shell */
  double prob[10];              /**< The probability for between 1 and 10 electrons being ejected */
  double I;                     /**< Ionization energy */
  double Ea;                    /**< Average electron energy */
} Inner_elec_yield, Inner_elec_yieldPtr;

extern Inner_elec_yield inner_elec_yield[N_INNER * NIONS];

/** This structure for the flourescent photon yield following inner shell ionization from Kaastra and Mewe*/
typedef struct inner_fluor_yield
{
  int nion;                     /**< Index to the ion which was the parent of the inner shell ionization */
  int z, istate;                /**< atomic number and state (astro notation) of the parent */
  int n, l;                     /**< shell and subshell of the parent */
  double freq;                  /**< the rest frequency of the photon emitted */
  double yield;                 /**< number of photons per ionization */
} Inner_fluor_yield, Inner_fluor_yieldPtr;

extern Inner_fluor_yield inner_fluor_yield[N_INNER * NIONS];



/** ground_fracs is the table of fractions of recombinations going straight to the ground state */
struct ground_fracs
{
  int z, istate;                /**< element and ion */
  double frac[20];              /**<  this is a table containing the fraction of recombinations going directly
                                   to the ground state as a function of temperature. ground_frac[0] is or t=5000
                                   and then we go in steps of 5000 to ground_frac[19] which is for t=1e5. these
                                   fractions must have been computed elsewhere */
};

extern struct ground_fracs ground_frac[NIONS];


#define MAX_DR_PARAMS 9         //This is the maximum number of c or e parameters.
#define DRTYPE_BADNELL	    0
#define DRTYPE_SHULL	    1
extern int ndrecomb;            //This is the actual number of DR parameters

/** This is the struture to store information about dielectronic recombination
  * 081115 nsh New structure and variables to hold the dielectronic recombination rate data
  *set up to accept the korista data from the university of strahclyde website.
  */

typedef struct dielectronic_recombination
{
  int nion;                     /**< Internal cross reference to the ion in the ion structure thsat it refers to
                                  */
  int nparam;                   /**< the number of parameters - it varies from ion to ion
                                  */
  double c[MAX_DR_PARAMS];      /**< c parameters */
  double e[MAX_DR_PARAMS];      /**< e parameters */
  double shull[4];              /**< schull DR parameters */
  int type;                     /**< defines wether we have a schull type or a badnell type */
} Drecomb, *Drecombptr;


extern Drecomb drecomb[NIONS];  //set up the actual structure

extern double dr_coeffs[NIONS]; //this will be an array to temprarily store the volumetric dielectronic recombination rate coefficients for the current cell under interest. We may want to make this 2D and store the coefficients for a range of temperatures to interpolate.


#define T_RR_PARAMS         6   //This is the number of parameters.
#define RRTYPE_BADNELL	    0
#define RRTYPE_SHULL	    1
extern int n_total_rr;
/** This is the toal recombination structure */
typedef struct total_rr
{
  int nion;                     /**< Internal cross reference to the ion that this refers to
                                  */
  double params[T_RR_PARAMS];   /**< There are up to six parameters. If the last two are zero, we still
                                   the data in the same way, but they have no effect - NB - 
                                   important to zero these! */
  /* NSH 23/7/2012 - This array will double up for Shull parameters */
  int type;                     /**< NSH 23/7/2012 - What type of parampeters we have for this ion */
} Total_rr, *total_rrptr;

extern Total_rr total_rr[NIONS];        //Set up the structure

#define BAD_GS_RR_PARAMS 19     //This is the number of points in the fit.
extern int n_bad_gs_rr;

/**  Badnell */
typedef struct badnell_gs_rr
{
  int nion;                     //Internal cross reference to the ion that this refers to
  double temps[BAD_GS_RR_PARAMS];       //temperatures at which the rate is tabulated
  double rates[BAD_GS_RR_PARAMS];       //rates corresponding to those temperatures
} Bad_gs_rr, *Bad_gs_rrptr;

extern Bad_gs_rr bad_gs_rr[NIONS];      //Set up the structure


#define DERE_DI_PARAMS 20       //This is the maximum number of points in the fit.
extern int n_dere_di_rate;
typedef struct dere_di_rate
{
  int nion;                     //Internal cross reference to the ion that this refers to
  int nspline;
  double temps[DERE_DI_PARAMS]; //temperatures at which the rate is tabulated
  double rates[DERE_DI_PARAMS]; //rates corresponding to those temperatures
  double xi;
  double min_temp;
} Dere_di_rate, *Dere_di_rateptr;

extern Dere_di_rate dere_di_rate[NIONS];        //Set up the structure

extern double di_coeffs[NIONS]; //This is an array to store the di_coeffs 
extern double qrecomb_coeffs[NIONS];    //JM 1508 analogous array for three body recombination 

#define MAX_GAUNT_N_GSQRD 100   //Space set aside for the number of parameters for scaled inverse temperature

extern int gaunt_n_gsqrd;       //The actual number of scaled temperatures

typedef struct gaunt_total
{
  float log_gsqrd;              //The scaled electron temperature squared for this array
  float gff;
  float s1, s2, s3;
} Gaunt_total, *Gaunt_totalptr;

extern Gaunt_total gaunt_total[MAX_GAUNT_N_GSQRD];      //Set up the structure



#define MAX_CHARGE_EXCHANGE 100 //Space set aside for charge exchange parameters

extern int n_charge_exchange;   //The actual number of scaled temperatures

typedef struct charge_exchange
{
  int nion1;                    //The ion which will be ionized - normally hydrogen
  int nion2;                    //The ion which will be recombining
  double a, b, c, d;            //The parameters for the fit
  double tmin;                  //The minimum temperature which the fit is valif for
  double tmax;                  //The maximum temperature which the fit is valif for
  double energy_defect;         //The energy defect for the reaction - used for heating
  double delta_e_ovr_k;         //The boltzman factor for charge exchange ionization (only a few records have this)

} Charge_exchange, *Charge_exchange_ptr;

extern Charge_exchange charge_exchange[MAX_CHARGE_EXCHANGE];    //Set up the structure

extern double charge_exchange_recomb_rates[NIONS];      //An array to store the actual recombination rates for a given temperature - 
//there is an estimated rate for ions without an actual rate, so we need to dimneions for ions.
extern double charge_exchange_ioniz_rates[MAX_CHARGE_EXCHANGE]; //An array to store the actual ionization rates for a given temperature



/* a variable which controls whether to save a summary of atomic data
   this is defined in atomic.h, rather than the modes structure */
extern int write_atomicdata;
