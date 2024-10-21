

/* Note that there are some limits on how large NMODS*NWAVES can be.  When you exceed this limit, you will get a segmentation
 * fault on startup */

/* Note that there is some overlap between these definitions and those in sirocco.h and so mixing these two
 * header files can be dangerous.  
 *
 * NDIM --> MOD_NDIM  -- 081026
 */

 /* 1405 JM -- Increased PDF array for use with disk14 models- also removed duplication of ncomps */
/* 1409 JM -- Increased LINELEN to 160 */

#define NWAVES  28000           // The maximum number of wavelength bins in models
#define NDIM	10              // The maximum number of free parameters
#define NCOMPS	10              //The maximum number of separate components
#define NPARS	10              //The maximum number of parameters in a component (not all variable)
#define NMODS   2000            //The total number of models read in in all components
#define LINELEN 160             //This should always be the same as LINELENGTH in sirocco.h!



extern int ncomps;              // The number of components that have been read.  
extern int nmods_tot;           // The total number of models that have been read in 


/** The Model structure that contais the  individual continuum model. 
  * 
  * mods is the set of all models that are read in; the various
  * sets of models are distinguished by the ModSum structure 
  */
typedef struct Model
{
  char name[LINELEN]; /**< file name for the model */
  double par[NPARS];  /**< Parameters, usually temperatur and gravity for the model */
  double w[NWAVES];  /**< Wavelengths (in A) */
  double f[NWAVES];  /**< Fluxes in ergs/cm**2/s/A */
  int nwaves;        /**< The number of wavelength flux pairs for this model */
} model_dummy, *ModelPtr;
extern ModelPtr mods;

/**
  * Modsum is the structure that is used to describe a set of models
  * used for various components of the system.
  *
  * When models are read in there is usual one set of models for each
  * component of the system.
  * There is one element of comp for each set of models of the same type, 
  * i.e. if
  *one reads in a list of WD atmosphers this will occupy one componenet here 
 */

struct ModSum
{
  char name[LINELEN];  /**< The name of the file that is read in to 
                         * describe each set of models
                         */
  int npars;            /**< Number of parameters that describe the model */
  int nmods, modstart, modstop; /**< number of models of this type and 
                                  * the index into the model struct
                                  */
  double min[NPARS];            /**< The minimum value for each 
                                  * "free" parameter of the model
                                  */
  double max[NPARS];            /**< The maxium for each "free" parameter
                                  * for this set of models
                                  */
  int nwaves;                   /**<  The number of wavelength flux pairs
                                  * for this set of moels
                                  * All models in each comp should have same wavelengths;
                                  */
  struct Model xmod;            /**< The current intepolated model of this type 
                                  */
  struct Cdf xcdf;              /**< The current cumulative distribution 
                                  * function for this component
                                  */
};

extern struct ModSum comp[NCOMPS];
