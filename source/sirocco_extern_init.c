
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "atomic.h"

#include "sirocco.h"


int np_mpi_global;              ///< Global variable which holds the number of MPI processes

int rank_global;                ///<  Rank of a particular thread

int verbosity;                  ///< verbosity level. 0 low, 10 is high 

int rel_mode;                   ///< How doppler effects and co-moving frames are treated

int run_xtest;                  ///< Variable if TRUE causes a special test mode to be run 

int NDIM2;                      ///< The total number of wind cells in wmain
int NPLASMA;                    ///< The number of cells with non-zero volume or the size of plasma structure

double DFUDGE;

double SMAX_FRAC;               /**< In translate_in_wind, a limit is placed on the maximum distance a
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
double DENSITY_PHOT_MIN;        /**< This constant is a minimum density for the purpose of calculating
                                   photoionization heating and recombination cooling.  It is important that heating and cooling
                                   be calculated self-consistently.  Program speed is somewhat sensitive 
                                   to this parameter, at the 10% level if raised from 1e-3 to 1.  There is a 
                                   trade-off since lower minima may give better results, especially for macro atoms. */

double PHOT_RANGE;              /**< When a variable number of photons are called in different ionization
                                   cycles this is the log of the difference between NPHOT_MAX
                                   and the value in the first cycle
                                 */
int NPHOT_MAX;                  /**< The maximum number of photon bundles created per cycle */
int NPHOT;                      /**< The number of photon bundles created, defined in setup.c */

int NWAVE_MAX;
int NWAVE_EXTRACT;              ///< The number of wavelength bins for spectra during the spectrum cycles
int NWAVE_NOW;                  ///< Either NWAVE_IONIZ or NWAVE_EXTRACT depending on whether in ionizaiton of spectrum cycles

plane_dummy plane_m2_near, plane_m2_far;

DomainPtr zdom;                 ///< This is the array pointer that contains the domains
int current_domain;             ///<  This integer is used by swind only

struct geometry geo;

struct xdisk disk, qdisk;   /**< disk defines zones in the disk which in a specified frequency band emit equal amounts
                                   of radiation. disk gets reinitialized whenever the frequency interval of interest
                                   is changed.  qdisk stores the amount of heating of the disk as a result of
                                   illumination by the star or wind. It's boundaries are fixed throughout a cycle */

struct blmodel blmod;

WindPtr wmain;

PlasmaPtr plasmamain;

PhotStorePtr photstoremain;

MatomPhotStorePtr matomphotstoremain;

MacroPtr macromain;

//int xxxpdfwind;                 ///< When 1, line luminosity calculates pdf

int size_Jbar_est, size_gamma_est, size_alpha_est;

PhotPtr photmain;               /**< A pointer to all of the photons that have been created in a subcycle. Added to ease 
                                   breaking the main routine of sirocco into separate rooutines for inputs and 
                                   running the program */

int nspectra;                   /**< After create_spectrum, the number of elements allocated for s, or
                                   alternatively the number of spectra one has to work with.  Note that
                                   general s[0],s[1] and s[2] are the escaping, scattered and absorbed photons,
                                   while elements higher than this will contain spectra as seen by different observers */

int nscat[MAXSCAT + 1], nres[MAXSCAT + 1], nstat[NSTAT];

SpecPtr xxspec;

int swind_min, swind_max, swind_delta, swind_project;
double *aaa;                    ///< A pointer to an array used by swind

struct Cdf cdf_ff;
struct Cdf cdf_fb;
struct Cdf cdf_vcos;
struct Cdf cdf_vdipole;
struct Cdf cdf_bb;
struct Cdf cdf_brem;

int itest, jtest;

char hubeny_list[132];          ///< Location of listing of files representing hubeny atmospheres

struct xbands xband;

double kap_bf[NLEVELS];

FILE *pstatptr;                 ///<  pointer to a diagnostic file that will contain photon data for given cells
int cell_phot_stats;            ///< 1=do  it, 0=dont do it
int ncstat;                     ///<  the actual number we are going to log
int ncell_stats[NCSTAT];        ///< the numbers of the cells we are going to log

int nerr_no_Jmodel;
int nerr_Jmodel_wrong_freq;

struct advanced_modes modes;

FILE *optr;                     ///< pointer to a diagnostic file that will contain dvds information

struct filenames files;

int xxxbound;

struct rdpar_choices zz_spec;

struct Import *imported_model;  ///<  MAX_DOM is defined in sirocco.h and as such import.h has to be included after


double velocity_electron[3];    // velocity of the electron when thermal effects are included
