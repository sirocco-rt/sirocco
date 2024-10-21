/** ************************************************************************* */
/**
 * @file     py_optical_depth.h
 * @author   Edward Parkinson
 * @date     February 2021
 *
 * @details
 *
 * This header file contains the main constants, macros and types used in
 * py_optical_depth.
 *
 * ************************************************************************** */

#define MAX_CUSTOM_ANGLES 10
#define NUM_FREQUENCY_BINS 10000
#define NAMELEN 32

// Error message macro, adds the file name and line to the start of the error
// message

#define errormsg(fmt, ...)                           \
{                                                    \
  fprintf(stderr, "(%s:%i): ", __FILE__, __LINE__);  \
  fprintf(stderr, fmt, ##__VA_ARGS__);               \
}

/** Structure to hold the angles to extract the optical depth/column density from
  */

typedef struct SightLines_s
{
  char name[NAMELEN];
  double angle;
  double lmn[3];
} SightLines_t;

/** Structure to hold the name and frequency of a photoionization edge to
  *evaluate the optical depth at
  */

typedef struct Edges_s
{
  char name[50];
  double freq;
} Edges_t;

/** Structure to save the angle and position of the electron scattering
  * photosphere surface
  */

typedef struct Positions_s
{
  double angle;
  double x, y, z;
} Positions_t;

/** Enumerator used to control the column density which is extracted, i.e. by
  * default mass density/N_H is extracted by the density of an ion can also
  * be extracted
  */

enum COLUMN_DENSITY
{
  COLUMN_MODE_RHO,
  COLUMN_MODE_ION,
};

extern int COLUMN_MODE;
extern int COLUMN_MODE_ION_NUMBER;

// Control which domain to initially send photons from

extern int N_DOMAIN;

// Control how the optical depth integration is done. This can be done in
// integrated tau mode, which gets the integrated tau along the path. The other
// mode aims to find the surface of the electron scattering photosphere

typedef enum RunModeEnum
{
  RUN_MODE_TAU_INTEGRATE = 0,
  RUN_MODE_ES_PHOTOSPHERE = 1,
  RUN_MODE_NO_ES_OPACITY = 2,
} RunMode_t;

extern RunMode_t RUN_MODE;

extern double TAU_DEPTH;

// External functions from other files

int create_photon (PhotPtr p_out, double freq, double *lmn);
SightLines_t *initialize_inclination_angles (int *n_angles, double *input_inclinations);
int integrate_tau_across_wind (PhotPtr photon, double *c_column_density, double *c_optical_depth);
void print_optical_depths (SightLines_t * inclinations, int n_inclinations, Edges_t edges[], int n_edges, double *optical_depth,
                           double *column_density);
void write_optical_depth_spectrum (SightLines_t * inclinations, int n_inclinations, double *tau_spectrum, double freq_min, double d_freq);
void write_photosphere_location_to_file (Positions_t * positions, int n_angles);
