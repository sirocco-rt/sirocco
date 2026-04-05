
/***********************************************************/
/** @file  windstruct2fits.c
 * @author ksl
 * @date   April, 2026
 *
 * @brief  Write wind struct (wmain) for all cells — including
 * lower-hemisphere transport cells — to a FITS binary table.
 *
 * ### Notes ###
 * windsave2fits only covers the plasma cells (upper hemisphere).
 * This program covers all 2*NDIM2 wind cells so that quantities
 * stored per wind cell (position, velocity, velocity gradient,
 * cell boundaries, dvds) can be inspected and compared between
 * the two hemispheres.
 *
 * Usage:
 *   windstruct2fits root
 *
 * Reads root.wind_save and writes root_windstruct.fits
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"
#include "fitsio.h"


/**********************************************************/
/**
 * @brief  Write one binary-table extension with one row per
 *         wind cell (all 2*NDIM2 cells).
 *
 * Columns written:
 *   nwind      J   - cell index in wmain
 *   ndom       J   - domain number
 *   nplasma    J   - plasma cell index (upper+lower share same)
 *   inwind     J   - W_ALL_INWIND / W_NOT_INWIND / ...
 *   x          E   - cylindrical R position of inner vertex (cm)
 *   z          E   - cylindrical z position of inner vertex (cm)
 *   xcen       E   - cylindrical R position of cell centre (cm)
 *   zcen       E   - cylindrical z position of cell centre (cm)
 *   i          J   - radial grid index
 *   j          J   - axial grid index
 *   x_3d      3E   - full 3-vector inner vertex position (cm)
 *   xcen_3d   3E   - full 3-vector cell centre position (cm)
 *   v         3E   - velocity at inner vertex (cm/s)
 *   v_grad    9E   - velocity gradient tensor (row-major, 1/s)
 *   cell_z_min E   - signed lower z boundary (cm)
 *   cell_z_max E   - signed upper z boundary (cm)
 *   dvds_ave   E   - direction-averaged |dv/ds| (1/s)
 *   dvds_max   E   - maximum |dv/ds| at vertex (1/s)
 *   nscat      J   - scatter count (hemisphere-specific)
 *   nrad       J   - emission count (hemisphere-specific)
 *   ntot       K   - wind_counts: total photon passages
 *   ntot_star  K   - wind_counts: passages from central star
 *   ntot_bl    K   - wind_counts: passages from boundary layer
 *   ntot_disk  K   - wind_counts: passages from disk
 *   ntot_wind  K   - wind_counts: passages from wind
 *   ntot_agn   K   - wind_counts: passages from AGN
 *   nioniz     K   - wind_counts: H-ionizing photon passages
 *
 **********************************************************/

int
write_wind_table (fitsfile *fptr)
{
  int status = 0;
  int num_rows, num_cols;
  int n;
  WindPtr cell;

  num_rows = 2 * NDIM2;
  num_cols = 27;

  /* Column definitions */
  char *ttype[] = {
    "nwind", "ndom", "nplasma", "inwind",
    "x", "z", "xcen", "zcen", "i", "j",
    "x_3d", "xcen_3d", "v", "v_grad",
    "cell_z_min", "cell_z_max",
    "dvds_ave", "dvds_max",
    "nscat", "nrad",
    "ntot", "ntot_star", "ntot_bl", "ntot_disk", "ntot_wind", "ntot_agn",
    "nioniz"
  };
  char *tform[] = {
    "J", "J", "J", "J",
    "E", "E", "E", "E", "J", "J",
    "3E", "3E", "3E", "9E",
    "E", "E",
    "E", "E",
    "J", "J",
    "K", "K", "K", "K", "K", "K",
    "K"
  };
  char *tunit[] = {
    "", "", "", "",
    "cm", "cm", "cm", "cm", "", "",
    "cm", "cm", "cm/s", "1/s",
    "cm", "cm",
    "1/s", "1/s",
    "", "",
    "", "", "", "", "", "",
    ""
  };

  if (fits_create_tbl (fptr, BINARY_TBL, num_rows, num_cols, ttype, tform, tunit, "windstruct", &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  /* Allocate column buffers */
  int *col_nwind = calloc (num_rows, sizeof (int));
  int *col_ndom = calloc (num_rows, sizeof (int));
  int *col_nplasma = calloc (num_rows, sizeof (int));
  int *col_inwind = calloc (num_rows, sizeof (int));
  float *col_x = calloc (num_rows, sizeof (float));
  float *col_z = calloc (num_rows, sizeof (float));
  float *col_xcen = calloc (num_rows, sizeof (float));
  float *col_zcen = calloc (num_rows, sizeof (float));
  int *col_i = calloc (num_rows, sizeof (int));
  int *col_j = calloc (num_rows, sizeof (int));
  float *col_x_3d = calloc (num_rows * 3, sizeof (float));
  float *col_xcen_3d = calloc (num_rows * 3, sizeof (float));
  float *col_v = calloc (num_rows * 3, sizeof (float));
  float *col_v_grad = calloc (num_rows * 9, sizeof (float));
  float *col_cell_z_min = calloc (num_rows, sizeof (float));
  float *col_cell_z_max = calloc (num_rows, sizeof (float));
  float *col_dvds_ave = calloc (num_rows, sizeof (float));
  float *col_dvds_max = calloc (num_rows, sizeof (float));
  int *col_nscat = calloc (num_rows, sizeof (int));
  int *col_nrad = calloc (num_rows, sizeof (int));
  long long *col_wc_ntot = calloc (num_rows, sizeof (long long));
  long long *col_wc_ntot_star = calloc (num_rows, sizeof (long long));
  long long *col_wc_ntot_bl = calloc (num_rows, sizeof (long long));
  long long *col_wc_ntot_disk = calloc (num_rows, sizeof (long long));
  long long *col_wc_ntot_wind = calloc (num_rows, sizeof (long long));
  long long *col_wc_ntot_agn = calloc (num_rows, sizeof (long long));
  long long *col_wc_nioniz = calloc (num_rows, sizeof (long long));

  if (!col_nwind || !col_ndom || !col_nplasma || !col_inwind ||
      !col_x || !col_z || !col_xcen || !col_zcen || !col_i || !col_j ||
      !col_x_3d || !col_xcen_3d || !col_v || !col_v_grad ||
      !col_cell_z_min || !col_cell_z_max || !col_dvds_ave || !col_dvds_max ||
      !col_nscat || !col_nrad ||
      !col_wc_ntot || !col_wc_ntot_star || !col_wc_ntot_bl || !col_wc_ntot_disk || !col_wc_ntot_wind || !col_wc_ntot_agn || !col_wc_nioniz)
  {
    fprintf (stderr, "windstruct2fits: memory allocation failed\n");
    return -1;
  }

  for (n = 0; n < num_rows; n++)
  {
    int ii, jj;
    cell = &wmain[n];

    col_nwind[n] = cell->nwind;
    col_ndom[n] = cell->ndom;
    col_nplasma[n] = cell->nplasma;
    col_inwind[n] = (int) cell->inwind;

    /* Scalar cylindrical coordinates (matching windsave2table convention) */
    col_x[n] = (float) cell->x[0];
    col_z[n] = (float) cell->x[2];
    col_xcen[n] = (float) cell->xcen[0];
    col_zcen[n] = (float) cell->xcen[2];

    /* Grid indices */
    wind_n_to_ij (cell->ndom, n, &ii, &jj);
    col_i[n] = ii;
    col_j[n] = jj;

    /* Full 3-vector positions */
    col_x_3d[n * 3 + 0] = (float) cell->x[0];
    col_x_3d[n * 3 + 1] = (float) cell->x[1];
    col_x_3d[n * 3 + 2] = (float) cell->x[2];

    col_xcen_3d[n * 3 + 0] = (float) cell->xcen[0];
    col_xcen_3d[n * 3 + 1] = (float) cell->xcen[1];
    col_xcen_3d[n * 3 + 2] = (float) cell->xcen[2];

    col_v[n * 3 + 0] = (float) cell->v[0];
    col_v[n * 3 + 1] = (float) cell->v[1];
    col_v[n * 3 + 2] = (float) cell->v[2];

    col_v_grad[n * 9 + 0] = (float) cell->v_grad[0][0];
    col_v_grad[n * 9 + 1] = (float) cell->v_grad[0][1];
    col_v_grad[n * 9 + 2] = (float) cell->v_grad[0][2];
    col_v_grad[n * 9 + 3] = (float) cell->v_grad[1][0];
    col_v_grad[n * 9 + 4] = (float) cell->v_grad[1][1];
    col_v_grad[n * 9 + 5] = (float) cell->v_grad[1][2];
    col_v_grad[n * 9 + 6] = (float) cell->v_grad[2][0];
    col_v_grad[n * 9 + 7] = (float) cell->v_grad[2][1];
    col_v_grad[n * 9 + 8] = (float) cell->v_grad[2][2];

    col_cell_z_min[n] = (float) cell->cell_z_min;
    col_cell_z_max[n] = (float) cell->cell_z_max;
    col_dvds_ave[n] = (float) cell->dvds_ave;
    col_dvds_max[n] = (float) cell->dvds_max;

    /* Scatter/emission counts: hemisphere-specific from plasma_derived */
    if (n < NDIM2)
    {
      col_nscat[n] = plasmamain[cell->nplasma].derived.nscat_upper;
      col_nrad[n] = plasmamain[cell->nplasma].derived.nrad_upper;
    }
    else
    {
      col_nscat[n] = plasmamain[cell->nplasma].derived.nscat_lower;
      col_nrad[n] = plasmamain[cell->nplasma].derived.nrad_lower;
    }

    /* Per-cell photon passage counters from wind_counts_main */
    col_wc_ntot[n] = (long long) wind_counts_main[n].ntot;
    col_wc_ntot_star[n] = (long long) wind_counts_main[n].ntot_star;
    col_wc_ntot_bl[n] = (long long) wind_counts_main[n].ntot_bl;
    col_wc_ntot_disk[n] = (long long) wind_counts_main[n].ntot_disk;
    col_wc_ntot_wind[n] = (long long) wind_counts_main[n].ntot_wind;
    col_wc_ntot_agn[n] = (long long) wind_counts_main[n].ntot_agn;
    col_wc_nioniz[n] = (long long) wind_counts_main[n].nioniz;
  }

  /* Write columns */
  fits_write_col (fptr, TINT, 1, 1, 1, num_rows, col_nwind, &status);
  fits_write_col (fptr, TINT, 2, 1, 1, num_rows, col_ndom, &status);
  fits_write_col (fptr, TINT, 3, 1, 1, num_rows, col_nplasma, &status);
  fits_write_col (fptr, TINT, 4, 1, 1, num_rows, col_inwind, &status);
  fits_write_col (fptr, TFLOAT, 5, 1, 1, num_rows, col_x, &status);
  fits_write_col (fptr, TFLOAT, 6, 1, 1, num_rows, col_z, &status);
  fits_write_col (fptr, TFLOAT, 7, 1, 1, num_rows, col_xcen, &status);
  fits_write_col (fptr, TFLOAT, 8, 1, 1, num_rows, col_zcen, &status);
  fits_write_col (fptr, TINT, 9, 1, 1, num_rows, col_i, &status);
  fits_write_col (fptr, TINT, 10, 1, 1, num_rows, col_j, &status);
  fits_write_col (fptr, TFLOAT, 11, 1, 1, num_rows * 3, col_x_3d, &status);
  fits_write_col (fptr, TFLOAT, 12, 1, 1, num_rows * 3, col_xcen_3d, &status);
  fits_write_col (fptr, TFLOAT, 13, 1, 1, num_rows * 3, col_v, &status);
  fits_write_col (fptr, TFLOAT, 14, 1, 1, num_rows * 9, col_v_grad, &status);
  fits_write_col (fptr, TFLOAT, 15, 1, 1, num_rows, col_cell_z_min, &status);
  fits_write_col (fptr, TFLOAT, 16, 1, 1, num_rows, col_cell_z_max, &status);
  fits_write_col (fptr, TFLOAT, 17, 1, 1, num_rows, col_dvds_ave, &status);
  fits_write_col (fptr, TFLOAT, 18, 1, 1, num_rows, col_dvds_max, &status);
  fits_write_col (fptr, TINT, 19, 1, 1, num_rows, col_nscat, &status);
  fits_write_col (fptr, TINT, 20, 1, 1, num_rows, col_nrad, &status);
  fits_write_col (fptr, TLONGLONG, 21, 1, 1, num_rows, col_wc_ntot, &status);
  fits_write_col (fptr, TLONGLONG, 22, 1, 1, num_rows, col_wc_ntot_star, &status);
  fits_write_col (fptr, TLONGLONG, 23, 1, 1, num_rows, col_wc_ntot_bl, &status);
  fits_write_col (fptr, TLONGLONG, 24, 1, 1, num_rows, col_wc_ntot_disk, &status);
  fits_write_col (fptr, TLONGLONG, 25, 1, 1, num_rows, col_wc_ntot_wind, &status);
  fits_write_col (fptr, TLONGLONG, 26, 1, 1, num_rows, col_wc_ntot_agn, &status);
  fits_write_col (fptr, TLONGLONG, 27, 1, 1, num_rows, col_wc_nioniz, &status);

  if (status)
    fits_report_error (stderr, status);

  free (col_nwind);
  free (col_ndom);
  free (col_nplasma);
  free (col_inwind);
  free (col_x);
  free (col_z);
  free (col_xcen);
  free (col_zcen);
  free (col_i);
  free (col_j);
  free (col_x_3d);
  free (col_xcen_3d);
  free (col_v);
  free (col_v_grad);
  free (col_cell_z_min);
  free (col_cell_z_max);
  free (col_dvds_ave);
  free (col_dvds_max);
  free (col_nscat);
  free (col_nrad);
  free (col_wc_ntot);
  free (col_wc_ntot_star);
  free (col_wc_ntot_bl);
  free (col_wc_ntot_disk);
  free (col_wc_ntot_wind);
  free (col_wc_ntot_agn);
  free (col_wc_nioniz);

  return status;
}


/**********************************************************/
/**
 * @brief  Write one binary-table extension with one row per plasma cell.
 *
 * Columns written (scalar unless noted):
 *
 *   Cross-references:
 *     nwind      J   - upper-hemisphere wind cell for this plasma cell
 *     nplasma    J   - self index
 *
 *   Grid location (matching windsave2table convention):
 *     x          E   - cylindrical R of inner vertex (cm)
 *     z          E   - cylindrical z of inner vertex (cm)
 *     xcen       E   - cylindrical R of cell centre (cm)
 *     zcen       E   - cylindrical z of cell centre (cm)
 *     i          J   - radial grid index
 *     jj         J   - axial grid index (named jj to avoid conflict with j mean intensity)
 *
 *   State (plasma_state):
 *     ne         E   - electron number density (cm^-3)
 *     rho        E   - mass density (g cm^-3)
 *     vol        E   - plasma cell volume (cm^3)
 *     t_r        E   - radiation temperature (K)
 *     t_e        E   - electron temperature (K)
 *     t_r_old    E   - radiation temperature previous cycle (K)
 *     t_e_old    E   - electron temperature previous cycle (K)
 *     w          E   - radiation dilution factor
 *
 *   Estimators (plasma_estimators):
 *     j          E   - mean intensity (erg s^-1 cm^-2 sr^-1)
 *     j_direct   E   - mean intensity, direct photons
 *     j_scatt    E   - mean intensity, scattered photons
 *     ave_freq   E   - intensity-averaged frequency (Hz)
 *     ip         E   - ionization parameter (Cloudy definition)
 *     heat_tot   E   - total heating rate (erg s^-1)
 *     heat_lines E
 *     heat_ff    E
 *     heat_photo E
 *     heat_auger E
 *     heat_comp  E
 *     cool_tot   E   - total cooling rate (erg s^-1)
 *     ntot       J   - total photon passages
 *     nioniz     J   - photon passages capable of ionizing H
 *     xj         20E - band-limited mean intensities
 *     xave_freq  20E - band-limited mean frequencies (Hz)
 *     rad_force_es 4E - radiative force, electron scattering (dyne)
 *     rad_force_bf 4E - radiative force, bound-free
 *     rad_force_ff 4E - radiative force, free-free
 *     F_vis      4E  - visible-band flux (erg s^-1 cm^-2)
 *     F_UV       4E  - UV-band flux
 *     F_Xray     4E  - X-ray-band flux
 *
 *   Derived (plasma_derived):
 *     lum_tot    E   - total radiative luminosity (erg s^-1)
 *     lum_lines  E
 *     lum_ff     E
 *     cool_rr    E   - radiative recombination cooling
 *     cool_comp  E   - Compton cooling
 *     cool_di    E   - direct-ionization cooling
 *     cool_dr    E   - dielectronic-recombination cooling
 *     dt_e       E   - change in t_e last iteration (K)
 *     gain       E   - iteration gain
 *     xi         E   - ionization parameter (TTS definition)
 *     converge_t_r E
 *     converge_t_e E
 *     converge_hc  E
 *     converge_whole J - sum of convergence check flags
 *     converging   J
 *     nscat_es   J   - electron scatters
 *     nscat_res  J   - resonant line scatters
 *     nscat_upper J  - scatters in upper-hemisphere transport cells
 *     nscat_lower J  - scatters in lower-hemisphere transport cells
 *
 **********************************************************/

int
write_plasma_table (fitsfile *fptr)
{
  int status = 0;
  int num_rows, num_cols;
  int n;
  PlasmaPtr cell;

  num_rows = NPLASMA;
  num_cols = 60;

  char *ttype[] = {
    "nwind", "nplasma",
    "x", "z", "xcen", "zcen", "i", "jj",
    "ne", "rho", "vol", "t_r", "t_e", "t_r_old", "t_e_old", "w",
    "j", "j_direct", "j_scatt", "ave_freq", "ip",
    "heat_tot", "heat_lines", "heat_ff", "heat_photo", "heat_auger", "heat_comp",
    "cool_tot", "ntot", "nioniz",
    "xj", "xave_freq",
    "rad_force_es", "rad_force_bf", "rad_force_ff",
    "F_vis", "F_UV", "F_Xray",
    "lum_tot", "lum_lines", "lum_ff", "cool_rr", "cool_comp", "cool_di", "cool_dr",
    "dt_e", "gain", "xi",
    "converge_t_r", "converge_t_e", "converge_hc",
    "converge_whole", "converging",
    "nscat_es", "nscat_res", "nscat_upper", "nscat_lower",
    "nrad", "nrad_upper", "nrad_lower"
  };

  char *tform[] = {
    "J", "J",
    "E", "E", "E", "E", "J", "J",
    "E", "E", "E", "E", "E", "E", "E", "E",
    "E", "E", "E", "E", "E",
    "E", "E", "E", "E", "E", "E",
    "E", "J", "J",
    "20E", "20E",
    "4E", "4E", "4E",
    "4E", "4E", "4E",
    "E", "E", "E", "E", "E", "E", "E",
    "E", "E", "E",
    "E", "E", "E",
    "J", "J",
    "J", "J", "J", "J",
    "J", "J", "J"
  };

  char *tunit[] = {
    "", "",
    "cm", "cm", "cm", "cm", "", "",
    "cm^-3", "g cm^-3", "cm^3", "K", "K", "K", "K", "",
    "erg s^-1 cm^-2 sr^-1", "erg s^-1 cm^-2 sr^-1", "erg s^-1 cm^-2 sr^-1", "Hz", "",
    "erg s^-1", "erg s^-1", "erg s^-1", "erg s^-1", "erg s^-1", "erg s^-1",
    "erg s^-1", "", "",
    "erg s^-1 cm^-2 sr^-1", "Hz",
    "dyne", "dyne", "dyne",
    "erg s^-1 cm^-2", "erg s^-1 cm^-2", "erg s^-1 cm^-2",
    "erg s^-1", "erg s^-1", "erg s^-1", "erg s^-1", "erg s^-1", "erg s^-1", "erg s^-1",
    "K", "", "",
    "", "", "",
    "", "",
    "", "", "", "",
    "", "", ""
  };

  if (fits_create_tbl (fptr, BINARY_TBL, num_rows, num_cols, ttype, tform, tunit, "plasma", &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  /* Allocate column buffers */
  int *col_nwind = calloc (num_rows, sizeof (int));
  int *col_nplasma = calloc (num_rows, sizeof (int));
  float *col_p_x = calloc (num_rows, sizeof (float));
  float *col_p_z = calloc (num_rows, sizeof (float));
  float *col_p_xcen = calloc (num_rows, sizeof (float));
  float *col_p_zcen = calloc (num_rows, sizeof (float));
  int *col_p_i = calloc (num_rows, sizeof (int));
  int *col_p_jj = calloc (num_rows, sizeof (int));
  float *col_ne = calloc (num_rows, sizeof (float));
  float *col_rho = calloc (num_rows, sizeof (float));
  float *col_vol = calloc (num_rows, sizeof (float));
  float *col_t_r = calloc (num_rows, sizeof (float));
  float *col_t_e = calloc (num_rows, sizeof (float));
  float *col_t_r_old = calloc (num_rows, sizeof (float));
  float *col_t_e_old = calloc (num_rows, sizeof (float));
  float *col_w = calloc (num_rows, sizeof (float));
  float *col_j = calloc (num_rows, sizeof (float));
  float *col_j_direct = calloc (num_rows, sizeof (float));
  float *col_j_scatt = calloc (num_rows, sizeof (float));
  float *col_ave_freq = calloc (num_rows, sizeof (float));
  float *col_ip = calloc (num_rows, sizeof (float));
  float *col_heat_tot = calloc (num_rows, sizeof (float));
  float *col_heat_lines = calloc (num_rows, sizeof (float));
  float *col_heat_ff = calloc (num_rows, sizeof (float));
  float *col_heat_photo = calloc (num_rows, sizeof (float));
  float *col_heat_auger = calloc (num_rows, sizeof (float));
  float *col_heat_comp = calloc (num_rows, sizeof (float));
  float *col_cool_tot = calloc (num_rows, sizeof (float));
  int *col_ntot = calloc (num_rows, sizeof (int));
  int *col_nioniz = calloc (num_rows, sizeof (int));
  float *col_xj = calloc (num_rows * NXBANDS, sizeof (float));
  float *col_xave_freq = calloc (num_rows * NXBANDS, sizeof (float));
  float *col_rad_force_es = calloc (num_rows * NFORCE_DIRECTIONS, sizeof (float));
  float *col_rad_force_bf = calloc (num_rows * NFORCE_DIRECTIONS, sizeof (float));
  float *col_rad_force_ff = calloc (num_rows * NFORCE_DIRECTIONS, sizeof (float));
  float *col_F_vis = calloc (num_rows * NFORCE_DIRECTIONS, sizeof (float));
  float *col_F_UV = calloc (num_rows * NFORCE_DIRECTIONS, sizeof (float));
  float *col_F_Xray = calloc (num_rows * NFORCE_DIRECTIONS, sizeof (float));
  float *col_lum_tot = calloc (num_rows, sizeof (float));
  float *col_lum_lines = calloc (num_rows, sizeof (float));
  float *col_lum_ff = calloc (num_rows, sizeof (float));
  float *col_cool_rr = calloc (num_rows, sizeof (float));
  float *col_cool_comp = calloc (num_rows, sizeof (float));
  float *col_cool_di = calloc (num_rows, sizeof (float));
  float *col_cool_dr = calloc (num_rows, sizeof (float));
  float *col_dt_e = calloc (num_rows, sizeof (float));
  float *col_gain = calloc (num_rows, sizeof (float));
  float *col_xi = calloc (num_rows, sizeof (float));
  float *col_converge_t_r = calloc (num_rows, sizeof (float));
  float *col_converge_t_e = calloc (num_rows, sizeof (float));
  float *col_converge_hc = calloc (num_rows, sizeof (float));
  int *col_converge_whole = calloc (num_rows, sizeof (int));
  int *col_converging = calloc (num_rows, sizeof (int));
  int *col_nscat_es = calloc (num_rows, sizeof (int));
  int *col_nscat_res = calloc (num_rows, sizeof (int));
  int *col_nscat_upper = calloc (num_rows, sizeof (int));
  int *col_nscat_lower = calloc (num_rows, sizeof (int));
  int *col_nrad = calloc (num_rows, sizeof (int));
  int *col_nrad_upper = calloc (num_rows, sizeof (int));
  int *col_nrad_lower = calloc (num_rows, sizeof (int));

  if (!col_nwind || !col_nplasma ||
      !col_p_x || !col_p_z || !col_p_xcen || !col_p_zcen || !col_p_i || !col_p_jj ||
      !col_ne || !col_rho || !col_vol ||
      !col_t_r || !col_t_e || !col_t_r_old || !col_t_e_old || !col_w ||
      !col_j || !col_j_direct || !col_j_scatt || !col_ave_freq || !col_ip ||
      !col_heat_tot || !col_heat_lines || !col_heat_ff || !col_heat_photo ||
      !col_heat_auger || !col_heat_comp || !col_cool_tot || !col_ntot || !col_nioniz ||
      !col_xj || !col_xave_freq ||
      !col_rad_force_es || !col_rad_force_bf || !col_rad_force_ff ||
      !col_F_vis || !col_F_UV || !col_F_Xray ||
      !col_lum_tot || !col_lum_lines || !col_lum_ff ||
      !col_cool_rr || !col_cool_comp || !col_cool_di || !col_cool_dr ||
      !col_dt_e || !col_gain || !col_xi ||
      !col_converge_t_r || !col_converge_t_e || !col_converge_hc ||
      !col_converge_whole || !col_converging || !col_nscat_es || !col_nscat_res || !col_nscat_upper || !col_nscat_lower
      || !col_nrad || !col_nrad_upper || !col_nrad_lower)
  {
    fprintf (stderr, "write_plasma_table: memory allocation failed\n");
    return -1;
  }

  for (n = 0; n < num_rows; n++)
  {
    int k, ii, jj;
    WindPtr wcell;
    cell = &plasmamain[n];

    col_nwind[n] = cell->nwind;
    col_nplasma[n] = cell->nplasma;

    /* Grid location from the corresponding upper-hemisphere wind cell */
    wcell = &wmain[cell->nwind];
    col_p_x[n] = (float) wcell->x[0];
    col_p_z[n] = (float) wcell->x[2];
    col_p_xcen[n] = (float) wcell->xcen[0];
    col_p_zcen[n] = (float) wcell->xcen[2];
    wind_n_to_ij (wcell->ndom, cell->nwind, &ii, &jj);
    col_p_i[n] = ii;
    col_p_jj[n] = jj;

    col_ne[n] = (float) cell->state.ne;
    col_rho[n] = (float) cell->state.rho;
    col_vol[n] = (float) cell->state.vol;
    col_t_r[n] = (float) cell->state.t_r;
    col_t_e[n] = (float) cell->state.t_e;
    col_t_r_old[n] = (float) cell->state.t_r_old;
    col_t_e_old[n] = (float) cell->state.t_e_old;
    col_w[n] = (float) cell->state.w;

    col_j[n] = (float) cell->est.j;
    col_j_direct[n] = (float) cell->est.j_direct;
    col_j_scatt[n] = (float) cell->est.j_scatt;
    col_ave_freq[n] = (float) cell->est.ave_freq;
    col_ip[n] = (float) cell->est.ip;
    col_heat_tot[n] = (float) cell->est.heat_tot;
    col_heat_lines[n] = (float) cell->est.heat_lines;
    col_heat_ff[n] = (float) cell->est.heat_ff;
    col_heat_photo[n] = (float) cell->est.heat_photo;
    col_heat_auger[n] = (float) cell->est.heat_auger;
    col_heat_comp[n] = (float) cell->est.heat_comp;
    col_cool_tot[n] = (float) cell->est.cool_tot;
    col_ntot[n] = cell->est.ntot;
    col_nioniz[n] = cell->est.nioniz;

    for (k = 0; k < NXBANDS; k++)
    {
      col_xj[n * NXBANDS + k] = (float) cell->est.xj[k];
      col_xave_freq[n * NXBANDS + k] = (float) cell->est.xave_freq[k];
    }
    for (k = 0; k < NFORCE_DIRECTIONS; k++)
    {
      col_rad_force_es[n * NFORCE_DIRECTIONS + k] = (float) cell->est.rad_force_es[k];
      col_rad_force_bf[n * NFORCE_DIRECTIONS + k] = (float) cell->est.rad_force_bf[k];
      col_rad_force_ff[n * NFORCE_DIRECTIONS + k] = (float) cell->est.rad_force_ff[k];
      col_F_vis[n * NFORCE_DIRECTIONS + k] = (float) cell->est.F_vis[k];
      col_F_UV[n * NFORCE_DIRECTIONS + k] = (float) cell->est.F_UV[k];
      col_F_Xray[n * NFORCE_DIRECTIONS + k] = (float) cell->est.F_Xray[k];
    }

    col_lum_tot[n] = (float) cell->derived.lum_tot;
    col_lum_lines[n] = (float) cell->derived.lum_lines;
    col_lum_ff[n] = (float) cell->derived.lum_ff;
    col_cool_rr[n] = (float) cell->derived.cool_rr;
    col_cool_comp[n] = (float) cell->derived.cool_comp;
    col_cool_di[n] = (float) cell->derived.cool_di;
    col_cool_dr[n] = (float) cell->derived.cool_dr;
    col_dt_e[n] = (float) cell->derived.dt_e;
    col_gain[n] = (float) cell->derived.gain;
    col_xi[n] = (float) cell->derived.xi;
    col_converge_t_r[n] = (float) cell->derived.converge_t_r;
    col_converge_t_e[n] = (float) cell->derived.converge_t_e;
    col_converge_hc[n] = (float) cell->derived.converge_hc;
    col_converge_whole[n] = cell->derived.converge_whole;
    col_converging[n] = cell->derived.converging;
    col_nscat_es[n] = cell->derived.nscat_es;
    col_nscat_res[n] = cell->derived.nscat_res;
    col_nscat_upper[n] = cell->derived.nscat_upper;
    col_nscat_lower[n] = cell->derived.nscat_lower;
    col_nrad[n] = cell->derived.nrad;
    col_nrad_upper[n] = cell->derived.nrad_upper;
    col_nrad_lower[n] = cell->derived.nrad_lower;
  }

  /* Write columns */
  int col = 1;
  fits_write_col (fptr, TINT, col++, 1, 1, num_rows, col_nwind, &status);
  fits_write_col (fptr, TINT, col++, 1, 1, num_rows, col_nplasma, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_p_x, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_p_z, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_p_xcen, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_p_zcen, &status);
  fits_write_col (fptr, TINT, col++, 1, 1, num_rows, col_p_i, &status);
  fits_write_col (fptr, TINT, col++, 1, 1, num_rows, col_p_jj, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_ne, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_rho, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_vol, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_t_r, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_t_e, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_t_r_old, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_t_e_old, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_w, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_j, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_j_direct, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_j_scatt, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_ave_freq, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_ip, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_heat_tot, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_heat_lines, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_heat_ff, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_heat_photo, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_heat_auger, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_heat_comp, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_cool_tot, &status);
  fits_write_col (fptr, TINT, col++, 1, 1, num_rows, col_ntot, &status);
  fits_write_col (fptr, TINT, col++, 1, 1, num_rows, col_nioniz, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows * NXBANDS, col_xj, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows * NXBANDS, col_xave_freq, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows * NFORCE_DIRECTIONS, col_rad_force_es, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows * NFORCE_DIRECTIONS, col_rad_force_bf, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows * NFORCE_DIRECTIONS, col_rad_force_ff, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows * NFORCE_DIRECTIONS, col_F_vis, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows * NFORCE_DIRECTIONS, col_F_UV, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows * NFORCE_DIRECTIONS, col_F_Xray, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_lum_tot, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_lum_lines, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_lum_ff, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_cool_rr, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_cool_comp, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_cool_di, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_cool_dr, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_dt_e, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_gain, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_xi, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_converge_t_r, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_converge_t_e, &status);
  fits_write_col (fptr, TFLOAT, col++, 1, 1, num_rows, col_converge_hc, &status);
  fits_write_col (fptr, TINT, col++, 1, 1, num_rows, col_converge_whole, &status);
  fits_write_col (fptr, TINT, col++, 1, 1, num_rows, col_converging, &status);
  fits_write_col (fptr, TINT, col++, 1, 1, num_rows, col_nscat_es, &status);
  fits_write_col (fptr, TINT, col++, 1, 1, num_rows, col_nscat_res, &status);
  fits_write_col (fptr, TINT, col++, 1, 1, num_rows, col_nscat_upper, &status);
  fits_write_col (fptr, TINT, col++, 1, 1, num_rows, col_nscat_lower, &status);
  fits_write_col (fptr, TINT, col++, 1, 1, num_rows, col_nrad, &status);
  fits_write_col (fptr, TINT, col++, 1, 1, num_rows, col_nrad_upper, &status);
  fits_write_col (fptr, TINT, col++, 1, 1, num_rows, col_nrad_lower, &status);

  if (status)
    fits_report_error (stderr, status);

  free (col_nwind);
  free (col_nplasma);
  free (col_p_x);
  free (col_p_z);
  free (col_p_xcen);
  free (col_p_zcen);
  free (col_p_i);
  free (col_p_jj);
  free (col_ne);
  free (col_rho);
  free (col_vol);
  free (col_t_r);
  free (col_t_e);
  free (col_t_r_old);
  free (col_t_e_old);
  free (col_w);
  free (col_j);
  free (col_j_direct);
  free (col_j_scatt);
  free (col_ave_freq);
  free (col_ip);
  free (col_heat_tot);
  free (col_heat_lines);
  free (col_heat_ff);
  free (col_heat_photo);
  free (col_heat_auger);
  free (col_heat_comp);
  free (col_cool_tot);
  free (col_ntot);
  free (col_nioniz);
  free (col_xj);
  free (col_xave_freq);
  free (col_rad_force_es);
  free (col_rad_force_bf);
  free (col_rad_force_ff);
  free (col_F_vis);
  free (col_F_UV);
  free (col_F_Xray);
  free (col_lum_tot);
  free (col_lum_lines);
  free (col_lum_ff);
  free (col_cool_rr);
  free (col_cool_comp);
  free (col_cool_di);
  free (col_cool_dr);
  free (col_dt_e);
  free (col_gain);
  free (col_xi);
  free (col_converge_t_r);
  free (col_converge_t_e);
  free (col_converge_hc);
  free (col_converge_whole);
  free (col_converging);
  free (col_nscat_es);
  free (col_nscat_res);
  free (col_nscat_upper);
  free (col_nscat_lower);
  free (col_nrad);
  free (col_nrad_upper);
  free (col_nrad_lower);

  return status;
}


/**********************************************************/
/**
 * @brief  main
 **********************************************************/

int
main (int argc, char *argv[])
{
  char inroot[LINELENGTH];
  char infile[LINELENGTH];
  char outfile[LINELENGTH];
  fitsfile *fptr;
  int status = 0;

  if (argc < 2)
  {
    fprintf (stderr, "Usage: windstruct2fits root\n");
    fprintf (stderr, "  Reads root.wind_save and writes root_windstruct.fits\n");
    return EXIT_FAILURE;
  }

  get_root (inroot, argv[argc - 1]);
  sprintf (infile, "%.150s.wind_save", inroot);
  sprintf (outfile, "!%.150s_windstruct.fits", inroot);

  printf ("Reading %s\n", infile);
  printf ("Writing %s\n", outfile + 1); /* skip the leading '!' in the message */

  zdom = calloc (MAX_DOM, sizeof (domain_dummy));
  if (zdom == NULL)
  {
    fprintf (stderr, "windstruct2fits: unable to allocate domain memory\n");
    return EXIT_FAILURE;
  }

  init_rand (1);                /* dvds calculations in make_transport_grid need the RNG */
  wind_read (infile);

  printf ("NDIM2 = %d  (total wind cells = %d)\n", NDIM2, 2 * NDIM2);

  if (fits_create_file (&fptr, outfile, &status))
  {
    fits_report_error (stderr, status);
    return EXIT_FAILURE;
  }

  /* Empty primary HDU required by FITS standard */
  if (fits_create_img (fptr, FLOAT_IMG, 0, NULL, &status))
  {
    fits_report_error (stderr, status);
  }

  write_wind_table (fptr);
  write_plasma_table (fptr);

  if (fits_close_file (fptr, &status))
    fits_report_error (stderr, status);

  return EXIT_SUCCESS;
}
