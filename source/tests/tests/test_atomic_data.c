/** ********************************************************************************************************************
 *
 *  @file test_atomic_data.c
 *  @date October 2026
 *
 *  @brief Unit tests for reading the atomic data: levels, lines, collision strengths, photoionization, inner shell
 *  and macro-atom Auger data
 *
 *  Each atomic data set in xdata/ is loaded, and every phot_top, inner_cross and auger_macro record is reduced to a
 *  64-bit digest which is compared with a reference value. For the cross sections the digest covers the integer
 *  fields of each record, its frequencies and cross sections, phot_freq_min, inner_freq_min and the frequency ordering
 *  in phot_top_ptr and inner_cross_ptr. For the Auger data it covers every field of each record and the iauger and
 *  nauger of each level. These values are exactly reproducible on any IEEE platform, since they come only from
 *  reading the files and from multiplication and division. The log versions of the cross sections depend on the
 *  maths library, so instead of being in the digest they are checked against the log of the linear values.
 *
 *  A second digest covers the levels (xconfig), lines (line, and their order in lin_ptr) and collision strengths
 *  (coll_stren): every field except log_g, which is instead checked against log (g), and the line power, which is
 *  set later.
 *
 *  The reference digests were made before these tables were sized to the data read (#1204 for the cross sections
 *  and Auger data, #1223 for the levels, lines and collision strengths), so the test shows that the changes left
 *  the data unaltered. The reference for fe_17to27.dat, which has more Auger records than the old fixed
 *  limit, was made with the old code compiled with LARGE_MACRO_ATOMS. If an atomic data file is changed, the digests
 *  must be remade.
 *
 * ****************************************************************************************************************** */

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <CUnit/CUnit.h>

#include "../../atomic.h"
#include "../../sirocco.h"
#include "../unit_test.h"

#define FNV_OFFSET 14695981039346656037ULL
#define FNV_PRIME 1099511628211ULL

static char ATOMIC_DATA_DEST[LINELENGTH];

static const struct
{
  const char *masterfile;
  int nphot_total;
  int n_inner_tot;
  int nauger_macro;
  uint64_t digest;
  int nlevels;
  int nlines;
  int n_coll_stren;
  uint64_t digest_levels_lines;
} REFERENCE[] = {
  {"standard80.dat", 374, 257, 0, 0x7d0e2417300c862eULL, 3328, 5459, 2294, 0xa466729af4920718ULL},
  {"standard80_reduced.dat", 232, 257, 0, 0x98f821552e352232ULL, 3328, 5459, 2294, 0xa52c2d0347c7ad8cULL},
  {"h20_hetop_standard80.dat", 290, 257, 0, 0x106c8ad1b47f569cULL, 3302, 5964, 2274, 0x9889297707270106ULL},
  {"master_cno.dat", 493, 242, 0, 0xfecb074b57433784ULL, 3027, 5660, 2573, 0x911032d9a6cee64aULL},
  {"fe_25fe_26_f10_fe_27_f1.dat", 234, 185, 0, 0x2bb6a1e785a19892ULL, 1785, 4275, 1745, 0x0708cee08e99baf5ULL},
  {"fe_17to27.dat", 452, 185, 1332, 0xdcc3c7467b157036ULL, 2172, 6699, 1745, 0x4e4eab0ca4c1479aULL},
};

/** *******************************************************************************************************************
 *
 * @brief Add a 64-bit value to an FNV-1a digest, least significant byte first
 *
 * ****************************************************************************************************************** */

static uint64_t
digest_add (uint64_t digest, uint64_t value)
{
  int i;

  for (i = 0; i < 8; i++)
  {
    digest ^= (value >> (8 * i)) & 0xff;
    digest *= FNV_PRIME;
  }
  return digest;
}

static uint64_t
digest_add_double (uint64_t digest, double value)
{
  uint64_t bits;

  memcpy (&bits, &value, sizeof (bits));
  return digest_add (digest, bits);
}

static uint64_t
digest_add_int (uint64_t digest, int value)
{
  return digest_add (digest, (uint64_t) (int64_t) value);
}

/** *******************************************************************************************************************
 *
 * @brief Add one cross-section record to a digest, and check its log arrays
 *
 * ****************************************************************************************************************** */

static uint64_t
digest_record (uint64_t digest, const Topbase_phot *xptr)
{
  int i;
  int nbad = 0;

  digest = digest_add_int (digest, xptr->nlev);
  digest = digest_add_int (digest, xptr->uplev);
  digest = digest_add_int (digest, xptr->nion);
  digest = digest_add_int (digest, xptr->z);
  digest = digest_add_int (digest, xptr->istate);
  digest = digest_add_int (digest, xptr->np);
  digest = digest_add_int (digest, xptr->n);
  digest = digest_add_int (digest, xptr->l);
  digest = digest_add_int (digest, xptr->nlast);
  digest = digest_add_int (digest, xptr->n_elec_yield);
  digest = digest_add_int (digest, xptr->macro_info);
  digest = digest_add_int (digest, xptr->down_index);
  digest = digest_add_int (digest, xptr->up_index);
  digest = digest_add_int (digest, xptr->use);

  for (i = 0; i < xptr->np; i++)
  {
    digest = digest_add_double (digest, xptr->freq[i]);
    digest = digest_add_double (digest, xptr->x[i]);
    if (xptr->log_freq[i] != log (xptr->freq[i]) || xptr->log_x[i] != log (xptr->x[i]))
      nbad++;
  }
  CU_ASSERT_EQUAL (nbad, 0);

  return digest;
}

/** *******************************************************************************************************************
 *
 * @brief Add one macro-atom Auger record to a digest
 *
 * ****************************************************************************************************************** */

static uint64_t
digest_auger (uint64_t digest, const auger_dummy *aptr)
{
  int i;

  digest = digest_add_int (digest, aptr->nion);
  digest = digest_add_int (digest, aptr->z);
  digest = digest_add_int (digest, aptr->istate);
  digest = digest_add_int (digest, aptr->nconfig);
  digest = digest_add_int (digest, aptr->iauger);
  digest = digest_add_int (digest, aptr->nauger);
  digest = digest_add_double (digest, aptr->Avalue_auger);
  for (i = 0; i < NAUGER_ELECTRONS; i++)
  {
    digest = digest_add_int (digest, aptr->nconfig_target[i]);
    digest = digest_add_double (digest, aptr->branching_ratio[i]);
  }
  return digest;
}

/** *******************************************************************************************************************
 *
 * @brief Add one level to a digest, and check its log_g
 *
 * ****************************************************************************************************************** */

static uint64_t
digest_level (uint64_t digest, const config_dummy *cptr)
{
  int i;

  digest = digest_add_int (digest, cptr->z);
  digest = digest_add_int (digest, cptr->istate);
  digest = digest_add_int (digest, cptr->nion);
  digest = digest_add_int (digest, cptr->nden);
  digest = digest_add_int (digest, cptr->isp);
  digest = digest_add_int (digest, cptr->ilv);
  digest = digest_add_int (digest, cptr->macro_info);
  digest = digest_add_int (digest, cptr->n_bbu_jump);
  digest = digest_add_int (digest, cptr->n_bbd_jump);
  digest = digest_add_int (digest, cptr->n_bfu_jump);
  digest = digest_add_int (digest, cptr->n_bfd_jump);
  digest = digest_add_double (digest, cptr->g);
  digest = digest_add_double (digest, cptr->q_num);
  digest = digest_add_double (digest, cptr->ex);
  digest = digest_add_double (digest, cptr->rad_rate);
  for (i = 0; i < cptr->n_bbu_jump; i++)
    digest = digest_add_int (digest, cptr->bbu_jump[i]);
  for (i = 0; i < cptr->n_bbd_jump; i++)
    digest = digest_add_int (digest, cptr->bbd_jump[i]);
  for (i = 0; i < cptr->n_bfu_jump; i++)
    digest = digest_add_int (digest, cptr->bfu_jump[i]);
  for (i = 0; i < cptr->n_bfd_jump; i++)
    digest = digest_add_int (digest, cptr->bfd_jump[i]);
  digest = digest_add_int (digest, cptr->bbu_indx_first);
  digest = digest_add_int (digest, cptr->bfu_indx_first);
  digest = digest_add_int (digest, cptr->bfd_indx_first);
  digest = digest_add_int (digest, cptr->iauger);
  digest = digest_add_int (digest, cptr->nauger);

  if (cptr->g > 0)
    CU_ASSERT (cptr->log_g == 0.0 || cptr->log_g == log (cptr->g));

  return digest;
}

/** *******************************************************************************************************************
 *
 * @brief Add one line to a digest
 *
 * ****************************************************************************************************************** */

static uint64_t
digest_line (uint64_t digest, const line_dummy *lptr)
{
  digest = digest_add_int (digest, lptr->nion);
  digest = digest_add_int (digest, lptr->z);
  digest = digest_add_int (digest, lptr->istate);
  digest = digest_add_double (digest, lptr->gl);
  digest = digest_add_double (digest, lptr->gu);
  digest = digest_add_int (digest, lptr->nconfigl);
  digest = digest_add_int (digest, lptr->nconfigu);
  digest = digest_add_int (digest, lptr->levl);
  digest = digest_add_int (digest, lptr->levu);
  digest = digest_add_int (digest, lptr->macro_info);
  digest = digest_add_double (digest, lptr->freq);
  digest = digest_add_double (digest, lptr->f);
  digest = digest_add_double (digest, lptr->el);
  digest = digest_add_double (digest, lptr->eu);
  digest = digest_add_int (digest, lptr->where_in_list);
  digest = digest_add_int (digest, lptr->down_index);
  digest = digest_add_int (digest, lptr->up_index);
  digest = digest_add_int (digest, lptr->coll_index);
  return digest;
}

/** *******************************************************************************************************************
 *
 * @brief Add one collision-strength record to a digest
 *
 * ****************************************************************************************************************** */

static uint64_t
digest_coll_stren (uint64_t digest, const Coll_stren *cptr)
{
  int i;

  digest = digest_add_int (digest, cptr->n);
  digest = digest_add_int (digest, cptr->lower);
  digest = digest_add_int (digest, cptr->upper);
  digest = digest_add_double (digest, cptr->energy);
  digest = digest_add_double (digest, cptr->gf);
  digest = digest_add_double (digest, cptr->hi_t_lim);
  digest = digest_add_int (digest, cptr->n_points);
  digest = digest_add_int (digest, cptr->type);
  digest = digest_add_double (digest, (double) cptr->scaling_param);
  for (i = 0; i < N_COLL_STREN_PTS; i++)
  {
    digest = digest_add_double (digest, cptr->sct[i]);
    digest = digest_add_double (digest, cptr->scups[i]);
  }
  return digest;
}

/** *******************************************************************************************************************
 *
 * @brief Load one atomic data set and compare its cross sections with the reference
 *
 * ****************************************************************************************************************** */

static void
check_atomic_data_set (int k)
{
  char masterfile[LINELENGTH];
  uint64_t digest = FNV_OFFSET;
  uint64_t digest2 = FNV_OFFSET;
  int n;

  snprintf (masterfile, LINELENGTH, "%s/%s", ATOMIC_DATA_DEST, REFERENCE[k].masterfile);
  get_atomic_data (masterfile);

  CU_ASSERT_EQUAL_FATAL (nphot_total, REFERENCE[k].nphot_total);
  CU_ASSERT_EQUAL_FATAL (n_inner_tot, REFERENCE[k].n_inner_tot);
  CU_ASSERT_EQUAL_FATAL (nauger_macro, REFERENCE[k].nauger_macro);

  digest = digest_add_double (digest, phot_freq_min);
  digest = digest_add_double (digest, inner_freq_min);
  for (n = 0; n < nphot_total; n++)
  {
    digest = digest_record (digest, &phot_top[n]);
    digest = digest_add_int (digest, (int) (phot_top_ptr[n] - phot_top));
  }
  for (n = 0; n < n_inner_tot; n++)
  {
    digest = digest_record (digest, &inner_cross[n]);
    digest = digest_add_int (digest, (int) (inner_cross_ptr[n] - inner_cross));
  }
  for (n = 0; n < nauger_macro; n++)
    digest = digest_auger (digest, &auger_macro[n]);
  for (n = 0; n < nlevels; n++)
  {
    digest = digest_add_int (digest, xconfig[n].iauger);
    digest = digest_add_int (digest, xconfig[n].nauger);
  }

  if (digest != REFERENCE[k].digest)
    fprintf (stderr, "\n%s: atomic data digest 0x%016llxULL, expected 0x%016llxULL\n", REFERENCE[k].masterfile,
             (unsigned long long) digest, (unsigned long long) REFERENCE[k].digest);
  CU_ASSERT_EQUAL (digest, REFERENCE[k].digest);

  for (n = 0; n < nlevels; n++)
    digest2 = digest_level (digest2, &xconfig[n]);
  for (n = 0; n < nlines; n++)
  {
    digest2 = digest_line (digest2, &line[n]);
    digest2 = digest_add_int (digest2, (int) (lin_ptr[n] - line));
  }
  for (n = 0; n < n_coll_stren; n++)
    digest2 = digest_coll_stren (digest2, &coll_stren[n]);

  if (nlevels != REFERENCE[k].nlevels || nlines != REFERENCE[k].nlines || n_coll_stren != REFERENCE[k].n_coll_stren
      || digest2 != REFERENCE[k].digest_levels_lines)
    fprintf (stderr, "\n%s: %d levels, %d lines, %d collision strengths, digest 0x%016llxULL; expected %d, %d, %d, 0x%016llxULL\n",
             REFERENCE[k].masterfile, nlevels, nlines, n_coll_stren, (unsigned long long) digest2, REFERENCE[k].nlevels,
             REFERENCE[k].nlines, REFERENCE[k].n_coll_stren, (unsigned long long) REFERENCE[k].digest_levels_lines);
  CU_ASSERT_EQUAL (nlevels, REFERENCE[k].nlevels);
  CU_ASSERT_EQUAL (nlines, REFERENCE[k].nlines);
  CU_ASSERT_EQUAL (n_coll_stren, REFERENCE[k].n_coll_stren);
  CU_ASSERT_EQUAL (digest2, REFERENCE[k].digest_levels_lines);
}

static void
test_standard80 (void)
{
  check_atomic_data_set (0);
}

static void
test_standard80_reduced (void)
{
  check_atomic_data_set (1);
}

static void
test_h20_hetop_standard80 (void)
{
  check_atomic_data_set (2);
}

static void
test_master_cno (void)
{
  check_atomic_data_set (3);
}

static void
test_fe_25fe_26 (void)
{
  check_atomic_data_set (4);
}

static void
test_fe_17to27 (void)
{
  check_atomic_data_set (5);
}

static int
suite_init (void)
{
  char target[LINELENGTH];
  struct stat sb;
  const char *sirocco_env = getenv ("SIROCCO");

  if (sirocco_env == NULL)
  {
    fprintf (stderr, "Failed to find SIROCCO environment variable\n");
    return EXIT_FAILURE;
  }
  snprintf (target, LINELENGTH, "%s/xdata", sirocco_env);
  if (!(stat (target, &sb) == EXIT_SUCCESS && S_ISDIR (sb.st_mode)))
  {
    perror ("Unable to find atomic data directory");
    return EXIT_FAILURE;
  }
  snprintf (ATOMIC_DATA_DEST, LINELENGTH, "data");
  if (symlink (target, ATOMIC_DATA_DEST) != EXIT_SUCCESS && errno != EEXIST)
  {
    perror ("Unable to created symbolic link for atomic data for test case");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

void
create_atomic_data_test_suite (void)
{
  CU_pSuite suite = CU_add_suite ("Atomic data", suite_init, NULL);

  if (suite == NULL)
  {
    fprintf (stderr, "Failed to create `Atomic data` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }

  if ((CU_add_test (suite, "standard80", test_standard80) == NULL) ||
      (CU_add_test (suite, "standard80_reduced", test_standard80_reduced) == NULL) ||
      (CU_add_test (suite, "h20_hetop_standard80", test_h20_hetop_standard80) == NULL) ||
      (CU_add_test (suite, "master_cno", test_master_cno) == NULL) ||
      (CU_add_test (suite, "fe_25fe_26_f10_fe_27_f1", test_fe_25fe_26) == NULL) ||
      (CU_add_test (suite, "fe_17to27", test_fe_17to27) == NULL))
  {
    fprintf (stderr, "Failed to add tests to `Atomic data` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
