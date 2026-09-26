
/***********************************************************/
/** @file  hydro2sirocco.c
 * @author ksl
 * @date   September, 2026
 *
 * @brief  Update a windsave file with the density, velocity and
 * temperature from a hydrodynamics model
 *
 * usage: hydro2sirocco -model_file model.txt [-out_root new] [-nonrel|-sr_doppler_only] root
 *
 * reads root.wind_save, applies the model and writes new.wind_save
 * (or out_root.wind_save).  The model file has the format used for
 * imported models.  The result can then be used with sirocco -r.
 *
 * ###Notes###
 *
 * This replaces the -model_file option of modify_wind, which did
 * not follow Sirocco's conventions for clumped winds (issue #1190).
 * The work is done by update_wind_from_model.
 *
 * The relativity switches must be the same as those used for the
 * sirocco runs, since they affect the gamma factors.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

char inroot[LINELENGTH], outroot[LINELENGTH], model_file[LINELENGTH];


/**********************************************************/
/**
 * @brief      parses the command line options
 *
 * @param [in]  int  argc   the number of command line arguments
 * @param [in]  char *  argv[]   The command line arguments
 *
 **********************************************************/

static void
h2s_parse_command_line (int argc, char *argv[])
{
  int i;
  char dummy[LINELENGTH];

  sprintf (outroot, "%s", "new");
  model_file[0] = '\0';
  rel_mode = REL_MODE_FULL;

  for (i = 1; i < argc - 1; i++)
  {
    if (strcmp (argv[i], "-model_file") == 0 && i + 1 < argc - 1)
    {
      sprintf (model_file, "%.*s", LINELENGTH - 1, argv[++i]);
    }
    else if (strcmp (argv[i], "-out_root") == 0 && i + 1 < argc - 1)
    {
      get_root (outroot, argv[++i]);
    }
    else if (strcmp (argv[i], "-gamma") == 0)
    {
      rel_mode = REL_MODE_FULL;
    }
    else if (strcmp (argv[i], "-sr_doppler_only") == 0)
    {
      rel_mode = REL_MODE_SR_FREQ;
    }
    else if (strcmp (argv[i], "-nonrel") == 0)
    {
      rel_mode = REL_MODE_LINEAR;
    }
    else
    {
      printf ("hydro2sirocco: Unknown or incomplete switch %s\n", argv[i]);
      exit (1);
    }
  }

  if (argc < 2 || strncmp (argv[argc - 1], "-", 1) == 0 || model_file[0] == '\0')
  {
    printf ("usage: hydro2sirocco -model_file model.txt [-out_root new] [-nonrel|-sr_doppler_only] root\n");
    exit (1);
  }

  strcpy (dummy, argv[argc - 1]);
  get_root (inroot, dummy);
}


/**********************************************************/
/**
 * @brief      read the windsave file, apply the model and write the result
 *
 * @param [in]  int  argc   the number of command line arguments
 * @param [in]  char *  argv[]   The command line arguments
 *
 * ###Notes###
 *
 * The model is applied to the one domain with an imported wind.
 *
 **********************************************************/

int
main (int argc, char *argv[])
{
  char infile[LINELENGTH], outfile[LINELENGTH];
  int ndom, n, nimport;

  Log_set_verbosity (3);
  h2s_parse_command_line (argc, argv);

  /* The average dv/ds in each cell is calculated with random directions */
  init_rand (1084515760);

  sprintf (infile, "%.150s.wind_save", inroot);
  sprintf (outfile, "%.150s.wind_save", outroot);

  zdom = calloc (MAX_DOM, sizeof (domain_dummy));
  if (zdom == NULL)
  {
    printf ("hydro2sirocco: Unable to allocate memory for domain\n");
    exit (1);
  }

  if (wind_read (infile) < 0)
  {
    printf ("hydro2sirocco: Unable to read %s\n", infile);
    exit (1);
  }

  ndom = -1;
  nimport = 0;
  for (n = 0; n < geo.ndomain; n++)
  {
    if (zdom[n].wind_type == IMPORT)
    {
      ndom = n;
      nimport++;
    }
  }
  if (nimport != 1)
  {
    printf ("hydro2sirocco: %s has %d imported domains; exactly one is needed\n", infile, nimport);
    exit (1);
  }

  printf ("hydro2sirocco: Applying %s to domain %d of %s and writing %s\n", model_file, ndom, infile, outfile);

  update_wind_from_model (ndom, model_file);

  wind_save (outfile);

  return (0);
}
