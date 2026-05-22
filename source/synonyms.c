/***********************************************************/
/** @file  synonyms.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  Translate the keywords used in previous versions of Sirocco to the current keywords
 *
 * The input variable used in Sirocco have evolved over time.  The routine provided hre
 * is intened to make it possible to associte a new keyword with one that was used in
 * an earlier version of Sirocco.  As long as the keyword has simply been renamed then one
 * can use a synomym to allow one to extact information from the old parameter file and translate
 * it so, one can use an old parmeter file with a new version of Sirocco.
 *
 * It is important to note that this process cannot continue indefinitely because , one may want
 * to change keyword/parameters so that they are not simple rplacements.
 ***********************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "log.h"

#define	LINELEN 256
char *old_names[] = { "mstar", "rstar", "Disk.illumination.treatment", "disk.type",
  "Disk_radiation", "Rad_type_for_disk", "disk.mdot", "T_profile_file",
  "disk.radmax",
  "stellar_wind_mdot", "stellar.wind.radmin", "stellar.wind_vbase",
  "stellar.wind.v_infinity", "stellar.wind.acceleration_exponent",
  "spectrum_wavemin", "spectrum_wavemax", "no_observers", "angle",
  "phase", "live.or.die", "spec.type", "mstar", "rstar", "Star_radiation",
  "tstar", "Rad_type_for_star", "Rad_type_for_star",
  "Rad_type_for_disk", "Rad_type_for_bl", "Boundary_layer_radiation",
  "Rad_type_for_bl", "t_bl", "lum_bl", "homologous_boundary_mdot", "msec",
  "period", "shell_wind_mdot", "Photon.sampling.approach", "Num.of.frequency.bands",
  "Lowest_energy_to_be_considered", "Highest_energy_to_be_considered", "Band.boundary",
  "Band.minimum_fraction",
  "agn_bremsstrahlung_temp", "agn_bremsstrahlung_alpha", "agn_blackbody_temp",
  "agn_power_law_cutoff", "geometry_for_pl_source", "lamp_post.height",
  "@Select_specific_no_of_scatters_in_spectra", "@Select_scatters", "@Select_photons_by_position",
  "@Select_location", "@rho", "@z", "@azimuth", "@r",
  "@save_cell_statistics", "@ispymode", "@keep_ioncycle_windsaves", "@make_ioncycle_tables",
  "@save_extract_photons", "@print_dvds_info", "@track_resonant_scatters",
  "@Use.standard.care.factors", "@Fractional.distance.photon.may.travel",
  "@Lowest.ion.density.contributing.to.photoabsorption", "@Keep.photoabs.during.final.spectrum",
  "@adjust_grid",
  "filling_factor", "Coord.system", "@write_atomicdata", "Fixed.concentrations.filename",
  "@Extra.diagnostics", "File.with.model2read", "Number.of.wind.components", "Old_windfile",
  "Model_file", "agn_power_law_index", "hydro_file", "Hydro_thetamax",
  "kn.acceleration_exponent", "kn.acceleration_length", "kn.d", "kn.mdot_r_exponent",
  "kn.rmax", "kn.rmin", "kn.v_infinity", "kn.v_zero", "QSO_BH_radiation", "lum_agn", "AGN.power_law_index",
  "AGN.blackbody_temp", "@AGN.power_law_cutoff", "AGN.geometry_for_pl_source", "Rad_type_for_agn", "Rad_type_for_agn",
  "wind.mdot", "Wind_ionization", "Wind_radiation", "Wind_type", "Thermal_balance_options",
  "wind.fixed_concentrations_file", "Disk-radiation",
  "AGN.bremsstrahlung_temp", "AGN.bremsstrahlung_alpha", "BH.blackbody_temp", "@BH.power_law_cutoff",
  "BH.geometry_for_pl_source", "BH.lamp_post_height", "BH.radiation", "BH.lum",
  "BH.rad_type_to_make_wind", "BH.rad_type_in_final_spectrum", "BH.power_law_index",
  "low_energy_break", "high_energy_break",
  "lum_agn", "AGN.power_law_index", "@AGN.power_law_cutoff",
  "AGN.lamp_post_height", "Reverb.filter_lines(0=off,-1=continuum,>0=count)",
  NULL
};

char *new_names[] = { "Central_object.mass", "Central_object.radius",
  "Surface.reflection.or.absorption", "Disk.type", "Disk.radiation",
  "Disk.rad_type_to_make_wind",
  "Disk.mdot", "Disk.T_profile_file", "Disk.radmax",
  "Stellar_wind.mdot", "Stellar_wind.radmin", "Stellar_wind.vbase",
  "Stellar_wind.v_infinity", "Stellar_wind.acceleration_exponent",
  "Spectrum.wavemin", "Spectrum.wavemax", "Spectrum.no_observers",
  "Spectrum.angle", "Spectrum.orbit_phase", "Spectrum.live_or_die",
  "Spectrum.type", "Central_object.mass", "Central_object.radius",
  "Central_object.radiation", "Central_object.temp", "Central_object.rad_type_to_make_wind",
  "Central_object.rad_type_in_final_spectrum",
  "Disk.rad_type_in_final_spectrum", "Boundary_layer.rad_type_in_final_spectrum",
  "Boundary_layer.radiation", "Boundary_layer.rad_type_to_make_wind", "Boundary_layer.temp",
  "Boundary_layer.luminosity", "Homologous.boundary_mdot", "Binary.mass_sec",
  "Binary.period", "Shell.wind_mdot", "Photon_sampling.approach", "Photon_sampling.nbands",
  "Photon_sampling.low_energy_limit", "Photon_sampling.high_energy_limit", "Photon_sampling.band_boundary",
  "Photon_sampling.band_min_frac",
  "Central_object.bremsstrahlung_temp", "Central_object.bremsstrahlung_alpha", "Central_object.blackbody_temp",
  "Central_object.power_law_cutoff", "Central_object.geometry_for_source", "Central_object.lamp_post_height",
  "@Spectrum.select_specific_no_of_scatters_in_spectra", "@Spectrum.select_scatters", "@Spectrum.select_photons_by_position",
  "@Spectrum.select_location", "@Spectrum.select_rho", "@Spectrum.select_z", "@Spectrum.select_azimuth", "@Spectrum.select_r",
  "@Diag.save_cell_statistics", "@Diag.ispymode", "@Diag.keep_ioncycle_windsaves", "@Diag.make_ioncycle_tables",
  "@Diag.save_extract_photons", "@Diag.print_dvds_info", "@Diag.track_resonant_scatters",
  "@Diag.use_standard_care_factors", "@Diag.fractional_distance_photon_may_travel",
  "@Diag.lowest_ion_density_for_photoabs", "@Diag.keep_photoabs_in_final_spectra",
  "@Diag.adjust_grid",
  "Wind.filling_factor", "Wind.coord_system", "@Diag.write_atomicdata", "Wind.fixed_concentrations_file",
  "@Diag.extra", "Wind.model2import", "Wind.number_of_components", "Wind.old_windfile",
  "Input_spectra.model_file", "AGN.power_law_index", "Hydro.file", "Hydro.thetamax",
  "KWD.acceleration_exponent", "KWD.acceleration_length", "KWD.d", "KWD.mdot_r_exponent",
  "KWD.rmax", "KWD.rmin", "KWD.v_infinity", "KWD.v_zero", "Central_object.radiation", "Central_object.luminosity",
  "Central_object.power_law_index",
  "Central_object.blackbody_temp", "@Central_object.power_law_cutoff", "Central_object.geometry_for_pl_source",
  "Central_object.rad_type_in_final_spectrum", "Central_object.rad_type_to_make_wind",
  "Wind.mdot", "Wind.ionization", "Wind.radiation", "Wind.type", "Wind_heating.extra_processes",
  "Wind.fixed_concentrations_file", "Disk.radiation",
  "Central_object.bremsstrahlung_temp", "Central_object.bremsstrahlung_alpha", "Central_object.blackbody_temp",
  "@Central_object.power_law_cutoff",
  "Central_object.geometry_for_source", "Central_object.lamp_post_height", "Central_object.radiation", "Central_object.luminosity",
  "Central_object.rad_type_to_make_wind", "Central_object.rad_type_in_final_spectrum", "Central_object.power_law_index",
  "Central_object.cloudy.low_energy_break", "Central_object.cloudy.high_energy_break",
  "Boundary_layer.luminosity", "Boundary_layer.power_law_index", "@Boundary_layer.power_law_cutoff",
  "Central_object.lamp_post_height", "Reverb.filter_lines(0=off,-1=continuum,-2=dumpall,>0=count)",
  NULL
};

int number_of_names = 122;
int synonyms_validated = 0;

#define MIN(a,b) ((a)<b ? a:b)


/**********************************************************/
/**
 * @brief  Given a question with variable name and info
 *   prompt, determines the length of the actual question name.
 *
 * @param [in] char  question[]   The full question e.g. xyz(cm/s)
 * @return number of characters in the string before (,
 *   whitespace or EOL.
 *
 * This simple helper function is intended to get string length
 *   so that full lines can be easily compared on only the
 *   subset that is actually the variable name. It can cope
 *   with both question-only strings e.g. "xyz(cm/s)" or question
 *   and answer strings e.g. "xyz(cm/s) 10.7e3"
 **********************************************************/
int
get_question_name_length (question)
     char question[];
{
  char *found_location;

  // strchr will return NULL if it can't find the character in the question string,
  // or a pointer to the memory location it finds the character in. The length of the
  // string is thus the difference in memory locations.
  if ((found_location = strchr (question, '(')))
  {
    return (int) (found_location - question);

  }
  else if ((found_location = strchr (question, ' ')))
  {
    return (int) (found_location - question);

  }
  else if ((found_location = strchr (question, '\t')))
  {
    return (int) (found_location - question);

  }
  else if ((found_location = strchr (question, '\n')))
  {
    return (int) (found_location - question);

  }
  else
  {
    return strlen (question);
  }
}


/**********************************************************/
/**
 * @brief  Validates the hard-coded lists of synonyms
 *
 * @return  1 if they're valid, else throws an Error.
 *
 * Tries to match the lengths of the lists of new and old names
 *   with the expected number of synonyms. On a failure, exits the
 *   run.
 *
 * ###Notes###
 *
 * Broken out from check_synonyms().
 *
 **********************************************************/

int
are_synonym_lists_valid ()
{
  int i, n;

  int n_old_names = -1;
  while (old_names[++n_old_names] != NULL)
  {                             /* do nothing */
  }

  int n_new_names = -1;
  while (new_names[++n_new_names] != NULL)
  {                             /* do nothing */
  }

  if (n_new_names != n_old_names || number_of_names != n_old_names)
  {
    Error ("check_synonyms: %d %d %d\n", number_of_names, n_old_names, n_new_names);
    n = MIN (n_new_names, n_old_names);
    for (i = 0; i < n; i++)
    {
      Log ("%3d %40s %40s\n", i, old_names[i], new_names[i]);
    }
    Exit (0);
  }
  return (1);
}


/**********************************************************/
/**
 * @brief  Checks to see if the question on an input line
 *   is an old version of the question being asked
 *
 * @param [in] char  new_question[]   The currnt keyword
 * @param [in] char  old_question[]   The keyword in an earlier version of Sirocco
 * @return  1 if the input line question is a synonym for
 *   the question we're trying to ask, else 0.
 *
 * The routine scans through the list of synonyms for a given question,
 * to check to see if the current line matches any of the possible old versions.
 *
 * ###Notes###
 *
 * It would be much more efficient to do this in file input
 * rather than checking the synonym status of every line for
 * every parameter.
 *
 **********************************************************/

int
is_input_line_synonym_for_question (question, input_line)
     char question[];
     char input_line[];
{
  int synonym_index;
  int question_name_length = get_question_name_length (question);

  // First, if we haven't already, let's see if the synonym lists are valid
  if (!synonyms_validated)
  {
    synonyms_validated = are_synonym_lists_valid ();
  }

  // We want to know what the possible 'old names' are for the question we're trying to ask.
  for (synonym_index = 0; synonym_index < number_of_names; synonym_index++)
  {

    // For each synonym, if the 'new name' refers to the question we're trying to ask...
    if (!strncmp (new_names[synonym_index], question, question_name_length))
    {
      // Does the 'old name' match the question on the input line?
      if (!strncmp (old_names[synonym_index], input_line, get_question_name_length (old_names[synonym_index])))
      {
        // If the question on the input line matches the 'old name' for this question (i.e. the difference is 0)
        return (1);
      }

    }
  }
  // If we've not found a match, then this line *isn't* a synonym
  return (0);
}
