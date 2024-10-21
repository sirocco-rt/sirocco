
/***********************************************************/
/** @file  setup_star_bh.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  Files that get input parameters for a star and 
 * for an agn and  * boundray layer
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"
#include "models.h"


/**********************************************************/
/** 
 * @brief      sets rstar, mstar, tstar as well
 *   as secondary parameters based on user inputs
 *
 * @return     The system returns the luminosity of
 * the central object using the temperature and radius
 *
 * @details
 * This obtains the inputs for the Central Object, whether
 * we are dealing with an ordianry star of an AGN)  If the
 * system is a binary system, the parameters of the secondary
 * are also obtained.
 *
 * ### Notes ###
 * @bug This routine is not completely self consistent, in the
 * sense that for an AGN we cannot estimate a luminosity. This
 * reflects the fact that we have not completely setlled on how
 * to homogenize inputs for diffent system types.  It's not 
 * clear that we need to rerun anything.
 *
 **********************************************************/

double
get_stellar_params ()
{
  char answer[LINELENGTH];

  rdpar_comment ("Parameters for the Central Object");

  /* Describe the basic binary star system */

  geo.mstar /= MSOL;            // Convert to MSOL for ease of data entry

  if (geo.system_type == SYSTEM_TYPE_STAR)
  {
    geo.mstar = 52.5;
    geo.rstar = 1.32e12;
    geo.tstar_init = 42000;
  }
  else if (geo.system_type == SYSTEM_TYPE_CV)
  {
    geo.mstar = 0.8;
    geo.tstar_init = 40000;
  }
  else if (geo.system_type == SYSTEM_TYPE_BH)
  {
    geo.mstar = 10;
  }
  else if (geo.system_type == SYSTEM_TYPE_AGN)
  {
    geo.mstar = 1e9;
  }

  rddoub ("Central_object.mass(msol)", &geo.mstar);
  geo.mstar *= MSOL;


  /* If a BH we want geo.rstar to be at least as large as the last stable orbit for
   * a non-rotating BH
   */

  if (geo.system_type == SYSTEM_TYPE_AGN || geo.system_type == SYSTEM_TYPE_BH)
  {
    geo.rstar = 6. * GRAV * geo.mstar / (VLIGHT * VLIGHT);      //Set value to ISCO, namely 6 x Rg 
  }
  else if (geo.system_type == SYSTEM_TYPE_CV)
  {
    geo.rstar = wdrad (geo.mstar);
  }

  rddoub ("Central_object.radius(cm)", &geo.rstar);
  if (geo.rstar < 6. * GRAV * geo.mstar / (VLIGHT * VLIGHT))
  {
    Log ("Warning: Central object size %.1e is less than ISCO %.1e\n", geo.rstar, 6. * GRAV * geo.mstar / (VLIGHT * VLIGHT));
    geo.disk_rad_min = 6. * GRAV * geo.mstar / (VLIGHT * VLIGHT);

  }
  else
  {
    geo.disk_rad_min = geo.rstar;       //Normally this is what we want
  }


  geo.rstar_sq = geo.rstar * geo.rstar;
  if (geo.system_type != SYSTEM_TYPE_AGN && geo.system_type != SYSTEM_TYPE_BH)
  {
    strcpy (answer, "yes");
    geo.star_radiation = rdchoice ("Central_object.radiation(yes,no)", "1,0", answer);
    get_spectype (geo.star_radiation, "Central_object.rad_type_to_make_wind(bb,models)", &geo.star_ion_spectype);
    
    if (geo.star_ion_spectype == SPECTYPE_BB_FCOL)
    {
      Error ("Colour corrected BB not implemented for star at this stage. Exiting.");
      Exit (0);
    }

    if (geo.star_radiation)
      rddoub ("Central_object.temp", &geo.tstar_init);
  }
  else
  {
    geo.star_radiation = 0;
    geo.tstar_init = 0;
  }

  /* tstar_init and lum_star_init refer to values without the effects of backscattering */

  geo.tstar = geo.tstar_init;

  geo.lum_star = geo.lum_star_init = 4 * PI * geo.rstar * geo.rstar * STEFAN_BOLTZMANN * pow (geo.tstar, 4.);


  /* Describe the secondary if that is required */

  if (geo.system_type == SYSTEM_TYPE_CV || geo.system_type == SYSTEM_TYPE_BH)   /* It's a binary system */
  {

    geo.m_sec /= MSOL;          // Convert units for ease of data entry
    rddoub ("Binary.mass_sec(msol)", &geo.m_sec);
    geo.m_sec *= MSOL;

    geo.period /= 3600.;        // Convert units to hours for easy of data entry
    rddoub ("Binary.period(hr)", &geo.period);
    geo.period *= 3600.;        // Put back to cgs immediately
  }

  return (geo.lum_star_init);
}






/**********************************************************/
/** 
 * @brief      sets up the boundary layer and compact object 
 *   power law parameters
 *   based on user input and system type
 *
 * @param [in] double  lstar   A luminosity used to initialize the luminosity
 * of the agn in the event that there is no disk
 * @return    Always returns 0
 *
 * @details
 * The routine reads a number of parameters that have to do
 * with defining the boundary layer or the BH/NS in an 
 * X-ray binary or an AGN.  
 *
 *
 * ### Notes ###
 * 
 * Although the routine was orignally written for CVs and AGN
 * it now is used for any system where either a boundary layer
 * or a compact object is involved.  The common element here
 * is that luminosity is not derived from a temperature and 
 * size, but rather is something that the user gives as an
 * input variable (e.g as a fraction of the accretion luminosity).
 *
 **********************************************************/

int
get_bl_and_agn_params (lstar)
     double lstar;
{
  double xbl;
  double temp_const_agn;
  char answer[LINELENGTH];


  rdpar_comment ("Parameters for Boundary Layer or the compact object in an X-ray Binary or AGN");

  if (geo.system_type == SYSTEM_TYPE_AGN || geo.system_type == SYSTEM_TYPE_BH)  /* If it is an AGN */
  {
    strcpy (answer, "yes");
    geo.agn_radiation = rdchoice ("Central_object.radiation(yes,no)", "1,0", answer);
    geo.star_radiation = 0;     // 70b - AGN do not have a star at the center */
    geo.bl_radiation = 0;
    if (geo.agn_radiation)
      get_spectype (geo.agn_radiation, "Central_object.rad_type_to_make_wind(bb,models,power,cloudy,brems,mono)", &geo.agn_ion_spectype);

    if (geo.agn_ion_spectype == SPECTYPE_BB_FCOL)
    {
      Error ("Colour corrected BB not implemented for AGN at this stage. Exiting.");
      Exit (0);
    }
  }
  else
  {
    strcpy (answer, "no");
    geo.bl_radiation = rdchoice ("Boundary_layer.radiation(yes,no)", "1,0", answer);
    geo.agn_radiation = 0;      // So far at least, our star systems don't have a BH

    if (geo.bl_radiation)
      get_spectype (geo.bl_radiation, "Boundary_layer.rad_type_to_make_wind(bb,models,power)", &geo.bl_ion_spectype);
  
    if (geo.bl_ion_spectype == SPECTYPE_BB_FCOL)
    {
      Error ("Colour corrected BB not implemented for boundary layer. Exiting.");
      Exit (0);
    }
  }


  if (geo.agn_radiation && geo.agn_ion_spectype >= 0 && comp[geo.agn_ion_spectype].nmods != 1)
  {
    Error ("get_bl_and_agn_params: When using models with an AGN, there should be exactly 1 model, we have %i for ion cycles\n",
           comp[geo.agn_ion_spectype].nmods);
    exit (0);
  }


  /* 130621 - ksl - This is a kluge to add a power law to stellar systems.  What is done
     is to remove the bl emission, which we always assume to some kind of temperature
     driven source, and replace it with a power law source

     Note that the next 3 or 4 lines just tell you that there is supposed to be a power
     law source.  They don't teel you what the parameters are.
   */

  if (geo.bl_ion_spectype == SPECTYPE_POW)
  {
    geo.agn_radiation = 1;
    geo.agn_ion_spectype = SPECTYPE_POW;
    geo.bl_radiation = 0;
    Log ("Trying to make a star with a power law boundary layer\n");
  }
  else
  {
    Log ("Not trying to make a star with a power law boundary layer %d\n", geo.bl_ion_spectype);
  }


  /* Describe the boundary layer */

  if (geo.bl_radiation && geo.bl_ion_spectype != SPECTYPE_POW)
  {
    xbl = geo.lum_bl = 0.5 * GRAV * geo.mstar * geo.disk_mdot / geo.rstar;

    rddoub ("Boundary_layer.luminosity(ergs/s)", &geo.lum_bl);
    Log ("OK, the boundary layer lum will be about %.2e the disk lum\n", geo.lum_bl / xbl);
    rddoub ("Boundary_layer.temp(K)", &geo.t_bl);
  }
  else
  {
    geo.lum_bl = 0;
    geo.t_bl = 0;
  }

  /* Describe the agn */

  if (geo.agn_radiation && (geo.system_type == SYSTEM_TYPE_AGN || geo.system_type == SYSTEM_TYPE_BH))   /* This peculiar line is to enamble us to add a star with a power law component */
  {
    xbl = geo.lum_agn = 0.5 * GRAV * geo.mstar * geo.disk_mdot / geo.rstar;

    /* If there is no disk, initilize geo.lum to the luminosity of a star */
    if (geo.disk_type == DISK_NONE)
    {
      geo.lum_agn = lstar;
    }

    // At present we have set geo.r_agn = geo.rstar, and encouraged the user
    // set the default for the radius of the BH to be 6 R_Schwartschild.
    // rddoub("R_agn(cm)",&geo.r_agn);

    /* if we have a "blackbody agn" the luminosity is set by Stefan Boltzmann law
       once the AGN blackbody temp is read in, otherwise set by user */
    if (geo.agn_ion_spectype != SPECTYPE_BB)
    {
      rddoub ("Central_object.luminosity(ergs/s)", &geo.lum_agn);
      Log ("OK, the black hole luminosity will be about %.2e the disk lum\n", geo.lum_agn / xbl);
    }

    if (geo.agn_ion_spectype == SPECTYPE_POW || geo.agn_ion_spectype == SPECTYPE_CL_TAB)
    {
      geo.alpha_agn = (-1.5);
      rddoub ("Central_object.power_law_index", &geo.alpha_agn);

      if (geo.alpha_agn == -1.0)        //deal with the pathological case
      {
        geo.const_agn = geo.lum_agn / (log (2.42e18) - log (4.84e17));
      }
      else
      {
        geo.const_agn = geo.lum_agn / (((pow (2.42e18, geo.alpha_agn + 1.)) - pow (4.84e17, geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));
      }
      Log ("BH Input parameters give a power law constant of %e\n", geo.const_agn);
    }
    else if (geo.agn_ion_spectype == SPECTYPE_BREM)
    {

      geo.brem_temp = 1.16e8;   //10kev
      geo.brem_alpha = -0.2;    //This is the cloudy form of bremstrahlung
      geo.const_agn = 1.0;
      rddoub ("Central_object.bremsstrahlung_temp(K)", &geo.brem_temp);
      rddoub ("Central_object.bremsstrahlung_alpha", &geo.brem_alpha);
      temp_const_agn = geo.lum_agn / num_int (integ_brem, 4.84e17, 2.42e18, 1e-4);
      geo.const_agn = temp_const_agn;
      Log ("AGN Input parameters give a Bremsstrahlung constant of %e\n", temp_const_agn);

    }
    else if (geo.agn_ion_spectype == SPECTYPE_BB)
    {
      /* note that alpha_agn holds the temperature in the case of "blackbody agn" */
      rddoub ("Central_object.blackbody_temp(K)", &geo.alpha_agn);
      geo.lum_agn = 4 * PI * geo.rstar * geo.rstar * STEFAN_BOLTZMANN * pow (geo.alpha_agn, 4.);
      Log ("OK, the black hole/AGN lum will be about %.2e the disk lum\n", geo.lum_agn / xbl);
    }

    /* JM 1502 -- lines to add a low frequency power law cutoff. accessible
       only in advanced mode and for non broken power law.
       default is zero which is checked before we call photo_gen_agn */
    geo.pl_low_cutoff = 0.0;
    if (modes.iadvanced && (geo.agn_ion_spectype == SPECTYPE_POW))
      rddoub ("@Central_object.power_law_cutoff", &geo.pl_low_cutoff);

    strcpy (answer, "sphere");
    geo.pl_geometry = rdchoice ("Central_object.geometry_for_source(sphere,lamp_post,bubble)", "0,1,2", answer);


    if (geo.pl_geometry == PL_GEOMETRY_LAMP_POST)
    {
      rddoub ("Central_object.lamp_post_height(r_g)", &geo.lamp_post_height);
      geo.lamp_post_height *= GRAV * geo.mstar / VLIGHT / VLIGHT;       //get it in CGS units
      Log ("lamp_post_height in cm is %g\n", geo.lamp_post_height);
    }
    else if (geo.pl_geometry == PL_GEOMETRY_BUBBLE)
    {
      rddoub ("Central_object.bubble_size(r_g)", &geo.bubble_size);
      geo.bubble_size *= GRAV * geo.mstar / VLIGHT / VLIGHT;    //get it in CGS units
      Log ("bubble size  in cm is %g\n", geo.bubble_size);
    }
    else if (geo.pl_geometry != PL_GEOMETRY_SPHERE)     // only three options at the moment
    {
      Error ("Did not understand power law geometry %i. Fatal.\n", geo.pl_geometry);
      Exit (0);
    }



    /* Computes the constant for the power law spectrum from the input alpha and 2-10 luminosity.
       This is only used in the sim correction factor for the first time through.
       Afterwards, the photons are used to compute the sim parameters. */



    if (geo.agn_ion_spectype == SPECTYPE_CL_TAB)        /*NSH 0412 - option added to allow direct comparison with cloudy power law table option */
    {
      geo.agn_cltab_low = 1.0;
      geo.agn_cltab_hi = 10000;
      rddoub ("Central_object.cloudy.low_energy_break(ev)", &geo.agn_cltab_low);        /*lo frequency break - in ev */
      rddoub ("Central_object.cloudy.high_energy_break(ev)", &geo.agn_cltab_hi);
      geo.agn_cltab_low_alpha = 2.5;    //this is the default value in cloudy
      geo.agn_cltab_hi_alpha = -2.0;    //this is the default value in cloudy
    }
  }
  else if (geo.agn_radiation)   /* We want to add a power law to something other than an AGN */
  {
    xbl = geo.lum_agn = 0.5 * GRAV * geo.mstar * geo.disk_mdot / geo.rstar;

    // At present we have set geo.rstar = geo.rstar, and encouraged the user
    // set the default for the radius of the BH to be 6 R_Schwartschild.
    // rddoub("R_agn(cm)",&geo.rstar);

    rddoub ("Boundary_layer.luminosity(ergs/s)", &geo.lum_agn);
    Log ("OK, the boundary layer lum will be about %.2e the disk lum\n", geo.lum_agn / xbl);
    geo.alpha_agn = (-1.5);
    rddoub ("Boundary_layer.power_law_index", &geo.alpha_agn);

    /* JM 1502 -- lines to add a low frequency power law cutoff. accessible
       only in advanced mode. default is zero which is checked before we call photo_gen_agn */
    geo.pl_low_cutoff = 0.0;
    if (modes.iadvanced)
      rddoub ("@Boundary_layer.power_law_cutoff", &geo.pl_low_cutoff);


    /* Computes the constant for the power law spectrum from the input alpha and 2-10 luminosity.
       This is only used in the sim correction factor for the first time through.
       Afterwards, the photons are used to compute the sim parameters. */


    if (geo.alpha_agn == -1.0)  //deal with the pathological case
    {
      geo.const_agn = geo.lum_agn / (log (2.42e18) - log (4.84e17));
    }
    else
    {
      geo.const_agn = geo.lum_agn / (((pow (2.42e18, geo.alpha_agn + 1.)) - pow (4.84e17, geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));
    }

    Log ("Boundary layer input parameters give a power law constant of %e\n", geo.const_agn);
  }

  else
  {
    geo.lum_agn = 0.0;
    geo.alpha_agn = 0.0;
    geo.const_agn = 0.0;
  }
  return (0);
}
