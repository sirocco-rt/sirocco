
/***********************************************************/
/** @file  compton.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  These are the routines to handle Compton scattering.
 *
 ***********************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "atomic.h"
#include "sirocco.h"

PlasmaPtr xplasma;              /// Pointer to current plasma cell



/**********************************************************/
/** 
 * @brief      carry out the process of Compton scattering a photon    
 *
 * @param [in] Photpr p  A photon                        
 * @return     0
 *
 * @details
 * 
 * ### Notes ###
 *
 * This routine oversees the compton scatterin of a single photon
 *
 **********************************************************/


int
compton_scatter (p)
     PhotPtr p;                 // Pointer to the current photon
{
  double t_e;
  double vel[3];
  double v;


  WindPtr one;
  PlasmaPtr xplasma;
  one = &wmain[p->grid];
  xplasma = &plasmamain[one->nplasma];




  t_e = xplasma->t_e;

  compton_get_thermal_velocity (t_e, velocity_electron);

  /* Now account for the fact that the photons sees fewer 
     electrons traveling in the direction of the initial
     photon than those which are moving towards it

     What we want to assure is that the fraction of photons
     that go in the direction of the photon is given
     by 

     0.5 (1-v/c)

     and the fraction that goes in the opposite direction
     of the photon is

     0.5 (1+v/c)

     where v is the component of the velocity vector of
     the electron going the same direction as the photon.

     If the electron were going near the speed of light
     in this direction, then an observer in the rest
     frame would see that the photon was always encontuering
     photons moving towards the photons, and none moveing
     with it.

     Note that because the velocities in the 3 directions
     are uncorrellated, we do not have to go into the frame
     aligned with the photon, and revese only that component
     and transform back

     This calculation is carried out non-relativistically,
     as is the thermal velocity calculation.
   */

  v = dot (p->lmn, velocity_electron);

  if (random_number (0.0, 1.0) < 0.5 * (1. + fabs (v / VLIGHT)))
  {
    /* Then we want the photon to be headed  towards the photon */
    if (v > 0)
    {
      velocity_electron[0] = (-velocity_electron[0]);
      velocity_electron[2] = (-velocity_electron[1]);
      velocity_electron[2] = (-velocity_electron[2]);
    }
  }
  else if (v < 0)
  {
    velocity_electron[0] = (-velocity_electron[0]);
    velocity_electron[2] = (-velocity_electron[1]);
    velocity_electron[2] = (-velocity_electron[2]);
  }





  lorentz_transform (p, p, velocity_electron);
  if (modes.save_extract_photons)
    save_photons (p, "BeforeC");
  compton_dir (p);
  if (modes.save_extract_photons)
    save_photons (p, "AfterC");
  rescale (velocity_electron, -1, vel);
  lorentz_transform (p, p, vel);


  return (0);
}




/**********************************************************/
/** 
 * @brief      computes the effective Compton opacity in the cell.
 *
 * @param [in] PlasmaPtr  xplasma   - pointer to the current plasma cell
 * @param [in] double  freq   - the frequency of the current photon packet being tracked
 * @return     kappa - the Compton opacity for the cell.
 *
 * @details
 * Calculate opacity of a cell due to Compton (Thompson) scattering.
 *
 * The routine gives e correct *heating* effect when high frequency photons
 * undergo Compton scattering.
 *
 * The implementation follows that given in Hazy 3 to calculate the frequency dependant
 * cross section. This cross section takes account of the energy exchange, i.e.
 * not only the cross section is frequeucy dependant, but so is the mean energy
 * transfer.
 *
 * ### Notes ###
 *
 * The treatment is necessary because the Klein-Nishina x-section
 * gives incorrect results in terms of the actual effect of heatin

 * This implements the equation kappa=(sigmaT_h*ne*freq)/(me*c^2)
 *             Initally it is called by radiation
 *
 **********************************************************/

double
kappa_comp (xplasma, freq)
     PlasmaPtr xplasma;
     double freq;
{
  double x;
  double sigma;                 /*The cross section, thompson, or KN if hnu/mec2 > 0.01 */
  int ndom;


  sigma = compton_alpha (freq) * THOMPSON;      //the energy exchange cross section

  x = (sigma * PLANCK) / (MELEC * VLIGHT * VLIGHT);
  x *= xplasma->ne * freq;

  ndom = wmain[xplasma->nwind].ndom;
  x *= zdom[ndom].fill;         // multiply by the filling factor-
  return (x);
}




/**********************************************************/
/** 
 * @brief      Computes the effective opacity in the cell due to *induced* Compton heating.
 *
 * @param [in] PlasmaPtr  xplasma   - pointer to the current plasma cell
 * @param [in] double  freq   - the photon frequency               
 * @return     kappa - the inducd Compton opacity for the cell.
 *
 * @details
 * This routine implements a rather unusual situation, where the mean intensity in
 * a cell exceeds the blackbody value and we get induced Compton heating.
 *
 * ### Notes ###
 * This implements the induced Compton term in equation 6.5 in Hazy vol 3 - this is
 * the last term in the equation. It will
 * rarely if ever occur in a real model, but was included to give better agreement
 * with Cloudy over a wide range of ionization parameters.
 *
 **********************************************************/

double
kappa_ind_comp (xplasma, freq)
     PlasmaPtr xplasma;
     double freq;
{
  double x;
  double sigma;                 /*The cross section, thompson, or KN if hnu/mec2 > 0.01 */
  double J;                     //The estimated intensity in the cell
  int ndom;


  /* Obtain a model for the mean intensity - we call this with mode=2, which means
     that if we have not yet completed a cycle, dont return a dilute blackbody
     estimate if we are in PL mode. */
  J = mean_intensity (xplasma, freq, MEAN_INTENSITY_ESTIMATOR_MODEL);


  sigma = THOMPSON * compton_alpha (freq);      //the energy exchange cross section

  x = (xplasma->ne) / (MELEC);
  x *= sigma * J;
  x *= 1 / (2 * freq * freq);

  ndom = wmain[xplasma->nplasma].ndom;
  x *= zdom[ndom].fill;         // multiply by the filling factor

  if (sane_check (x))
  {
    Error ("kappa_ind_comp:sane_check - undefined value for Kappa_ind_comp - setting to zero\n");
    return (0.0);
  }

  return (x);
}


/**********************************************************/
/** 
 * @brief      computes the cooling in the cell due to inverse Compton scattering.
 *
 * @param [in] WindPtr  one   - pointer to the current wind cell
 * @param [in] double  t_e   - electron temperature
 * @return     the Compton luminosity for the cell.
 *
 * @details
 * This routine computes the temperature dependant Compton cooling rate in ergs/s.
 * It is mainly called as part of the search for a new temperature balancing heating
 * and cooling rates after an ionization cycle.
 *
 *
 * ### Notes ##
 * Compton cooling or 'inverse Compton scattering' is the process of energy transfer
 * from hot electrons to 'cool' photons, that is when 4kT > hnu.
 *
 * For each interaction, the fractional energy loss from the electrons is given
 * by \frac{}\delta e}{e} = 4k_BT/mc^2.
 * where we have a model for the mean intensity of radiation in a cell. 
 *
 * We 
 * compute the cooling rate by doing an integration over the mean intensity
 * multiplied by the frequency dependant cooling cross section.
 * which looks like C=4k_BT/mc^2 * 4/pi \int{\sigma_T * beta * J_{\nu} d\nu}.
 * Beta is the effective cooling cross section from fits in Hazy 3.
 *
 * If we are in the low frequency limit, we dont need to integrate, because the
 * cooling cross section is not frequency dependant in the Thompson limit so we
 * can just multiply by the band limited mean intensity.
 *
 * Similarly, if we dont have a model - we just do the best we can which is
 * to multiply by the total integrated mean indensity times the Thormpson
 * cross section.
 *
 * Like other cooling routines, this function uses quantities in the CMF
 *
 **********************************************************/

double
total_comp (one, t_e)
     WindPtr one;
     double t_e;
{
  double x, f1, f2;
  int nplasma, j;


  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  x = 0.0;

  //Since J_nu is constant for a given cycle - we only need to compute the integral once when searching for a thermal balance
  if (xplasma->comp_nujnu < 0.0)
  {
    if (geo.spec_mod)           //Check to see if we have generated a spectral model
    {
      for (j = 0; j < geo.nxfreq; j++)
      {
        if (xplasma->spec_mod_type[j] != SPEC_MOD_FAIL) //Only bother doing the integrals if we have a model in this band
        {
          f1 = xplasma->fmin_mod[j];
          f2 = xplasma->fmax_mod[j];
          if (f1 > 1e18)
            x += num_int (comp_cool_integrand, f1, f2, 1e-6);

          else
            x += THOMPSON * xplasma->xj[j];     //If in the Thompson limit, we just multiply the band limited frequency integrated mean intensity by the Thompson cross section
        }
      }
    }

    else                        //If no spectral model - we do the best we can, and multply the mean intensity by the Thompson cross section.
    {
      x = THOMPSON * xplasma->j;
    }
    xplasma->comp_nujnu = x;
  }
  else
    x = xplasma->comp_nujnu;

  //Multply by the other terms - including temperature - this gives the temperature dependance of this cooling term.
  x *= (16. * PI * BOLTZMANN * t_e * xplasma->ne) / (MELEC * VLIGHT * VLIGHT) * xplasma->vol;


  return (x);
}


/**********************************************************/
/** 
 * @brief      computes the (total) Klein Nishina cross section for a photon.
 *
 * @param [in] double  nu   photon frequency
 * @return     the KN cross section.
 *
 * @details
 * This routine is used to compute the frequency dependant cross section
 * of a photon to electron scattering. 
 * 
 * ### Notes ###
 * This implements equation 7.5 in Rybicki and Lightman
 *
 * For a frequency of less than
 * 1e18Hz this is just the Thompson scattering cross section. It is currently
 * only used to compute the chance of a photon interacting in this way. The
 * actual eneergy loss when an interaction occurs requires a cross section
 * that takes account of energy balance (provided by the alpha function.)
 *
 **********************************************************/

double
klein_nishina (nu)
     double nu;                 //The frequency of the photon packet
{
  double x;                     //h nu / kt
  double x1, x2, x3, x4;        //variables to store intermediate results.
  double kn;                    // the final cross section

  kn = THOMPSON;
  x1 = x2 = x3 = x4 = 0.0;

  /* x is the photon energey relative to the electron mass */
  x = (PLANCK * nu) / (MELEC * VLIGHT * VLIGHT);

  /* Use the full KN formula if the photon energy is high enough.  
   * Otherwise just use the Thompson x-section.
   */
  if (x > 0.0001)
  {
    x1 = 1. + x;
    x2 = 1. + (2. * x);
    x3 = log (x2);
    x4 = ((2. * x * x1) / x2) - x3;
    x4 *= x1 / (x * x * x);
    x4 = x4 + x3 / (2. * x);
    x4 = x4 - (1 + 3. * x) / (x2 * x2);
    kn *= 0.75 * x4;
  }

  return (kn);
}

//External variables to allow zero_find to search for the correct fractional energy change

double sigma_rand;              //The randomised cross section that our photon will see
double sigma_max;               //The cross section for the maxmimum energy loss
double x1;                      //The ratio of photon eneergy to electron energy

/** ****************************************************************************
 *
 * @brief Set the value for parameters used in comp_func
 *
 * @param [in] rand_cs The random cross section to set
 * @param [in] max_cs The maximum cross section to set
 * @param [in] energy_ratio The photon to electron energy ratio to set
 *
 * @details
 *
 * `sigma_rand`, `sigma_max` and `x1` are used in `comp_func` to calculate the
 * fractional energy change for a given frequency and photon. These variables
 * are "external" (i.e. global in this file), as they are not passed as
 * parameters to `comp_func` due to how a root finding algorithm is setup.
 *
 * The existence of this function is justified by some usages of comp_func
 * outside of this file, e.g. for unit testing. The other option would be to
 * write the unit tests in this file, or to have a compton.h header where the
 * globals are stored.
 *
 * ************************************************************************** */

void
set_comp_func_values (double rand_cs, double max_cs, double energy_ratio)
{
  sigma_rand = rand_cs;
  sigma_max = max_cs;
  x1 = energy_ratio;
}

/**********************************************************/
/** 
 * @brief      find a new (random) direction for a photon undergoing Compton scattering.
 *
 * @param [in,out] PhotPtr  p   - the photon currently being scattered
 * @return     0 if successful - the new direction is returned as part of the photon structure
 *
 * @details
 * This routine computes a new direction (and energy change) for a photon undergoing Compton
 * scattering. 
 *
 * If the frequency is low with respect to the electron rest mass then it is
 * totally elastic and one obtains isotropic Thompson scattering.
 *
 * However if the frequency is higher one needs to do actually find the direction via a zero 
 * finding routine 
 * 
 * ### Notes ###
 * Everything is calculated in the rest frame of the electron.
 *
 * In this frame, the photon energy change is E/E'=1+(hnu/mc^2)(1+cos\theta)
 * where \theta is the angle through which the photon is deflected.
 * This is a maximum for \theta=180 - E/E'=1+2hnumc^2
 * and minimum for \theta=0 - E/E'=1
 * 
 * We first draw a random cross section that our photon packet will see.
 * This cross section represents an energy change and hence a direction. The
 * distribnution of angles is taken care of by using a differential cross section
 * vs energy change function. So we work out what energy change and hence direction
 * is implied by our random cross section.
 *
 **********************************************************/

int
compton_dir (p)
     PhotPtr p;                 // Pointer to the current photon

{
  double f_min, f_max, f;       //Fractional energy changes - E_old/E_new - minimum possible, maximum possible, actual as implied by random cross section
  double n, l, m, phi, len;     //The direction cosines of the new photon direction in the frame of reference with q along the photon path
  struct basis nbasis;          //The basis function which transforms between the photon frame and the observer frame
  double lmn[3];                /* the individual direction cosines in the rotated frame */
  double x[3];                  /*photon direction in the frame of reference of the original photon */
  double dummy[3], c[3];
  int ierr = FALSE;

  x1 = PLANCK * p->freq / MELEC / VLIGHT / VLIGHT;      //compute the ratio of photon energy to electron energy. In the electron rest frame this is just the electron rest mass energy

  n = l = m = 0.0;

  if (x1 < 0.0001)              //If the photon energy is low, we use the diple approximation            
  {
    randvdipole (lmn, p->lmn);
    stuff_v (lmn, p->lmn);
  }
  else
  {

    /* The process is as follows:
       Generate a random number between 0 and 1 - this represents a randomised cross section (normalised to the maximum which out photon packet sees
       Calculate the mininumum and max energy change, corresponding to scattering angles of 0 and 180 degrees
       Calculate sigma_max, a variable that is communicated externally to the zero_find routine.  This is used to scale the K_N function to lie between 
       0 and 1. This is essentually the chance of a photon scattering through 180 degrees - or the angle giving the maximum energy loss
       Find the zero that represents our randomised fractional energy loss z_rand.
     */
    sigma_rand = random_number (0.0, 1.0);
    f_min = 1.;
    f_max = 1. + (2. * x1);
    sigma_max = sigma_compton_partial (f_max, x1);
    f = zero_find (compton_func, f_min, f_max, 1e-8, &ierr);
    if (ierr)
    {
      Error ("compton_dir: zero_find failed\n");
    }
/*We now have the fractional energy change f - we use the 'normal' equation for Compton scattering 
  to obtain the angle cosine n=cos(\theta)	for the scattering direction*/

    n = (1. - ((f - 1.) / x1)); //This is the angle cosine of the new direction in the frame of reference of the photon - this gives a 2D scattering angle

    if (isfinite (len = sqrt (1. - (n * n))) == 0)      //Compute the length of the other angle cosines - the isfinite is to take care of the very rare occasion where n=1!
      len = 0.0;                // If n=1, then the photon has either been undeflected or bounced straight back. Both are vanishingly unlikely but need to be treated.
    phi = 0.0;                  //no need to randomise phi, the random rotation of the vector generating the basis function takes care of this

    l = len * cos (phi);
    m = len * sin (phi);

    randvec (dummy, 1.0);
    cross (dummy, p->lmn, c);   //c will be perpendicular to p->lmn (the original photon direction) *and* the random vector just computed
    create_basis (p->lmn, c, &nbasis);  //create a basis with the first axis in the direction of the original photon direction, c will be perpendicular to the photon direction. Orientaion of the y/z axes are will give  arandomization of the phi axis

    x[0] = n;                   //This is the cosine direction of the new photon direction, in the frame of reference where the first axis is in the original photon direction
    x[1] = l;
    x[2] = m;

    project_from (&nbasis, x, lmn);     /* Project the vector from the FOR of the original photon into the observer frame */
    renorm (lmn, 1.0);          //Make sure the length of the direction vector is equal to 1
    stuff_v (lmn, p->lmn);      //Put the new photon direction into the photon structure

    p->freq = p->freq / f;      //reduce the photon frequency by the fractional energy change
    p->w = p->w / f;            //reduce the photon weight by the same ammount to conserve photon numbers
  }
  return (0);
}

/**********************************************************/
/** 
 * @brief      return the (unormalized) proability density Maxwell-Boltzman distribution
 *
 * @param [in] double x  
 * @param [in] void *params  dummy variable needed for gsl compatability
 * @return     the unormalized proablitly density for a Maxwell-Boltzmann distribution         
 *
 * @details
 * 
 * ### Notes ###
 *
 **********************************************************/

struct Cdf cdf_thermal;
int init_cdf_thermal = TRUE;


double
pdf_thermal (double x, void *params)
{
  return (x * x * exp (-x * x));
}

/**********************************************************/
/** 
 * @brief      get a random thermal velocity of an electron in a plasma
 *
 * @param [in] double t   temperature of the plasma
 * @param [out] double *v the velocity of a random electron
 * @return     0
 *
 * @details
 * 
 * ### Notes ###
 *
 * The PDF for the speed of thermalized particles is given by the Maxwell-Boltzmann
 * distribution.  Given a temperature, this routine gets a random
 * speed and then converts this to a velocity assuming that the 
 * directions are isotropic.
 *
 **********************************************************/


int
compton_get_thermal_velocity (t, v)
     double t, *v;
{
  double vel;


  if (init_cdf_thermal)
  {
    double dummy[2] = { 0, 1 };
    cdf_gen_from_func (&cdf_thermal, &pdf_thermal, 0, 5, 0, dummy);
    init_cdf_thermal = FALSE;
  }

  vel = cdf_get_rand (&cdf_thermal);

  vel *= sqrt ((2. * BOLTZMANN / MELEC) * t);

  if (vel > 0.5 * VLIGHT)
  {
    Error ("compton_get_thermal_velocity: v (%e) > 0.5 C\n", vel);
    vel = 0.5 * VLIGHT;
  }


  randvec (v, vel);



  return (0);



}


/**********************************************************/
/** 
 * @brief      Function used to find the fractional energy change for a given cross-section in the Compton scattering limit.
 *
 * @param [in] double  f   The fractional energy change of a photon undergoing scattering
 * @return     the difference between sigma(f) and the randomised cross section we are searching for
 *
 * @details
 * This function calculates the difference between the KN cross section for a given fractional
 * energy loss f/x1 and the randomly generated cross section for our current interaction. A zero
 * finding function is used to work out where this function is equal to zero - we need to
 * do this because we cannot invert the function sigma_compton_partial.
 *
 * ### Notes ###
 * Since this function is called by zbrent, many of the variables have to be communicated
 * externally. These are
 * x1, the ratio of incoming photon energy to electron rest mass,
 * sigma_max, the KN cross section for the maximum possible energy loss
 * sigma_rand - a randomised number between 0 and 1 representing the normalised cross section we want to find
 *
 **********************************************************/

double
compton_func (double f, void *params)
{
  double ans;
  ans = (sigma_compton_partial (f, x1) / sigma_max) - sigma_rand;
  return (ans);
}

/**********************************************************/
/** 
 * @brief      is the Klein Nishina cross section as a function of the fractional energy change
 *
 * @param [in] double  f   - the fractional energy change
 * @param [in] double  x   the energy of the incoming photon divided by the rest mass energy of an electron
 * @return     the cross section for the scattering angle that gives this fractional energy change.
 *
 * @details
 * This is a functiuon which returns the Klein-Nishina cross section implied by a supplied
 * fractional energy change for an incoming photon.
 *
 * ### Notes ###
 * Stuart Sim supplied this formula, but coundn't recall where he had found it, although
 * 			he had rederived it and thinks it is correct!
 *
 **********************************************************/

double
sigma_compton_partial (f, x)
     double f;                  //This is the fractional energy change, nu/nu'
     double x;                  //h nu/mec**2 - the energy of the photon divided by the rest energy of an electron
{
  double term1, term2, term3, tot;

  term1 = ((x * x) - (2 * x) - 2) * log (f) / x / x;
  term2 = (((f * f) - 1) / (f * f)) / 2;
  term3 = ((f - 1) / x) * ((1 / x) + (2 / f) + (1 / (x * f)));

  tot = 3 * THOMPSON * (term1 + term2 + term3) / (8 * x);

  return (tot);

}

/**********************************************************/
/** 
 * @brief      The 'heating' cross section correction to Thompson scattering
 *
 * @param [in] double  nu   - frequency
 * @return     alpha- the effective heating cross section relative to sigma_T
 *
 * @details
 * We calculate Compton heating due to photons passing through
 * a parcel of gas using an opacity. In order to work out the opacity, we
 * need the cross section that the electrons present to the photons, and
 * in order to take proper account of the energy transfer during a scatter
 * we use a 'heating crosss ection'. This is equal to 'alpha' times
 * the Thompson cross section - we calculate alpha here.
 *
 * ### Notes ###
 * The fit used here is given in Hazy 3, eq 6.6
 *
 **********************************************************/

double
compton_alpha (nu)
     double nu;
{
  double alpha;
  if (nu < 1e17)
    alpha = 1.0;
  else
    alpha = 1. / (1. + nu * HRYD * (1.1792e-4 + 7.084e-10 * nu * HRYD));
  return (alpha);
}


/**********************************************************/
/** 
 * @brief      The 'cooling' cross section correction to Thompson
 *
 * @param [in] double  nu   - frequency
 * @return     beta - the effective cooling cross section relative to sigma_T
 *
 * @details
 * The inverse Compton cooling rate in a cell is computed by ingegrating
 * the mean intensity in that cell by the frequency dependant cross section.
 * In order to take account of the energy transfer in an interaction, the
 * cross section must be the effective cooling cross section (not just the
 * Klein Nishina cross section)
 *
 * ### Notes ###
 * The fit used here is given in Hazy 3, eq 6.6
 *
 **********************************************************/

double
compton_beta (nu)
     double nu;
{
  double alp, beta;
  if (nu < 1e17)
    beta = 1.0;
  else
  {
    alp = compton_alpha (nu);
    beta = (1. - alp * nu * HRYD * (1.1792e-4 + (2 * 7.084e-10 * nu * HRYD)) / 4.);
  }
  return (beta);
}



/**********************************************************/
/** 
 * @brief      The integrand in the integral to obtain the total cooling rate due to inverse Compton scattering.
 *
 * @param [in] double  nu   - frequency
 * @param [in] void  params   An extra (unused) variable to make it paletable for the gsl integrator
 * @return     frequwncy dependant cross section multiplied by the mean intensity at frequency nu
 *
 * @details
 * This is is the integrand sigma x J_nu that is integrated
 * to obtain the Compton cooling rate in a cell. The sigma in question is the
 * effective energy exchange cross section, which for coolnig is /alpha/beta/sigma_T.
 * These alpha and beta terms are computed in ther subroutines, but here we only
 * need to use \beta, because the \alpha term is taken account of when we compute
 * J_nu during the photon transport phase of the code.
 * If we are calling this function, we know we have a model so there is no need to
 * protect against the problem of not having a model
 *
 * ### Notes ###
 * See Hazy 3 eq 6.6 for more details of this approach
 *
 **********************************************************/

double
comp_cool_integrand (double nu, void *params)
{
  double value;
  value = THOMPSON * compton_beta (nu) * mean_intensity (xplasma, nu, MEAN_INTENSITY_ESTIMATOR_MODEL);
  return (value);
}



/**********************************************************/
/** 
 * @brief      Calculate the renormalization needed for reweighting photons due to 
 *  the anisotropic nature of Compton scattering 
 *
 * @param [in] double  nu   - frequency of the unscattered photon
 *
 * @return     the value to use for the renormaliztion
 *
 * @details
 *
 * The normalization required depends only on the input photon
 *
 * ### Notes ###
 *
 * At very low values of frequency, provide the answer at
 * in the Thompson limit to avoid spurious results
 *
 **********************************************************/


double
compton_reweight_norm (nu)
     double nu;
{
  double v;
  double x1, x2, x3, x4, x5, xx;
  double r;

  r = PLANCK * nu / (MELEC * VLIGHT * VLIGHT);

  if (r < 0.00001)
  {
    v = 16. * PI / 3.;
    return v;
  }

  x1 = (2 / (r * r));
  xx = 1. + 2. * r;
  x2 = log (xx);
  x3 = 1 / r - 2 / (r * r) - 2 / (r * r * r);

  x4 = 1. / (2 * r) * (1 / (xx * xx) - 1.);

  x5 = 1 / xx * (-2 / (r * r) - 1. / (r * r * r)) + (2 / (r * r) + 1. / (r * r * r));


  v = x1 + x2 * x3 - x4 + x5;

  v *= 2. * PI;
  return v;
}


/**********************************************************/
/** 
 * @brief      Reweight a phton in order to extract it in a new direction due to 
 *  the anisotropic nature of Compton scattering 
 *
 * @param [in] double  nu   - frequency of the unscattered photon
 *
 * @details
 *
 * ### Notes ###
 *
 * An explanation of the philosophy behind the reweighting
 * can be found in docs/notes/compton_scattering/Compton_reweighting.ipynb
 *
 **********************************************************/


int
compton_reweight (p_in, p_out)
     PhotPtr p_in, p_out;


{
  double nu_in, nu_out, xr, reweight;
  double theta, ctheta;

  nu_in = p_in->freq;

  ctheta = dot (p_in->lmn, p_out->lmn);

  p_out->freq = nu_out = nu_in / (1 + PLANCK * nu_in / (MELEC * VLIGHT * VLIGHT) * (1 - ctheta));

  xr = nu_out / nu_in;
  theta = acos (ctheta);

  reweight = xr * xr * (xr + 1. / xr - sin (theta) * sin (theta));

  reweight *= 4. * PI / compton_reweight_norm (p_in->freq);



  /* This accounts for the reweighting that is due to
   * The anistorpy, but we must also change the
   * weight due simply to the change in freqency
   */

  p_out->w *= nu_out / nu_in * reweight;


  return (0);
}
