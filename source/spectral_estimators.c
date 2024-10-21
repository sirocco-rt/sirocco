
/***********************************************************/
/** @file  spectral_estimators.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  The routines in this file are all associated with 
 * the process of modelling the spectrum in a cell based
 * on information that was accrued during an ionization
 * cycle.  
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

/* external variables set up so zbrent can solve for various variables  */

/// Minimum, maximum and mean frequency in a band
double spec_numin, spec_numax, spec_numean;
/// Log versions of numin and numax - the band ends
double lspec_numin, lspec_numax;


/**********************************************************/
/**
 * @brief      Calculate the parameters of a model for J_nu in a cell
 *
 * @param [in] PlasmaPtr  xplasma   Pointer to the a specific element in the plasma structure
 * @return     0 if successful
 *
 * @details
 * This routine uses the mean frequency (calculated during photon
 * transport) find a power law model and an exponential model.
 * It then uses the standard deviation to decide which is best.
 * The routine is called from ionization.c
 * The results are stored in parameters in the PlasmaPtr structure
 * in variables, such as spec_mod_type, pl_alpah, pl_log,exp_temp, etc.
 *
 * ### Notes ###
 * TODO this function only returns 0, the documentation claims otherwise however
 *
 **********************************************************/

int
spectral_estimators (xplasma)
     PlasmaPtr xplasma;
{
  double pl_alpha_min, pl_alpha_max, pl_alpha_temp, pl_w_temp, j;
  double exp_temp_min, exp_temp_max, exp_temp_store;    /* The 'temperature' range we are going to search for an effective temperature for the exponential model */
  double exp_temp_temp, exp_w_temp;     /* The temporary values for temperature and weight of the exponential model */
  int n, n1;
  double pl_sd, exp_sd;         /* Computed standard deviations for two models for comparison with true value */
  int plflag, expflag;          /* Two flags to say if we have a reasonable PL or EXP model,
                                   set to 1 initially, -1 means there has been some failure that means
                                   we must not use this model, +1 means it is OK */
  double genmin, genmax;        /*The min and max frequencies over which we have made photons originally (actually band ends) */
  double dfreq;                 /* A number to help work out if we have fully filled a band */
  int ierr = FALSE;

  /* This call is after a photon flight, so we *should* have access to j and ave_freq,
     and so we can calculate proper values for W and alpha
     To avoid problems with solving, we need to find a reasonable range of values within
     which to search for a solution. A reasonable guess is that it is around the current value....
   */
  genmin = xband.f1[0];
  genmax = xband.f2[xband.nbands - 1];

  /* NSH 131108 The next few lines just work out which bands actually should have photons in -
     i.e. if we dont generate any photons in a band, we will not worry if we didnt see any photons in that band */

  for (n1 = 0; n1 < xband.nbands; n1++)
  {
    if (xband.nphot[n1] > 0)
    {
      genmin = xband.f1[n1];    //genmin will get set to the lower frecuency band of the first band with any photons generated
      break;
    }
  }

  for (n1 = xband.nbands - 1; n1 > -1; n1--)
  {
    if (xband.nphot[n1] > 0)
    {
      genmax = xband.f2[n1];    //genmax will get set to the upper frecuency band of the last band with any photons generated
      break;
    }
  }

  /* We loop over all of the bands, the first band is band number 0, and the last is band nxfreq-1 */
  for (n = 0; n < geo.nxfreq; n++)

  {
    Log_silent
      ("Starting out band %i in cell %i. mean=%e, sd=%e, minfreq=%e, maxfreq=%e, nphot=%i\n",
       n, xplasma->nplasma, xplasma->xave_freq[n], xplasma->xsd_freq[n], xplasma->fmin[n], xplasma->fmax[n], xplasma->nxtot[n]);

    plflag = expflag = 1;       //Both potential models are in the running

    if (xplasma->nxtot[n] <= 1) /* Catch the situation where there are only 1 or 0 photons in a band -
                                   we cannot reasonably try to model this situation */
    {
      if (geo.xfreq[n] >= genmax || geo.xfreq[n + 1] <= genmin)
      {
        /* The band is outside where photons were generated, so not very
           surprising that there are no photons - just generate a log */
        Log_silent ("spectral_estimators: too few photons (1 or 0) in cell %d band %d but we weren't expecting any \n", xplasma->nplasma, n);   // This is just a log, clearly it is not a problem
      }

      else
      {
        Error ("spectral_estimators: too few photons (1 or 0) in cell %d (%d) band %d to produce a model\n", xplasma->nplasma,
               wmain[xplasma->nwind].inwind, n);
      }

      /* We also want to make sure that the weight will be zero, this way we make
         sure there is no contribution to the ionization balance from this frequency. */
      xplasma->pl_log_w[n] = -999;      //A very tiny weight
      xplasma->pl_alpha[n] = 999.9;     //Give alpha a value that will show up as an error
      xplasma->exp_w[n] = 0.0;  //Make sure that w is zero, so no chance of mucking up ionization balance even if for some reason we end up integrating
      xplasma->exp_temp[n] = 1e99;      //Give the effective temperature a value that will show up as an error
      xplasma->spec_mod_type[n] = SPEC_MOD_FAIL;        //This tells the code that we have failed to model the spectrum in this band/cell/
    }


    /*  If all the photons in the cell are concentrated in a tiny range then we will also not
       expect to make a sensible model - this check could be reviewed later if lots of warning are produced */

    else if (xplasma->fmax[n] == xplasma->fmin[n])
    {
      Error ("spectral_estimators: multiple photons but only one frequency seen in %d band %d\n", xplasma->nplasma, n); /* Flag as a warning, so one can see if it is an issue */

      xplasma->pl_log_w[n] = -999;      //A very tiny weight
      xplasma->pl_alpha[n] = 999.9;     //Give alpha a value that will show up as an error
      xplasma->exp_w[n] = 0.0;  //Make sure that w is zero, s no chance of mucking up ionization balance
      xplasma->exp_temp[n] = 1e99;      //Give temp a value that will show up as an error
      xplasma->spec_mod_type[n] = SPEC_MOD_FAIL;        //This tells the code that we have failed to model the spectrum in this band/cell/
    }


    else                        //All is well, lets move on to try and model the photon distribution
    {

      /* The next lines check and assign band limits.
         If the min/max frequency in a cell is within (fmax-fmin)/(sqrt(nphot)) of the end of
         a band, we say that the photons fill the band to that end - i.e. the fact we didnt
         see the minimum frequency is just because of photon numbers. If on the other hand, one
         end of the band is 'surprisingly' empty then we wasume this is because absolutely
         no photons are here - its probably an edge so we should modify the model bands. */

      dfreq = (geo.xfreq[n + 1] - geo.xfreq[n]) / sqrt (xplasma->nxtot[n]);     //This is a measure of the spacing between photons on average
      if ((xplasma->fmin[n] - geo.xfreq[n]) < dfreq)    //If true, this check suggests that there are no edges
      {
        spec_numin = geo.xfreq[n];      //Use the photon generation band edge to set the lower frequency band for the model
      }
      else
      {
        spec_numin = xplasma->fmin[n];  //There may be an edge, use the lowest observed photon frequency for the lower nu band in the model
      }
      if ((geo.xfreq[n + 1] - xplasma->fmax[n]) < dfreq)        //Repeat above but for upper band edge
      {
        spec_numax = geo.xfreq[n + 1];
      }
      else
      {
        spec_numax = xplasma->fmax[n];
      }

      xplasma->fmin_mod[n] = spec_numin;        //This is the low frequency limit of any model we might make
      xplasma->fmax_mod[n] = spec_numax;        //This is the high frequency limit of any model we might make
      lspec_numax = log10 (spec_numax);
      lspec_numin = log10 (spec_numin);
      spec_numean = xplasma->xave_freq[n];
      j = xplasma->xj[n];

      /* Try to find the exponent of a power law model that fits the cell spectrum */

      pl_alpha_min = -0.1;      /*Lets just start the search around zero */
      pl_alpha_max = +0.1;

      while (pl_alpha_func_log (pl_alpha_min) * pl_alpha_func_log (pl_alpha_max) > 0.0) //Bracket the correct value
      {
        pl_alpha_min = pl_alpha_min - 1.0;
        pl_alpha_max = pl_alpha_max + 1.0;
      }

      if (isfinite (pl_alpha_func_log (pl_alpha_min)) == 0 || isfinite (pl_alpha_func_log (pl_alpha_max)) == 0) //We can't backet alpha!!!!
      {
        Error ("spectral_estimators: Alpha cannot be bracketed (%e %e)in band %i cell %i- setting w to zero\n", pl_alpha_min, pl_alpha_max,
               n, xplasma->nplasma);
        xplasma->pl_log_w[n] = -999.0;
        xplasma->pl_alpha[n] = -999.0;  //Set this to a value that might let us diagnose the problem
        plflag = -1;            //Discount a PL model
      }

      else                      //We have bracketed alpha
      {
        pl_alpha_temp = zero_find (pl_alpha_func_log2, pl_alpha_min, pl_alpha_max, 0.00001, &ierr);     //find the actual value of alpha that matches our mean frequency
        if (ierr)
        {
          Error ("spectral_estimators: alpha not bracketed\n");
        }

        /* This next line computes the PL weight using an external function. */

        pl_w_temp = pl_log_w (j, pl_alpha_temp, lspec_numin, lspec_numax);

        if ((isfinite (pl_w_temp)) == 0)        //Catch problems
        {
          Error
            ("spectral_estimators: New PL parameters (%e) unreasonable, using existing parameters. Check number of photons in this cell\n",
             pl_w_temp);

          plflag = -1;          // Dont use this model
          xplasma->pl_log_w[n] = -999.0;
          xplasma->pl_alpha[n] = -999.0;
        }
        else                    //All is well, assign model parameters to the plasma structure - we still need to work out if this is the *best* model
        {
          xplasma->pl_alpha[n] = pl_alpha_temp;
          xplasma->pl_log_w[n] = pl_w_temp;
        }
      }

      /* We will start the search around the temperature that we know will yield a sensible answer  */
      exp_temp_min = ((PLANCK * spec_numax) / (BOLTZMANN)) * 0.9;
      exp_temp_max = ((PLANCK * spec_numax) / (BOLTZMANN)) / 0.9;       /* NSH 131107 - and the same for the maximum temp */

      /* Bracket the temperature of an exponential model */

      while ((exp_temp_func (exp_temp_min) * exp_temp_func (exp_temp_max) > 0.0)
             && ((exp_temp_func (-1.0 * exp_temp_min) * exp_temp_func (-1.0 * exp_temp_max) > 0.0)))
      {
        /* In this case we are going to get errors since the temperature is too low
           give a result in the exponential, and we will divide by zero */
        if ((PLANCK * spec_numax) < (100.0 * BOLTZMANN * exp_temp_min * 0.9))
        {
          exp_temp_min = exp_temp_min * 0.9;    // Reduce the mininmum temperature, only if we will not end up with problems
        }
        exp_temp_max = exp_temp_max * 1.1;      // The maximum temperature can go up forever with no fear of numerical problems
      }



      if (isfinite (exp_temp_func (exp_temp_min)) == 0 || isfinite (exp_temp_func (exp_temp_max)) == 0)
      {
        Error ("spectral_estimators: Exponential temperature cannot be bracketed (%e %e) in band %i - setting w to zero\n", exp_temp_min,
               exp_temp_max, n);
        xplasma->exp_w[n] = 0.0;
        xplasma->exp_temp[n] = -1e99;
        expflag = -1;           //Discount an exponential model
      }

      else
      {
        /* But first see if we have a positive or negative solution. The temperatures are positive at the moment,
           if it was the negatives that worked, change the sign of the temperatures - we also need to swap min and max */
        if (exp_temp_func (-1.0 * exp_temp_min) * exp_temp_func (-1.0 * exp_temp_max) < 0.0)
        {
          exp_temp_store = -1.0 * exp_temp_min;
          exp_temp_min = -1.0 * exp_temp_max;
          exp_temp_max = exp_temp_store;
        }
        /* Solve for the effective temperature */

        exp_temp_temp = zero_find (exp_temp_func2, exp_temp_min, exp_temp_max, 0.00001, &ierr);
        if (ierr)
        {
          Error ("spectral_estimators: temp not bracketed\n");
        }

        /* Calculate the weight */
        exp_w_temp = exp_w (j, exp_temp_temp, spec_numin, spec_numax);


        if ((isfinite (exp_w_temp)) == 0)
        {
          Error ("spectral_estimators: New exponential parameters (%e) unreasonable, using existing parameters. Check number of photons in this cell\n", exp_w_temp);   //NSH 131108 - now a warning, this should no longer happen

          expflag = -1;         //discount an exponential model
          xplasma->exp_w[n] = 0.0;
          xplasma->exp_temp[n] = -1e99;
        }
        else                    //We have a reasonable exponential function model
        {
          xplasma->exp_temp[n] = exp_temp_temp;
          xplasma->exp_w[n] = exp_w_temp;
        }
      }

      /* compute standard deviations for exponential and power law models - these will be used to check the models */
      exp_sd = exp_stddev (xplasma->exp_temp[n], spec_numin, spec_numax);

      pl_sd = pl_log_stddev (xplasma->pl_alpha[n], lspec_numin, lspec_numax);

      Log_silent ("NSH in this cell %i band %i PL estimators are log(w)=%10.2e, alpha=%5.3f giving sd=%e compared to %e\n",
                  xplasma->nplasma, n, xplasma->pl_log_w[n], xplasma->pl_alpha[n], pl_sd, xplasma->xsd_freq[n]);

      Log_silent ("NSH in this cell %i band %i exp estimators are w=%10.2e, temp=%10.2e giving sd=%e compared to %e\n",
                  xplasma->nplasma, n, xplasma->exp_w[n], xplasma->exp_temp[n], exp_sd, xplasma->xsd_freq[n]);

      /*Compute the fractionasl errors in standard dev */

      exp_sd = fabs ((exp_sd - xplasma->xsd_freq[n]) / xplasma->xsd_freq[n]);
      pl_sd = fabs ((pl_sd - xplasma->xsd_freq[n]) / xplasma->xsd_freq[n]);

      /* These commands decide upon the best model,
         based upon how well the models predict the standard deviation */
      if (expflag > 0 && plflag > 0)    //Both models are in the running - see which has the lowest error in stdev
      {
        if (exp_sd < pl_sd)
          xplasma->spec_mod_type[n] = SPEC_MOD_EXP;
        else
          xplasma->spec_mod_type[n] = SPEC_MOD_PL;
      }

      else if (plflag > 0)      //Only PL model in running, no point in testing for STDEV
      {
        xplasma->spec_mod_type[n] = SPEC_MOD_PL;
      }

      else if (expflag > 0)     //Only EXP model in running, no point in testing for STDEV
      {
        xplasma->spec_mod_type[n] = SPEC_MOD_EXP;
      }

      else
      {
        xplasma->spec_mod_type[n] = SPEC_MOD_FAIL;      //Oh dear, there is no suitable model - this should be an error

        Error ("No suitable model in band %i cell %i (nphot=%i fmin=%e fmax=%e)\n",
               n, xplasma->nplasma, xplasma->nxtot[n], xplasma->fmin[n], xplasma->fmax[n]);

        /* We will set the applicable frequency bands for the model to values that will cause errors if the model is used */
        xplasma->fmin_mod[n] = spec_numax;
        xplasma->fmax_mod[n] = spec_numin;
      }

      Log_silent ("NSH In cell %i, band %i, the best model is %i\n", xplasma->nplasma, n, xplasma->spec_mod_type[n]);
    }                           //End of loop that does things if there are more than zero photons in the band.

  }                             //End of loop over bands

  geo.spec_mod = TRUE;          //Tell the code we have a model

  return (0);

}

/**********************************************************/
/**
 * @brief      Function used in zbrent to work out best value of alpha for a power law model
 *
 *
 * @param [in] double  alpha  - the exponent of a power law model for J_nu
 * @return     Difference between the mean of a modelled photon spectrum and the observed one
 *
 * @details
 * The function is just subtracts the 'observed' mean frequency of the photons
 * which passed thruogh as cell from thr mean frequency of a power law spectral model
 * (which is computed externally). This is used in a zero finding function to work out
 * the best value of alpha.
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
pl_alpha_func_log (double alpha)
{
  double answer;
  answer = pl_logmean (alpha, lspec_numin, lspec_numax) - spec_numean;
  return (answer);
}

double
pl_alpha_func_log2 (double alpha, void *params)
{
  return (pl_alpha_func_log (alpha));
}

/**********************************************************/
/**
 * @brief      Computes the mean frequency of a power law with exponent alpha
 *
 * @param [in] double  alpha   The exponent of the power law
 * @param [in] double  lnumin   Miniumum frequency
 * @param [in] double  lnumax   Maximum frequency
 * @return     the mean frequency
 *
 * @details
 * This function computes the mean frequency of a power law disturbution.
 * It is essentially a ratio of two integrals \int(nuJnu)  / \int(Jnu)
 * For a power law this can be done algegraically, so we do that
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
pl_logmean (alpha, lnumin, lnumax)
     double alpha;
     double lnumin, lnumax;
{
  double k, answer, numerator, denominator, a, b, c, d;

  k = -1.0 * (lnumax * (alpha + 2.0));  //This is a prefector that scales everything to sensible numbers.

  a = pow (10.0, (k + (lnumax * (alpha + 2.0))));
  b = pow (10.0, (k + (lnumin * (alpha + 2.0))));
  c = pow (10.0, (k + (lnumax * (alpha + 1.0))));
  d = pow (10.0, (k + (lnumin * (alpha + 1.0))));

  numerator = (a - b) / (alpha + 2.0);  //The integral of nu J_nu
  denominator = (c - d) / (alpha + 1.0);        //The integral of J_nu
  answer = numerator / denominator;


  return (answer);
}

/**********************************************************/
/**
 * @brief      works out the scaling factor for a power law representation of a cell spectrum
 *
 * @param [in] double  j   The 'measured' mean intensity for a cell/band
 * @param [in] double  alpha   The exponent of the power law representation
 * @param [in] double  lnumin   Band lower frequency limit
 * @param [in] double  lnumax   Band upper frequency limit
 * @return     logw - the log of the scaling factor
 *
 * @details
 * This routine works out the scaling factor for a power law representation of a spectrum in
 * a cell - the W in J_nu=Wnu**alpha. Assuming alpha is known, it computes w for J=W x integral
 * from numin to numax of nu^alpha. For a power law model, the integral is easily done
 * numercially. We work in log space to allow huge numbers to work and not give errors.
 *
 * ### Notes ###
 *
 **********************************************************/

double
pl_log_w (j, alpha, lnumin, lnumax)
     double j;                  //the band limited spectral density
     double alpha;              //Computed spectral index for the cell
     double lnumin, lnumax;     //Range of frequencies we are considering
{
  double logw;                  //the answer
  double logk;                  //scaling prefactor to permit huge numbers to be dealt with
  double log_integral;          /*This will hold the unscaled integral of jnu from numin to numax */

  logk = -1.0 * lnumax * (alpha + 1.0); //Compute a temporary scaling factor to avoid numerical problems

  log_integral = log10 ((pow (10, (logk + (alpha + 1.0) * lnumax)) - pow (10, (logk + (alpha + 1.0) * lnumin))) / (alpha + 1.0));

  logw = log10 (j) + logk - log_integral;

  return (logw);
}

/**********************************************************/
/**
 * @brief      Work out the standard deviation of a power law spectral model
 *
 * @param [in] double  alpha   The exponent of the power law representation
 * @param [in] double  lnumin   Band lower frequency limit
 * @param [in] double  lnumax   Band upper frequency limit
 * @return     The standard deviation
 *
 * @details
 * This subroutine works out the standard deviation of a photon distribution modelled
 * 	by a power law. First we work out the integral of nu^2 jnu / int jnu.
 * 	Then we subtract the (mean/j)^2 and then take the square root.
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
pl_log_stddev (alpha, lnumin, lnumax)
     double alpha;              //Computed spectral index for the cell
     double lnumin, lnumax;     //Range of frequencies we are considering
{
  double answer;                //the answer

  double numerator;             /*This will hold the unscaled integral of nu^2 jnu from numin to numax */
  double denominator;           /*This is the unscaled integral of jnu from numin to numax */
  double mean;
  double k, a, b, c, d;

  k = -1.0 * (lnumax * (alpha + 3.0));  //This is a prefector that scales everything to sensible numbers.

  a = pow (10.0, (k + (lnumax * (alpha + 3.0))));       //Integral of nu^2 jnu
  b = pow (10.0, (k + (lnumin * (alpha + 3.0))));
  c = pow (10.0, (k + (lnumax * (alpha + 1.0))));       //Integral of jnu
  d = pow (10.0, (k + (lnumin * (alpha + 1.0))));

//printf ("a=%e,b=%e,c=%e,d=%e\n",a,b,c,d);

  numerator = (a - b) / (alpha + 3.0);
  denominator = (c - d) / (alpha + 1.0);

  answer = numerator / denominator;

  mean = pl_logmean (alpha, lnumin, lnumax);    //Get the mean
  answer = sqrt ((numerator / denominator) - (mean * mean));


  return (answer);
}

/**********************************************************/
/**
 * @brief      Function used in zbrent to work out best value of alpha for a power law model

 *
 * @param [in] double  exp_temp   The effective temperature of an expontential model for J_nu
 * @return     Difference between the mean of a modelled photon spectrum and the observed one
 *
 * @details
 * This function is used in zbrent to search for the value of 'temperature' which makes the
 * 	mean frequewncy from an expoential model (coputed in a subroutine) function equal to
 * the mean frequency computed for the cell from photon
 * 	packets. So we subtract the mean frequency at the end to make the value of the
 * 	function equal to zero at the correct temperature.
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
exp_temp_func (double exp_temp)
{
  double answer;
  answer = exp_mean (exp_temp, spec_numin, spec_numax) - spec_numean;
  return (answer);
}

double
exp_temp_func2 (double exp_temp, void *params)
{
  return (exp_temp_func (exp_temp));
}

/**********************************************************/
/**
 * @brief      Calculate the mean frequency of an exponential distribution
 *
 * @param [in] double  exp_temp   The effective temperature of the distribution
 * @param [in] double  lnumin   Miniumum frequency
 * @param [in] double  lnumax   Maximum frequency
 * @return     The mean frequency
 *
 * @details
 *  The mean of a frequency depedant mean inteisntiy is the integral of
 * nu j(nu) divided by the integral of j(nu).
 * 	The prefactor w is lost in the division, so this function should equal the
 * 	mean frequency. Since this is a simple algrebraic frequency distribution we
 * can integrate algebraically
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
exp_mean (exp_temp, numin, numax)
     double exp_temp;
     double numin, numax;
{
  double answer, numerator, denominator;
  double exp1;                  /* We supply a temperature, but actually we expect the correct function to be of the form e^-hnu/kt, so this will hold -1*h/kt */
  double emin, emax;            /*The exponential evaluated at numin and numax */
  double t1, t2;

  exp1 = (-1.0 * PLANCK) / (BOLTZMANN * exp_temp);

  t1 = exp1 * numin;
  t2 = exp1 * numax;

  emin = exp (t1);
  emax = exp (t2);

  numerator = ((emax / exp1) * ((exp1 * numax) - 1.0)) - ((emin / exp1) * ((exp1 * numin) - 1.0));      /*THe integral of nu e^(-hnu/kt) from numin to numax */
  denominator = emax - emin;    /*NSH 120817 This is the integral of (J(nu)=e^(-hnu/jkt)) */
  answer = numerator / denominator;     /*NB there is a factor of exp1 top and bottom which I have cancelled out */
  return (answer);
}



/**********************************************************/
/**
 * @brief      This program works out the scaling factor for an exponential representation of a spectrum
 *
 * @param [in] double  j   The 'measured' mean intensity for a cell/band
 * @param [in out] double  exp_temp  - the temperature of the exponential model of J_nu
 * @param [in] double  numin   Band lower frequency limit
 * @param [in] double  numax   Band upper frequency limit
 * @return     w - the scaling factor
 *
 * @details
 * We compute J/W= integral from numin to numax of e^(-nhu/kt). We then
 * obtain W by comparing this number to the measured mean intensity in the band
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
exp_w (j, exp_temp, numin, numax)
     double j;                  //the band limited spectral density
     double exp_temp;           //Computed effective temperature for the cell
     double numin, numax;       //Range of frequencies we are considering
{
  double w;                     //the answer

  double integral;              /*This will hold the unscaled integral of jnu from numin to numax */
  double exp1;                  /* We supply a temperature, but actually we expect the correct function to be of the form e^-hnu/kt, so this will hold -1*h/kt */
  double emin, emax;            /*The exponential evaluated at numin and numax */

  exp1 = (-1.0 * PLANCK) / (BOLTZMANN * exp_temp);


  emin = exp (exp1 * numin);
  emax = exp (exp1 * numax);
  integral = (emax - emin) / exp1;      /* This is the integral */
  w = j / integral;             /*THe scaling factor is the actual band limited J / the integral. */

  return (w);
}


/**********************************************************/
/**
 * @brief      This program works out the standard deviation
 *
 * @param [in out] double  exp_temp  - the temperature of the exponential model of J_nu
 * @param [in] double  numin   Band lower frequency limit
 * @param [in] double  numax   Band upper frequency limit
 * @return     The standard deviation of the mean intensity model
 *
 * @details
 * First we work out the integral of nu^2 jnu / int jnu.
 * 	Then we subtract the (mean/j)^2 and then take the square root.
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
exp_stddev (exp_temp, numin, numax)
     double exp_temp;           //Computed spectral index for the cell
     double numin, numax;       //Range of frequencies we are considering
{
  double answer;                //the answer
  double exp1;                  /* We supply a temperature, but actually we expect the correct function to be of the form e^-hnu/kt, so this will hold -1*h/kt */
  double emin, emax;            /*The exponential evaluated at numin and numax */
  double nmax, nmin;            /*The integral is a bit complicated, so we have split it up into value for numin and numax */
  double numerator, denominator;
  double mean;

  exp1 = (-1.0 * PLANCK) / (BOLTZMANN * exp_temp);
  emin = exp (exp1 * numin);
  emax = exp (exp1 * numax);

  nmax = emax * (((numax * numax) / (exp1)) - ((2 * numax) / (exp1 * exp1)) + ((2.0) / (exp1 * exp1 * exp1)));  /* Is the integral for nu^2 bnu exaluated for numax */
  nmin = emin * (((numin * numin) / (exp1)) - ((2 * numin) / (exp1 * exp1)) + ((2.0) / (exp1 * exp1 * exp1)));  /* Is the integral for nu^2 bnu exaluated for numin */
  numerator = nmax - nmin;      /*The definite integral */
  denominator = (emax - emin) / exp1;   /* This is the integral for nbu */
  mean = exp_mean (exp_temp, numin, numax);
  answer = sqrt ((numerator / denominator) - (mean * mean));


  return (answer);
}
