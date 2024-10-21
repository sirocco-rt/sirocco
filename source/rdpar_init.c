
/***********************************************************/
/** @file   rdpar_init.c
 * @author  ksl
 * @date   January, 2019
 * @brief  Routines for setting up choices for a rdchoice 
 * command in the situation where only some of the options
 * that are potentially available are actually made available
 *
 * The routines are currently only used in connection with
 * choosing what type of spectra (bb, power law, etc) to be
 * used for different sources (e.g the disk, the central object).
 *
 * In principle, one could use the mechanisms here to make
 * sure that a rdchoice command was set up properly.  
***********************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "log.h"
#include "atomic.h"
#include "sirocco.h"

/* the form of the structure which contains
 * the map from words to values
    struct  xchoices {
        char choices[10][LINELENGTH];
        int  vals[10];
        int  n;
    } zz_spec; 
*/



/* this is a routine which is intended to initialize a structure 
 * used in the special instance when we are going to use the rdchoice
 * command, but when we want different options to be present in different
 * calls.  init_choices only needs to be called once
 */

int xinit_choices = 0;


/**********************************************************/
/**
 * @brief   initilize a structure that provides a connection
 * between the string names for an input accessed by rdchoice
 * to an integer.
*
 * @return Alwsys returns 0
 *
 *
 * @details     
 *
 * This routine is associated with reading data in via the rdchoice
 * commaind, which lets the user type in a string to specify an
 * option in the program.  The structue contains all of the 
 * possible rdchoice options for a given command.
 *
 *
 * ### Notes ###
 * 
 * Currently, one one structure is intialized, one which is 
 * associated with spectral types.
 *
 * However, it would be straightforward to add additional
 * structures for different rdchoice variables, by adding 
 * addional structures.
 *
**********************************************************/


int
init_choices ()
{
  /* Initialize the structure that contains all of the types of possible radiation types */

  char *xchoices[] = { "bb", "uniform", "power", "cloudy", "brems", "none", "models", "mod_bb", "mono" };
  int xvals[] =
    { SPECTYPE_BB, SPECTYPE_UNIFORM, SPECTYPE_POW, SPECTYPE_CL_TAB, SPECTYPE_BREM, SPECTYPE_NONE, SPECTYPE_MODEL, SPECTYPE_BB_FCOL,
    SPECTYPE_MONO
  };
  int num_choices = 9;          //Should match the length of xchoices and xvals above. must be <= MAX_RDPAR_CHOICES in sirocco.h 

  if (xinit_choices)
    return (0);

  int i = 0;
  for (i = 0; i < num_choices; i++)
  {
    strcpy (zz_spec.choices[i], xchoices[i]);
    zz_spec.vals[i] = xvals[i];
  }
  zz_spec.n = num_choices;

  /* To initialize more than one sets of choices one needs essentially to repeat the lines above
   * making sure that the names of the input arrays changed.  Alternatively, one could create
   * a routine for each structure one wanted to initialize
   *
   * If one had large numbers of variables to intialize by this method, one could put all oft he
   * relevant data into a file
   */

  xinit_choices = 1;


  return (0);
}



/**********************************************************/
/**
 * @brief   construct a string that contains the allowed set
 * of choices posed in a rdchoice command
 * @param  [in] char *question  A string containing the question
 * @param  [out] char *choices  The choices that correspond to the 
 * options in the question
 * @param  struct rdpar_choices *qstruct  The structure that contains
 * the full set of choices
*
 * @return Always returns 0
 * 
 *
 * @details     
 *
 * This routine is general, but was implemented specifically to deal
 * with spectral types, where the spectral choices for the disk may
 * be different from those for a BH as a central object.
 *
 * A typical question to be parsed will be of the form: 
 *
 * disk.spec_type(bb,models)
 *
 * wheres there are other types of spectra, power, brems, that 
 * Python can produce.
 *
 * This routine reads the quesion, and returns a string that
 * contains the numbers that correspond to the choices bb, and models
 *
 * Following this, rdchoice can be called with the allowed choices.
 *
 *
 * ### Notes ###
 *
 * This routine does not do any kind of minimum match.  The choices,
 * items within the parenthese must exactly match those in the
 * structure.
 *
 * Init_choices must have been called, prior to get_choices.   
 * If not, Python exits with an error.
 * 
**********************************************************/
int
get_choices (question, choices, qstruct)
     char *question;
     char *choices;
     struct rdpar_choices *qstruct;
{
  char cur_choices[MAX_RDPAR_CHOICES][LINELENGTH];
  int cur_values[MAX_RDPAR_CHOICES] = { -999 };
  int cur_num;
  char cur_string[LINELENGTH];
  char xcur_string[LINELENGTH];


  if (xinit_choices == 0)
  {
    Error ("get_choices: init_choices needs to have been called before get_choices.  Programming error\n");
    exit (0);
  }

  /*
   * A typical question will be disk.spec_type(bb,power), that is to
   * say * it will not include all of the possibileities
   */


  /* Now match the current possibilities to the integer choices */

  int n, nstart, nstop;
  char dummy[LINELENGTH];
  nstart = nstop = 0;
  for (n = 0; n < strlen (question); n++)
  {
    if (question[n] == '(')
    {
      nstart = n;
    }
    if (question[n] == ')')
    {
      nstop = n;
    }
  }


  snprintf (dummy, nstop - nstart, "%s", &question[nstart + 1]);

  strcpy (dummy, " ");
  strcpy (dummy, question);
  int i = 0;
  int j = 0;
  for (i = 0; i < strlen (dummy); i++)
  {
    if (dummy[i] == '(')
    {
      dummy[i] = ' ';
    }
    if (dummy[i] == ')')
    {
      dummy[i] = ' ';
    }
    if (dummy[i] == ',')
    {
      dummy[i] = ' ';
    }
  }

  cur_num = sscanf (dummy, "%*s %s %s %s %s %s %s %s %s %s %s",
                    cur_choices[0], cur_choices[1],
                    cur_choices[2], cur_choices[3],
                    cur_choices[4], cur_choices[5], cur_choices[6], cur_choices[7], cur_choices[8], cur_choices[9]);


  /* Now we have to find out the values associated with these choices */


  for (i = 0; i < cur_num; i++)
  {
    for (j = 0; j < qstruct->n; j++)
    {
      if (strncmp (cur_choices[i], qstruct->choices[j], strlen (qstruct->choices[j])) == 0)
      {
        cur_values[i] = qstruct->vals[j];
        break;
      }
    }
  }
  /*
   * If we have gotten to this point, we can construct the integer
   * string that represents the choices that are currently available
   */

  sprintf (cur_string, "%d", cur_values[0]);
  for (i = 1; i < cur_num; i++)
  {
    sprintf (xcur_string, "%.100s,%d", cur_string, cur_values[i]);
    sprintf (cur_string, "%.200s", xcur_string);
  }

  strcpy (choices, cur_string);

  return (0);
}
