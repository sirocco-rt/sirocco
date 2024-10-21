
/***********************************************************/
/** @file  rdpar.c
 * @author ksl
 * @date   February, 2018
 *
 * @brief  A general purpose interface for obtaining input data from the command line 
 * and/or from a file
 *
 *  This file contains various routines which implement a relatively straightforward
 *  keyword-based interface to get data either from the user or a command line.
 *  
 *  The basic routines are as follows:
 *  	- rdstr(question,answer)	gets a single contiguous string 
 *  	- rdchar(question,answer)	gets a single character
 *  	- rdchoice(question,answer) gets one of a number of allowed string answers
 *  	- rdint(question,answer)	gets a single integer
 *  	- rdflo(question,answer)	gets a single precision floating point number
 *  	- rddoub(question,answer)	gets a double precision floating point number
 *  	- rdline(question,answer)	gets an entire line
 *  	- rdpar_comment(string)   Adds a comment to any output.pf files
 *  
 *  
 *  
 *  The file also contains routines to control how rdpar operates:
 *  
 *  - opar(infilename)   --- allows one to read from a file instead of from the
 *  	command line.  If opar is not called then interactive mode is 
 *  	assumed.  If the file does not exist, the routines drop into interactive
 *  	mode. 
 *  - cpar(outfilename)  ---- causes the temporary file being created by the various
 *  	routines to be copied to outfilename.  If there was a file with
 *  	that name previously it will be copied to filename.old.   If cpar
 *  	is not called, then the datat which is aquired with the readpar
 *  	routines will be found in "tmp.rdpar"
 *  - int get_root(root,total) -- strips off the trailing characters in a .pf file
 *  	name.
 *  - int rdpar_set_mpi_rank(rank) -- passes the rank of the parallel process to rdpar
 *  - int rdpar_set_verbosity  -- Control how much diagnostic information is printed out
 *  
 *  The file also contains a subroutine to put out a message of
 *  one or more lines.
 *  
 *  message(string)
 *  
 *  Arguments:
 *  
 *  	For the io routines the arguments are
 *  
 *  	char question[];	The string printed to the screen to prompt the user
 *  					or read from the file and matched
 *  	*answer			The value returned; the type depends on the call.
 *  		
 *  
 *  Returns:
 *  	The individual io routines return 
 *  		1 or more for a successful read
 *  		EOF if it receives done or EOF 
 *   
 *  Description:	
 *  	
 *  	In interactive mode, the programs print out the question and
 *  	expects an answer.  Only one integer (or whatever) is read at a time.
 *  	The calls are generally of the form
 *  		- rdxxx(question,whatever)
 *  		- char question[]
 *  		- whatevertype *whatever
 *
 *  	One can answere in several ways:
 *  	-To provide a new value, enter it
 *  	-To keep the current value, respond with \n
 *  	-To issue a system level command, respond with ! and the command
 *  	-To exist the program entirely, type ^D
 *  	
 *  	In noninteractive mode, the program reads in a file into the input structure.  
 *  	As inputs are requested, it reads this structure  to match the requested
 *  	parameter to a line in the structure from which it gets the parameter value.
 *  	The keywords can be in any order but once a line has been used, it is not
 *  	reused.  There can be multiple lines with the same keyword, as for example
 *  	is the case in Python when one wants to extract spectra from multiple angles
 *  	
 *  	In noninteractive mode, if the second string in the file is "$" (or 
 *  	there is no second string, it will
 *  	ask the user that question interactively and then continue to read
 *  	other answers from the file.  (Thus if one wants to run the program
 *  	a bunch of times, changind only one or two of the values in the
 *  	program, one can include most of the values in the file, but leave
 *  	the values which will change as $)
 *  
 *  	In noninteractive mode, if there is no matching keyword, the routine will
 *  	switch to interactive mode for that keyword, and try to get the value of
 *  	the parameter from the user.  It will then switch back to nonintercative mode.
 *  	interactive mode and try that.
 *
 *  	Since keywords occaassionaly change their names, the rdpar routines have a
 *  	synonyms facility that makes it possible (at least for a time) to use 
 *  	old parameter files.
 *  
 *  	The routine keeps track of all of the input variables as they are entered
 *  	or modified, and when cpar is called and at various other times can write
 *  	the answers to a file, including a new .pf file, generally called rootname.out.pf.
 *
 *      The output files use current parameter file names.
 *
 *  		
 *  ###Notes###
 *  
 *  	The original rdpar was written at SAO for handling the reduction programs
 *  	assocated with Einstein data.  C and Fortran versions of the code exist.
 *  	The SAO version of rdpar always works from parameter (".pf" files".  
 *  	When I began working with c programs on Macs (and Atari's), I wrote
 *  	a set of routines which I called "irdpar" which were intended to provide
 *  	a simple interactive version of rdpar, i.e. a very readable way to get
 *  	single input values into a program.  Eventually, it became clear that I
 *  	was writing programs that would benefit from a file i/o interface as
 *  	well and so this version of rdpar has been written.  
 *  	
 *  	This version of rdpar does not have all of the advantages of the old
 *  	rdpar.  It does allow one to use $ in file to switch to interactive mode
 *  	for that particular input.  But this version of rdpar cannot read more 
 *  	than one variable from a line.  Finally, unlike rdpar this routine expects 
 *  	variables to be in the file in the prescribed manner.  It bails out
 *  	if the variables are not in the prescibed order.
 *  	
 *  	Note that you can have extra lines in the input file.  These will simply
 *  	be copied to the output file.  But if you are thinking of using several
 *  	programs which have some of the same keywords, be sure to put them
 *  	in the same order in both programs.  It is OK to have 
 *  		- key1
 *  		- newkey
 *  		- key2
 *  	
 *  	where newkey is used by the second program but not the first.  But
 *  	it is not a good idea to look for key2 before key1 in the second program.
 *  	
 *  	If you are writing programs using rdpar, and you seem to be having
 *  	problems with the inputs, the usual problem is that you did not pass
 *  	a pointer to rdxxx.
 *  	
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#include "log.h"
//OLD #include "strict.h"

#define LINELEN		256
#define	OLD		100
#define	NORMAL		1
#define REISSUE		-199

FILE *rdin_ptr, *rdout_ptr;     /* Pointers to the input and output files */
int rdpar_stat = 0;             /*rdpar_stat=0 initially, 
                                   rd_stat=1 while reading from stdin,
                                   rdpar_stat=2 if reading from file */
int verbose = 1;                /* verbose = 1 -> printout all errors
                                   verbose=2 -->suppresses some error messages */

char current_filename[LINELEN];

int rd_rank = 0;                // rank of mpi process, set to zero

/* Array to record the actual values of rdpar accepted by rdpar.  All values are 
 * stored as strings */
#define  MAX_RECORDS   500

int rdpar_nrec = 0;

int rdpar_warning = 0;          // A variable used to allow one to print out a warning the first time oen moves to interactive mode

struct rdpar
{
  char name[LINELEN];
  char value[LINELEN];
}
rdpar_record[MAX_RECORDS];

/* Struture to hold the input rdpar file */

#define UNUSED 0
#define USED  1

int rdpar_cursor = 0;
int rdpar_ntot = 0;             // total number of raw input lines
int rdpar_choice = 0;           // 0 if this is not rdchoice; 1 if so
struct rdpar_raw
{
  char line[LINELEN];           // Raw input line 
  int icheck;                   // Check off to indicate a line has been used
}
input[MAX_RECORDS];


int strict = 0;                 // Initialize to a value that indicates everyting is OK



/**********************************************************/
/** 
 * @brief      Open a parameter file for reading
 *
 * @param [in] char  filename[]   string giving the name of the parameter file
 * @return     A status indicating whether one is 
 * reading an existing parameter file or proceegin interactivley
 *
 *  
 * Open and read (or try to read) a parameter file and the information
 * in the input structure
 *
 * ###Notes###
 *
 * If a parameter file is currently opened, the routine prints an
 * error message, but continues to read from the current file. If one
 * wants to read from a new file, the existing one must be closed
 * with cpar.
 *
 **********************************************************/

int
opar (filename)
     char filename[];
{
  FILE *fopen (), *tmp_ptr;
  int rdpar_init ();

  /* Check that an input file is not currently open */
  if (rdpar_stat == 2)
  {
    printf ("Error: opar: A file for reading input is already open.\n");
    return (rdpar_stat);
  }

  /*Open a temporary output file */
  rdpar_init ();
  if ((tmp_ptr = fopen (filename, "r")) == NULL)
  {
    printf ("Error: opar: Could not open filename %s\n", filename);
    printf ("               Proceeding in interactive mode\n");
    return (rdpar_stat);
  }
  else
  {
    rdin_ptr = tmp_ptr;
    rdpar_stat = 2;             /* implies we are now trying to read from a file */
    while (fgets (input[rdpar_ntot].line, LINELEN, rdin_ptr) != NULL)
    {
      input[rdpar_ntot].icheck = UNUSED;
      rdpar_ntot++;
    }
    fclose (rdin_ptr);
  }
  strcpy (current_filename, filename);


  return (rdpar_stat);
}






/**********************************************************/
/** 
 * @brief      open another parameter file and append information from that parameter 
 * file to a previous one
 *
 * @param [in] char  filename[]   name of the file
 * @return     A status indicating whether the file was read
 *
 * This routines reads an auxiliary .pf file and appends data to the primary one,
 * or more correctly appends data to the input data structure
 *
 * ###Notes###
 *
 *  This is intended for the case where one would like to use a second parameter
 * 	file for extra parameters or to add a parameter to an old file. 
 *
 * 	As written one needs to open and read the primary file first
 *
 * 	Note that cpar, will write all of the parameters out to the .out.pf file
 * 	in the order that they were used.  Paramerters that were not used are not
 * 	written out.
 * 
 *
 **********************************************************/

int
add_par (filename)
     char filename[];
{
  int rdpar_init ();

  /* Check that an input file is not currently open */
  if (rdpar_stat != 2)
  {
    printf ("Error: add_par: opar needs to be called before add_par\n");
    return (rdpar_stat);
  }

  if ((rdin_ptr = fopen (filename, "r")) == NULL)
  {
    printf ("Error: add_par: Could not additional file %s\n", filename);
    return (rdpar_stat);
  }

  while (fgets (input[rdpar_ntot].line, LINELEN, rdin_ptr) != NULL)
  {
    input[rdpar_ntot].icheck = UNUSED;
    rdpar_ntot++;
  }
  fclose (rdin_ptr);


  return (rdpar_stat);
}




/**********************************************************/
/** 
 * @brief      This routine finishing writing out a new parameter file with 
 * all of the parameters that were actually used.
 *
 * @param [in] char  filename[]   Name of the new parameter file
 * @return     NORMAL
 *
 * While input parameters are being processed they are written to 
 * a file tmp.rdpar. This routine usually just closes tmp.rdpar and
 * renames it to rootname.out.pf
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
cpar (filename)
     char filename[];
{
  char old_filename[LINELEN];

  if (rd_rank != 0)
    return 0;

  fclose (rdout_ptr);
  rdpar_stat = 0;


  strcpy (old_filename, filename);
  strcat (old_filename, ".old");

  if (remove (old_filename) != 0 && verbose)
    printf ("Old file %s did not exist\n", old_filename);

  if (rename (filename, old_filename) != 0 && verbose)
    printf ("Could not rename %s to %s\n", filename, old_filename);

  if (rename ("tmp.rdpar", filename) != 0 && verbose)
    printf ("Could not rename %s to %s", "tmp.rdpar\n", filename);

  return (NORMAL);
}

/* Despite its name, this routine is internal to the other rdpar parts of rdpar. It
   does however set put rdpar into a known mode and set up the output file */


/**********************************************************/
/** 
 * @brief      open a temporary file to store inputs as they are processed.
 *
 * @return     NORMAL, unles one cannot get a pointer to the file
 *
 * Simply open a ptr to tmp.rdapr
 *
 * ###Notes###
 *
 * Storing information in a temporaray file prevents overwriting 
 * a permannt one before one is ready.
 *
 **********************************************************/

int
rdpar_init ()
{
  FILE *fopen ();
  rdin_ptr = stdin;             /* Initialize rdin_ptr to standard input */
  if ((rdout_ptr = fopen ("tmp.rdpar", "w")) == NULL)
  {
    printf ("Error: rdpar_init: Problem opening tmp.rdpar\n");
    exit (1);
  }
  rdpar_stat = 1;
  strcpy (current_filename, "tmp.rdpar.out");
  return (NORMAL);
}



/**********************************************************/
/** 
 * @brief      process a line of input regardless of whether 
 * it comes from a file or the command line
 *
 * is the routine called by all of the individual
 * 	routines, like rdpar for processing a line of input whether
 * 	it comes from the command line or from a .pf file
 *
 * @param [in] char  question[]    The keyword and or question 
 * @param [in, out] char  dummy[]   The suggested answer
 * @return     When the routine returns NORMAL or OLD the input has been 
 * 	successfully processed
 *
 *
 * ###Notes###
 *
 * At present this is merely a steering routine, which calls
 * the routines that process an input variable from either 
 * the command line or from the input file
 *
 **********************************************************/

int
string_process (question, dummy)
     char question[], dummy[];
{

  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */


  if (rdpar_stat == 1)
  {                             /* Then we are reading from stdin */
    return (string_process_from_command_line (question, dummy));
  }


  else
  {                             /*Then read from the file */

    return (string_process_from_file (question, dummy));

  }
}




/**********************************************************/
/** 
 * @brief      Ask for and process one input from the command line 
 *
 * @param [in ] char  question[]   parameter for which an answer is sought
 * @param [in, out] char  dummy[]   current value of parameter which the user may change
 * @return     When the routine returns NORMAL or OLD the input has been 
 * 	successfully processed
 *
 * The routine writes out a question and possible answer to the command
 * line and records the answer in the charcater string dummy
 *
 *  - If the user eneters a carriage return, dummy is unchanged
 *  - If the user enters a new value, that value is returned as string
 *  - If the user enters an ! followed by one or more words, those words will
 *  be issued as a system command, and another request for an answer generated
 *  - If the user enters ^D the program will exit
 *
 *
 *
 * ###Notes###
 *
 * This routine is a low level routine called by string_process
 *
 **********************************************************/

int
string_process_from_command_line (question, dummy)
     char question[], dummy[];
{
  char tdummy[LINELEN];
  fprintf (stderr, "%s (%s) :", question, dummy);
  fflush (stderr);
  strcpy (tdummy, "");
  if (fgets (tdummy, LINELEN, stdin) == NULL)
  {
    printf ("Exiting since rdpar got EOF in interactive mode\n");
    exit (1);
  }
  else if (tdummy[0] == '\n')
  {                             //Use the current value

    printf ("%s %s\n", question, dummy);
    fprintf (rdout_ptr, "%-30s %20s\n", question, dummy);
    rdpar_store_record (question, dummy);
    return (OLD);
  }
  else if (tdummy[0] == '!')
  {
    if (system (&tdummy[1]) == -1)      /* Send a command to the system */
      Error ("string_process_from_command_line: '%s' returned error code", tdummy[1]);

    return (REISSUE);
  }
  else if (strncmp (tdummy, "done", 4) == 0)
    return (EOF);
  else
  {                             //Use the value input from the command line

    strcpy (dummy, tdummy);
    if (rdpar_choice == 0)
    {
      fprintf (rdout_ptr, "%-30s %20s", question, dummy);
    }
    rdpar_store_record (question, dummy);
    return (NORMAL);
  }
}




/**********************************************************/
/** 
 * @brief      locate a keyword in the parameter file and retrieve the string associated with it
 *
 *  @param [in ] char  question[]   parameter for which an answer is sought
 *  @param [in, out] char  dummy[]   current value of parameter which the user may change
 *
 * @return     When the routine returns NORMAL or OLD the input has been 
 * 	successfully processed
 *
 *  When a parameter file is opened the data in the parameter file
 *  are read into a structure called input.   This routine 
 *  locates the appropriate line in the input structure 
 *  and returns the value as a string 
 *
 * ###Notes###
 *
 *  When a keyword is missing the routine checks if there would be
 *  a match if both strings are converted to lower case  (This change
 *  was made to facilitate standardization of keywords).  If that
 *  does not work it checks to see if 
 *  the keyword corresponds to an name that was changed by calling
 *  the routine synonyms.
 *
 *  If no match if foundl the rotine invokes the interactive version
 * 	of rdpar, but it does not force all future inputs to come from
 * 	the command line version, as was the case for version of rdpar
 * 	prior to 16sept.
 *
 * 	If the value for the keyword in the input list is a $, then
 * 	the user will be queried for that value interactively.
 *
 **********************************************************/

int
string_process_from_file (question, dummy)
     char question[], dummy[];
{

  char firstword[LINELEN], secondword[LINELEN];
  char *line, *fgets_rc;
  char *ccc, *index (), *fgets ();
  int nwords = 0;               // Initialise to avoid warning
  int wordlength;
  char xfirstword[LINELEN], xquestion[LINELEN];
  int i;

  for (rdpar_cursor = 0; rdpar_cursor < rdpar_ntot; rdpar_cursor++)
  {
    if (input[rdpar_cursor].icheck == USED)
    {
      continue;
    }

    line = input[rdpar_cursor].line;

    // Parse the first two words of the input line

    strcpy (firstword, "");
    strcpy (secondword, "");
    nwords = sscanf (line, "%s %s", firstword, secondword);

    wordlength = strlen (firstword);
    if (nwords < 2 || wordlength == 0)
    {
      continue;                 // Read the next line in the file
    }

    // Strip off everthing in paren
    if ((ccc = index (firstword, '(')) != NULL)
    {
      wordlength = (int) (ccc - firstword);
      if (wordlength == 0)
      {
        continue;
      }
    }


    // If there is a match go to the section that handles the various possibilities

    if (strncmp (question, firstword, wordlength) == 0)
    {
      break;                    // We have matched the keywords
    }



    /* Check for synonyms - That is look for keywords in the file we are reading 
     * that have been replaced by a new keyword
     */

    if (is_input_line_synonym_for_question (question, line))
    {
      strict = 1;
      Error ("Had to parse a synonym. Program will stop after writing out a new parameter file\n");
      break;
    }

    strncpy (xquestion, question, wordlength);
    strncpy (xfirstword, firstword, wordlength);
    for (i = 0; i < wordlength; i++)
    {
      xquestion[i] = tolower (xquestion[i]);
      xfirstword[i] = tolower (xfirstword[i]);
    }


    if (strncmp (xquestion, xfirstword, wordlength) == 0)
    {
      break;                    // We have matched the keywords
    }




    // Print a warning if the lines do not match
    if (verbose > 1 && rd_rank == 0)
      printf ("Warning:  Question (%s) does not match word (%s) in file. Continuing!\n", question, firstword);
  }

  // Handle the EOF since we were not successful in identifying the keyword
  if (rdpar_cursor == rdpar_ntot)
  {
    // printf ("Error: string_proces: Unexpectedly reached EOF\n");
    if (rdpar_warning == 0)
    {
      printf
        ("Switching to interactive mode for this variable.\nAdditional inputs from the terminal may be required to complete the .pf file\n");
      rdpar_warning = 1;
    }
    return (string_process_from_command_line (question, dummy));
  }
  else
  {
    input[rdpar_cursor].icheck = USED;
    rdpar_cursor++;             // This is needed because we have already processed the earlier line
  }


  // At this point we have the correct line in the input list

  if (strncmp (secondword, "$", 1) == 0 || nwords == 1)
  {                             // This if for the case where one wants to read this value only from the command line
    fprintf (stderr, "%s (%s) :", question, dummy);
    fflush (stderr);

    strcpy (secondword, dummy);

    fgets_rc = fgets (dummy, LINELEN, stdin);
    if (!fgets_rc)
    {
      Error ("Input value is NULL or invalid\n");
      exit (1);
    }

    if (strcmp (dummy, "\n") == 0)
    {                           /* Store the provided value since \n */
      rdpar_store_record (question, secondword);
      if (rdpar_choice == 0)
      {
        fprintf (rdout_ptr, "%-40s   %s\n", question, secondword);
      }
    }
    else
    {                           /* Store the value received via the command line */

      rdpar_store_record (question, dummy);
      if (rdpar_choice == 0)
      {
        fprintf (rdout_ptr, "%-40s   %s\n", question, dummy);
      }
    }

    return (NORMAL);
  }
  else                          // This handles the situation where the variable is actually read from the rdpar file
    strcpy (dummy, secondword);
  rdpar_store_record (question, secondword);
  if (rdpar_choice == 0)
  {
    fprintf (rdout_ptr, "%-40s   %s\n", question, secondword);
  }
  return (NORMAL);
}




/**********************************************************/
/** 
 * @brief      store the parsed version of an input line into 
 * 	a structure so that it can be written out to a file
 * 	, generally
 * 	the header for one of the ascii files produced by sirocco
 *
 * @param [in] char *  name   name of the parameter
 * @param [in] char *  value   string containing the value for the parameter
 * @return     Always returns 0
 *
 * rdpar_store_record(name,value) adds a line to the structrue rdpar_record,
 * which records the accepted values for all of the variables.
 *
 * ###Notes###
 *
 * rdpar_save is the routine that writes out the parameters 
 * to a file.  
 *
 **********************************************************/

int
rdpar_store_record (name, value)
     char *name, *value;
{
  strcpy (rdpar_record[rdpar_nrec].name, name);
  strcpy (rdpar_record[rdpar_nrec].value, value);
  rdpar_nrec++;

  return (0);
}




/**********************************************************/
/** 
 * @brief      save the actual values input to rdpar to a file
 *
 * @param [in] FILE *  file_ptr   An open file
 * @return     Always returns 0
 *
 * saves the stored parameter file vectors to an already open 
 * 	file, e.g. a spectrum file.  It allows one to attach the
 * 	inputs variables to the file for future reference.
 *
 * ###Notes###
 *
 *  The input is a file_ptr because generally this is intended simply
 * 	to document the inputs that were used to produce a model
 * 	
 * 	In some cases additional information, e.g the spectrum is
 * 	also writen to the file
 *
 **********************************************************/

int
rdpar_save (file_ptr)
     FILE *file_ptr;
{
  int i;


  for (i = 0; i < rdpar_nrec; i++)
  {
    fprintf (file_ptr, "# Var	 %-40s	%s\n", rdpar_record[i].name, rdpar_record[i].value);
  }
  return (0);
}




/**********************************************************/
/** 
 * @brief      Add a comment to the header of the output file
 *
 * @param [in] char *  format   format string for what is to be writtne out
 * @param [in]   ...   Any number of additional variables that are parsed by the stirng
 * @return     The number of characters succesfully written out out
 *
 * This routine allows one to add comments to the root.out.pf file
 * which represents a cleaned up version of the input.pf file
 *
 * ###Notes###
 *
 * 
 *
 **********************************************************/

int
rdpar_comment (char *format, ...)
{
  va_list ap, ap2;
  int result = 0;

  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */


  va_start (ap, format);
  va_copy (ap2, ap);            /* ap is not necessarily preserved by vprintf */

  if (rd_rank == 0)
  {
    fprintf (rdout_ptr, "\n### ");
    result = vfprintf (rdout_ptr, format, ap2);
    fprintf (rdout_ptr, "\n");
  }
  va_end (ap);
  return (result);
}


/**********************************************************/
/** 
 * @brief      Send a string to stderr
 *
 * @param [in] char  string[]   The message to be seent
 * @return     NORMAL           
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
message (string)
     char string[];
{
  fprintf (stderr, "%s\n", string);
  fflush (stderr);
  return (NORMAL);
}




/**********************************************************/
/** 
 * @brief      process a input line where the answer is expected to be string
 * 
 * @param [in] char  question[]   the name of the parameter for which an answer is sought
 * @param [in, out] char  answer[]   The current/final value of the string
 * @return    A status, indicating whether the answering string was successuly captured.
 *
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
rdstr (question, answer)
     char question[], answer[];
{
  int query;
  char dummy[LINELEN];
  query = REISSUE;
  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */
  while (query == REISSUE)
  {
    sprintf (dummy, "%s", answer);
    query = string_process (question, dummy);
    if (query == NORMAL)
    {
      if (sscanf (dummy, "%s", answer) != 1)
      {
        printf ("Could not convert input (%s) to string. Try again\n", dummy);
        query = REISSUE;
      }
      if (rd_rank == 0 && verbose == 1)
        printf ("%s   %s\n", question, answer);
    }
  }
  return (query);
}



/**********************************************************/
/** 
 * @brief      process a input line where the answer is expected to be a single character
 *
 * @param [in] char  question[]   the name of the parameter for which an answer is sought
 * @param [in, out] char  answer[]   The current/final value of the parameter
 * @return    A status, indicating whether the answering string was successuly captured.
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
rdchar (question, answer)
     char question[];
     char *answer;
{
  int query;
  char dummy[LINELEN];
  query = REISSUE;
  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */
  while (query == REISSUE)
  {
    sprintf (dummy, "%c", *answer);
    query = string_process (question, dummy);
    if (query == NORMAL)
    {
      if (sscanf (dummy, " %1c", answer) != 1)
      {
        printf ("Could not convert input (%s) to character. Try again\n", dummy);
        query = REISSUE;
      }

      if (rd_rank == 0 && verbose == 1)
        printf ("%s   %1c\n", question, *answer);
    }
  }
  return (query);
}


/**********************************************************/
/** 
 * @brief      process a input line where the answer is expected to be an integer
 *
 * @param [in] char  question[]   the name of the parameter for which an answer is sought
 * @param [in, out] char  answer[]   The current/final value of the parameter
 * @return    A status, indicating whether the answering string was successuly captured.
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
rdint (question, answer)
     char question[];
     int *answer;
{
  int query;
  char dummy[LINELEN];
  query = REISSUE;
  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */
  while (query == REISSUE)
  {
    sprintf (dummy, "%d", *answer);
    query = string_process (question, dummy);
    if (query == NORMAL)
    {
      if (sscanf (dummy, "%d", answer) != 1)
      {
        printf ("Could not convert input (%s) to integer. Try again\n", dummy);
        query = REISSUE;
      }
      if (rd_rank == 0 && verbose == 1)
        printf ("%s	  %d\n", question, *answer);
    }
  }
  return (query);
}



/**********************************************************/
/** 
 * @brief      process a input line where the answer is expected to be a float
 *
 * @param [in] char  question[]   the name of the parameter for which an answer is sought
 * @param [in, out] char  answer[]   The current/final value of the string
 * @return    A status, indicating whether the answering string was successuly captured.
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
rdflo (question, answer)
     char question[];
     float *answer;
{
  int query;
  char dummy[LINELEN];
  query = REISSUE;
  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */
  while (query == REISSUE)
  {
    sprintf (dummy, "%g", *answer);
    query = string_process (question, dummy);
    if (query == NORMAL)
    {
      if (sscanf (dummy, "%e", answer) != 1)
      {
        printf ("Could not convert input (%s) to float. Try again\n", dummy);
        query = REISSUE;
      }
      if (rd_rank == 0 && verbose == 1)
        printf ("%s	  %e\n", question, *answer);
    }
  }
  return (query);
}


/**********************************************************/
/** 
 * @brief      process a input line where the answer is expected to be double precision number
 *
 * @param [in] char  question[]   the name of the parameter for which an answer is sought
 * @param [in, out] char  answer[]   The value of the parameter
 * @return    A status, indicating whether the answering string was successuly captured.
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
rddoub (question, answer)
     char question[];
     double *answer;
{
  int query;
  char dummy[LINELEN];
  query = REISSUE;
  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */
  while (query == REISSUE)
  {
    sprintf (dummy, "%g", *answer);
    query = string_process (question, dummy);
    if (query == NORMAL)
    {
      if (sscanf (dummy, "%le", answer) != 1)
      {
        printf ("Could not convert input (%s) to double. Try again\n", dummy);
        query = REISSUE;
      };
      if (rd_rank == 0 && verbose == 1)
        printf ("%s	  %e\n", question, *answer);
    }
  }
  return (query);
}


/**********************************************************/
/** 
 * @brief      process a input line where the answer is expected be an entire line with multipe words
 *
 * @param [in] char  question[]   the name of the parameter for which an answer is sought
 * @param [in, out] char  answer[]   The current/final value of the line   
 * @return    A status, indicating whether the answering string was successuly captured.
 *
 *
 * ###Notes###
 *
 * This is used by setup_reverb
 *
 **********************************************************/

int
rdline (question, answer)
     char question[];
     char answer[];
{
  int query;
  char dummy[LINELEN];
  query = REISSUE;
  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */
  while (query == REISSUE)
  {
    sprintf (dummy, "%s", answer);
    query = string_process (question, dummy);
    if (query == NORMAL)
    {
      strcpy (answer, dummy);
      if (rd_rank == 0 && verbose == 1)
        printf ("%s	  %s\n", question, answer);
    }
  }
  return (query);
}



/**********************************************************/
/** 
 * @brief      find the integer value that corresponds to an input word
 *
 * @param [in] char  *word  The input string that we want to match                       
 * @param [in] char  *string_choices A comma separated string containing the possible choices
 * @param [in] char  *string_values    A comma separated string containing integers that correspond to the input strign
 * @param [out] char  *string_answer   The complete string version of the answer                                        
 * @return    The integer represents the choice indicated by the input string.           
 *
 * ###Notes###
 *
 * This routine called by rdchoice expects an input keyword that should
 * match one of several possible string choices.  string_choices and string_values
 * are two strings that both contain comma separated words or integers.  
 *
 * And example of the two strings are
 *
 * * string_cohoices "raw,medium,well"
 * * string_values  "1,8,16"
 *
 * As indicated above the two strings have to be in some sense parallel.
 *
 * We want to match *word to one of these values.  
 *
 * The routine splits the two comma separated strings into two parallel arrays.
 * It then performs a mimimum match between the input word and the version of the chooices.
 *
 * It returns the integer that corresponds to the one which is matched.
 *
 *
 *
 **********************************************************/


#define MAX_CHOICES 10
int
string2int (word, string_choices, string_values, string_answer)
     char *word;
     char *string_choices;
     char *string_values;
     char *string_answer;
{
  int i;
  int nchoices;
  char xs[MAX_CHOICES][LINELEN];
  int xv[MAX_CHOICES];
  char choices[LINELEN];
  char values[LINELEN];
  int ivalue, matched, ibest;



  /*Blank out the arrays we will be using here */

  for (i = 0; i < LINELEN; i++)
  {
    choices[i] = ' ';
    values[i] = ' ';
  }


  for (i = 0; i < strlen (word); i++)
  {
    word[i] = tolower (word[i]);
  }


  for (i = 0; i < strlen (string_choices); i++)
  {
    choices[i] = tolower (string_choices[i]);
  }


  for (i = 0; i < strlen (string_values); i++)
  {
    values[i] = tolower (string_values[i]);
  }


  for (i = 0; i < strlen (string_choices); i++)
  {
    if (choices[i] == ',')
    {
      choices[i] = ' ';
    }
  }



  for (i = 0; i < strlen (string_values); i++)
  {
    if (values[i] == ',')
    {
      values[i] = ' ';
    }
  }




  nchoices = sscanf (choices, "%s %s %s %s %s %s %s %s %s %s", xs[0], xs[1], xs[2], xs[3], xs[4], xs[5], xs[6], xs[7], xs[8], xs[9]);
  nchoices =
    sscanf (values, "%d %d %d %d %d %d %d %d %d %d", &xv[0], &xv[1], &xv[2], &xv[3], &xv[4], &xv[5], &xv[6], &xv[7], &xv[8], &xv[9]);


  /* Check that one of the values is not one of the values that is an error return value, -9998 or -9999.  If
   * xv is one of these values then the program will remain confused and so we exit 
   */

  for (i = 0; i < nchoices; i++)
  {
    if (xv[i] == -9998 || xv[i] == -9999)
    {
      Error ("string2int: Internal programming error: value for rdchoice is an error retrun value\n");
      exit (1);
    }
  }

  /* Perform a minimum match on the answer */

  matched = 0;
  ivalue = -9998;
  ibest = -1;                   //Set this to a sensible initial value
  for (i = 0; i < nchoices; i++)
  {
    if (strncmp (word, xs[i], strlen (word)) == 0)
    {
      ivalue = xv[i];
      ibest = i;
      matched += 1;
    }
  }

  strcpy (string_answer, "none");

  /* Check for more than one match to answer */

  if (matched > 1)
  {
    ivalue = -9999;
    return (ivalue);
  }



  if (ibest >= 0)
  {
    strcpy (string_answer, xs[ibest]);
  }


  return (ivalue);
}



/**********************************************************/
/** 
 * @brief      process a input line where the answer is expected to be one of several choices
 *
 * @param [in] char  question[]   the name of the parameter and the possible choices for which an answer is sought
 * @param [in] char  answers[]   The possible integer answers that are expected   
 * @param [in, out] char  answer[]   The current/final value of the parameter
 * @return    An integer indicatcating what choice was made                             
 *
 *
 * ###Notes###
 *
 * This routine reads a command line of the form
 *
 * "cooked(rare,medium,well)
 *
 * and returns an integer that corresponds to one of the choices of the words in parenthesis.
 *
 * The integer that is returned is that in the string answers,  for example
 *
 * "1,4,8"
 *
 * It gets the choice from the command line or the parameter file in the standard manner.  Unlike
 * most of the other routines in the package, this routine returns the integer value of the choice that was made.  
 * answer contains the string version of the answer.
 *
 * ###Programming Comment###
 *
 * For backward compatibility, if the user enters an integer rather than a string.  The integer is returned. No
 * error checking is done to see that the integer is one of the allowed choices.   A comment is written to the
 * whatever.out.pf file to indicate that one should change integers to strings as ultimately we plan to remove 
 * the code that provides backward compatibility.
 *
 * Very important:  in rdchoice, the string answer is changed. c is very picky about when a string can be updated.  One cannot
 * change a string that has simply been set to a fixed value,  e. g.
 *
 * char whatever="whatever"
 *
 * This is fixed in memory.  Instead, one nees to create an array with a certain length, e.g
 *
 * char answer[LINELENGTH];
 *
 * and use strcpy or some other routine to initialize it.  Then one can call rdchoice.
 **********************************************************/

int
rdchoice (question, answers, answer)
     char question[];
     char answers[];
     char *answer;
{
  char dummy[LINELEN];
  char string_answer[LINELEN];
  int n, nstart, nstop;
  int ianswer;
  int query;
  char full_answer[LINELEN];
  rdpar_choice = 1;
  strcpy (string_answer, answer);
  query = REISSUE;
  while (query == REISSUE)
  {
    query = rdstr (question, string_answer);
    /* First check to see if we have returned an integer.  The fact taht we attempt to find a string
     * after the integer is to make it possible for a string answer to be something like "2dcoords" */
    if (sscanf (string_answer, "%d%s", &ianswer, dummy) == 1)
    {
      strcpy (answer, string_answer);
      rdpar_comment ("Deprecated use of rdchoice. NO ERROR CHECKS! For %s replace answer %s in %s with its string equivalent",
                     question, string_answer, answers);
      fprintf (rdout_ptr, "%-30s %20s\n", question, string_answer);
      Error ("rdchoice: Deprecated use of rdchoice. NO ERROR CHECKS! For %s replace answer %s in %s with its string equivalent \n",
             question, string_answer, answers);
      strict = 1;
      rdpar_choice = 0;
      return (ianswer);
    }

    /* Otherwise we assume it is a new style input */

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

    strcpy (dummy, " ");
    strncpy (dummy, &question[nstart + 1], nstop - nstart - 1);
    strcpy (dummy, &question[nstart + 1]);
    dummy[strlen (dummy) - 1] = ' ';
    ianswer = string2int (string_answer, dummy, answers, full_answer);
    if (ianswer == -9998)
    {
      Error ("rdchoice: Could not match %s input to one of answers: %s\nTry again\n", string_answer, dummy);
      query = REISSUE;
    }
    if (ianswer == -9999)
    {
      Error ("rdchoice: Multiple matches of  %s input to answers: %s\nTry again\n", string_answer, dummy);
      query = REISSUE;
    }
  }

  if (query != OLD)
  {
    fprintf (rdout_ptr, "%-30s %20s\n", question, full_answer);
  }
  strcpy (answer, string_answer);
  rdpar_choice = 0;
  return (ianswer);
}

/* This is the end of the various routines which parse different kinds of inputs */





/**********************************************************/
/** 
 * @brief      get the rootnmae of a .pf file
 *
 * @param [out] char  root[]   The root name of a parameter file
 * @param [in] char  total[]   A character string containing a file name, assumed to be .pf file
 * @return  0  
 *
 * Bascially the point of this is just to get a base or root name
 * for various files produced by Python.
 *
 * ###Notes###
 *
 * The routine just strips off the .pf if it exists
 * and returns the rest of the name
 *
 * The order of the arguments may be some what odd, but it resembles
 * routines like strcpy in which the output is the first string
 *
 *
 *
 **********************************************************/

int
get_root (root, total)
     char root[], total[];
{
  int j;
  char *pf;
  int position;
  /* Check whether total is an empty string */
  j = strcspn (total, "\n");
  if (j == 0)
  {
    strcpy (root, "rdpar");
    return (0);
  }

  /* Check for .pf at the end of the string 
   * Note that there is no easy way to actually
   * check for the last occurence.  Here we
   * assume there is only one occurrcne of .pf
   */

  pf = strstr (total, ".pf");
  if (pf != NULL)
  {
    position = pf - total;
    strncpy (root, total, position);
    root[position] = '\0';
    return (0);
  }

  strncpy (root, total, j);
  root[j] = '\0';
  return (0);
}




/**********************************************************/
/** 
 * @brief      Comminicate the mpi_rank to rdpar  
 *
 * @param [in, out] int  rank   an integer that givves the process number rdapr should use
 * @return     0
 *
 * The next routine simply sets the rank of the process 
 * if not in parallel mode then we set rd_rank to zero
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
rdpar_set_mpi_rank (rank)
     int rank;
{
  rd_rank = rank;
  return (0);
}



/**********************************************************/
/** 
 * @brief      Set the level of verbosity for logging within the rdpar routines
 *
 * @param [in, out] int  vlevel   
 * @return     Always returns 0
 *
 * if vlevel is 2 or greater then more diagnostic information
 * will be printed out from the various rdpar toutines
 *
 * ###Notes###
 *
 * In writing to log files, Python allows one to control the
 * amount of information which is wrtten out by defining a verbosity
 * level.  The rdpar routines (because they are intended to be
 * useful in other applications), have their own verbosity 
 * control, and these routine is used to set it, based on the
 * verbosity for the entire program..
 *
 **********************************************************/

int
rdpar_set_verbose (vlevel)
     int vlevel;
{
  if (vlevel < 2)
    verbose = 0;
  return (0);
}


/* Check that inputs were correctly updated */

int
rdpar_check ()
{
  return (strict);
}
