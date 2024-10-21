
/***********************************************************/
/** @file  xlog.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  These are a series of routines designed to store comments and errors
 * in a diagnostic file or files.  The routines also provide a mechanism for tracking the numbers of errors
 * of each type
 *
 * Instead of using printf and fprintf statements throughout varius subroutines, the set of routines
 * here are intended to provide a standard interface to various diagnostic files, and to manage the interaction of 
 * logging with the various mpi threads that can exist when running Python in multiprocessor mode.  The routines
 * also contain a verbosity mechanism to allow one to control how much information is written to the screen, and a 
 * mechanism to keep track of the number of times a particular error message has been generated.  The overall goal 
 * is to keep the log files produced by Python a manageable size.
 *
 * Messages in Python runs are sent to the screen and to diagnostic files.  During multiprocessing runs, a diagnostic 
 * file is opened for each thread.  With some exceptions, messages to the screen arise from thread 0.  
 * 		
 *  The routines that control what files are open and closed for logging are as follows.
 * - Log_init(filename)			Open a logfile
 * - Log_append(filename)			Reopen an existing log file (so that one can continue to log to it)
 * - Log_close()				Close the current logfile
 *
 *   The routines that write to the log files (and optionally the screen) for informational reasons are:
 *  - Log ( char *format, ...)			Send a message to the screen and logfile
 *  - Log_silent ( char *format, ...)		Send a message to the logfile
 *
 * The routines that are designed to report errors are as follows: 
 * - Shout( char *format, ...)	   Send a message to the screen and logfile that cannot be suppressed.
 * - Error ( char *format, ...)			Send a message prefaced by the word "Error:" to the screen and to the logfile.
 * - Error_silent ( char *format, ...)		Send a message prefaced by the word "Error:" to
 * 										the logfile
 *
 *  One can control how much information is printed to the screen and how many times a specific error message is logged 
 *  to a file with several routines
 *  - Log_set_verbosity(vlevel)			Set the verbosity of the what is printed to the screen and the log file
 *  - Log_print_max(print_max)			Set the number of times a single error will be 
 * 	    output to the screen and the log file
 *
 *  The amount printed out is controlled by a command line switch -v (in Python).   The default value causes Shout, Log, 
 *  and Error to be printed and logged, but not Log_silent Error_silent.  The default value is 5. If the value is raised 
 *  more data will be printed or logged.  
 *
 *  In most cases, it is sufficient to log to the screen only from thread 0.  there are a few times, one might want to 
 *  send a message to the screen from any thread. For this purpose there is a specific command:
 *  - Log_parallel(message)			Log statement for parallel reporting
 *
 *
 *  There are several specific commands that have been included for debugging problems:
 *  - Log_flush()					simply flushes the logfile to disk (before the program crashes).
 *  - Debug( char *format, ...) 			Log an statement to the screen and to a file.  This is essentially a 
 *						        intended to replace a printf statement in situations where
 *							one is debugging code.  The use of Debug instead of log
 *							means that a future developer should be free to remove the
 *							Debug statement from the code.  Because the command prepends the word Debug to 
 *							the write staatement, it should be easy to grep debug statements out of the log files
 *
 *
 *  For errors, the logging routines keep track of how many times a particular error has been reported, using the format 
 *  statement as a proxy for the error.  After a certin number of times a particular error has been roported the error is
 *  no longer written out to the diag file, but the routens still keep trank of the nubmer of times the error is reported.
 *  If this number becomes too large then, the progam exits, and when it does it indicates why. There are a nubmer of 
 *  routines associated with this.
 *
 *
 *  - error_summary(char *format)	Summarize all of the erors that have been
 * 					logged to this point in time
 *
 *
 *  In addition there are several routines that largely internal
 *
 *  - Log_set_mpi_rank(rank, n_mpi)		Tells  the rank of the parallel process in parallel mode,
 * 						and divides max errors by n_mpi
 *     
 *
 *
 *	
 *
 ***********************************************************/

#ifdef MPI_ON
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "log.h"

#define LINELENGTH 256
#define NERROR_MAX 500          // Number of different errors that are recorded

/* definitions of what is logged at what verbosity level */

#define SHOW_PARALLEL	    1
#define SHOW_ERROR	    2
#define SHOW_LOG  	    3
#define SHOW_DEBUG	    4
#define SHOW_LOG_SILENT     5
#define SHOW_ERROR_SILENT   5


int my_rank = 0;                // rank of mpi process, set to zero
int n_mpi_procs = 0;            // the number of mpi processes

int log_print_max = 100;        // Maximum number of times a single error will be reported.  
                                // Note that it will still be counted.
int max_errors = 1000000;       // Maximum number of times an error can occur before giving up

typedef struct error_log
{
  char description[LINELENGTH];
  int n;
} error_dummy, *ErrorPtr;

ErrorPtr errorlog;


int nerrors;

FILE *diagptr;
int init_log = 0;
int log_verbosity = 5;          // A parameter which can be used to suppress what would normally be logged or printed


/**********************************************************/
/** 
 * @brief      Open a log file 
 * 		
 * @param [in] char *  filename   The name of the file where logging will occur
 * @return     Always returns 0
 *
 * ###Notes###
 *
 * The routine will exit if log the file can not be opened. If the log file can
 * be opened, then atexit is used to ensure that upon normal termination that
 * the log file buffer is flushed and the log file closed.
 *
 **********************************************************/

int
Log_init (filename)
     char *filename;
{
  if ((diagptr = fopen (filename, "w")) == NULL)
  {
    printf ("Yikes: could not even open log file %s\n", filename);
    Exit (0);
  }
  init_log = 1;

  nerrors = 0;
  errorlog = (ErrorPtr) calloc (sizeof (error_dummy), NERROR_MAX);

  if (errorlog == NULL)
  {
    printf ("There is a problem in allocating memory for the errorlog structure\n");
    Exit (0);
  }

  return (0);
}



/**********************************************************/
/** 
 * @brief      Opens an existing log file so that diagnostic message can be added to it.
 *
 * @param [in] char *  filename   Name of the existing file to reopen
 * @return     0 unless the log file cannot be reopened in which case the program terminates
 *
 * The routine opens a logfile that should have existed previously so that one 
 * can continue an early run of sirocco 
 *
 * ###Notes###
 *
 * This routine is called on a restart of a run, if for example one
 * needed to checkpoint a run because a limitation on the time for running
 * an individual program 
 *
 **********************************************************/

int
Log_append (filename)
     char *filename;
{
  if ((diagptr = fopen (filename, "a")) == NULL)
  {
    printf ("Yikes: could not even open log file %s\n", filename);
    Exit (0);
  }
  init_log = 1;

  nerrors = 0;
  errorlog = (ErrorPtr) calloc (sizeof (error_dummy), NERROR_MAX);

  if (errorlog == NULL)
  {
    printf ("There is a problem in allocating memory for the errorlog structure\n");
    Exit (0);
  }

  return (0);
}


/**********************************************************/
/** 
 * @brief      Close the log file 
 *
 * @return     Always returns 0
 *
 * Close the current log file
 *
 * ###Notes###
 *
 * One should not need to worry about when closing the log file and missing some
 * output, as fflush is used to flush the buffer and write to file.
 *
 **********************************************************/

void
Log_close ()
{
  init_log = 0;
  fflush (diagptr);
  fclose (diagptr);
  free (errorlog);              // Release the error summary structure
}

/* The next routine allows the user to change the amount of reporting that is
 * carried out through the logging subroutines
 */


/**********************************************************/
/** 
 * @brief      Control the range of messages which are logged to the command line
 *
 * @param [in] int  vlevel   An integer which controls what level of 
 *    logging goes to the screen or to a file
 * @return     Always returns 0
 *
 *
 * ###Notes###
 * The default value of the verbosity is 5.  Setting the
 * verbosity to a higher valued will cause more data to
 * be printed or logged.  Setting it to a lower value
 * will cause less.
 *
 * This routine is accessed via the command line with the 
 * swithc -v
 *
 *
 **********************************************************/

int
Log_set_verbosity (vlevel)
     int vlevel;
{
  log_verbosity = vlevel;
  rdpar_set_verbose (vlevel);
  return (0);
}



/**********************************************************/
/** 
 * @brief      Limit the number of times an error is printed out
 *
 * @param [in] int  print_max   a number which specifies the number of times an error message is logged
 * @return     Always returns 0
 *
 * This routines allows the user to change the number of times an error message is printed out
 *
 *
 * ###Notes###
 * Even if the error is not printed out, occurences are counted
 *
 * The number (which be default is 100) can be modified ont he
 * commandline with the switch -e_write
 *
 *
 **********************************************************/

int
Log_print_max (print_max)
     int print_max;
{
  log_print_max = print_max;
  return (0);
}



/**********************************************************/
/** 
 * @brief      Set the number of errors of a specific type which will cause the program to exit
 *
 * @param [in] int  n   the maximum number of errors
 * @return     0               
 *
 * When running the program keeps track of how many errors of any type have been issues, the assumption
 * being that if there are two many errors something is drastically wrong and the code should
 * exit gracefully.  This routine resets that number from the default
 *
 * ###Notes###
 *
 * The current default is 1e5.
 *
 *
 **********************************************************/

int
Log_quit_after_n_errors (n)
     int n;
{
  max_errors = n;
  return (0);
}


/**********************************************************/
/** 
 * @brief      Print/write an informational message
 *
 * @param [in] char *  format   The format string for the message
 * @param [in]   ...   The various values which fill out the format screen
 * @return     The number of characters sucessfully written
 *
 * This is the standard way of printing a message to screen and to the diag file
 *
 * ###Notes###
 *
 * Printing to the screen can be suppressed by setting the verbosity level
 *
 **********************************************************/

int
Log (char *format, ...)
{
  va_list ap, ap2;
  int result;

  if (init_log == 0)
    Log_init ("logfile");

  if (log_verbosity < SHOW_LOG)
    return (0);

  va_start (ap, format);
  va_copy (ap2, ap);            /* ap is not necessarily preserved by vprintf */

  if (my_rank == 0)
    result = vprintf (format, ap);
  result = vfprintf (diagptr, format, ap2);
  va_end (ap);
  return (result);
}


/**********************************************************/
/** 
 * @brief      Write a message to the diagnostic file
 *
 * @param [in] char *  format   The format string for the message
 * @param [in]   ...   The various values which fill out the format screen
 * @return     The number of characters sucessfully written
 *
 * This is the standard way of writing a message to the diag file (without also writing to the screen)
 *
 * ###Notes###
 *
 * Printing to the screen can be enabled if the level of verbosity is set to high enough a level
 *
 **********************************************************/

int
Log_silent (char *format, ...)
{
  va_list ap, ap2;
  int result;

  if (init_log == 0)
    Log_init ("logfile");

  if (log_verbosity < SHOW_LOG_SILENT)
    return (0);
  va_start (ap, format);
  va_copy (ap2, ap);            /* ap is not necessarily preserved by vprintf */

  result = vfprintf (diagptr, format, ap);
  va_end (ap);
  return (result);
}


/**********************************************************/
/** 
 * @brief      Print/write out an error message
 *
 * @param [in] char *  format   The format string for the message
 * @param [in]   ...   The various values which fill out the format statement
 * @return     The number of characters sucessfully written
 *
 *
 * This is the standard way of writing an error message to the screen and to
 * a file
 *
 * ###Notes###
 *
 * Writing to the screen can be suppressed depending on the level of verbosity
 *
 **********************************************************/

int
Error (char *format, ...)
{
  va_list ap, ap2;
  int result;

  if (init_log == 0)
    Log_init ("logfile");

  if (error_count (format) > log_print_max || log_verbosity < SHOW_ERROR)
    return (0);

  va_start (ap, format);
  va_copy (ap2, ap);            /*NSH 121212 - Line added to allow error logging to work */
  if (my_rank == 0)             // only want to print errors if master thread
    result = vprintf (format, ap);

  fprintf (diagptr, "Error: ");
  result = vfprintf (diagptr, format, ap2);
  va_end (ap);
  return (result);
}




/**********************************************************/
/** 
 * @brief      Write an error message to the diag file
 *
 * @param [in] char *  format   The format string for the message
 * @param [in]   ...   The various values which fill out the format statement
 * @return     The number of characters sucessfully written
 *
 *
 * These routine normally only writes an error message to the diag file
 *
 * ###Notes###
 *
 * Writing the message to the screen can be enabled using the verbosity settings
 *
 **********************************************************/

int
Error_silent (char *format, ...)
{
  va_list ap, ap2;
  int result;
  if (init_log == 0)
    Log_init ("logfile");

  if (error_count (format) > log_print_max || log_verbosity < SHOW_ERROR_SILENT)
    return (0);

  va_start (ap, format);
  va_copy (ap2, ap);            /* ap is not necessarily preserved by vprintf */

  if (my_rank == 0)             // only want to print errors if master thread
    result = vprintf (format, ap);
  fprintf (diagptr, "Error: ");
  result = vfprintf (diagptr, format, ap2);
  va_end (ap);
  return (result);
}




/**********************************************************/
/** 
 * @brief      Write an error message to the screen and to a file
 *
 * @param [in] char *  format   The format string for the message
 * @param [in]   ...   The various values which fill out the format statement
 * @return     The number of characters sucessfully written
 *
 *
 * This routine writes and error message to the screen and to the diagnostic file.
 * This message is not suppressed by any of the verbosity flags and is intended
 * to be for the most important errors
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
Shout (char *format, ...)
{
  va_list ap, ap2;
  int result;

  if (init_log == 0)
    Log_init ("logfile");

  if (error_count (format) > log_print_max)
    return (0);

  printf ("Error: ");
  va_start (ap, format);
  va_copy (ap2, ap);            /* ap is not necessarily preserved by vprintf */
  result = vprintf (format, ap);
  fprintf (diagptr, "Error: ");
  result = vfprintf (diagptr, format, ap);
  va_end (ap);
  return (result);
}


/**********************************************************/
/** 
 * @brief      Check that a variable is valid double precision number
 *
 * @param [in] double  x   the variable to check
 * @return     0 if the number is valid, -1 if it is a NAN, or infinity
 *
 * This is a diagnostic routine to allow one to check whether NANs or 
 * infinities have crept into a calculation.  
 *
 * ###Notes###
 *
 * The routine was intended to make it easy to set break points when 
 * something untoward is happening.
 *
 **********************************************************/

int
sane_check (x)
     double x;
{
  int i;
  if ((i = isfinite (x)) == 0)
  {
    Error ("sane_check: %d %e\n", i, x);
    return (-1);
  }
  return (0);
}



/**********************************************************/
/** 
 * @brief      track the number of errors of each type
 *
 * @param [in, out] char *  format   A format statement for an error message
 * @return     The number of errors wwith associated with a certain format
 *
 * When an error message is received, the format statement is used to identify the
 * error.  Here the format statement is compared to the format statement associated
 * with previous errors and a counter is kept of the number of times a particular error
 * has occured.
 *
 * Once the error count for a particular error has been reached then the error is no
 * longer printed out; once the error count has reached a much larger number the
 * program terminates on the assumption that something is drastically wrong.
 *
 * ###Notes###
 *
 * The number for stopping the print out is contolled by NERROR_MAX and is hardcoded
 *
 * The number for stopping  the program is controled by max_errors and can be altered, see
 * log_set_max_errors 
 *
 **********************************************************/

int
error_count (char *format)
{
  int n;
  n = 0;
  while (n < nerrors)
  {
    if (strcmp (errorlog[n].description, (format)) == 0)
      break;
    n++;
  }

  if (n == nerrors)
  {
    strcpy (errorlog[nerrors].description, format);
    errorlog[n].n = 1;
    if (nerrors < NERROR_MAX)
    {
      nerrors++;
    }
    else
    {
      printf ("Exceeded number of different errors that can be stored\n");
      error_summary ("Quitting because there are too many differnt types of errors\n");
      Exit (0);
    }
  }
  else
  {
    n = errorlog[n].n++;
    if (n == log_print_max)
      Error ("error_count: This error will no longer be logged: %s\n", format);
    if (n == max_errors)
    {
      error_summary ("Something is drastically wrong for any error to occur so much!\n");
      Exit (0);
    }
  }
  return (n + 1);
}




/**********************************************************/
/** 
 * @brief      Summarize the errors that have occurred 
 *
 * @param [in] char *  message   A message that can accompany the error summary
 * @return     Always returns 0
 *
 * Print out the errors that have occured and the number of times each 
 * error has occured
 *
 * ###Notes###
 *
 * This is printed out at the end of all Python runs.  When running in multiprocessor mode, 
 * the number of 
 * errors referes only to the errors that have occurred in that particular theread.
 *
 **********************************************************/

int
error_summary (message)
     char *message;
{
  int n;

  if (nerrors == 0)
    return 0;

  Log ("\nError summary: %s\n", message);
  Log ("Recurrences --  Description\n");
  for (n = 0; n < nerrors; n++)
  {
    Log ("%9d -- %s", errorlog[n].n, errorlog[n].description);
  }

  Log_flush ();
  return (0);
}




/**********************************************************/
/**
 * @brief      Summarize the errors that have occurred for an mpi process
 *
 * @param [in] char *  message   A message that can accompany the error summary
 * @return    void
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
error_summary_parallel (char *msg)
{
  int i;

  if (nerrors == 0)
    return 0;

  Log_parallel ("\nError summary for thread %i: %s\n", my_rank, msg);
  Log_parallel ("Recurrences --  Description\n");

  for (i = 0; i < nerrors; i++)
    Log_parallel ("%9d -- %s", errorlog[i].n, errorlog[i].description);

  Log_flush ();

  return 0;
}




/**********************************************************/
/** 
 * @brief      Flush the diagnostic file to assure that one has an up-to-date version of the log file
 *
 * @return     Always returns 0
 *
 * This routine is intended to assure that the log file is complete before
 * a possilbe program crash
 *
 * ###Notes###
 *
 *  
 *
 **********************************************************/

int
Log_flush ()
{
  if (init_log == 0)
    Log_init ("logfile");

  fflush (diagptr);
  return (0);
}





/**********************************************************/
/** 
 * @brief      Store the rank of this particular thread
 *
 * @param [in] int  rank   The rank of this thread
 * @param [in] int  n_mpi   The total number of threads 
 * @return     Always returns 0
 *
 * The routine simply sets the rank of the process 
 * if not in parallel mode then we set my_rank to zero
 *
 * ###Notes###
 *
 * The number of thread is not used by the program even though it
 * is passed.
 *
 **********************************************************/

int
Log_set_mpi_rank (rank, n_mpi)
     int rank, n_mpi;
{
  my_rank = rank;
  n_mpi_procs = n_mpi;
  rdpar_set_mpi_rank (rank);    //this just communicates the rank to rdpar      

  return (0);
}


/* log statement that prints parallel statements to screen */


/**********************************************************/
/** 
 * @brief      In addition to writing to the diag file, print a message to the screen from all threads
 *
 * @param [in] char *  format   The format string for the message
 * @param [in]   ...   The various values which fill out the format statement
 * @return     The number of characters sucessfully written
 *
 *
 * The normal logging routines write a message to the diagnositc file associated
 * with a paticular thread, and write messages to the screen only for thread 
 * zero.  This routine writes a message to the screen from all threads.
 *
 * ###Notes###
 *
 * 26/11/18 EP: removed if statement so Log_parallel really does print a message
 * to stdout for all MPI processes.
 *
 **********************************************************/

int
Log_parallel (char *format, ...)
{
  va_list ap, ap2;
  int result;

  if (init_log == 0)
    Log_init ("logfile");

  if (log_verbosity < SHOW_PARALLEL)
    return (0);

  va_start (ap, format);
  va_copy (ap2, ap);

  result = vprintf (format, ap);

  fprintf (diagptr, "Para: ");
  result = vfprintf (diagptr, format, ap2);

  return (result);
}




/**********************************************************/
/** 
 * @brief      Print/write to a file a special message for debugging
 *
 * @param [in] char *  format   The format string for the message
 * @param [in]   ...   The various values which fill out the format statement
 * @return     The number of characters sucessfully written
 *
 * Straight  printf and fprintf states are strongly discouraged in Python except
 * as issued withing the routines here.  This routine preface an fprintf statement
 * with Dobug so such statements are issue to grep out of a log file, and in the code
 *
 * The intention that this routine would be used to log messages when one
 * is trying to debug a problme and that these lines would be removed from
 * the code once the debugging was completed. 
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
Debug (char *format, ...)
{
  va_list ap, ap2;
  int result;

  if (log_verbosity < SHOW_DEBUG)
    return (0);

  if (init_log == 0)
    Log_init ("logfile");

  va_start (ap, format);
  va_copy (ap2, ap);
  if (my_rank == 0)
    vprintf ("Debug: ", ap);
  result = vprintf (format, ap);
  fprintf (diagptr, "Debug: ");
  result = vfprintf (diagptr, format, ap2);
  va_end (ap);
  return (result);
}


/**********************************************************/
/**
 *  @brief Wrapper function to exit MPI/Python.
 *
 *  @param[in] int error_code. An integer classifying the error. Note that this
 *  number should be a non-zero integer.
 *
 *  @details
 *  When MPI is in use, using the standard C library exit function can be somewhat
 *  dangerous, as if a process exits with return code 0, MPI assumes that the
 *  process has exited correctly and the program will continue. This results in
 *  a deadlock and without any safety mechanism can result in an MPI run never
 *  finishing. Hence, in multiprocessor mode, one should use use MPI_Abort
 *  (which isn't a very graceful) to ensure that if a process does need to exit,
 *  then all processes will exit with it.
 *
 *  Before exiting, an error summary is printed to screen and and log is flushed
 *  to ensure everything is up to date to aid with error diagnosis.
 *
 *  ### Notes ###
 *
 *  Exit codes should be non-zero. In the past, we used to warn the user that
 *  a non-zero exit code was sent, but we found that this confused non-expert
 *  users (and us sometimes). Hence we silently now change the exit code to
 *  EXIT_FAILURE and move on with our lives.
 *
 **********************************************************/

void
Exit (int error_code)
{
  Log_flush ();

#ifdef MPI_ON
  Log_parallel ("--------------------------------------------------------------------------\n"
                "Aborting rank %04i: exiting all processes with error %i\n", my_rank, error_code);
  error_summary_parallel ("summary prior to abort");

  if (n_mpi_procs > 1)
  {
    MPI_Abort (MPI_COMM_WORLD, error_code);
  }
  else
  {
    exit (error_code);
  }
#else
  Log ("--------------------------------------------------------------------------\n" "Aborting: exiting with error %i\n", error_code);
  error_summary ("summary prior to abort");
  exit (error_code);
#endif
}
