/***********************************************************/
/** @file  time.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  A few simple routines relating to estimate
 * wallclock time for the program to run and to get the current
 * date an time
 *
 ***********************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include "log.h"

/*
Return the time in seconds since the timer was initiated
*/

double t_init = 0.0;


/**********************************************************/
/** 
 * @brief      Return the time in seconds since the timer was initiated
 *
 * @return     The time since t_init
 *
 * The first time the routine is called it gets the time in seconds
 * sets t_init to this * and returns 0,  Subsequently it returns 
 * the differecne
 *
 *
 * ###Notes###
 *
 * Uses gettimeofday
 *
 *
 **********************************************************/

double
timer ()
{
  struct timeval tv;
  struct timezone tz;
  double t;

  gettimeofday (&tv, &tz);
  t = tv.tv_sec + 1.e-6 * tv.tv_usec;
  if (t_init == 0.0)
  {
    t_init = t;
    return (0.0);
  }
  else
  {
    return (t - t_init);
  }
}




/**********************************************************/
/** 
 * @brief      Get the current time as an ascii string
 *
 * @param [out] char  curtime[]   The current time and datae
 * @return     Always returns 0  
 *
 * This simply gets the data and time
 *
 * ###Notes###
 *
 * Uses ctime.  ctime returns a string with a trailing \n, which we need
 * to strip off.  This accounts for the rather bizarre handling
 * of the string.
 *
 *
 **********************************************************/

int
get_time (char curtime[])
{
  time_t tloc;
  time (&tloc);
  strcpy (curtime, "");
  strncpy (curtime, ctime (&tloc), 24);
  curtime[24] = '\0';           // We need to end the string properly
  return (0);
}




/**********************************************************/
/**
 * @brief initialise a timeval structure with the current time
 *
 * param[in] void
 *
 * @return struct timeval timer_t0 The current CPU time as a timespec structure
 *
 * ###Notes###
 *
 * Note that NULL is passed as the second argument to gettimeofday as we most
 * likely do not require to know the time zone.
 *
 **********************************************************/


struct timeval
init_timer_t0 (void)
{
  struct timeval timer_t0;

  gettimeofday (&timer_t0, NULL);

  return timer_t0;
}




/**********************************************************/
/**
 * @brief Calculate the time difference between a previous timer and print the
 *        duration and a message to the screen
 *
 * @param[in] char *msg    A descriptive measure to print with the time difference
 * @param[in] struct timespec timer_t0  An initialised timer to calculate the time difference from
 *
 * @return void
 *
 * The time difference between this function being called and when the timer
 * time_t0 was initialised will be printed to screen down to nanosecond resolution.
 * The time difference will be preceded by the message provided by *msg.
 *
 * ###Notes###
 *
 * Note that NULL is passed as the second argument to gettimeofday as we most
 * likely do not require to know the time zone.
 *
 **********************************************************/

void
print_timer_duration (char *msg, struct timeval timer_t0)
{
  double td;
  struct timeval timer_t1;

  gettimeofday (&timer_t1, NULL);
  td = (timer_t1.tv_sec - timer_t0.tv_sec) + (timer_t1.tv_usec - timer_t0.tv_usec) * 1e-6;
  Log ("%s: %f seconds\n", msg, td);
}
