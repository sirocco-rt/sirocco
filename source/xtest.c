

/***********************************************************/
/** @file  xtest.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  A test routine to check a part of the program 
 * for any special test.  This is not intended to be 
 * used for anything permanent.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/** 
 * @brief      Conduct a diagnostic test for some purpose
 *
 *
 * @details
 *  This is a totally uncommented routine for use 
 *  with Python to perform some specific diagnostic 
 *  test
 *
 * ### Notes ###
 *
 * xtest is called based on a commmand line switch
 * -xtest.  Depending on whether an exit command is
 * included here (the normal situation), the code
 * will continue to run after that.
 *
 * The routine is called after all of the inputs
 * of a run have been read in.  In many cases
 * (namely when the -r (restart option is used,
 * a previous wind file.
 *
 * xtest is intended to be temporary.  Developers should
 * not expect to be the same the next time they want to
 * debug something.  As a result, developers may want
 * to save their version of xtest.c outside of the
 * git repository for later use.
 *
 * Ideally, developers would indicate what issue one trying
 * to solve and recored the commit assoicated with this 
 * along with the issue, in case one wants to locate
 * what was done
 *
 * This version is assocated with issue 900
 *
 *
 **********************************************************/

int
xtest ()
{
  double pos[3];
  double x, xmin, xmax;
  double y, z, delta;
  double rzero, theta, vel[3];
  FILE *fopen (), *fptr;
  modes.run_xtest_diagnostics = TRUE;
  double v_escape;

  z = 1.e15;
// z = 1e+14;
  y = 0;
  xmin = 0;
  xmax = 5.67e17;
  xmax = 2e16;
  xmax = 10. * z;

  pos[1] = y;
  pos[2] = z;

  fptr = fopen ("xtest.diag", "w");

  delta = (xmax - xmin) / 1000.;

  for (x = xmin; x <= xmax; x += delta)
  {
    pos[0] = x;
    rzero = sv_find_wind_rzero (0, pos);
    theta = sv_theta_wind (0, rzero);

    v_escape = zdom[0].sv_v_infinity * sqrt (2. * GRAV * geo.mstar / rzero);
    sv_velocity (pos, vel, 0);
    fprintf (fptr, "%.3e %.3e %.3e   %.3e %.3e  %.3e %.3e %.3e %.3e\n", pos[0], pos[1], pos[2], rzero, theta, vel[0], vel[1], vel[2],
             v_escape);

  }


  wind_save (files.windsave);
  do_windsave2table (files.root, 0, FALSE);




//OLD  double x[3];
//OLD  x[0] = 4.02e+10;
//OLD  x[1] = 0;
//OLD  x[2] = 1.24e+09;
//OLD  double v[3];


//OLD  Log ("Start  %.5e %.5e %.5e\n", x[0], x[1], x[2]);
//OLD  kn_velocity (0, x, v);

//OLD  Log ("Result %.5e %.5e %.5e\n", v[0], v[1], v[2]);



  exit (0);

}
