
/***********************************************************/
/** @file  bilinear.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  As the name indicates the routines solve for 
 * the fractional position of a 2d vector in a 2d cell with
 * 4 corners.
 *
 ***********************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "log.h"


/**********************************************************/
/** 
 * @brief      calculates the fractional position in the coordinate
 *  	grid of a 2d vector x.
 *
 * @param [in] double  x[]   The 2d position within the cell
 * @param [in] double  x00[] The 2d position of the lower-left edge of cell  
 * @param [in] double  x01[] The 2d position of the upper-left edge of the cell
 * @param [in] double  x10[] The 2d position of the lower-right edge of the cell
 * @param [in] double  x11[] The 2d position of the upper-right edge of the cell
 * @param [out] double *  f  The fractional position in the "x" direction
 * @param [out] double *  g  The fractional postion in the "z" direction
 * @return     
 *     The routine itself returns 0 if the postion is inside the cell, and 
 *     -1 otherwise.
 *
 * @details
 * The input 2d input vectors are assumed to have the relationship
 *   shown in real space, that is x00 and x01 are more or less along
 *   the z axis, while x00 and x10 are more or less along the +x acis
 * 
 * * x01 ....... x11
 * * .            .
 * * .            .
 * * x00 ....... x10
 *
 * ### Notes ###
 *
 * Basically
 *         
 * X= (1-g) [ (1-f) X00 + f X10 ] * g [ (1-f)*X01+ f* X11 ]
 * 
 * where the capitalized values of X,etc denote vectors, and
 * f and g are the fractional weightings.  
 * 
 * X = (1-g) A + g B 
 * 
 * It is easy to solve for g in terms of f
 * 
 * g = (X-A)/(B-A)
 * 
 * where A = [ (1-f) X00 + f X10 ] and B =  [ (1-f)*X01+ f* X11 ]
 * 
 * But this gives you two expressions for g depending on whether you
 * consider the 0 or 1 component of the vector, and therefore you
 * can solve for f.  More specifically, 
 * 
 * g= [(X-XOO) + f (XOO-X10)] / [(X01-X10)+f (X11+X00-X01-X10)]
 * 
 * or 
 * 
 * g = (Q + f R) / (S+f T) 
 * 
 * where Q, R, S, and T are defined from the equation above.
 * 
 * So expanding in terms of the component veciots what we have is
 * 
 * (Q[0]+f R[0]) (S[2]+f T[2]) = (Q[2]+f R[2]) (S[0]+f T[0])
 * 
 * which can now be put in the form of a quadratic
 * 
 * a f**2 + b f +c = 0
 * 
 * where 
 * 
 * * a = r[0]*t[2] - r[2]*t[0]
 * * b=  r[0]*s[2] + q[0]*t[2] - ( r[2]*s[0] + q[2]*t[0] )
 * * c = q[0]*s[2] - q[2]*s[0]
 * 
 * 
 * A problem with this approach is that does not work in the case 
 * where the coordinates are the x axis positions of X01 and X01
 * are the same and the X10 and X11 are the same.  This is because
 * the denominator becomes 0 in the equation for g.  This 
 * particular case is important because it is the case we are
 * trying to solve for.
 * 
 * In this case we need to solve for f by itself. Similarly there 
 * is a degnereate case when
 *               (x00[2] == x10[2] && x01[2] == x11[2])
 *
 **********************************************************/

int
bilin (double x[], double x00[], double x01[], double x10[], double x11[], double *f, double *g)
{
  double z;
  double root[2];
  double a, b, c, d;
  double q[3], r[3], s[3], t[3];
  double zz[3];
  int i;
  int quadratic (double a, double b, double c, double r[]);
  void Exit (int error_code);


  z = 0;                        /* Initialize z to prevent warning on compilation */

  if (x00[0] == x01[0] && x10[0] == x11[0])
  {
    *f = z = (x[0] - x00[0]) / (x10[0] - x00[0]);
    a = (1. - z) * x00[2] + z * x10[2];
    b = (1. - z) * x01[2] + z * x11[2];
    *g = (x[2] - a) / (b - a);

  }
  else if (x00[2] == x10[2] && x01[2] == x11[2])
  {
    *g = z = (x[2] - x00[2]) / (x01[2] - x00[2]);
    a = (1. - z) * x00[0] + z * x01[0];
    b = (1. - z) * x10[0] + z * x11[0];
    *f = (x[0] - a) / (b - a);
  }
  else
  {
    q[0] = x[0] - x00[0];
    r[0] = x00[0] - x10[0];
    s[0] = x01[0] - x00[0];
    t[0] = x11[0] - x01[0] - x10[0] + x00[0];

    q[2] = x[2] - x00[2];
    r[2] = x00[2] - x10[2];
    s[2] = x01[2] - x00[2];
    t[2] = x11[2] - x01[2] - x10[2] + x00[2];


    a = r[0] * t[2] - r[2] * t[0];
    b = r[0] * s[2] + q[0] * t[2] - (r[2] * s[0] + q[2] * t[0]);
    c = q[0] * s[2] - q[2] * s[0];


    if (c == 0)
    {
      *f = z = 0;
    }
    else
    {
      i = quadratic (a, b, c, root);    /* root[i] is the smallest positive root unless i is
                                           negative in which case either both roots were negative or both roots 
                                           were imaginary */
      if (i != (-1))
      {                         // Find the root that closest to 0.5
        zz[0] = fabs (root[0] - 0.5);
        zz[1] = fabs (root[1] - 0.5);
        if (zz[0] < zz[1])
          z = root[0];
        else
          z = root[1];
        *f = z;

      }
      else
      {
        Error ("Bilinear -- imaginary roots from quadratic\n");
      }

    }


    if ((d = (s[2] + z * t[2])) != 0)
    {
      *g = (q[2] + r[2] * z) / (d);
    }
    else if ((d = (s[0] + z * t[0])) != 0)
    {
      *g = (q[0] + r[0] * z) / (d);
    }
    else
    {
      Error ("bilin: Denominator zero\n");
      Exit (0);
    }
  }

  if (*f < 0. || *f > 1. || *g < 0 || *g > 1.)
  {
    return (-1);
  }

  return (0);

}
