
/***********************************************************/
/** @file  vvector.c
 * @author ksl
 * @date   January, 2018
 *
 *
 * @brief      A series of subroutines designed to carry out common 3 vector operations
 * 		
 * These routines should be kept SEPARATE from routines that require the Python specific
 * structures in sirocco.h
 *
 * Many of the simpler routines are intened to simplify the code, and are not complicated.
 * Others especially those that transform between cylindrical and cartesian coordiantes
 * reflect the fact that photons travel in cartesian spece but the models are azimuthally
 * (sometimes) spherically symmetric.
 *
 * 	These are simple utilites
 * 	- dot(a,b)		returns the  product of the 3 vectors a and b
 * 	- double length(x)	returns the length of the 3 vector x
 * 	- int renorm(a,scalar)		renormalizes the 3-vector a to a length of
 * 							scalaar; retuns -1 if a has length 0.
 * 	- int cross(a,b,c)		c = a x b  (where x is the cross-product)
 * 	- vmove(u,lmn,s,result)
 * 	- vsub(u,v,result)
 * 	- vadd(u,v,result)
 * 	- stuff_v(vin,vout)
 *
 * 	These are transformations between cylindricatl and cartesian coordiantes
 * 
 * 	- int project_from_xyz_cyl(a,b,result)
 * 	- int project_from_cyl_xyz(a,b,result)
 * 					
 * 	The next few routines are designed to allow one to handle relatively arbitrary basis manipulations
 * 	
 *  -int create_basis(u,v,basis_new)
 * 	-project_from(basis_from, v_in,v_out)
 * 	-int project_to(basis_from, v_in,v_out)
 * 	-reorient(basis_from,basis_to,v_from,v_to)
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "log.h"
#include "math_struc.h"
#include "math_proto.h"
#define EPS 1.e-10

/* A basis is defined such that if x is a 3 vector as expressed an unprimed cartesian coordinate
   frame, and if y is the same vector in some rotated frame, then
   x[i] = a[i][j] y[j]

History:
	02jan	ksl	Encountered some weird Nan errors with py_ray and so rewrote to
			include a sane_check
 */



/**********************************************************/
/** 
 * @brief      The dot product of two vectors
 * 		
 * @param [in, out] double  a[]   the first vector
 * @param [in, out] double  b[]   the second vector 
 * @return     The resulting dot product
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

double
dot (a, b)
     double a[], b[];
{
  double x;

  x = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

  if (sane_check (x))
  {
    Error ("dot:sane_check of x: a %f %f %f b %f %f %f\n", a[0], a[1], a[2], b[0], b[1], b[2]);
  }
  return (x);
}



/**********************************************************/
/** 
 * @brief      The length of a 3 vector
 *
 * @param [in, out] double  a[]   a vector
 * @return     The leencth of the vector
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

double
length (a)
     double a[];
{
  double x, y;
  double sqrt ();
  y = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  if (sane_check (y))
    Error ("length:sane_check of y: a %f %f %f \n", a[0], a[1], a[2]);
  x = sqrt (y);
  if (sane_check (x))
    Error ("length:sane_check of x=sqrt(y): y= %f\n", y);
  return (x);
}




/**********************************************************/
/** 
 * @brief      renormalize a vector to a length
 *
 * @param [in, out] double  a[]  the vector
 * @param [in] double  scalar   the desired length
 * @return     Always returns 0 
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
renorm (a, scalar)
     double a[], scalar;
{
  double x;

  x = dot (a, a);
  if (x < EPS)
  {
    Log ("renorm: Cannot renormalize a vector [%e, %e, %e] of length 0\n", a[0], a[1], a[2]);
    return -EXIT_FAILURE;
  }

  x = scalar / sqrt (x);
  a[0] *= x;
  a[1] *= x;
  a[2] *= x;

  return (0);
}




/**********************************************************/
/** 
 * @brief      multiply  a vector by a scalar
 *
 * @param [in] double  a[]  the input vector
 * @param [in] double  scalar   the number to scale the vector by
 * @param [out]double  b[]  the output vector      
 * @return     Always returns 0 
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
rescale (a, scalar, b)
     double a[], scalar, b[];
{

  b[0] = a[0] * scalar;
  b[1] = a[1] * scalar;
  b[2] = a[2] * scalar;
  return (0);
}

/**********************************************************/
/** 
 * @brief      the cross product of two 3 vectors
 *
 * @param [in] double  a[]   the first vector
 * @param [inout] double  b[]   the second vector
 * @param [out] double  c[]   the cross product
 * @return     Always returns 0 
 *
 * c = a x b 
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
cross (a, b, c)
     double a[], b[], c[];
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];

  return (0);
}



/**********************************************************/
/** 
 * @brief      move a 3 vector a certain distance in a given direction
 *
 * @param [in] double  u[]   the starting position
 * @param [in] double  lmn[]   the directon
 * @param [in] double  s   the distance
 * @param [out] double  result[]   the resulting vector
 * @return     Always returns 0
 *
 *
 * result = u + s * lmn
 *
 * ###Notes###
 *
 * lmn must be a normalized direction vector.
 *
 **********************************************************/

int
vmove (u, lmn, s, result)
     double u[], lmn[], s, result[];
{
  result[0] = lmn[0] * s + u[0];
  result[1] = lmn[1] * s + u[1];
  result[2] = lmn[2] * s + u[2];
  return (0);
}


/**********************************************************/
/** 
 * @brief      Subtract two 3 vectors from one another
 *
 * @param [in] double  u[]   The initial 3 vector (or position)
 * @param [in] double  v[]   The vector to subtact
 * @param [out] double  result[]   The answeer
 * @return     Always returns 0
 *
 * result = u - v
 * ###Notes###
 *
 *
 **********************************************************/

int
vsub (u, v, result)
     double u[], v[], result[];
{
  result[0] = u[0] - v[0];
  result[1] = u[1] - v[1];
  result[2] = u[2] - v[2];
  return (0);

}


/**********************************************************/
/** 
 * @brief      Add two 3 vectors 
 *
 * @param [in] double  u[]   The first vector 
 * @param [in] double  v[]   The second vector
 * @param [out] double  result[]   the answere
 * @return     Always returns 0
 *
 * result = u + v
 * ###Notes###
 *
 *
 **********************************************************/

int
vadd (u, v, result)
     double u[], v[], result[];
{
  result[0] = u[0] + v[0];
  result[1] = u[1] + v[1];
  result[2] = u[2] + v[2];
  return (0);

}


/**********************************************************/
/** 
 * @brief      Simple stuff one 3 vector into another one
 *
 * @param [in] double  vin[]   The vector to be copied
 * @param [out] double  vout[]   The place where the vector was copied
 * @return     0                
 *
 *  vout = vin
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
stuff_v (vin, vout)
     double vin[], vout[];
{
  vout[0] = vin[0];
  vout[1] = vin[1];
  vout[2] = vin[2];
  return (0);

}



/**********************************************************/
/** 
 * @brief      Dot product of a 2d tensor and a vector to produce
 * a vector
 *
 * @param [in] double  tensor[3][3]   The tensor
 * @param [inout] double  vin[3]   The input vector
 * @param [out] double  vout[3]   The resulting vector
 * @return     The dot product of vin and vout
 *
 *
 * ###Notes###
 *
 * The program returns the dot product of the 
 * resultant, so if for example the tensor is the velocity
 * gradient tensor, and if vin is a direction (of unit length)
 * then the program will return the change in speed along
 * the direction of vin. 
 *
 *
 **********************************************************/

double
dot_tensor_vec (tensor, vin, vout)
     double tensor[3][3], vin[3], vout[3];
{
  double dot ();
  vout[0] = dot (tensor[0], vin);
  vout[1] = dot (tensor[1], vin);
  vout[2] = dot (tensor[2], vin);
  return (dot (vin, vout));
}



/**********************************************************/
/** 
 * @brief       Project a vector b in xyz coords from position a in xyz coords
 * into cylindrical coordinates
 *
 * @param [in, out] double  a[]   The position of the vector in cartesiona coords
 * @param [in, out] double  b[]   The vector (also in cartesiona coordinates)
 * @param [in, out] double  result[]   The result in cylindrical coordinates
 * @return     Always returns 0    
 *
 * Start with two vectors a and b in cartesiona coordinates, where a is a position
 * in space and v is a vector (for example a velocity at that position), and 
 * express v in cylindrical coordinates 
 * velocity at some position in a flow.
 *
 * ###Notes###
 *
 **********************************************************/

int
project_from_xyz_cyl (a, b, result)
     double a[], b[], result[];
{

  double n_rho[3], n_phi[3], n_z[3];

  n_rho[0] = a[0];
  n_rho[1] = a[1];
  n_rho[2] = 0;
  if (renorm (n_rho, 1.0) == -1)
  {
    Log ("Position on z axis; conversion cylintrical coords indeterminate\n");
    return (-1);
  }
  n_z[0] = n_z[1] = 0;
  n_z[2] = 1;
  cross (n_z, n_rho, n_phi);

  result[0] = dot (b, n_rho);
  result[1] = dot (b, n_phi);
  result[2] = b[2];


  return (0);
}




/**********************************************************/
/** 
 * @brief      Transform a vector in cylindrical coordinates at a position in xyz coordinates and convert the vector to
 * cartesian coordinates
 *
 * @param [in] double  a[]   The position at which one wishes to project the vector 
 * @param [in] double  b[]   The vector one wishes to project
 * @param [out] double  result[]   The resulting vector in cartesian coordiantes
 * @return     0
 *
 * 
 * Start with a vector b in cylindrical coordinates, and a position a, expressed in
 * cartesian coordinates, and express b in cartesian coordinates.
 *
 * ###Notes###
 * The naming and nomeclature are confusing.  a is a position in cartesian coordinates.
 * b is a 3 vector whose direction relative to a cylindrically symmetric coordianate
 * system is fixed.  The routine rotates b so it is in the correct direction in at the
 * postion given by a.  When a is in the xz plane no rotation is required.
 *
 *
 *
 * Note this is clearly different
 * from the situation where both a and b are in cylindrical coordinates.
 *
 *
 **********************************************************/

int
project_from_cyl_xyz (a, b, result)
     double a[], b[], result[];
{
  double x, ctheta, stheta;

  if ((x = sqrt (a[0] * a[0] + a[1] * a[1])) == 0)
  {
    Error ("project_from_cyl_xyz: indeterminate a =%f %f %f\n", a[0], a[1], a[2]);
  }
  ctheta = a[0] / x;
  stheta = a[1] / x;

  result[0] = b[0] * ctheta - b[1] * stheta;
  result[1] = b[0] * stheta + b[1] * ctheta;
  result[2] = b[2];

  return (0);
}



/*  The next few routine have to do with various basis manipulations */


/**********************************************************/
/** 
 * @brief      created an orthonomal basis from two vectors
 *
 * @param [in] double  u[]   the first of two vectors requiried to created a basis
 * @param [in] double  v[]   the second of two vectors needed to create a set of basis vectors
 * @param [out] struct basis *  basis_new   the new bais
 * @return   0  
 *
 *  create_basis creates an orthnormal basis with the first axis in the u direction, the
 *  y axis in the uv plane, and the the third axis in the perpendicular direction.  The
 *  system is right-handed 
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
create_basis (u, v, basis_new)
     double u[], v[];
     struct basis *basis_new;
{
  int i;
  double x[3], y[3], z[3];
  double mu_x;
  double dot ();

  for (i = 0; i < 3; i++)
  {
    x[i] = u[i];
    y[i] = v[i];
  }
  if (renorm (x, 1.) || renorm (y, 1.))
  {
    Log ("Problem creating basis: Either u or v had length 0\n");
    Log ("create_basis: u %e %e %e\n", u[0], u[1], u[2]);
    Log ("create_basis: v %e %e %e\n", v[0], v[1], v[2]);
    return (-1);
  }
  if ((mu_x = dot (x, y)) > 1. - EPS)
  {
    Log ("Problem creating basis u,v parallel\n");
    Log ("create_basis: u %e %e %e\n", u[0], u[1], u[2]);
    Log ("create_basis: v %e %e %e\n", v[0], v[1], v[2]);
    return (-1);
  }
  for (i = 0; i < 3; i++)
    y[i] -= mu_x * x[i];        /* y=y-dot(x,y)*x is orthogonal to x */

  renorm (y, 1.);


  cross (x, y, z);


  for (i = 0; i < 3; i++)
  {
    basis_new->a[i][0] = x[i];
    basis_new->a[i][1] = y[i];
    basis_new->a[i][2] = z[i];
  }

  return (0);
}



/**********************************************************/
/** 
 * @brief      Project a a vector in the locally rotated frame onto the unrotated fram
 *
 * @param [in] struct basis *  basis_from   basis which goes from the rotated to unrotated frame
 * @param [in] double  v_in[]   the vector in the rotate frame
 * @param [out] double  v_out[]   the vector in the unrotate frame
 * @return     0                   
 *
 *
 * ###Notes###
 *
 * basis from contains the direction cosines to go from rotated to unrotated fram
 *
 **********************************************************/

int
project_from (basis_from, v_in, v_out)
     struct basis *basis_from;  /* direction cosines to go from rotated to unrotated frame */
     double v_in[], v_out[];    /*v_in here is in rotated frame, v_out in unrotated frame */

{
  int i, j;
  for (i = 0; i < 3; i++)
  {
    v_out[i] = 0;
    for (j = 0; j < 3; j++)
    {
      v_out[i] += basis_from->a[i][j] * v_in[j];
    }
  }
  return (0);
}


/**********************************************************/
/** 
 * @brief       projects a vector in the unrotated frame onto the rotated frame
 *
 * @param [in] struct basis *  basis_from  basis which goes from the rotated to unrotated frame
 * @param [in] double  v_in[]   vector in the unrotated rotated frame
 * @param [out] double  v_out[]   vector in the rotated frame
 * @return     Always returns 0
 *
 *
 * ###Notes###
 *
 * basis_from contains the directon cosines to go from rotated to unrotated frame
 *
 **********************************************************/

int
project_to (basis_from, v_in, v_out)
     struct basis *basis_from;  /* direction cosines to go from rotated to unrotated frame */
     double v_in[], v_out[];    /*v_in here is in unrotated frame, v_out in rotated frame */

{
  int i, j;
  for (i = 0; i < 3; i++)
  {
    v_out[i] = 0;
    for (j = 0; j < 3; j++)
    {
      v_out[i] += basis_from->a[j][i] * v_in[j];
    }
  }
  return (0);
}




/**********************************************************/
/** 
 * @brief      Project a vector in one rotated frame to another rotated frame
 *
 * @param [in] struct basis *  basis_from   The basis which goes from the rotated to the unrotated frame
 * @param [in] struct basis *  basis_to   The basis which goes from the unrotated to another rotated frame
 * @param [in] double  v_from[]  A vector in the first frame
 * @param [out] double  v_to[]   The same vector in the second frame
 * @return     0               
 *
 * 
 *  projects a vector in one rotated frame to another rotated frame.  
 *
 * ###Notes###
 *
 *  a[i][j] is always the project from a rotated frame to the unrotated (Master) 
 *  frame.  The projection in the other direction is just the transpose of this
 *
 **********************************************************/

int
reorient (basis_from, basis_to, v_from, v_to)
     struct basis *basis_from, *basis_to;
     double v_from[], v_to[];
{
  double a[3][3];
  int i, j, k;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      a[i][j] = 0;
      for (k = 0; k < 3; k++)
      {
        a[i][j] += basis_to->a[k][i] * basis_from->a[k][j];
      }
    }
  }
  for (i = 0; i < 3; i++)
  {
    v_to[i] = 0;
    for (j = 0; j < 3; j++)
    {
      v_to[i] += a[i][j] * v_from[j];
    }
  }
  return (0);
}
