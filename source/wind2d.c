
/***********************************************************/
/** @file  wind2d.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  Routines that are used to define the wind initially
 *
 * The file also contains a number of utility routines, several of
 * which are steering routines to coordinate system specific
 * routines.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

int wig_n;
double wig_x, wig_y, wig_z;

/**********************************************************/
/**
 * @brief      locates the element in wmain associated with a position 
 *
 * @param [in] int  ndom   The domain number for the search
 * @param [in] double  x[]   The position
 * @return     where_in_grid normally  returns the element of wmain 
 * associated that position.  
 * 
 * If the positions is in the grid (of the domains) this will be 
 * a positive integer.  If the position is not in the grid of the
 * domain, a negative number will be returned.  The negative number
 * will be -1 if the positionis inside the grid, or -2 if it is outside
 * the grid.,
 *
 * @details

 * where_in_grid is mainly a steering routine that calls other various
 * coordinate system specific routines.
 *
 * ### Notes ###
 *
 * Where_in grid does not tell you whether the position is in the wind!
 *
 * What one means by inside or outside the grid may well be different
 * for different coordinate systems.
 *
 * This is one of the routines in Python that has a long history, in that
 * it existed prior to the implementation of domains.  This is why it
 * returns the position in wmain, rather than the position in one of
 * the plasma domains.  It would make sense to revise this.
 *
 **********************************************************/

int
where_in_grid (ndom, x)
     int ndom;
     double x[];
{
  int n;
  double fx, fz;


  if (wig_x != x[0] || wig_y != x[1] || wig_z != x[2])  // Calculate if new position
  {

    if (zdom[ndom].coord_type == CYLIND)
    {
      n = cylind_where_in_grid (ndom, x);
    }
    else if (zdom[ndom].coord_type == RTHETA)
    {
      n = rtheta_where_in_grid (ndom, x);
    }
    else if (zdom[ndom].coord_type == SPHERICAL)
    {
      n = spherical_where_in_grid (ndom, x);
    }
    else if (zdom[ndom].coord_type == CYLVAR)
    {
      n = cylvar_where_in_grid (ndom, x, 0, &fx, &fz);
    }
    else
    {
      Error ("where_in_grid: Unknown coord_type %d for domain %d\n", zdom[ndom].coord_type, ndom);
      Exit (0);
      return (0);
    }


    /* Store old positions to short-circuit calculation if asked for same position more
       than once */
    wig_x = x[0];
    wig_y = x[1];
    wig_z = x[2];
    wig_n = n;
  }

  return (wig_n);
}

int ierr_vwind = 0;


/**********************************************************/
/**
 * @brief      finds the velocity vector v for the wind in cartesian
 * 	coordinates at the position of the photon p in the Observer frame.
 *
 * @param [in] int  ndom   The domain of interest
 * @param [in] PhotPtr  p   A photon
 * @param [out] double  v[]   The velocity at the position given by the photon
 * @return     0 	on successful completion, or a negative number if
 * photon position is not in the in the grid for the domain of interest.
 *
 * @details
 * This routine carries out a bilinear interpolation of the wind velocity.  The velocities
 * at the edges of the grid cells must have been defined elsewhere (e.g in wind_define).
 *
 * If the position is outside the wind grid, then the velocities will be returned as 0.
 *
 * ### Notes ###
 * For all of the 2d coordinate systems, the positions in the WindPtr array are
 * calculated in the xz plane.  As a result, the velocities in that array are
 * calcuated there as well.  Hence, to obtain the velocity in cartesion coordinates
 * one needs to do a simple rotation of the velocity vector.
 *
 * For spherical (1d) coordinates the situation is complex, the philosopy that
 * has been adopted is to let the program run as if the wind model is intrinsically
 * 2d.   The xyz positions of the grid are defined are defined so that
 * are in the xz plane at an angle of 45 degrees to the x axis.  This was done
 * so that if one used, say an sv wind, in spherical coordinates, one would select
 * the wind velocities at a plausible location.   But it implies that the velocities
 * that have been accepted are not purely radial as one might expect.  The choice
 * that has been made is to assume that v[1] is preserved as a non zero value in
 * making the projection.  An alternative might have been to preserve the angular
 * rotation rate.
 *
 * If we ever implement a real 3d coordinate system, one will need to look at
 * this routine again.
 *
 * The routine checks to see whether the position for which the velocity is needed
 * and if so short-circuits the calculation returning a stored value
 *
 * vwind_xyz expects the photon position to be in the Observer frame
 * It returns the velocity in the Observer frame.
 *
 **********************************************************/

#define NVWIND  3
int nvwind = 0;
int nvwind_last = 0;
struct vwind
{
  int iorder;
  double v[3], pos[3];
} xvwind[NVWIND];

int
vwind_xyz (ndom, p, v)
     int ndom;
     PhotPtr p;
     double v[];
{
  int i;
  double rho, r;
  double vv[3];
  double ctheta, stheta;
  double x, frac[4];
  int nn, nnn[4], nelem;
  int n;

  /* Check if the velocity for this position is in the buffer, and if so return that */
  for (n = 0; n < nvwind; n++)
  {
    if (xvwind[n].pos[0] == p->x[0] && xvwind[n].pos[1] == p->x[1] && xvwind[n].pos[2] == p->x[2])
    {
      v[0] = xvwind[n].v[0];
      v[1] = xvwind[n].v[1];
      v[2] = xvwind[n].v[2];
      return (0);
    }
  }

  if (ndom < 0 || ndom >= geo.ndomain)
  {
    Error ("vwind_xyz: Received invalid domain  %d\n", ndom);
  }

  coord_fraction (ndom, 0, p->x, nnn, frac, &nelem);

  for (i = 0; i < 3; i++)
  {

    x = 0;
    for (nn = 0; nn < nelem; nn++)
      x += wmain[nnn[nn]].v[i] * frac[nn];

    vv[i] = x;
  }

  rho = sqrt (p->x[0] * p->x[0] + p->x[1] * p->x[1]);

  if (zdom[ndom].coord_type == SPHERICAL)
  {                             // put the velocity into cylindrical coordinates on the xz axis
    x = sqrt (vv[0] * vv[0] + vv[2] * vv[2]);   //v_r
    r = length (p->x);
    vv[0] = rho / r * x;
    vv[2] = p->x[2] / r * x;
  }
  else if (p->x[2] < 0)         // For 2d coord syatems, velocity is reversed if the photon is in the lower hemisphere.
    vv[2] *= -1;

  if (rho == 0)
  {                             // Then we will not be able to project from a cylindrical ot a cartesian system
    Error ("vwind_xyz: Cannot determine an xyz velocity on z axis. Returnin 0,0,v[2]\n");
    v[0] = v[1] = 0;
    v[2] = vv[2];
    return (0);
  }

  /* Now we have v in cylindrical coordinates, but we would like it in cartesian coordinates.
     Note that could use project_from_cyl_xyz(p->x,vv,v) ?? */

  ctheta = p->x[0] / rho;
  stheta = p->x[1] / rho;
  v[0] = vv[0] * ctheta - vv[1] * stheta;
  v[1] = vv[0] * stheta + vv[1] * ctheta;
  v[2] = vv[2];

  /* Now populate the buffer */


  if (nvwind < NVWIND)
  {
    nvwind_last = nvwind;
    nvwind++;
  }
  else
  {
    nvwind_last = (nvwind_last + 1) % NVWIND;
  }


  xvwind[nvwind_last].pos[0] = p->x[0];
  xvwind[nvwind_last].pos[1] = p->x[1];
  xvwind[nvwind_last].pos[2] = p->x[2];
  xvwind[nvwind_last].v[0] = v[0];
  xvwind[nvwind_last].v[1] = v[1];
  xvwind[nvwind_last].v[2] = v[2];


  return (0);
}




/**********************************************************/
/**
 * @brief      calculates the divergence of the velocity at the center of all the grid cells.
 *
 * @param [in] int ndom  The domain of the cell
 * @param [in,out] WindPtr  w   An individual wind cell
 * @return     Always returns 0
 *
 *
 * @details
 * This is one of the initialization routines for the wind.
 * The divergence of velocity of the wind is used to calculated 
 * PdV cooling of a grid cell
 *
 * The results are stored in wmain[n_wind].div_v
 *
 * ### Notes ###
 *
 * The divergence, like other scalar quantities, does not
 * need to be "rotated" differently in different coordinate
 * systems
 *
 **********************************************************/

int
wind_div_v (int ndom, WindPtr cell)
{
  static int wind_div_err = -3; /* Stop reporting the error after 4 occurrences for cells not in the wind */

  /* Compute the divergence for at the center of the cell */
  cell->div_v = get_div_v_in_cmf_frame (ndom, cell->xcen);

  if (cell->div_v < 0 && (wind_div_err < 0 || cell->inwind == W_ALL_INWIND))
  {
    Error ("wind_div_v: div v %e negative in cell %d Domain %d. Major problem if inwind (%d) == 0\n", div, cell->nwind, ndom, cell->inwind);
    wind_div_err++;
  }

  return (0);
}

/**********************************************************/
/**
 * @brief      find the density of the wind at x
 *
 * @param [in] WindPtr  w   The entire wind
 * @param [in] double  x[]   A position
 * @return
 * The density at x, if the position is in the active region of the wind.
 * If the position is not in the active region of one of the domains, then
 * 0 is returned.  No error is reported for this
 *
 * @details
 * The routine first determines what domain the position is located in,
 * and then interpolates to find the density at a specific position
 *
 * ### Notes ###
 *
 *
 **********************************************************/

double
rho (w, x)
     WindPtr w;
     double x[];
{
  int n;
  double dd;
  double frac[4];
  int nn, nnn[4], nelem;
  int nplasma;
  int ndom, ndomain;

  if (where_in_wind (x, &ndomain) < 0)  //note that where_in_wind is independent of grid.
    return (0.0);

  if ((ndom = ndomain) < 0)
  {
    Error ("rho: Domain not fournd %d\n", ndomain);
    return (0.0);
  }

  n = coord_fraction (ndom, 1, x, nnn, frac, &nelem);

  if (n < 0)
  {
    dd = 0;
  }
  else
  {

    dd = 0;
    for (nn = 0; nn < nelem; nn++)
    {
      nplasma = w[nnn[nn]].nplasma;
      dd += plasmamain[nplasma].rho * frac[nn];
    }

  }

  return (dd);
}


#define NSTEPS 100

/**********************************************************/
/**
 * @brief      The routine calculates and then logs the mass loss rate in two ways
 *
 * @param [in] WindPtr  w   The entire wind domain
 * @param [in] double  z   A zheight above the disk to calculate the mass loss rate
 * @param [in] double  rmax   A radius at which to calculate the mass loss rate
 * @return     Always retruns 0
 *
 * @details
 * The routine calculates the mass loss rate in two ways, one in a plane with
 * a zheight of 0 above the disk (out to a distance rmax, and one at a spherical
 * radius rmax
 *
 * ### Notes ###
 *
 * ### Programming Comment ###
 * This routine has errors.  It is only correct if there is
 * a single domain (although it does provide a warning to this effect.
 * In particular rho(w,x) gives the correct answer
 * for rho regardless of domains, but ndom is set to 0 for vind(ndom,&p,v).
 * This should be fixed.  This is now #395
 *
 **********************************************************/

int
mdot_wind (w, z, rmax)
     WindPtr w;
     double z;                  // The height (usually small) above the disk at which mdot will be calculated
     double rmax;               // The radius at which mdot will be calculated
{
  struct photon p;
  double r, dr, rmin;
  double theta, dtheta;
  double den, rho ();
  double mdot, mplane, msphere;
  double x[3], v[3], q[3], dot ();
  int ndom;

  ndom = 0;

  Log ("For simplicity, mdot wind checks only carried out for domain 0\n");

/* Calculate the mass loss rate immediately above the disk */

  rmin = geo.rstar;

  if (rmax <= rmin)
  {
    Error ("mdot_wind: rmax %g is less than %g, so returning\n", rmax, rmin);
    return (0.0);
  }

  dr = (rmax - rmin) / NSTEPS;
  mdot = 0;
  p.x[1] = x[1] = 0.0;
  p.x[2] = x[2] = z;
  for (r = dr / 2; r < rmax; r += dr)
  {
    p.x[0] = x[0] = r;
    den = rho (w, x);
    vwind_xyz (ndom, &p, v);
    mdot += 2 * PI * r * dr * den * v[2];
  }

  mplane = 2. * mdot;

/* Calculate the mass loss rate in a sphere */

  dtheta = PI / (2 * NSTEPS);
  mdot = 0;
  p.x[1] = x[1] = 0.0;
  q[1] = 0;
  for (theta = dtheta / 2; theta < (PI / 2.); theta += dtheta)
  {
    p.x[0] = x[0] = r * sin (theta);
    p.x[2] = x[2] = r * cos (theta);
    q[0] = sin (theta);
    q[2] = cos (theta);
    den = rho (w, x);
    vwind_xyz (ndom, &p, v);
    mdot += 2 * PI * rmax * rmax * sin (theta) * dtheta * den * dot (v, q);
  }
  msphere = 2 * mdot;           // Factor of two because wind is in both hemispheres

  Log ("Wind Mdot Desired %e Plane %e Sphere(%e) %e\n", zdom[ndom].wind_mdot, mplane, rmax, msphere);

  return (0);
}



/**********************************************************/
/**
 * @brief      is simply will produce a postion at a random place in
 * 	a cell n
 *
 * @param [in] int  n   An element number in the wind structure
 * @param [out] double  x[]   The random location that is calculated
 * @return     Always returns 0
 *
 * @details
 * This is really just a driver routine for coordinate system specific
 * routines
 *
 * ### Notes ###
 *
 **********************************************************/

int
get_random_location (n, x)
     int n;                     // Cell in which to create position
     double x[];                // Returned position
{
  int ndom;

  ndom = wmain[n].ndom;

  if (zdom[ndom].coord_type == CYLIND)
  {
    cylind_get_random_location (n, x);
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {
    rtheta_get_random_location (n, x);
  }
  else if (zdom[ndom].coord_type == SPHERICAL)
  {
    spherical_get_random_location (n, x);
  }
  else if (zdom[ndom].coord_type == CYLVAR)
  {
    cylvar_get_random_location (n, x);
  }
  else
  {
    Error ("get_random_location: Don't know this coord_type %d\n", zdom[ndom].coord_type);
    Exit (0);
  }

  return (0);
}



/**********************************************************/
/**
 * @brief      zero out the portion of plasmamain that records
 * the number of scatters in each cell
 *
 * @return     Always returns 0
 *
 * @details
 *
 * ### Notes ###
 * The routine is called at the point where we begin to
 * calculate the detailed spectrum.  It's a little unclear
 * why it is simply not incorporated into spectrum_init
 *
 **********************************************************/

int
zero_scatters ()
{
  int n, j;

  for (n = 0; n < NPLASMA; n++) /*NSH 1107 - changed loop to only run over nions to avoid running over the
                                   end of the array after the arry was dynamically allocated in 78a */
  {
    for (j = 0; j < nions; j++)
    {
      plasmamain[n].scatters[j] = 0;
    }
  }

  return (0);
}


/**********************************************************/
/**
 * @brief      The routine basically just calls where_in_wind for each of 
 ' the 4 corners of a cell in one of the 2d coordinate systems.
 *
 * @param [in] int  n   the cell number in wmain
 * @return     the number of corners of the cell that
 * are in the wind
 *
 * @details
 * This is really just a driver routine which calls coordinate system specific
 * routines
 *
 * It is used to  standardize this check for calculations of 
 * the volume of a cell that is in the wind in a domain.
 *
 * ### Notes ###
 * This is one of the routines that makes explicity assumptions
 * about guard cells, and it is not entirely clear why
 *
 * This rutine is not called for a spherical coordinate system
 *
 **********************************************************/

int
check_corners_inwind (n)
     int n;
{
  int n_inwind;
  int i, j;
  DomainPtr one_dom;
  int ndom, ndomain;

  /* find the domain */
  ndom = wmain[n].ndom;
  one_dom = &zdom[ndom];

  wind_n_to_ij (ndom, n, &i, &j);

  n_inwind = 0;

  if (i < (one_dom->ndim - 2) && j < (one_dom->mdim - 2))
  {
    if (where_in_wind (wmain[n].x, &ndomain) == W_ALL_INWIND && ndom == ndomain)
      n_inwind++;
    if (where_in_wind (wmain[n + 1].x, &ndomain) == W_ALL_INWIND && ndom == ndomain)
      n_inwind++;
    if (where_in_wind (wmain[n + one_dom->mdim].x, &ndomain) == W_ALL_INWIND && ndom == ndomain)
      n_inwind++;
    if (where_in_wind (wmain[n + one_dom->mdim + 1].x, &ndomain) == W_ALL_INWIND && ndom == ndomain)
      n_inwind++;
  }

  return (n_inwind);
}




/**********************************************************/
/**
 * @brief check whether some quantities in a grid cell are
 * changing so rapidly in the grid cell that one might
 * wish to have smaller grid cells
 *
 * @return     Always returns 0
 *
 * The results are logged.
 *
 * @details
 * loops over the entire wind points and checks quantities
 * like velocity gradients and optical depths across cells for large changes.
 * If there are large changes in cells, it tells the user how many cells and
 * suggests that the grid may need to be modified.
 *
 *
 * ### Notes ###
 * Must be called after define wind so NPLASMA and densities are set up.
 *
 * The information that is checked is very basic and it is clear that
 * someone was thinking about doing somewhat better checks.
 *
 * Errors are not generated, or rather have been commented out, and
 * so one needs to consider whether this is useful or not.  If so there
 * should be errors.
 *
 **********************************************************/

int
check_grid ()
{
  int ndom, n;
  double lambda_t, nh;
  double delta_r, delta_x, delta_z, delta_vz, delta_vx;
  double v1[3], v2[3];
  WindPtr one;
  PlasmaPtr xplasma;
  int n_dv, n_tau;

  n_dv = n_tau = 0;

  for (n = 0; n < NPLASMA; n++)
  {
    /* get the relevant pointers for this cell */
    xplasma = &plasmamain[n];
    one = &wmain[xplasma->nplasma];
    ndom = one->ndom;

    /* Hydrogen density, ne should be roughly this */
    nh = xplasma->rho * rho2nh;

    /* thermal speed */
//OLD    vth = sqrt (1.5 * BOLTZMANN * xplasma->t_e / MPROT);

    /* sobolev length -- this could be used for a check but isn't yet */
//OLD    l_sob = vth / one->dvds_ave;

    /* get approximate cell extents */
    delta_z = 2.0 * (one->xcen[2] - one->x[2]);
    delta_x = 2.0 * (one->xcen[0] - one->x[0]);

    delta_r = sqrt (delta_x * delta_x + delta_z * delta_z);

    /* Thomson mean free path 1/sigma*nh */
    lambda_t = 1.0 / THOMPSON * nh;

    /* get the velocity at cell corner and cell edge in x and z */
    model_velocity (ndom, one->x, v1);
    model_velocity (ndom, one->xcen, v2);

    delta_vx = fabs (v1[0] - v2[0]);
    delta_vz = fabs (v1[1] - v2[1]);

    /* if we have errors then increment the error counters */
    if (delta_r / lambda_t > 1)
    {
      n_tau++;
      //Error("check_grid: optical depth may be high in cell %i domain %i\n",
      //       n, ndom);
    }

    if (delta_vx > 1e8 || delta_vz > 1e8)
    {
      n_dv++;
      //Error("check_grid: velocity changes by >1,000 km/s across cell %i domain %i\n",
      //       n, ndom);
    }

  }

  /* report problems to user in a summary error message */
  if (n_dv > 0)
    Error ("check_grid: velocity changes by >1,000 km/s in %i cells\n", n_dv);

  if (n_tau > 0)
    Error ("check_grid: optical depth may be high in %i\n", n_tau);

  if (n_dv > 0 || n_tau > 0)
    Error ("check_grid: some cells have large changes. Consider modifying zlog_scale or grid dims\n");

  return (0);
}
