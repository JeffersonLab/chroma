// $Id: srcfil.cc,v 1.1 2003-02-15 04:08:02 edwards Exp $
/*! \file
 *  \brief Point source construction
 */

#include "qdp.h"
#include "szin.h"

using namespace QDP;

//! Fill a specific color and spin index with 1.0
/*!
 * This routine is specific to Wilson fermions! 
 *
 * \param a      Source fermion
 * \param coord  Lattice coordinate
 * \param color_index  Color index
 * \param spin_index   Spin index
 */

void srcfil(LatticeFermion& a, multi1d<int>& coord, int color_index, int spin_index)
{
  if (color_index >= Nc || color_index < 0)
    QDP_error_exit("invalid color index", color_index);

  if (spin_index >= Ns || spin_index < 0)
    QDP_error_exit("invalid spin index", spin_index);

  Real one = 1;
  Complex sitecomp = cmplx(one,0);
  ColorVector sitecolor = zero;
  Fermion sitefield = zero;

  // Put [1,0] into the fermion a 
  pokeSite(a, 
	   pokeSpin(sitefield,
		    pokeColor(sitecolor,sitecomp,color_index),
		    spin_index), 
	   coord);
}


//! Fill a specific color index with 1.0
/*!
 * This routine is specific to Wilson fermions! 
 *
 * \param a      Source lattice Color Vector
 * \param coord  Lattice coordinate
 * \param color_index  Color index
 */

void srcfil(LatticeColorVector& a, multi1d<int>& coord, int color_index)
{
  if (color_index >= Nc || color_index < 0)
    QDP_error_exit("invalid color index", color_index);

  Real one = 1;
  Complex sitecomp = cmplx(one,0);
  ColorVector sitecolor = zero;

  // Put [1,0] into the fermion a 
  pokeSite(a, 
	   pokeColor(sitecolor,sitecomp,color_index),
	   coord);
}


