// $Id: srcfil_s.cc,v 1.1 2004-11-20 15:10:38 mcneile Exp $
/*! \file
 *  \brief Point source construction
 */

#include "chromabase.h"
#include "meas/hadron/srcfil.h"

using namespace QDP;

//! Fill a specific color and spin index with 1.0
/*!
 * \ingroup hadron
 *
 * This routine is specific to Staggered fermions! 
 *
 * \param a      Source fermion
 * \param coord  Lattice coordinate
 * \param color_index  Color index
 */

void srcfil(LatticeStaggeredFermion& a, multi1d<int>& coord, int color_index)
{
  if (color_index >= Nc || color_index < 0)
    QDP_error_exit("invalid color index", color_index);

  Real one = 1;
  Complex sitecomp = cmplx(one,0);
  ColorVector sitecolor = zero;
  StaggeredFermion sitefield = zero;

  const int spin_index = 0 ; 

  // Put [1,0] into the fermion a 
  pokeSite(a, 
	   pokeSpin(sitefield,
		    pokeColor(sitecolor,sitecomp,color_index),
		    spin_index), 
	   coord);
}
