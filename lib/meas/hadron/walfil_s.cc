// $Id: walfil_s.cc,v 1.4 2003-12-30 17:27:15 bjoo Exp $
/*! \file
 *  \brief Wall source construction
 */

#include "chromabase.h"
#include "meas/hadron/walfil_s.h"

using namespace QDP;

//! Fill a specific color and spin index with 1.0 on a wall
/*!
 * \ingroup hadron
 *
 * This routine is specific to Staggered fermions! 
 *
 * Fill a specific color index with 1.0, on sites in a slice
 * where everything has even coordinates.
 *
 * \param a            Source fermion (write)
 * \param slice        time slice
 * \param mu           direction of slice
 * \param color_index  Color index
 */

void walfil(LatticeFermion& a, int slice, int mu, int color_index)
{
  START_CODE("walfil");

  if ( color_index >= Nc )
    QDP_error_exit("Color index out of bounds", color_index, Nc);

  if ( (slice % 2) != 0 )
    QDP_error_exit("Require an even valued slice", slice);

  /* 
   * Write ONE to the even coordinates of the appropriate slice. 
   * Note: staggered fermions are represented as "Ns=1" spinors! 
   */
  // Write ONE to all field
  int spin_index = 0;
  Real one = 1;
  Complex sitecomp = cmplx(one,0);
  ColorVector sitecolor = zero;
  Fermion sitefield = zero;

  pokeSpin(sitefield,
	   pokeColor(sitecolor,sitecomp,color_index),
	   spin_index);


  // Narrow the context to the desired slice and to even coordinates
  LatticeBoolean ltest = (Layout::latticeCoordinate(mu) == slice);

  for(int m = 0; m < Nd; ++m)
    if ( m != mu )
      ltest &= (  (Layout::latticeCoordinate(m) % 2) == 0);

  // Write onto the appropriate slice
  LatticeFermion tmp;
  tmp = sitefield;  // QDP (not installed version) now supports   construct OLattice = OScalar

  a = where(ltest, tmp, LatticeFermion(zero));

  END_CODE("walfil");
}
