// $Id: walfil_s.cc,v 1.1 2003-12-15 04:23:32 edwards Exp $
/*! \file
 *  \brief Wall source construction
 */

#include "chromabase.h"
#include "meas/hadron/walfil.h"

using namespace QDP;

//! Fill a specific color and spin index with 1.0 on a wall
/*!
 * \ingroup hadron
 *
 * Fill a specific color index with 1.0, on sites in a slice
 * where everything has even coordinates.
 *
 * \param a            Source fermion (write)
 * \param slice        time slice
 * \param mu           direction of slice
 * \param color_index  Color index
 * \param spin_index   Spin index
 */

void walfil(LatticeFermion& a, int slice, int mu, int color_index, int spin_index)
{
  START_CODE("walfil");

  if ( color_index >= Nc )
    QDP_error_exit("Color index out of bounds", color_index, Nc);

  if ( spin_index != 0 )
    QDP_error_exit("Spin index out of bounds", spin_index, Ns);

  if ( (slice % 2) != 0 )
    QDP_error_exit("Require an even valued slice", slice);

  /* 
   * Write ONE to the even coordinates of the appropriate slice. 
   * Note: staggered fermions are represented as "Ns=1" spinors! 
   */
  // Write ONE to all field
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
      ltest &= (mod(Layout::latticeCoordinate(m), 2) == 0);

  // Write onto the appropriate slice
  a = where(ltest, LatticeFermion(sitefield), zero);

  END_CODE("walfil");
}
