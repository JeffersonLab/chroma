// $Id: walfil_s.cc,v 3.0 2006-04-03 04:59:06 edwards Exp $
/*! \file
 *  \brief Wall source construction
 */

#include "chromabase.h"
#include "meas/sources/walfil_s.h"

namespace Chroma {

//! Fill a specific color and spin index with 1.0 on a wall
/*!
 * \ingroup sources
 *
 * This routine is specific to Staggered fermions! 
 *
 * Fill a specific color index with 1.0, on sites at the corners 
 * of the cubes in a slice
 * 
 * 
 * \param a            Source fermion (write)
 * \param slice        time slice
 * \param mu           direction of slice
 * \param color_index  Color index
 * \param src_index    Index which defines which corner of a cube on
 *                     the source time slice you want your source to 
 *		       be on. The mapping from src_index to site is 
 *                     lexicographic, i.e: 0 is (0,0,0), 1 is (1,0,0),
 *		       2 is (0,1,0), 3 is (1,1,0), 4 is (0,0,1),
 *                     5 is (1,0,1), 6 is (0,1,1) and 7 is (1,1,1).
 * 
 * This is probably not the cleverest way to do this and realistically
 * you are not interested in all the sources at once so you have to be
 * careful to call this routine with the right index.   
 */

void walfil(LatticeStaggeredFermion& a, int slice, int mu, int color_index, int src_index)
{
  START_CODE();

  if ( color_index >= Nc )
    QDP_error_exit("Color index out of bounds", color_index, Nc);

  if ( (slice % 2) != 0 )
    QDP_error_exit("Require an even valued slice", slice);

  /* 
   * Write ONE to the required coordinates of the appropriate slice. 
   * Note: staggered fermions are represented as "Ns=1" spinors! 
   */
  // Write ONE to all field
  int spin_index = 0;
  Real one = 1;
  Complex sitecomp = cmplx(one,0);
  ColorVector sitecolor = zero;
  StaggeredFermion sitefield = zero;

  pokeSpin(sitefield,
	   pokeColor(sitecolor,sitecomp,color_index),
	   spin_index);


  // Narrow the context to the desired slice and to even coordinates
  LatticeBoolean ltest = (Layout::latticeCoordinate(mu) == slice);

  switch(src_index){
  case 0:
    for(int m = 0; m < Nd; ++m)
      if ( m != mu )
        ltest &= (  (Layout::latticeCoordinate(m) % 2) == 0);
    break;
  
  case 1:
    for(int m = 0; m < Nd-1; ++m){
      if ( m == 0 )
        ltest &= (  (Layout::latticeCoordinate(m) % 2) == 1);
      else
        ltest &= (  (Layout::latticeCoordinate(m) % 2) == 0);
    }
  break;

  case 2:
    for(int m = 0; m < Nd-1; ++m){
      if ( m == 1 )
        ltest &= (  (Layout::latticeCoordinate(m) % 2) == 1);
      else
        ltest &= (  (Layout::latticeCoordinate(m) % 2) == 0);
    }
  break;

  case 3:
    for(int m = 0; m < Nd-1; ++m){
      if ( m == 2 )
        ltest &= (  (Layout::latticeCoordinate(m) % 2) == 0);
      else
        ltest &= (  (Layout::latticeCoordinate(m) % 2) == 1);
    }
  break;

  case 4:
    for(int m = 0; m < Nd-1; ++m){
      if ( m == 2 )
        ltest &= (  (Layout::latticeCoordinate(m) % 2) == 1);
      else
        ltest &= (  (Layout::latticeCoordinate(m) % 2) == 0);
    }
  break;

  case 5:
    for(int m = 0; m < Nd-1; ++m){
      if ( m == 1 )
        ltest &= (  (Layout::latticeCoordinate(m) % 2) == 0);
      else
        ltest &= (  (Layout::latticeCoordinate(m) % 2) == 1);
    }
  break;

  case 6:
    for(int m = 0; m < Nd-1; ++m){
      if ( m == 0 )
        ltest &= (  (Layout::latticeCoordinate(m) % 2) == 0);
      else
        ltest &= (  (Layout::latticeCoordinate(m) % 2) == 1);
    }
  break;

  case 7:
    for(int m = 0; m < Nd-1; ++m)
        ltest &= (  (Layout::latticeCoordinate(m) % 2) == 1);
  break;

  default:
    QDPIO::cerr << "walfil_s: There are only 8 corners of a cube! " << endl;  
    QDP_abort(1);
  };

  // Write onto the appropriate slice
  LatticeStaggeredFermion tmp;
  tmp = sitefield;  // QDP (not installed version) now supports   construct OLattice = OScalar

  a = where(ltest, tmp, LatticeStaggeredFermion(zero));

  END_CODE();
}

}  // end namespace Chroma
