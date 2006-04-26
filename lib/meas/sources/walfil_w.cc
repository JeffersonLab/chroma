// $Id: walfil_w.cc,v 3.1 2006-04-26 03:42:18 edwards Exp $
/*! \file
 *  \brief Wall source construction
 */

#include "chromabase.h"
#include "meas/sources/walfil_w.h"

namespace Chroma 
{

  //! Fill a specific color and spin index with 1.0 on a wall
  /*!
   * \ingroup sources
   *
   * This routine is specific to Wilson fermions! 
   *
   * \param a            Source fermion
   * \param slice        time slice
   * \param mu           direction of slice
   * \param color_index  Color index
   * \param spin_index   Spin index
   */

  void walfil(LatticeFermion& a, int slice, int mu, int color_index, int spin_index)
  {
    START_CODE();

    if (color_index >= Nc || color_index < 0)
      QDP_error_exit("invalid color index", color_index);

    if (spin_index >= Ns || spin_index < 0)
      QDP_error_exit("invalid spin index", spin_index);

    // Write ONE to all field
    Real one = 1;
    Complex sitecomp = cmplx(one,0);
    ColorVector sitecolor = zero;
    Fermion sitefield = zero;

    pokeSpin(sitefield,
	     pokeColor(sitecolor,sitecomp,color_index),
	     spin_index);

    // Narrow the context to the desired slice.
    LatticeFermion tmp;
    tmp = sitefield;  // QDP (not installed version) now supports   construct OLattice = OScalar

    a = where(Layout::latticeCoordinate(mu) == slice, tmp, LatticeFermion(zero));

    END_CODE();
  }

}  // end namespace Chroma
