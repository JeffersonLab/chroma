//  $Id: su3proj.cc,v 1.5 2003-12-06 20:59:21 edwards Exp $
/*! \file
 *  \brief Project a complex Nc x Nc matrix W onto SU(Nc) by maximizing Tr(VW)
 */

#include "chromabase.h"
#include "util/gauge/su3proj.h"
#include "util/gauge/su2extract.h"
#include "util/gauge/sunfill.h"

using namespace QDP;

//! Project a GL(3,C) color matrix onto SU(3)
/*!
 * \ingroup gauge
 *
 * Arguments:
 *
 *  \param u            the projected SU(3) Matrix (Modify)
 *  \param w            matrix against which to maximize (Read)
 *  \param su2_index    SU(2) subgroup index (Read)
 */

void su3proj(LatticeColorMatrix& u, const LatticeColorMatrix& w, int su2_index)
{
  START_CODE("su3proj");

  // V = U*W
  LatticeColorMatrix v = u * w;

  /*
   * Extract components r_k proportional to SU(2) submatrix su2_index
   * from the "SU(3)" matrix V. The SU(2) matrix is parameterized in the
   * sigma matrix basis.
   */
  multi1d<LatticeReal> r(4);
  su2Extract(r, v, su2_index, all);

  /*
   * Now project onto SU(2)
   */
  LatticeReal r_l = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + r[3]*r[3]);

  // Normalize
  LatticeBoolean lbtmp = r_l > fuzz;
  LatticeReal lftmp = 1.0 / where(lbtmp, r_l, LatticeReal(1));

  // Fill   (r[0]/r_l, -r[1]/r_l, -r[2]/r_l, -r[3]/r_l) for r_l > fuzz
  //  and   (1,0,0,0)  for sites with r_l < fuzz
  multi1d<LatticeReal> a(4);
  a[0] = where(lbtmp, r[0] * lftmp, LatticeReal(1));
  a[1] = where(lbtmp, -(r[1] * lftmp), LatticeReal(0));
  a[2] = where(lbtmp, -(r[2] * lftmp), LatticeReal(0));
  a[3] = where(lbtmp, -(r[3] * lftmp), LatticeReal(0));

  /*
   * Now fill an SU(3) matrix V with the SU(2) submatrix su2_index
   * paramtrized by a_k in the sigma matrix basis.
   */
  sunFill(v, a, su2_index, all);

  // U = V*U
  LatticeColorMatrix tmp = v * u;
  u = tmp;

  END_CODE("su3proj");
}
