// -*- C++ -*-
//  $Id: su3proj.cc,v 3.1 2007-02-22 21:11:50 bjoo Exp $
/*! \file
 *  \ingroup gauge
 *  \author Subsetting added by A. Hart
 *  \param[in,out] u        the projected SU(3) Matrix (Modify)
 *  \param[in] w            matrix against which to maximize (Read)
 *  \param[in] su2_index    SU(2) subgroup index (Read)
 *  \param[in] mstag        An (un)ordered subset of sites (Read)
 *  \brief Project a GL(3,C) color matrix onto SU(3)
 *
 *  Project a complex Nc x Nc matrix W onto SU(Nc) by maximizing Tr(VW)
 */

#include "chromabase.h"
#include "util/gauge/su2extract.h"
#include "util/gauge/sunfill.h"
#include "util/gauge/su3proj.h"

namespace Chroma 
{

  template<typename S>
  void su3proj_t(LatticeColorMatrix& u, 
		 const LatticeColorMatrix& w, 
		 int su2_index,
		 const S& mstag)
  {
    START_CODE();

    // V = U*W
    LatticeColorMatrix v;
    v[mstag] = u * w;

    /*
     * Extract components r_k proportional to SU(2) submatrix su2_index
     * from the "SU(3)" matrix V. The SU(2) matrix is parameterized in the
     * sigma matrix basis.
     */
    multi1d<LatticeReal> r(4);
    su2Extract(r, v, su2_index, mstag);

    /*
     * Now project onto SU(2)
     */
    LatticeReal r_l = 0;
    r_l[mstag] = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + r[3]*r[3]);

    // Normalize
    LatticeBoolean lbtmp = false;
    lbtmp[mstag] = r_l > fuzz;
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
    sunFill(v, a, su2_index, mstag);

    // U = V*U
    LatticeColorMatrix tmp;
    tmp[mstag] = v * u;
    u[mstag] = tmp;

    END_CODE();
  }

  void su3proj(LatticeColorMatrix& u, 
	       const LatticeColorMatrix& w, 
	       int su2_index)
  {
    su3proj_t(u, w, su2_index, all);
  }

  void su3proj(LatticeColorMatrix& u, 
	       const LatticeColorMatrix& w, 
	       int su2_index,
	       const Subset& mstag)
  {
    su3proj_t(u, w, su2_index, mstag);
  }


}  // end namespace Chroma

