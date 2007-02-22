// $Id: su3over.cc,v 3.1 2007-02-22 21:11:49 bjoo Exp $
/*! \file
 *  \brief Do one SU(2) subgroup microcanonical overrelaxation update of SU(Nc)
 */

#include "chromabase.h"
#include "update/heatbath/su3over.h"
#include "util/gauge/su2extract.h"
#include "util/gauge/sunfill.h"

namespace Chroma 
{

  //! Do one SU(2) subgroup microcanonical overrelaxation update of SU(Nc) matrix
  /*!
   * \ingroup heatbath
   *
   * Do one SU(2) subgroup microcanonical overrelaxation update of SU(Nc)
   * matrix U keeping action tr(U*W) constant.
   *
   * \param u            field to be updated ( Modify )
   * \param w            "staple" field in the action ( Read )
   * \param su2_index    SU(2) subgroup index ( Read ) 
   * \param sub          Subset for operations ( Read ) 
   */

  void su3over(LatticeColorMatrix& u,
	       const LatticeColorMatrix& w,
	       int su2_index,
	       const Subset& sub)
  {
    START_CODE();
                      
    /* V = U*W */
    LatticeColorMatrix v;
    v[sub] = u * w;
  
    /* Extract components r_k proportional to SU(2) submatrix su2_index */
    /* from the "SU(Nc)" matrix V. The SU(2) matrix is parametrized in the */
    /* sigma matrix basis. */
    multi1d<LatticeReal> r(4);
    su2Extract(r, v, su2_index, sub);
  
    /*
     * Now project onto SU(2)
     */
    LatticeReal r_l;
    r_l[sub] = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + r[3]*r[3]);

    // Normalize
    LatticeBoolean lbtmp;
    lbtmp[sub] = r_l > fuzz;
    LatticeReal lftmp;
    lftmp[sub] = 1.0 / where(lbtmp, r_l, LatticeReal(1));

    // Fill   (r[0]/r_l, -r[1]/r_l, -r[2]/r_l, -r[3]/r_l) for r_l > fuzz
    //  and   (1,0,0,0)  for sites with r_l < fuzz
    multi1d<LatticeReal> a(4);
    a[0][sub] = where(lbtmp, r[0] * lftmp, LatticeReal(1));
    a[1][sub] = where(lbtmp, -(r[1] * lftmp), LatticeReal(0));
    a[2][sub] = where(lbtmp, -(r[2] * lftmp), LatticeReal(0));
    a[3][sub] = where(lbtmp, -(r[3] * lftmp), LatticeReal(0));

    /* Microcanonical updating matrix is the square of this */
    r[0][sub] = a[0]*a[0] - a[1]*a[1] - a[2]*a[2] - a[3]*a[3];
  
    a[0][sub] *= 2;
    r[1][sub] = a[0]*a[1];
    r[2][sub] = a[0]*a[2];
    r[3][sub] = a[0]*a[3];
  
    /*
     * Now fill an SU(3) matrix V with the SU(2) submatrix su2_index
     * paramtrized by a_k in the sigma matrix basis.
     */
    sunFill(v, r, su2_index, sub);

    // U = V*U
    LatticeColorMatrix tmp;
    tmp[sub] = v * u;
    u[sub] = tmp;
  
    END_CODE();
  }

}  // end namespace Chroma
