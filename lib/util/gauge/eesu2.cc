// $Id: eesu2.cc,v 3.2 2009/11/14 20:01:46 eneil Exp $
/*! \file
 *  \brief Exactly exponentiate a SU(2) lie algebra element
 */

#include "chromabase.h"
#include "util/gauge/eesu2.h"
#include "util/gauge/su2extract.h"
#include "util/gauge/sunfill.h"

namespace Chroma 
{
  //! Exactly exponentiate a SU(2) lie algebra element
  /*!
   * \ingroup gauge
   *
   *  \param m        LatticeColorMatrix          (Modify)
   */
  void eesu2(LatticeColorMatrix& m)
  {
    START_CODE();

    if ( Nc != 2 )
    {
      QDPIO::cerr << __func__ << ": can only handle SU(2)" << endl;
      QDP_abort(1);
    }

    /* Extract components r_k proportional to SU(2) submatrix su2_index */
    /* from the "SU(Nc)" matrix V. The SU(2) matrix is parametrized in the */
    /* sigma matrix basis. */
    /* NOTE: r_0 should be 0 */
    multi1d<LatticeReal> r(4);
    su2Extract(r, m, 0, all);

    r[0]=Real(0.5)*r[0];
    r[1]=Real(0.5)*r[1];
    r[2]=Real(0.5)*r[2];
    r[3]=Real(0.5)*r[3];

    LatticeReal r_l = sqrt(r[1]*r[1] + r[2]*r[2] + r[3]*r[3]);
  
    // Normalize
    LatticeBoolean lbtmp = r_l > fuzz;
    LatticeReal lftmp = sin(r_l) / where(lbtmp, r_l, LatticeReal(1));

    /* Euler expand SU(2) lie algebra element to SU(2) lie group element */
    /* Be careful to check for elements below fuzz, then set to identity */
        
    // Fill   (r[0]/r_l, -r[1]/r_l, -r[2]/r_l, -r[3]/r_l) for r_l > fuzz
    //  and   (1,0,0,0)  for sites with r_l < fuzz
    multi1d<LatticeReal> a(4);
    a[0] = where(lbtmp, cos(r_l), LatticeReal(1));
    a[1] = where(lbtmp, r[1] * lftmp, LatticeReal(0));
    a[2] = where(lbtmp, r[2] * lftmp, LatticeReal(0));
    a[3] = where(lbtmp, r[3] * lftmp, LatticeReal(0));

    /*
     * Now fill an SU(3) matrix V with the SU(2) submatrix su2_index
     * paramtrized by a_k in the sigma matrix basis.
     */
    sunFill(m, a, 0, all);

    END_CODE();
  }

}  // end namespace Chroma
