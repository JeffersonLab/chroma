// -*- C++ -*-
// $Id: remez_stub.h,v 1.2 2005-02-02 16:27:55 edwards Exp $
/*! \file
 *  \brief Remez algorithm for finding nth roots
 */

#ifndef __remez_stub_h__
#define __remez_stub_h__

#include "chromabase.h"

namespace Chroma
{

  //! Dummy class for case when gmp is not present

  class RemezStub
  {
  public:
    RemezStub(const Real& lower, const Real& upper, long prec) 
    {
      QDPIO::cerr << "RemezStub: Remez algorithm not supported without at least QMP" << endl;
      QDP_abort(1);
    }
    ~RemezStub() {}
    void setBounds(const Real& lower, const Real& upper) {}
    const Real generateApprox(int num_degree, int den_degree, 
			      unsigned long power_num, unsigned long power_den) {return 0;}
    const Real generateApprox(int degree, 
			      unsigned long power_num, unsigned long power_den) {return 0;}
    int getPFE(multi1d<Real>& res, multi1d<Real>& pole, Real& norm) {return 0;}
    int getIPFE(multi1d<Real>& res, multi1d<Real>& pole, Real& norm) {return 0;}

    //! Given a partial fraction expansion, evaluate it at x
    Real evalPFE(const Real& x, const multi1d<Real>& res, const multi1d<Real>& pole)
      {QDP_error_exit("%s: not implemented",__func__);}

  };
}  // end namespace Chroma

#endif  // Include guard



