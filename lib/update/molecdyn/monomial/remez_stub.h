// -*- C++ -*-
/*! \file
 *  \brief Remez algorithm for finding nth roots
 */

#ifndef __remez_stub_h__
#define __remez_stub_h__

#include "chromabase.h"
#include "update/molecdyn/monomial/remez_coeff.h"

namespace Chroma
{

  //! Dummy class for case when gmp is not present
  /*! @ingroup monomial
   *
   */
  class RemezStub
  {
  public:
    RemezStub(const Real& lower, const Real& upper, long prec) 
    {
      QDPIO::cerr << "RemezStub: Remez algorithm not supported without at least GMP (Gnu Muli-Precision library)" << std::endl;
      QDP_abort(1);
    }
    ~RemezStub() {}
    void setBounds(const Real& lower, const Real& upper) {}
    const Real generateApprox(int num_degree, int den_degree, 
			      unsigned long power_num, unsigned long power_den) {return 0;}
    const Real generateApprox(int degree, 
			      unsigned long power_num, unsigned long power_den) {return 0;}

    //! Return the partial fraction expansion of the approximation x^(pnum/pden)
    RemezCoeff_t getPFE()
    {
      QDP_error_exit("RemezStub not implemented");
      return RemezCoeff_t();
    }

    //! Return the partial fraction expansion of the approximation x^(-pnum/pden)
    RemezCoeff_t getIPFE()
    {
      QDP_error_exit("RemezStub not implemented");
      return RemezCoeff_t();
    }

    //! Given a partial fraction expansion, evaluate it at x
    Real evalPFE(const Real& x,  const RemezCoeff_t& coeff)
    {QDP_error_exit("RemezStub not implemented"); return Real(0);}

  };
}  // end namespace Chroma

#endif  // Include guard



