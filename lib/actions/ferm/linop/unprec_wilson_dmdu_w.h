#ifndef unprec_wilson_dmdu_h
#define unprec_wilson_dmdu_h

//! Linear Operator to evaluate dM/dU where M is the Unpreconditioned
//  Wilson Linear Operator
#include "chromabase.h"
#include "linearop.h"


using namespace QDP;
using namespace Chroma;

namespace Chroma {

  //! Unpreconditioned Wilson dM/dU
  /*! \ingroup Linear Operator
   *
   * This routine is specific to Wilson Fermions
   *
   * This subroutine applies the derivative of the Unprecontitioned
   * Wilson Operator with respect to the gauge fields to a vector 
   * Psi,
   * 
   *  isign = PLUS: 
   *
   *  Chi = dM/dU 
   *      = -(1/2) [  ( 1 - gamma_mu ) psi(x+mu) 
   *                - ( 1 + gamma_mu ) psi(x-mu)  ]
   *
   *  isign = MINUS:
   * 
   *  Chi = dM^dagger/dU 
   *      = -(1/2) [  ( 1 + gamma_mu ) psi(x+mu)
   *                - ( 1 - gamma_mu ) psi(x-mu)  ]
   *
   * Amusingly enough this is a linearop that does not depend 
   * on the gauge fields....
   *
   */
  class UnprecWilsondMdU : public LinearOperator<LatticeFermion>
  {
  public:

    //! Full constructor
    UnprecWilsondMdU(const int mu_) : mu(mu_) {}

    //! Destructor is automatic
    ~UnprecWilsondMdU() {}

    //! Unpreconditioned so defined on the whole lattice
    const OrderedSubset& subset() const { return all; }

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

  private:
    int mu;
    
  }; // End Class Definition

}; // End namespace Chroma


#endif
