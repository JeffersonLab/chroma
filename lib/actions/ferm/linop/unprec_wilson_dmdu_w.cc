#include "chromabase.h"
#include "actions/ferm/linop/unprec_wilson_dmdu_w.h"

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
   *      = -(1/2) [  ( 1 - gamma_mu ) psi(x+mu) ]
   *
   *  isign = MINUS:
   * 
   *  Chi = dM^dagger/dU 
   *      = -(1/2) [  ( 1 + gamma_mu ) psi(x+mu) ]
   *
   * Amusingly enough this is a linearop that does not depend 
   * on the gauge fields....
   *
   * The mu is done at creation...
   *
   */
  //! Apply the operator onto a source vector
  void 
  UnprecWilsondMdU::operator() (LatticeFermion& chi, 
				const LatticeFermion& psi, 
				enum PlusMinus isign) const
  {
    START_CODE();

    LatticeFermion tmp;
    LatticeFermion gamma_mu_psi;
    

    // Get gamma_mu * psi
    gamma_mu_psi = Gamma(1 << mu)*psi;

    // Construct the right derivative term...
    if( isign == PLUS ) {
      
      // Undaggered:
      // ( 1 - gamma_mu ) psi
      tmp = psi -  gamma_mu_psi;
    }
    else { 
      
      // Daggered:
      // ( 1 + gamma_mu) psi
      tmp = psi + gamma_mu_psi;
    }
    
    // Do the shift here for delta_y, x+mu
    chi = -Real(0.5)*shift(tmp, FORWARD, mu);
    
    
    END_CODE();
  }

}; // End Namespace Chroma
