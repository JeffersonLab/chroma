// $Id: unprec_dwf4d_linop_w.cc,v 3.1 2006-07-03 15:26:09 edwards Exp $
/*! \file
 *  \brief Unpreconditioned projected DWF operator to 4D
 */

#include "chromabase.h"
#include "actions/ferm/invert/invcg2_array.h"
#include "actions/ferm/linop/unprec_dwf4d_linop_w.h"
#include "actions/ferm/linop/dwffld_w.h"

using namespace QDP::Hints;

namespace Chroma
{

  // Apply unpreconditioned linear operator
  template<>
  void 
  UnprecDWF4DLinOp<LatticeFermion>::operator()(LatticeFermion& chi, 
					       const LatticeFermion& psi, 
					       enum PlusMinus isign) const
  {
    START_CODE();

    SystemSolverResults_t res;
    const int  N5 = size();   // array size better match
  
    // Initialize the 5D fields
    multi1d<LatticeFermion> psi5(N5); moveToFastMemoryHint(psi5);
    //  psi5 = (psi,0,0,0,..,0)^T
    psi5 = zero;
    psi5[0] = psi;

    // tmp5 = P . psi5
    multi1d<LatticeFermion> tmp5(N5); moveToFastMemoryHint(tmp5);
    DwfFld(tmp5, psi5, PLUS);

    multi1d<LatticeFermion> chi5(N5); moveToFastMemoryHint(chi5);

    switch(isign)
    {
    case PLUS:
    {
      // chi5 = D5(m_q) . tmp5 =  D5(m_q) . P . (chi,0,0,..,0)^T 
      (*D)(chi5, tmp5, PLUS);

      // Solve  D5(1) . psi5 = chi5
      psi5 = chi5;
  
//      case CG_INVERTER: 
      {
	/* tmp5 = M_dag(u) * chi5 */
	(*PV)(tmp5, chi5, MINUS);
	  
	/* psi5 = (M^dag * M)^(-1) chi5 */
	res = InvCG2(*PV, tmp5, psi5, invParam.RsdCG, invParam.MaxCG);
      }
    }
    break;

    case MINUS:
    {
      // Solve  D5(1) . psi5 = tmp5
      psi5 = tmp5;
  
//      case CG_INVERTER: 
      {
	/* psi5 = (M^dag * M)^(-1) tmp5 */
	res = InvCG2(*PV, tmp5, psi5, invParam.RsdCG, invParam.MaxCG);

	/* tmp5 = M_dag(u) * psi5 */
	(*PV)(tmp5, psi5, PLUS);
      }

      // psi5 = D5(m_q)^dag . tmp5 
      (*D)(psi5, tmp5, MINUS);

    }
    break;

    default:
      QDP_error_exit("dwf4d: unknown isign");
    }
    
    if ( res.n_count == invParam.MaxCG )
      QDP_error_exit("no convergence in the inverter", res.n_count);
  
    // Project out first slice after  chi <- chi5 <- P^(-1) . psi5
    DwfFld(chi5, psi5, MINUS);
    chi = chi5[0];

    END_CODE();
  }

}
