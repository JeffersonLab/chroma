// $Id: unprec_dwf_dsdu_array_w.cc,v 1.1 2003-11-12 22:16:22 edwards Exp $
/*! \file
 *  \brief dS/dU_f for unpreconditioned domain-wall fermions
 */

#include "chromabase.h"
#include "fermact.h"
#include "actions/ferm/fermacts/unprec_dwf_fermact_array_w.h"



//! Computes the derivative of the fermionic action respect to the link field
/*! \ingroup fermacts
 *
 *  u -- gauge field ( Read )
 *
 *         |  dS      dS_f
 * ds_u -- | ----   + -----   ( Modify )
 *         |  dU       dU
 *
 * psi -- [1./(M_dag*M)]*chi_  ( read ) 
 */

multi1d<LatticeColorMatrix> UnprecDWFermAct::dsdu(const multi1d<LatticeColorMatrix>& u, 
						  const multi1d<LatticeFermion>& psi) const
{
  multi1d<LatticeColorMatrix> ds_u(Nd);

  // hack
  ds_u = 0;

  QDPIO::cerr << "UnprecDWFermAct::dsdu not implemented" << endl;
  QDP_abort(1);

  END_CODE("UnprecDWFermAct::dsdu");

  return ds_u;
}
