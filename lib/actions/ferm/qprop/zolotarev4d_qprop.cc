#include "chromabase.h"
#include "fermact.h"
#include "linearop.h"
#include "actions/ferm/invert/invcg1.h"
#include "actions/ferm/invert/invcg2.h"
#include "actions/ferm/fermacts/zolotarev4d_fermact_bj_w.h"
#include "meas/eig/eig_w.h"

void
Zolotarev4DFermActBj::qprop(LatticeFermion& psi,
			    Handle< const ConnectState > state,
			    const LatticeFermion& chi,
			    enum InvType invType,
			    const Real& RsdCG,
			    int MaxCG,
			    int& ncg_had) const
{

  Handle< const LinearOperator<LatticeFermion> > M(linOp(state));
  int n_count;

  switch( invType ) {
  case CG_INVERTER:
    {
      // We are doing CGNE
      //  M^{dag} M psi = M^{dag} chi
      
      // First make M^{dag} chi
      LatticeFermion tmp;
      (*M)(tmp, chi, MINUS);
      
      // Check whether the source is chiral.
      Chirality ichiral = isChiralVector(chi);
      if( ichiral == CH_NONE ) { 
	
	// Source is not chiral. In this case we should use,
	// InvCG2 with M
	InvCG2(*M, tmp, psi, RsdCG, MaxCG, n_count);
      }
      else {
	
	// Source is chiral. In this case we should use InvCG1
	// with the special MdagM
	Handle< const LinearOperator<LatticeFermion> > MM(lMdagM(state, ichiral));
	InvCG1(*MM, tmp, psi, RsdCG, MaxCG, n_count);
      }
    }
    break;
  default:
    QDP_error_exit("Zolotarev4DFermActBj::qprop Solver Type not implemented\n");
    break;
  };

  if ( n_count == MaxCG ) { 
    QDP_error_exit("Zolotarev4DFermAct::qprop: No convergence in solver: n_count = %d\n", n_count);
  }

  // Update the ncg_had counter
  ncg_had += n_count;
  
  // Normalize and remove contact term 
  Real ftmp = Real(1) / ( Real(1) - m_q );
  
  psi -= chi;
  psi *= ftmp;
}
