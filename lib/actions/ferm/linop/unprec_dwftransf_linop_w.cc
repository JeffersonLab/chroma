// $Id: unprec_dwftransf_linop_w.cc,v 3.2 2007-08-27 18:24:05 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"
#include "actions/ferm/linop/unprec_dwftransf_linop_w.h"
#include "actions/ferm/invert/invcg2.h"

using namespace QDP::Hints;

namespace Chroma 
{ 
  void UnprecDWFTransfLinOp::create(Handle< FermState<T,P,Q> > fs,
				    const Real& Mass_,
				    const Real& b5_,
				    const Real& c5_,
				    const SysSolverCGParams& invParam_)
  {
    Mass = Mass_;
    b5 = b5_;
    c5 = c5_;
    invParam = invParam_;

    // Need to create a handle for a wilson linop
    QDPIO::cout << "Creating UnprecDWFTransfLinOp with ";
    QDPIO::cout << " b5=" << b5 << " c5=" << c5 << " Mass=" << Mass;
    QDPIO::cout << " RsdCG=" << invParam.RsdCG << endl;


    Real b5_minus_c5 = b5 - c5;
    D_w = new UnprecWilsonLinOp(fs, Mass_);
    D_denum = new UnprecDWFTransfDenLinOp(b5_minus_c5, D_w);
    fbc = fs->getFermBC();
  }



  void UnprecDWFTransfLinOp::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
					 enum PlusMinus isign) const
  {
    START_CODE();

    SystemSolverResults_t res;

    // OK Copy from SZIN Code. DOn't DO hermitian op... just do the D's
    switch(isign) 
    {
    case PLUS:
    {
      // Apply chi = (b5+c5)* gamma_5 * D_w * [ D_denum ]^{-1} * psi
      (*D_denum)(chi, psi, MINUS);
      LatticeFermion tmp;  moveToFastMemoryHint(tmp);

      tmp = psi; 
      
      res = InvCG2(*D_denum, 
		   chi, 
		   tmp,
		   invParam.RsdCG,
		   invParam.MaxCG);
      (*D_w)(chi, tmp, PLUS);
   
      break;
    }
    case MINUS:
    {
      // Apply  chi = (b5+c5) * gamma_5 * D_w * [D_denum]^(-1) * gamma_5 * psi
      LatticeFermion tmp;      moveToFastMemoryHint(tmp);
      chi = GammaConst<Ns,Ns*Ns-1>()*psi;      
      (*D_denum)(tmp, chi, MINUS);

      res = InvCG2(*D_denum, 
		   tmp, 
		   chi,
		   invParam.RsdCG,
		   invParam.MaxCG);
      
      (*D_w)(tmp, chi, PLUS);
      chi  = GammaConst<Ns,Ns*Ns-1>()*tmp;
      break;
    }
    default:
    {
      QDPIO::cerr << "Bad option " << isign << endl;
      QDP_abort(1);
      break;
    }
    }
    
    chi *= (b5 + c5);
    QDPIO::cout << "UnprecDWFTransfLinOp: ncount= " << res.n_count << endl;
    END_CODE();
  }



  //------------------------------------------------------------------------

  void UnprecDWFTransfMdagMLinOp::create(Handle< FermState<T,P,Q> > fs,
					 const Real& Mass_,
					 const Real& b5_,
					 const Real& c5_,
					 const SysSolverCGParams& invParam_)
  {
    Mass = Mass_;
    b5 = b5_;
    c5 = c5_;
    invParam = invParam_;

    // Need to create a handle for a wilson linop
    QDPIO::cout << "Creating UnprecDWFTransfMdagMLinOp with ";
    QDPIO::cout << " b5=" << b5 << " c5=" << c5 << " Mass=" << Mass;
    QDPIO::cout << " RsdCG=" << invParam.RsdCG << endl;

    Real b5_minus_c5 = b5 - c5;
    D_w = new UnprecWilsonLinOp(fs, Mass_);
    D_denum = new UnprecDWFTransfDenLinOp(b5_minus_c5, D_w);
    fbc = fs->getFermBC();
  }



  void UnprecDWFTransfMdagMLinOp::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
					      enum PlusMinus isign) const
  {
    START_CODE();

    int n_count;

    // Apply chi = (b5+c5)^2 * gamma_5 * D_w * [ D_denum^dag * D_denum ]^{-1} * gamma_5 * D_w * psi
    LatticeFermion tmp;       moveToFastMemoryHint(tmp);
    (*D_w)(chi, psi, PLUS);
    tmp = GammaConst<Ns,Ns*Ns-1>()*chi;

    chi = tmp;
    InvCG2(*D_denum, 
	   tmp,
	   chi, 
	   invParam.RsdCG,
	   invParam.MaxCG);

    (*D_w)(tmp, chi, PLUS);
    chi = GammaConst<Ns,Ns*Ns-1>()*tmp;
    chi *= Real((b5 + c5)*(b5 + c5));

    QDPIO::cout << "UnprecDWFTransfMdagMLinOp: ncount= " << n_count << endl;
    END_CODE();
  }

} // End Namespace Chroma

