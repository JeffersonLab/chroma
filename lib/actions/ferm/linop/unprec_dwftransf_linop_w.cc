// $Id: unprec_dwftransf_linop_w.cc,v 1.3 2004-11-13 17:32:11 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"
#include "actions/ferm/linop/unprec_dwftransf_linop_w.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "invtype.h"
#include "actions/ferm/invert/invcg2.h"

using namespace QDP;
using namespace Chroma;

void UnprecDWFTransfLinOp::create(const multi1d<LatticeColorMatrix>& u_, 
				  const Real& Mass_,
				  const Real& b5_,
				  const Real& c5_,
				  const InvertParam_t& invParam_)
{
  Mass = Mass_;
  b5 = b5_;
  c5 = c5_;
  inv_param = invParam_;

  // Need to create a handle for a wilson linop
  QDPIO::cout << "Creating UnprecDWFTransfLinOp with ";
  QDPIO::cout << " b5=" << b5 << " c5=" << c5 << " Mass=" << Mass;
  QDPIO::cout << " RsdCG=" << inv_param.RsdCG << endl;


  Real b5_minus_c5 = b5 - c5;
  D_w = new UnprecWilsonLinOp(u_, Mass_);
  D_denum = new UnprecDWFTransfDenLinOp(u, b5_minus_c5, D_w);
}

void UnprecDWFTransfMdagMLinOp::create(const multi1d<LatticeColorMatrix>& u_, 
				  const Real& Mass_,
				  const Real& b5_,
				  const Real& c5_,
				  const InvertParam_t& invParam_) 
{
  Mass = Mass_;
  b5 = b5_;
  c5 = c5_;
  inv_param = invParam_;

  // Need to create a handle for a wilson linop
  // Drop into handle
  QDPIO::cout << "Creating UnprecDWFTransfMdagMLinOp with ";
  QDPIO::cout << " b5=" << b5 << " c5=" << c5 << " Mass=" << Mass;
  QDPIO::cout << " RsdCG=" << inv_param.RsdCG << endl;

  // For the denominator
  Real b5_minus_c5 = b5 - c5;
  D_w = new UnprecWilsonLinOp(u_, Mass_);
  D_denum = new UnprecDWFTransfDenLinOp(u, b5_minus_c5, D_w);

}


void UnprecDWFTransfLinOp::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
				    enum PlusMinus isign) const
{
  START_CODE();

  // Apply   (b5+c5) D / ( 2 + (b5-c5) D )

  // First do tmp = (b5 + c5) D psi
  LatticeFermion tmp,tmp2;
  int n_count;

  // tmp2 = H_w psi
  (*D_w)(tmp, psi, PLUS);
  tmp2=GammaConst<Ns,Ns*Ns-1>()*tmp;
    
  (*D_denum)(tmp,tmp2 , MINUS);
      
  // ii) Solve  tmp = [ D_denum^{dag} D_denum ]^{-1} tmp2
  //               = D_denum^{-1} D_denum^{-dag} D_denum^{dag} psi
  //               = D_denum^{-1} psi
  InvCG2<LatticeFermion>(*D_denum, 
			 tmp, 
			 chi, 
			 inv_param.RsdCG, 
			 inv_param.MaxCG,
			 n_count);
  chi *= (b5 + c5);
  QDPIO::cout << "DWFTransf Denominator invert n_count = " << n_count << endl;
  
  END_CODE();
}

void UnprecDWFTransfMdagMLinOp::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
				    enum PlusMinus isign) const
{
  START_CODE();
  int n_count;

  // First do tmp = H_w psi
  LatticeFermion tmp, tmp2;
  (*D_w)(tmp, psi, PLUS);
  tmp2 = GammaConst<Ns, Ns*Ns-1>()*tmp;


  // Now apply 1 /[ ( 2 + (b5-c5)D^{dag} ) ( 2 + (b5-c5)D ) ]
  //
  // solve D_denum^{dag} D_denum = tmp
  InvCG2<LatticeFermion>(*D_denum, 
			 tmp2, 
			 tmp, 
			 inv_param.RsdCG, 
			 inv_param.MaxCG,
			 n_count);


  // Now multiply in the second H_w on top
  (*D_w)(tmp2, tmp, PLUS );
  chi = GammaConst<Ns,Ns*Ns-1>()*tmp2;

  // Now multiply in the b5_c5 factor
  chi *= ((b5 + c5)*(b5 + c5));
  

  QDPIO::cout << "DWFTransfMdagM Denominator invert n_count = " << n_count << endl;
  
  END_CODE();
}
