// $Id: unprec_dwftransf_linop_w.cc,v 1.1 2004-11-02 10:33:50 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"
#include "actions/ferm/linop/unprec_dwftransf_linop_w.h"
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
  // Drop into handle
  QDPIO::cout << "Creating UnprecDWFTransfLinOp with ";
  QDPIO::cout << " b5=" << b5 << " c5=" << c5 << " Mass=" << Mass;
  QDPIO::cout << " RsdCG=" << inv_param.RsdCG << endl;
  D_w = new UnprecWilsonLinOp(u_, Mass_);

  Real b5_minus_c5 = b5 - c5;

  D_denum = new UnprecDWFTransfDenLinOp(u, b5_minus_c5, D_w);
  D_sq_denum = new UnprecDWFTransfMdagMDenLinOp(u, b5_minus_c5, D_w);

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
  D_w = new UnprecWilsonLinOp(u_, Mass_);

  // Create D_denum = ( 2 + (b5-c5) D )
  Real b5_minus_c5 = b5 - c5;
  D_sq_denum = new UnprecDWFTransfMdagMDenLinOp(u, b5_minus_c5, D_w);

}


void UnprecDWFTransfLinOp::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
				    enum PlusMinus isign) const
{
  START_CODE();

  // Apply   (b5+c5) D / ( 2 + (b5-c5) D )

  // First do tmp = (b5 + c5) D psi
  LatticeFermion tmp,tmp2
;
  enum PlusMinus isign_dag = ( isign == PLUS ? MINUS : PLUS );
  (*D_denum)(tmp2, psi, isign_dag);

  int n_count;
  InvCG2<LatticeFermion>(*D_sq_denum, 
			 tmp2, 
			 tmp, 
			 inv_param.RsdCG, 
			 inv_param.MaxCG,
			 n_count);

  (*D_w)(chi, tmp, isign);
  chi *= (b5 + c5);


  QDPIO::cout << "DWFTransf Denominator invert n_count = " << n_count << endl;
  
  END_CODE();
}

void UnprecDWFTransfMdagMLinOp::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
				    enum PlusMinus isign) const
{
  START_CODE();

  enum PlusMinus isign_dag = ( isign == PLUS ? MINUS : PLUS );

  // Apply   (b5+c5)^2 D^{dag}D  / (2 + (b5 - c5) D^{dag} ) ( 2 + (b5-c5) D )

  // First do tmp = (b5 + c5) D psi
  LatticeFermion tmp, tmp2;
  int n_count;
  InvCG2<LatticeFermion>(*D_sq_denum, 
			 psi, 
			 tmp, 
			 inv_param.RsdCG, 
			 inv_param.MaxCG,
			 n_count);

  (*D_w)(tmp2, tmp, isign);
  (*D_w)(chi, tmp2, isign_dag);

  chi *= ((b5 + c5)*(b5 + c5));
  

  QDPIO::cout << "DWFTransfMdagM Denominator invert n_count = " << n_count << endl;
  
  END_CODE();
}
