// $Id: unprec_dwftransf_linop_w.cc,v 1.2 2004-11-02 12:21:54 bjoo Exp $
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

  // FOr the numerator
  H_w = new lgherm<LatticeFermion>(new UnprecWilsonLinOp(u_,Mass_));

  // For the denominator
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

  // For the numearator
  H_w = new lgherm<LatticeFermion>(new UnprecWilsonLinOp(u_,Mass_));

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

  switch( isign ) {
  case PLUS :
    {
      // DO CGNE: 
      //  i) tmp_2 = D_denum^{dag} psi (purely for CGNE)
      
      (*D_denum)(tmp2, psi, MINUS);
      
      // ii) Solve  tmp = [ D_denum^{dag} D_denum ]^{-1} tmp2
      //               = D_denum^{-1} D_denum^{-dag} D_denum^{dag} psi
      //               = D_denum^{-1} psi
      InvCG2<LatticeFermion>(*D_denum, 
			     tmp2, 
			     tmp, 
			     inv_param.RsdCG, 
			     inv_param.MaxCG,
			     n_count);
      
      // Now multiply in Numerator: 
      (*H_w)(chi, tmp, isign);

    }
    break;
  case MINUS:
    {
      int n_count;
      // THis bit is hermitian already just order changes
      (*H_w)(tmp, psi, PLUS);

      // Now solve D^{dag} D tmp2 = tmp
      // => tmp2 = D^{-1} D^{-dag} tmp
      InvCG2<LatticeFermion>(*D_denum,
			     tmp,
			     tmp2,
			     inv_param.RsdCG,
			     inv_param.MaxCG,
			     n_count);

      // Now multiply by D
      // chi = D * D^{-1} D^{-dag} tmp
      //     = D^{-dag} tmp = D^{-dag} H psi
      (*D_denum)(chi, tmp2, PLUS);
      break;
    }
  default: 
    {
      QDPIO::cerr << "SHould never get here " << endl;
      QDP_abort(1);
      break;
    }
  };

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

  // First do tmp = H_w psi
  LatticeFermion tmp, tmp2;
  (*H_w)(tmp, psi, PLUS);

  int n_count;
  // Now apply 1 /[ ( 2 + (b5-c5)D^{dag} ) ( 2 + (b5-c5)D ) ]
  //
  // solve D_denum^{dag} D_denum = tmp
  InvCG2<LatticeFermion>(*D_denum, 
			 tmp, 
			 tmp2, 
			 inv_param.RsdCG, 
			 inv_param.MaxCG,
			 n_count);


  // Now multiply in the second H_w on top
  (*H_w)(chi, tmp2, PLUS );

  // Now multiply in the b5_c5 factor
  chi *= ((b5 + c5)*(b5 + c5));
  

  QDPIO::cout << "DWFTransfMdagM Denominator invert n_count = " << n_count << endl;
  
  END_CODE();
}
