// $Id: unprec_parwilson_linop_w.cc,v 1.2 2004-01-12 04:54:04 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator with parity breaking term
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_parwilson_linop_w.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param u_ 	   gauge field     	       (Read)
 * \param Mass_    fermion kappa   	       (Read)
 * \param H__      parity breaking term	       (Read)
 */
void UnprecParWilsonLinOp::create(const multi1d<LatticeColorMatrix>& u_, 
				  const Real& Mass_, const Real& H_)
{
  Mass = Mass_;
  H = H_;
  u = u_;
  D.create(u);

//    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;
}


//! Apply unpreconditioned Wilson fermion linear operator with parity breaking term
/*!
 * \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param chi 	  Pseudofermion field     	       (Read)
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void UnprecParWilsonLinOp::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
				       enum PlusMinus isign) const
{
  START_CODE("UnprecParWilsonLinOp");

  //
  //  Chi   =  (Nd+Mass)*Psi  -  (1/2) * D' Psi
  //
  LatticeFermion tmp;
  Real fact1 = Nd + Mass;
  Real fact2 = -0.5;

  // D is a Dslash - must apply to both CB-s
  D(tmp, psi, isign);

  chi = fact1*psi + Gamma(Ns*Ns-1)*(H*timesI(psi)) + fact2*tmp;
  
  END_CODE("UnprecParWilsonLinOp");
}
