// $Id: unprec_wilson_linop_w.cc,v 1.9 2004-07-28 02:38:02 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param u_ 	    gauge field     	       (Read)
 * \param Mass_   fermion kappa   	       (Read)
 */
void UnprecWilsonLinOp::create(const multi1d<LatticeColorMatrix>& u_, const Real& Mass_)
{
  Mass = Mass_;
  u = u_;
  D.create(u);

//    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;
}


//! Apply unpreconditioned Wilson fermion linear operator
/*!
 * \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param chi 	  Pseudofermion field     	       (Read)
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void UnprecWilsonLinOp::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
				    enum PlusMinus isign) const
{
  START_CODE();

  //
  //  Chi   =  (Nd+Mass)*Psi  -  (1/2) * D' Psi
  //
  LatticeFermion tmp;
  Real fact1 = Nd + Mass;
  Real fact2 = -0.5;

  // D is a Dslash - must apply to both CB-s
  D(tmp, psi, isign);

  chi = fact1*psi + fact2*tmp;
  
  END_CODE();
}
