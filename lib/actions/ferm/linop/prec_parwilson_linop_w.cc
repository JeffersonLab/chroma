// $Id: prec_parwilson_linop_w.cc,v 1.2 2004-07-28 02:38:02 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion linear operator with parity breaking term
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_parwilson_linop_w.h"

//! Creation routine
/*!
 * \param u_ 	   gauge field     	       (Read)
 * \param Mass_    fermion mass   	       (Read)
 * \param H_       parity breaking term	       (Read)
 */
void EvenOddPrecParWilsonLinOp::create(const multi1d<LatticeColorMatrix>& u_, 
				       const Real& Mass_, const Real& H_)
{
  Mass = Mass_;
  H = H_;
  u = u_;
  D.create(u);

  fact = Nd + Mass;
  Real tmp = 1.0 / (fact*fact + H*H);
  invfact1 = fact * tmp;
  invfact2 = H * tmp;
}


//! Apply even-odd linop component
/*!
 * The operator acts on the entire even sublattice
 *
 * \param chi 	  Pseudofermion field     	       (Write)
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void 
EvenOddPrecParWilsonLinOp::evenOddLinOp(LatticeFermion& chi, 
					const LatticeFermion& psi, 
					enum PlusMinus isign) const
{
  START_CODE();

  Real mhalf = -0.5;

  D.apply(chi, psi, isign, 0);
  chi[rb[0]] *= mhalf;
  
  END_CODE();
}

//! Apply odd-even linop component
/*!
 * The operator acts on the entire odd sublattice
 *
 * \param chi 	  Pseudofermion field     	       (Write)
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void 
EvenOddPrecParWilsonLinOp::oddEvenLinOp(LatticeFermion& chi, 
					const LatticeFermion& psi, 
					enum PlusMinus isign) const
{
  START_CODE();

  Real mhalf = -0.5;

  D.apply(chi, psi, isign, 1);
  chi[rb[1]] *= mhalf;
  
  END_CODE();
}

