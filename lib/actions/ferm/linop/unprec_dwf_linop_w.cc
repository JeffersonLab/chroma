// $Id: unprec_dwf_linop_w.cc,v 1.1 2003-10-20 20:31:50 edwards Exp $
/*! \file
 *  \brief Unpreconditioned domain-wall linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_dwf_linop_w.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param u_            gauge field   (Read)
 * \param WilsonMass_   DWF height    (Read)
 * \param m_q_          quark mass    (Read)
 */
void UnprecDWLinOp::create(const multi1d<LatticeColorMatrix>& u_, const Real& WilsonMass_, const Real& m_q_)
{
  u = u_;
  WilsonMass = WilsonMass_;
  m_q = m_q_;

  D.create(u,WilsonMass);

//    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;
}


//! Apply unpreconditioned domain-wall fermion linear operator
/*!
 * \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
LatticeDWFermion UnprecDWLinOp::operator() (const LatticeDWFermion& psi, enum LinOpSign isign) const
{
  LatticeDWFermion chi;

  START_CODE("UnprecDWLinOp");

  //
  //  Chi   =  Psi  -  k   D' Psi
  //
  chi = m_q*psi - D(psi, isign);
  
  END_CODE("UnprecDWLinOp");

  return chi;
}
