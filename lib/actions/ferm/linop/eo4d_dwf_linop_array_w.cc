// $Id: eo4d_dwf_linop_array_w.cc,v 1.1 2003-11-22 19:37:23 kostas Exp $
/*! \file
 *  \brief  4D Even Odd preconditioned domain-wall linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/eo4d_dwf_linop_array_w.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param u_            gauge field   (Read)
 * \param WilsonMass_   DWF height    (Read)
 * \param m_q_          quark mass    (Read)
 */
void EvenOdd4dDWLinOpArray::create(const multi1d<LatticeColorMatrix>& u_, 
				const Real& WilsonMass_, const Real& m_q_, int N5_)
{
  WilsonMass = WilsonMass_;
  m_q = m_q_;
  a5  = 1.0;
  N5  = N5_;

  Qslash.create(u_,N5);
  Qdiag.create(WilsonMass,m_q,N5);
}



//-----------------------------------------------------------------------------



//! Apply the 4D even odd preconditioned domain-wall fermion linear operator
/*!
 * \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void EvenOdd4dDWLinOpArray::operator()(multi1d<LatticeFermion>& chi, 
				       const multi1d<LatticeFermion>& psi, 
				       enum PlusMinus isign) const
{
//  multi1d<LatticeFermion> chi(N5);   // probably should check psi.size() == N5

  START_CODE("EvenOdd4dDWLinOpArray");

  //
  //  Chi   =  D' Psi
  //
  multi1d<LatticeFermion>  tmp;

  Qslash.apply  (chi,psi,isign,0);
  Qdiag.applyInv(tmp,chi,isign,0);
  Qslash.apply  (chi,tmp,isign,1);
  Qdiag.applyInv(tmp,chi,isign,1);
  for(int s(0);s<N5;s++)
    chi[s][odd] = psi[s] - tmp[s] ;
 
  END_CODE("EvenOdd4dDWLinOpArray");
}
