// -*- C++ -*-
// $Id: unprec_ovdwf_linop_array_w.h,v 1.5 2004-09-03 14:24:36 kostas Exp $
/*! \file
 *  \brief Unpreconditioned Overlap-DWF (Borici) linear operator
 */

#ifndef __unprec_ovdwf_linop_array_w_h__
#define __unprec_ovdwf_linop_array_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"

using namespace QDP;

//! Unpreconditioned Overlap-DWF (Borici) linear operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 */

class UnprecOvDWLinOpArray : public LinearOperator< multi1d<LatticeFermion> >
{
public:
  //! Partial constructor
  UnprecOvDWLinOpArray() {}

  //! Full constructor
  UnprecOvDWLinOpArray(const multi1d<LatticeColorMatrix>& u_, const Real& WilsonMass_, const Real& m_q, int N5_)
    {create(u_,WilsonMass_,m_q,N5_);}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& u_, const Real& WilsonMass_, const Real& m_q_, int N5_);

  //! Length of DW flavor index/space
  int size() const {return N5;}

  //! Destructor is automatic
  ~UnprecOvDWLinOpArray() {}

  //! Only defined on the entire lattice
  const OrderedSubset& subset() const {return all;}

  //! Apply the operator onto a source vector
  void operator() (multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
		   enum PlusMinus isign) const;

  //! Apply the Dminus operator on a vector in Ls. See my notes ;-)
  inline
  void Dminus(multi1d<LatticeFermion>& chi,
	      const multi1d<LatticeFermion>& psi,
	      enum PlusMinus isign)
  {
    Real c5InvTwoKappa =  1.0 - (Nd-WilsonMass) ;
    multi1d<LatticeFermion> tt(N5) ;
    for(int s(0);s<N5;s++){
      D.apply(tt[s],psi[s],isign,0);
      D.apply(tt[s],psi[s],isign,1);
      chi[s] = c5InvTwoKappa*psi[s] - tt[s] ;
    }
  }
  
  //! Apply the Dminus operator on a lattice fermion. See my notes ;-)
  inline
  void Dminus(LatticeFermion chi,
	      const LatticeFermion& psi,
	      enum PlusMinus isign)
  {
    LatticeFermion tt ;
    D.apply(tt,psi,isign,0);
    D.apply(tt,psi,isign,1);
    chi = (1.0 - (Nd-WilsonMass))*psi - tt ;
  }
  

private:
  Real WilsonMass;
  Real m_q;
  Real a5;
  int  N5;
  WilsonDslash  D;
};

#endif
