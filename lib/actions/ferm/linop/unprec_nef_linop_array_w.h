// -*- C++ -*-
// $Id: unprec_nef_linop_array_w.h,v 1.4 2004-09-02 22:39:25 kostas Exp $
/*! \file
 *  \brief Unpreconditioned NEF domain-wall fermion linear operator
 */

#ifndef __unprec_nef_linop_array_w_h__
#define __unprec_nef_linop_array_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"

using namespace QDP;

//! Unpreconditioned domain-wall Dirac operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 */

class UnprecNEFDWLinOpArray : public LinearOperator< multi1d<LatticeFermion> >
{
public:
  //! Partial constructor
  UnprecNEFDWLinOpArray() {}

  /*!
    Full constructor
    Set b5 = 1.0 and c5=0.0 to get Shamir DWF with a5=1.
    Set b5 = 1.0 and c5=1.0 to get Borichi DWF.
   */
  UnprecNEFDWLinOpArray(const multi1d<LatticeColorMatrix>& u_, 
		     const Real& WilsonMass_, const Real& b5_,
		     const Real& c5_, const Real& m_q, int N5_)
    {create(u_,WilsonMass_, b5_,c5_,m_q,N5_);}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& u_, 
	      const Real& WilsonMass_, const Real& b5_, const Real& c5_,
	      const Real& m_q_, int N5_);

  //! Length of DW flavor index/space
  inline int size() const {return N5;}

  //! Destructor is automatic
  ~UnprecNEFDWLinOpArray() {}

  //! Only defined on the entire lattice
  const OrderedSubset& subset() const {return all;}

  //! Apply the operator onto a source vector
  void operator() (multi1d<LatticeFermion>& chi, 
		   const multi1d<LatticeFermion>& psi, 
		   enum PlusMinus isign) const;

  //! Apply the Dminus operator on a vector in Ls. See my notes ;-)
  inline
  void Dminus(multi1d<LatticeFermion>& chi,
	      const multi1d<LatticeFermion>& psi,
	      enum PlusMinus isign)
  {
    multi1d<LatticeFermion> tt(N5) ;
    for(int s(0);s<N5;s++){
      D.apply(tt[s],psi[s],isign,0);
      D.apply(tt[s],psi[s],isign,1);
      chi[s] = psi[s] - c5*tt[s] ;
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
    chi = psi - c5*tt ;
  }


private:
  Real WilsonMass;
  Real b5;
  Real c5;
  Real m_q;
  int  N5;
  WilsonDslash  D;
};

#endif
