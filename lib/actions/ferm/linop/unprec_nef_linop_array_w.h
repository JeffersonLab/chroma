// -*- C++ -*-
// $Id: unprec_nef_linop_array_w.h,v 1.2 2004-08-20 20:29:17 kostas Exp $
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

  void CompbineHoppingVectors(multi1d<LatticeFermion>& tmp,
			      const multi1d<LatticeFermion>& chi) const
  {
    Real fb5(-.5*b5);
    Real fc5(-.5*c5);
    for(int s(1);s<N5;s++)
      tmp[s] = fb5*chi[s] + 
	       fc5*(chi[s+1] + chi[s-1] + 
		    GammaConst<Ns,Ns*Ns-1>()*(chi[s-1] - chi[s+1])) ;
    //Now the boundary terms
    tmp[0] = fb5*chi[0] +
             fc5*(chi[1] - m_q*chi[N5-1] +
		  GammaConst<Ns,Ns*Ns-1>()*( m_q*chi[N5-1] -chi[1])) ;
    tmp[N5-1] = fb5*chi[N5-1] + 
                fc5*(chi[N5-2]  - m_q*chi[0] + 
		     GammaConst<Ns,Ns*Ns-1>()*(chi[N5-2] + m_q*chi[0]));
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
