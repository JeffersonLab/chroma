// -*- C++ -*-
// $Id: eo4d_dwf_factors_array_w.h,v 1.1 2003-11-22 19:37:23 kostas Exp $
/*! \file
 *  \brief 4D Even Odd preconditioned domain-wall fermion linear operator
 *   factors 
 */

// Check Convensions... Currently I am using Blum et.al.

#ifndef __eo4d_dwf_factors_array_w_h__
#define __eo4d_dwf_factors_array_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"

using namespace QDP;

//! The even-odd coupling piece of the domain wall operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 */

class QDPEvenOdd4dDWDslashArray : public LinearOperator< multi1d<LatticeFermion> >
{
public:
  //! Partial constructor
  QDPEvenOdd4dDWDslashArray() {}

  //! Full constructor
  QDPEvenOdd4dDWDslashArray(const multi1d<LatticeColorMatrix>& u_, const int N5_)
    {create(u_,N5_);}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& u_, const int N5_);

  //! Destructor is automatic
  ~QDPEvenOdd4dDWDslashArray() {}

  //! Only defined on the entire lattice
  const OrderedSubset& subset() const {return all;}

  //! Apply the even-odd coupling piece  of the  domain-wall fermion operator
  /*!
   * \ingroup linop
   *
   * The operator acts on the entire lattice
   *
   * \param psi 	  Pseudofermion field          (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void operator() (multi1d<LatticeFermion>& chi, 
		   const multi1d<LatticeFermion>& psi, 
		   enum PlusMinus isign) const
  {
    for(int s(0);s<N5;s++){
      D(chi[s],psi[s],isign);
      chi[s] *= 0.5 ;
    }
  }

  //! Apply the even-odd coupling piece  of the  domain-wall fermion operator
  /*!
   * \ingroup linop
   *
   * The operator acts on the entire lattice
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	               (Read)
   * \param cb      the reslult leaves on the cb checkerboard  (Read)
   */
  void apply(multi1d<LatticeFermion>& chi, 
	     const multi1d<LatticeFermion>& psi, 
	     enum PlusMinus isign, const int cb) const
  {
    for(int s(0);s<N5;s++){
      D.apply(chi[s],psi[s],isign,cb);
      chi[s] *= 0.5 ;
    }
  }

private:
  Real a5;
  int  N5;
  WilsonDslash  D;
};


//! The even-even (and odd-odd)  coupling piece of the domain wall operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 */

class QDPEvenOdd4dDWDiagArray : public LinearOperator< multi1d<LatticeFermion> >
{
public:
  //! Partial constructor
  QDPEvenOdd4dDWDiagArray() {}

  //! Full constructor
  QDPEvenOdd4dDWDiagArray(const Real& WilsonMass_, const Real& m_q, const int N5_)
  {create(WilsonMass_,m_q,N5_);}

  //! Creation routine
  void create(const Real& WilsonMass_, const Real& m_q_, const int N5_);

  //! Destructor is automatic
  ~QDPEvenOdd4dDWDiagArray() {}

  //! Only defined on the entire lattice
  const OrderedSubset& subset() const {return all;}

  //! Apply the operator onto a source vector
  void operator() (multi1d<LatticeFermion>& chi, 
		   const multi1d<LatticeFermion>& psi, 
		   enum PlusMinus isign) const;
  
  //! Apply the operator onto a source vector on a given checkerboard
  void apply(multi1d<LatticeFermion>& chi, 
	     const multi1d<LatticeFermion>& psi, 
	     enum PlusMinus isign, const int cb) const;

  //! Apply the inverse operator onto a source vector on a given checkerboard
  void applyInv(multi1d<LatticeFermion>& chi, 
		const multi1d<LatticeFermion>& psi, 
		enum PlusMinus isign, const int cb) const;
private:
  Real InvTwoKappa ;
  Real TwoKappa ;
  Real Kappa;
  Real WilsonMass;
  Real m_f ;
  Real invDfactor ;
  Real a5;
  int  N5;
};

#endif
