// -*- C++ -*-
// $Id: prec_nef_linop_array_w.h,v 1.4 2004-09-01 03:32:59 kostas Exp $
/*! \file
 *  \brief 4D Even Odd preconditioned NEF domain-wall fermion linear operator
 */

#ifndef __prec_nef_linop_array_w_h__
#define __prec_nef_linop_array_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"

using namespace QDP;

//! 4D Even Odd preconditioned NEF domain-wall Dirac operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 */

class EvenOddPrecNEFDWLinOpArray : public EvenOddPrecLinearOperator< multi1d<LatticeFermion> >
{
public:
  //! Partial constructor
  EvenOddPrecNEFDWLinOpArray() {}

  /*!
    Full constructor 
    Set b5 = 1.0 and c5=0.0 to get Shamir DWF with a5=1.
    Set b5 = 1.0 and c5=1.0 to get Borichi DWF.
   */
  EvenOddPrecNEFDWLinOpArray(const multi1d<LatticeColorMatrix>& u_, 
			     const Real& WilsonMass_, const Real &b5_, 
			     const Real &c5_, const Real& m_q, int N5_)
  {create(u_,WilsonMass_,b5_,c5_,m_q,N5_);}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& u_, 
	      const Real& WilsonMass_, const Real &b5_, 
	      const Real &c5_, const Real& m_q_, int N5_);


  //! set b5 and c5 given kappa and a5
  //void set_b5_c5(const Real &kappa_, const Real &a5_);


  //! Destructor is automatic
  ~EvenOddPrecNEFDWLinOpArray() {}

  //! Length of DW flavor index/space
  int size() const {return N5;}

  //! Apply the even-even block onto a source vector
  inline
  void evenEvenLinOp(multi1d<LatticeFermion>& chi, 
		     const multi1d<LatticeFermion>& psi, 
		     enum PlusMinus isign) const
    {
      applyDiag(chi, psi, isign, 0);
    }
  
  //! Apply the inverse of the even-even block onto a source vector
  inline
  void evenEvenInvLinOp(multi1d<LatticeFermion>& chi, 
			const multi1d<LatticeFermion>& psi, 
			enum PlusMinus isign) const
    {
      applyDiagInv(chi, psi, isign, 0);
    }
  
  void CompbineHoppingVectors(multi1d<LatticeFermion>& tmp,
			      const multi1d<LatticeFermion>& chi,
			      const int cb) const
  {
    for(int s(1);s<N5-1;s++)
      tmp[s][rb[cb]] = b5*chi[s] + 
	               c5*(chi[s+1] + chi[s-1] + 
	                   GammaConst<Ns,Ns*Ns-1>()*(chi[s-1] - chi[s+1])) ;
    //Now the boundary terms
    tmp[0][rb[cb]] = b5*chi[0] +
	             c5*(chi[1] - m_q*chi[N5-1] +
			 GammaConst<Ns,Ns*Ns-1>()*( m_q*chi[N5-1] -chi[1])) ;
    tmp[N5-1][rb[cb]] = b5*chi[N5-1] + 
	                c5*(chi[N5-2]  - m_q*chi[0] + 
			    GammaConst<Ns,Ns*Ns-1>()*(chi[N5-2] + m_q*chi[0]));
  }

  //! Apply the the even-odd block onto a source vector
  void evenOddLinOp(multi1d<LatticeFermion>& chi, 
		    const multi1d<LatticeFermion>& psi, 
		    enum PlusMinus isign) const
    {
      multi1d<LatticeFermion> tmp(N5) ;
      CompbineHoppingVectors(tmp,psi,0) ;
      for(int s(0);s<N5;s++)
	{
	  D.apply(chi[s],tmp[s],isign,0);
	  chi[s][rb[0]] *= (-0.5);
	}
    }

  //! Apply the the odd-even block onto a source vector
  void oddEvenLinOp(multi1d<LatticeFermion>& chi, 
		    const multi1d<LatticeFermion>& psi, 
		    enum PlusMinus isign) const
    {
      multi1d<LatticeFermion> tmp(N5) ;
      CompbineHoppingVectors(tmp,psi,1) ;
      for(int s(0);s<N5;s++)
      {
	D.apply(chi[s],tmp[s],isign,1);
	chi[s][rb[1]] *= (-0.5);
      }
    }

  //! Apply the the odd-odd block onto a source vector
  inline
  void oddOddLinOp(multi1d<LatticeFermion>& chi, 
		   const multi1d<LatticeFermion>& psi, 
		   enum PlusMinus isign) const
    {
      applyDiag(chi, psi, isign, 1);
    }

  //! Apply the inverse of the odd-odd block onto a source vector
  inline
  void oddOddInvLinOp(multi1d<LatticeFermion>& chi, 
			const multi1d<LatticeFermion>& psi, 
			enum PlusMinus isign) const
    {
      applyDiagInv(chi, psi, isign, 1);
    }

protected:

  //! Apply the even-even (odd-odd) coupling piece of the domain-wall fermion operator
  /*!
   * \param chi     result     	                   (Modify)
   * \param psi     source     	                   (Read)
   * \param isign   Flag ( PLUS | MINUS )          (Read)
   * \param cb      checkerboard ( 0 | 1 )         (Read)
   */
  void applyDiag(multi1d<LatticeFermion>& chi, 
	     const multi1d<LatticeFermion>& psi, 
	     enum PlusMinus isign,
	     const int cb) const;

  //! Apply the inverse even-even (odd-odd) coupling piece of the domain-wall fermion operator
  /*!
   * \param chi     result     	                   (Modify)
   * \param psi     source     	                   (Read)
   * \param isign   Flag ( PLUS | MINUS )   	   (Read)
   * \param cb      checkerboard ( 0 | 1 )         (Read)
   */
  void applyDiagInv(multi1d<LatticeFermion>& chi, 
		    const multi1d<LatticeFermion>& psi, 
		    enum PlusMinus isign,
		    const int cb) const;
    

private:
  Real WilsonMass;
  Real c5;
  Real b5;
  Real m_q;
  int  N5;

  Real InvTwoKappa ;
  Real TwoKappa ;
  Real Kappa;
  Real invDfactor ;

  WilsonDslash  D;
};

#endif
