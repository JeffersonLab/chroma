// -*- C++ -*-
// $Id: prec_contfrac_linop_array_w.h,v 1.2 2004-09-07 15:53:42 bjoo Exp $
/*! \file
 *  \brief 4D Even Odd preconditioned domain-wall fermion linear operator
 */

#ifndef __prec_contfrac_linop_array_w_h__
#define __prec_contfrac_linop_array_w_h__

#include "linearop.h"
#include "state.h"

#include "actions/ferm/linop/dslash_w.h"

using namespace QDP;

//! 4D Even Odd preconditioned domain-wall Dirac operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 */

class EvenOddPrecContFracLinOpArray : public EvenOddPrecLinearOperator< multi1d<LatticeFermion> >
{
public:
  //! Partial constructor
  EvenOddPrecContFracLinOpArray() {}

  //! Full constructor
  //
  // state_ the connect state to fix in the LinOp
  //
  // N5odd_ the length of the 5th dimension. Must be odd for now
  //        and perhaps forever
  //
  // m_q_    the quark mass
  //
  // WilsonMass_ -- the mass of the Auxiliary Wilson Op
  //
  // k_          -- the parameterisation of the continued fraction
  //
  // c_          -- the parameterisation of the equivalence transformation
  //
  EvenOddPrecContFracLinOpArray(Handle< ConnectState > state_,
				const int N5_,
				const Real& m_q_,
				const Real& WilsonMass_,
				const multi1d<Real>& k_,
				const multi1d<Real>& c_) {
    create(state_, N5, m_q_, WilsonMass_, k_, c_) ;
  }

  //! Creation routine
  //
  // state_ the connect state to fix in the LinOp
  //
  // N5_ the length of the 5th dimension. Must be odd for now
  //        and perhaps forever
  //
  // m_q_    the quark mass
  //
  // WilsonMass_ -- the mass of the Auxiliary Wilson Op
  //
  // k_          -- the parameterisation of the continued fraction
  //
  // c_          -- the parameterisation of the equivalence transformation
  //

  void create(Handle< ConnectState > state_,
	      const int N5_,
	      const Real& m_q_,
	      const Real& WilsonMass_,
	      const multi1d<Real>& k_,
	      const multi1d<Real>& c_) ;

  //! Destructor is automatic
  ~EvenOddPrecContFracLinOpArray() {}

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
    START_CODE();
    applyDiagInv(chi, psi, isign, 0);
    END_CODE();
  }
  
  //! Apply the the even-odd block onto a source vector
  void evenOddLinOp(multi1d<LatticeFermion>& chi, 
		    const multi1d<LatticeFermion>& psi, 
		    enum PlusMinus isign) const
  {
    START_CODE();
    // Gamma(5) index. What are the GammaConst thingies?
    const int G5=Ns*Ns-1;
    
    LatticeFermion tmp;
    
    // Our mass term is in the bottom component, so we 
    // start there and work up.
    // signH=1 for bottom component and wiggles as we come up.
    for(int signH=1, s=N5-1, index=0; 
	s >=0 ; 
	s--, signH *=-1, index++) {
      
      // tmp = D_{eo} psi_s
      //
      // Note -- possible optimisation -- apply all Dslash-es at once
      D.apply(tmp,psi[s],isign,0);

      // Tmp2 = gamma_5 * tmp
      chi[s][rb[0]] = Gamma(G5)*tmp;
      
      // Now multiply in the right coefficieng
      Real coeff;
      if( index == 0 ) {
	// s = 0 case.
	coeff = -0.5*k[index];
      }
      else {
	coeff = -0.5*c[index-1]*c[index-1]*k[index]*signH;
      }
      
      chi[s][rb[0]] *= coeff;
    }
    END_CODE();
  }

  //! Apply the the odd-even block onto a source vector
  void oddEvenLinOp(multi1d<LatticeFermion>& chi, 
		    const multi1d<LatticeFermion>& psi, 
		    enum PlusMinus isign) const
  {
    START_CODE();
    // Gamma(5) index. What are the GammaConst thingies?
    const int G5=Ns*Ns-1;
    
    LatticeFermion tmp;
    
    // Our mass term is in the bottom component, so we 
    // start there and work up.
    // signH=1 for bottom component and wiggles as we come up.
    for(int signH=1, s=N5-1, index=0; 
	s >=0 ; 
	s++, signH *=-1, index++) {
      
      // tmp = D_{eo} psi_s
      //
      // Note -- possible optimisation -- apply all Dslash-es at once
      D.apply(tmp,psi[s],isign,1);

      // Tmp2 = gamma_5 * tmp
      chi[s][rb[1]] = Gamma(G5)*tmp;
      
      // Now multiply in the right coefficieng
      Real coeff;
      if( index == 0 ) {
	// s = 0 case.
	coeff = -0.5*k[index];
      }
      else {
	coeff = -0.5*c[index-1]*c[index-1]*k[index]*signH;
      }
      
      chi[s][rb[1]] *= coeff;
    }
    END_CODE();
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
  int N5;
  Real m_q;
  Real WilsonMass;
  multi1d<Real> c;
  multi1d<Real> k;

  Real M_factor;              // (Nd + M_wilson)
  multi1d<Real> a;            // Diag Coeffs of tridiag matrix
  multi1d<Real> b;            // Off Diag Coeffs of tridiag matrix

  multi1d<Real> d;            // Coeffs of Diag Matrix
  multi1d<Real> u_l;          // Coeffs of upper/lower decomposition
 
  WilsonDslash D;
};

#endif
