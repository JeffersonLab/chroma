// -*- C++ -*-
// $Id: prec_ovdwf_linop_array_w.h,v 1.5 2004-09-03 14:24:35 kostas Exp $
/*! \file
 *  \brief 4D Even Odd preconditioned Overlap-DWF (Borici) linear operator
 */

#ifndef __prec_ovdwf_linop_array_w_h__
#define __prec_ovdwf_linop_array_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"

using namespace QDP;

//! 4D Even Odd preconditioned Overlap-DWF (Borici) linear operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 */

class EvenOddPrecOvDWLinOpArray : public EvenOddPrecLinearOperator< multi1d<LatticeFermion> >
{
public:
  //! Partial constructor
  EvenOddPrecOvDWLinOpArray() {}

  //! Full constructor
  EvenOddPrecOvDWLinOpArray(const multi1d<LatticeColorMatrix>& u_, 
			    const Real& WilsonMass_, const Real& m_q, int N5_)
    {create(u_,WilsonMass_,m_q,N5_);}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& u_, 
	      const Real& WilsonMass_, const Real& m_q_, int N5_);

  //! Destructor is automatic
  ~EvenOddPrecOvDWLinOpArray() {}

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
  
  //! Apply the the even-odd block onto a source vector
  void evenOddLinOp(multi1d<LatticeFermion>& chi, 
		    const multi1d<LatticeFermion>& psi, 
		    enum PlusMinus isign) const
    {
      applyOffDiag(chi, psi, isign, 0);
    }

  //! Apply the the odd-even block onto a source vector
  void oddEvenLinOp(multi1d<LatticeFermion>& chi, 
		    const multi1d<LatticeFermion>& psi, 
		    enum PlusMinus isign) const
    {
      applyOffDiag(chi, psi, isign, 1);
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
    chi = c5InvTwoKappa*psi - tt ;
  }
  
protected:

  //! Apply the even-even (odd-odd) coupling piece of the Borici fermion operator
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

  //! Apply the inverse even-even (odd-odd) coupling piece of the Borici fermion operator
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
    

  //! Apply the even-odd (odd-even) coupling piece of the Borici operator
  /*!
   * \param chi     result     	                   (Modify)
   * \param psi     source     	                   (Read)
   * \param isign   Flag ( PLUS | MINUS )          (Read)
   * \param cb      checkerboard ( 0 | 1 )         (Read)
   */
  void applyOffDiag(multi1d<LatticeFermion>& chi, 
		    const multi1d<LatticeFermion>& psi, 
		    enum PlusMinus isign,
		    const int cb) const;

private:
  Real WilsonMass;
  Real m_q;
  Real a5;
  int  N5;

  Real c5TwoKappa ;
  Real c5InvTwoKappa ;
  Real b5TwoKappa ;
  Real b5InvTwoKappa ;

  Real TwoKappa ;
  Real Kappa;
  Real invDfactor ;

  WilsonDslash  D;
};

#endif
