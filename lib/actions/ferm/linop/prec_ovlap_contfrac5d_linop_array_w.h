// -*- C++ -*-
// $Id: prec_ovlap_contfrac5d_linop_array_w.h,v 1.1 2004-09-30 14:52:47 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) linear operator
 */

#ifndef __prec_ovlap_contfrac5d_linop_array_w_h__
#define __prec_ovlap_contfrac5d_linop_array_w_h__

#include "linearop.h"
#include "fermact.h"
#include "state.h"

using namespace QDP;

//! Unpreconditioned Extended-Overlap (N&N) linear operator
/*!
 * \ingroup linop
 *
 * This operator applies the extended version of the hermitian overlap operator
 *   Chi  =   ((1+m_q)/(1-m_q)*gamma_5 + B) . Psi
 *  where  B  is the continued fraction of the pole approx. to eps(H(m))
 *
 * This operator implements  hep-lat/0005004
 */

class EvenOddPrecOvlapContFrac5DLinOpArray : public EvenOddPrecLinearOperator< multi1d<LatticeFermion> >
{
public:

  //! Full constructor
  /*! Pretty darn the same as for the unprec case
    except that the auxiliary linop M is no longer supplied, 
    but is created here 
  */
  EvenOddPrecOvlapContFrac5DLinOpArray(Handle<const ConnectState> state,
				       const Real& _m_q,
				       const Real& _OverMass,
				       int _N5,
				       const Real& _scale_fac,
				       const multi1d<Real>& _alpha,
				       const multi1d<Real>& _beta,
				       const bool _isLastZeroP );

  
  //! Length of DW flavor index/space
  int size() const {return N5;}

  //! Destructor is automatic
  ~EvenOddPrecOvlapContFrac5DLinOpArray() {}

  //! Only defined on the entire lattice
  // INHERIT THIS?
  // const OrderedSubset& subset() const {return all;}

  //! Apply the even-even block onto a source vector
  inline
  void evenEvenLinOp(multi1d<LatticeFermion>& chi, 
		     const multi1d<LatticeFermion>& psi, 
		     enum PlusMinus isign) const
  {
    START_CODE();
    applyDiag(chi, psi, isign, 0);
    END_CODE();
  }

  //! Apply the the odd-odd block onto a source vector
  inline
  void oddOddLinOp(multi1d<LatticeFermion>& chi, 
		   const multi1d<LatticeFermion>& psi, 
		   enum PlusMinus isign) const
    {
      START_CODE();
      applyDiag(chi, psi, isign, 1);
      END_CODE();
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

  //! Apply the inverse of the odd-odd block onto a source vector
  inline
  void oddOddInvLinOp(multi1d<LatticeFermion>& chi, 
			const multi1d<LatticeFermion>& psi, 
			enum PlusMinus isign) const
  {
    START_CODE();
    applyDiagInv(chi, psi, isign, 1);
    END_CODE();
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
    
  //! Apply the off diagonal block
  /*!
   * \param chi     result     	                   (Modify)
   * \param psi     source     	                   (Read)
   * \param isign   Flag ( PLUS | MINUS )   	   (Read)
   * \param cb      checkerboard ( 0 | 1 )         (Read)
   */
  void applyOffDiag(multi1d<LatticeFermion>& chi, 
		    const multi1d<LatticeFermion>& psi,
		    enum PlusMinus isign,
		    const int cb) const;

private:
  Handle< const DslashLinearOperator<LatticeFermion> > Dslash; //Dslash Op
  const Real m_q;
  const Real OverMass;
  const int  N5;    // Size of the 5th dimension
  const Real scale_fac;
  const multi1d<Real> alpha;
  const multi1d<Real> beta;
  const bool isLastZeroP;
  multi1d<Real> beta_tilde; // The beta_tilde_i
  multi1d<Real> a;  // The a_i
  multi1d<Real> d;  // The d_i
  multi1d<Real> u;  // The u_i = l_i

};

#endif
