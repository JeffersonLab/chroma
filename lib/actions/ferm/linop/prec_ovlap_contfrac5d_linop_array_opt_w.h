// -*- C++ -*-
// $Id: prec_ovlap_contfrac5d_linop_array_opt_w.h,v 1.1 2005-05-27 18:39:38 edwards Exp $
/*! \file
 *  \brief Optimized Even-odd prec. 5D continued fraction linop
 */

#ifndef __prec_ovlap_contfrac5d_linop_array_opt_w_h__
#define __prec_ovlap_contfrac5d_linop_array_opt_w_h__

#include "linearop.h"
#include "fermact.h"
#include "state.h"


namespace Chroma 
{ 
  //! Optimized Even-odd prec. 5D continued fraction linop
  /*!
   * \ingroup linop
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+m_q)/(1-m_q)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the pole approx. to eps(H(m))
   *
   * The details of this operator are in some lattice proceedings 
   * by Joo,Kennedy,Wenger
   */

  class OptEvenOddPrecOvlapContFrac5DLinOpArray : public EvenOddPrecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >
  {
  public:

    //! Full constructor
    /*! Pretty darn the same as for the unprec case
      except that the auxiliary linop M is no longer supplied, 
      but is created here 
    */
    OptEvenOddPrecOvlapContFrac5DLinOpArray(Handle<const ConnectState> state,
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
    ~OptEvenOddPrecOvlapContFrac5DLinOpArray() {}

    //! Only defined on the entire lattice
    // INHERIT THIS?
    // const OrderedSubset& subset() const {return all;}

    //! Apply the even-even block onto a source vector
    inline
    void evenEvenLinOp(multi1d<LatticeFermion>& chi, 
		       const multi1d<LatticeFermion>& psi, 
		       enum PlusMinus isign) const
    {
      applyDiag(chi, psi, isign, 0);
    }

    //! Apply the the odd-odd block onto a source vector
    inline
    void oddOddLinOp(multi1d<LatticeFermion>& chi, 
		     const multi1d<LatticeFermion>& psi, 
		     enum PlusMinus isign) const
    {
      applyDiag(chi, psi, isign, 1);
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
      applyDiagInv(chi, psi, isign, 0);
    }

    //! Apply the inverse of the odd-odd block onto a source vector
    inline
    void oddOddInvLinOp(multi1d<LatticeFermion>& chi, 
			const multi1d<LatticeFermion>& psi, 
			enum PlusMinus isign) const
    {
      applyDiagInv(chi, psi, isign, 1);
    }
  
    //! Apply the even-even block onto a source vector
    void derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			    const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
			    enum PlusMinus isign) const
    {
      ds_u.resize(Nd);
      ds_u = zero;
    }

    //! Apply the the odd-odd block onto a source vector
    void derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			  const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
			  enum PlusMinus isign) const
    {
      ds_u.resize(Nd);
      ds_u = zero;
    }

    //! Apply the the even-odd block onto a source vector
    void derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			   const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
			   enum PlusMinus isign) const
    {
      applyDerivOffDiag(ds_u, chi, psi, isign, 0);
    }
 
    //! Apply the the odd-even block onto a source vector
    void derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			   const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
			   enum PlusMinus isign) const
    {
      applyDerivOffDiag(ds_u, chi, psi, isign, 1);
    }

    // Override virtual function for efficiency.
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
	       enum PlusMinus isign) const;

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

    //! Apply the even-odd (odd-even) coupling piece of the NEF operator
    /*!
     * \param ds_u    conjugate momenta	               (Read)
     * \param psi     left pseudofermion field	       (Read)
     * \param psi     right pseudofermion field        (Read)
     * \param isign   Flag ( PLUS | MINUS )   	       (Read)
     * \param cb      checkerboard ( 0 | 1 )           (Read)
     */
    void applyDerivOffDiag(multi1d<LatticeColorMatrix>& ds_u, 
			   const multi1d<LatticeFermion>& chi, 
			   const multi1d<LatticeFermion>& psi, 
			   enum PlusMinus isign,
			   int cb) const;

  private:
    Handle< const DslashLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> > > Dslash; //Dslash Op
    const Real m_q;
    const Real OverMass;
    const int  N5;    // Size of the 5th dimension
    const Real scale_fac;
    const multi1d<Real> alpha;
    const multi1d<Real> beta;
    const bool isLastZeroP;
    multi1d<Real> beta_tilde; // The beta_tilde_i
    multi1d<Real> a;      // The a_i
    multi1d<Real> invd;   // The 1/d_i
    multi1d<Real> u;      // The u_i = l_i
    multi1d<Real> off_diag_coeff; // -0.5*beta_tilde[i]
  };

}; // End Namespace Chroma


#endif
