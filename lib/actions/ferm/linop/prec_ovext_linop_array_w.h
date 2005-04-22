// -*- C++ -*-
// $Id: prec_ovext_linop_array_w.h,v 1.1 2005-04-22 13:27:42 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) linear operator
 */

#ifndef __prec_ovext_linop_array_w_h__
#define __prec_ovext_linop_array_w_h__

#include "chromabase.h"
#include "handle.h"
#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"

using namespace QDP;

namespace Chroma 
{ 
  //! EvenOddPreconditioned Extended-Overlap (N&N) linear operator
  /*!
   * \ingroup linop
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+m_q)/(1-m_q)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the pole approx. to eps(H(m))
   *
   * This operator implements  hep-lat/0005004
   */

  class EvenOddPrecOvExtLinOpArray : public EvenOddPrecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Partial constructor
    EvenOddPrecOvExtLinOpArray() {}

    //! Full constructor
    EvenOddPrecOvExtLinOpArray(const multi1d<LatticeColorMatrix>& u_,
			  const int Npoles_,
			  const Real& coeffP_,
			  const multi1d<Real>& resP_,
			  const multi1d<Real>& rootQ_,
			  const multi1d<Real>& beta_,
			  const Real& OverMass_,
			  const Real& Mass_,
			  const Real& b5_,
			  const Real& c5_)

    {create(u_,Npoles_, coeffP_, resP_, rootQ_, beta_, 
	    OverMass_,Mass_,b5_,c5_);}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, 
		const int Npoles_,
		const Real& coeffP_,
		const multi1d<Real>& resP_,
		const multi1d<Real>& rootQ_,
		const multi1d<Real>& beta_,
		const Real& OverMass_,
		const Real& m_q_,
		const Real& b5_,
		const Real& c5_);

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Destructor is automatic
    ~EvenOddPrecOvExtLinOpArray() {}

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

    //! Apply the inverse of the odd-odd block onto a source vector
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
    int Npoles;
    int N5;
 
    // Needed for applying Even-Even
    Real Q; // (Nd-m)
    Real A; // -alpha Q
    multi1d<Real> B;  // B_p = sqrt(q_p beta_p) [ 2 + a5 Q ];
    multi1d<Real> D;  // D_p = sqrt(p_p beta_p) [ 2 + a5 Q ];
    multi1d<Real> C;  // C_p = alpha beta_p Q
    Real E; // (2R + (Ra5 + alpha p0)Q)
    Handle< const DslashLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> > > Dslash; //Dslash Op

    // Needed for applying Off Diag
    Real Aprime;  // -alpha/2
    Real Eprime;  // -(1/2)( Ra5 + p0 alpha )
    multi1d<Real> Bprime;  // -(1/2) a5 sqrt( q_p beta_p )
    multi1d<Real> Dprime;  // -(1/2) a5 sqrt( p_p beta_p )
    multi1d<Real> Cprime;  // -alpha beta_p / 2

    // Needed for applying diagInv
    multi1d<Real> Atilde; // (1/det(Block_p) * A )
    multi1d<Real> Btilde; // (1/det(Block_p) * B )
    multi1d<Real> Ctilde; // (1/det(Block_p) * C )
    multi1d<Real> D_bd_inv; // [ 0, D_N,...,0,D_0 ] Block Diag Inv
    Real S;  // Schur's complement: E + sum_p D_p^2 A tilde_p
  };

}; // End Namespace Chroma


#endif
