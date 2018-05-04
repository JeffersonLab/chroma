// -*- C++ -*-
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) linear operator
 */

#ifndef __prec_ht_contfrac5d_linop_array_w_h__
#define __prec_ht_contfrac5d_linop_array_w_h__

#include "eoprec_constdet_linop.h"
#include "state.h"
#include "dslash_array_w.h"

namespace Chroma 
{ 
  //! Even odd preconditioned Continued Fraction with H=H_t (Shamir Kernel)
  /*!
   * \ingroup linop
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+m_q)/(1-m_q)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the pole approx. to eps(H(m))
   *
   * and H(m) is H_t = (2 + (b5-c5)H_w gamma_5 )^{-1} [ (b5+c5) H_w ]
   *
   *                 = (2 + (b5-c5) D^{dagger} )^{-1} [ (b5+c5) g5 D_w ]
   * 
   */

  class EvenOddPrecHtContFrac5DLinOpArray : public EvenOddPrecConstDetLinearOperatorArray<LatticeFermion,
					    multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Full constructor
    /*! Pretty darn the same as for the unprec case
      except that the auxiliary linop M is no longer supplied, 
      but is created here 
    */
    EvenOddPrecHtContFrac5DLinOpArray(Handle< FermState<T,P,Q> > fs,
				      const Real& _m_q,
				      const Real& _OverMass,
				      int _N5,
				      const Real& _scale_fac,
				      const multi1d<Real>& _alpha,
				      const multi1d<Real>& _beta,
				      const Real& _b5,
				      const Real& _c5,
				      const bool _isLastZeroP );

  
    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Destructor is automatic
    ~EvenOddPrecHtContFrac5DLinOpArray() {}

    //! Only defined on the entire lattice
    // INHERIT THIS?
    // const Subset& subset() const {return all;}

    //! Apply the even-even block onto a source std::vector
    inline
    void evenEvenLinOp(multi1d<LatticeFermion>& chi, 
		       const multi1d<LatticeFermion>& psi, 
		       enum PlusMinus isign) const
      {
	applyDiag(chi, psi, isign, 0);
      }

    //! Apply the the odd-odd block onto a source std::vector
    inline
    void oddOddLinOp(multi1d<LatticeFermion>& chi, 
		     const multi1d<LatticeFermion>& psi, 
		     enum PlusMinus isign) const
      {
	applyDiag(chi, psi, isign, 1);
      }
  
  
    //! Apply the the even-odd block onto a source std::vector
    void evenOddLinOp(multi1d<LatticeFermion>& chi, 
		      const multi1d<LatticeFermion>& psi, 
		      enum PlusMinus isign) const
      {
	applyOffDiag(chi, psi, isign, 0);
      }

    //! Apply the the odd-even block onto a source std::vector
    void oddEvenLinOp(multi1d<LatticeFermion>& chi, 
		      const multi1d<LatticeFermion>& psi, 
		      enum PlusMinus isign) const
      {
	applyOffDiag(chi, psi, isign, 1);
      }

    //! Apply the inverse of the even-even block onto a source std::vector
    inline
    void evenEvenInvLinOp(multi1d<LatticeFermion>& chi, 
			  const multi1d<LatticeFermion>& psi, 
			  enum PlusMinus isign) const
      {
	applyDiagInv(chi, psi, isign, 0);
      }

    //! Apply the inverse of the odd-odd block onto a source std::vector
    inline
    void oddOddInvLinOp(multi1d<LatticeFermion>& chi, 
			const multi1d<LatticeFermion>& psi, 
			enum PlusMinus isign) const
      {
	applyDiagInv(chi, psi, isign, 1);
      }
  
    //! Apply the even-even block onto a source std::vector
    void derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			    const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
			    enum PlusMinus isign) const
      {
	ds_u.resize(Nd);
	ds_u = zero;
      }

    //! Apply the the odd-odd block onto a source std::vector
    void derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			  const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
			  enum PlusMinus isign) const
      {
	ds_u.resize(Nd);
	ds_u = zero;
      }

    //! Apply the the even-odd block onto a source std::vector
    void derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			   const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
			   enum PlusMinus isign) const
      {
	applyDerivOffDiag(ds_u, chi, psi, isign, 0);
      }
 
    //! Apply the the odd-even block onto a source std::vector
    void derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			   const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
			   enum PlusMinus isign) const
      {
	applyDerivOffDiag(ds_u, chi, psi, isign, 1);
      }

    //! Return flops performed by the evenEvenLinOp
    unsigned long evenEvenNFlops(void) const { 
      return diagNFlops();
    }

    //! Return flops performed by the oddOddLinOp
    unsigned long oddOddNFlops(void) const { 
      return diagNFlops();
    }

    //! Return flops performed by the evenOddLinOp
    unsigned long evenOddNFlops(void) const { 
      return offDiagNFlops();
    }

    //! Return flops performed by the oddEvenLinOp
    unsigned long oddEvenNFlops(void) const { 
      return offDiagNFlops();
    }

    //! Return flops performed by the evenEvenInvLinOp
    unsigned long evenEvenInvNFlops(void) const { 
      return diagInvNFlops();
    }

    //! Return flops performed by the operator()
    unsigned long nFlops() const { 
      // Flopcount is the oddEven EvenEvenInv evenOdd
      // the oddOdd
      // and the subtraction OddOdd - (  )
      unsigned long flops=oddEvenNFlops() + evenEvenInvNFlops() + evenOddNFlops() + oddOddNFlops() + (2*Nc*Ns*N5*(Layout::sitesOnNode()/2));

      return flops;
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

    //! Return flops performed by the diagonal part
    unsigned long diagNFlops(void) const { 
      unsigned long cbsite_flops = (10*N5-8)*Nc*Ns;
      return cbsite_flops*(Layout::sitesOnNode()/2);
    }

    //! Return flops performed by the diagonal part
    unsigned long offDiagNFlops(void) const { 
      unsigned long cbsite_flops = (10*N5-8)*Nc*Ns + N5*1320;
      return cbsite_flops*(Layout::sitesOnNode()/2);
    }

    //! Return flops performed by the diagonal part
    unsigned long diagInvNFlops(void) const { 
      unsigned long cbsite_flops = (10*N5-14)*Nc*Ns;
      return cbsite_flops*(Layout::sitesOnNode()/2);
    }

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return Dslash.getFermBC();}

  protected:
    //! Hide partial constructor
//    EvenOddPrecHtContFrac5DLinOpArray() {}
    //! Hide =
    void operator=(const EvenOddPrecHtContFrac5DLinOpArray&) {}

  private:
    WilsonDslashArray  Dslash;
    const Real m_q;
    const Real OverMass;
    const int  N5;    // Size of the 5th dimension
    const Real scale_fac;
    const multi1d<Real> alpha;
    const multi1d<Real> beta;
    const bool isLastZeroP;
    multi1d<Real> beta_tilde;  // The beta_tilde_i
    multi1d<Real> alpha_tilde; // The alpha_tilde_i

    multi1d<Real> a;  // The a_i
    multi1d<Real> d;  // The d_i
    multi1d<Real> u;  // The u_i = l_i
    multi1d<Real> invd;
    const Real b5;
    const Real c5;
    Real f_plus;
    Real f_minus;

  };

} // End Namespace Chroma


#endif
