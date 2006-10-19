// -*- C++ -*-
// $Id: eoprec_dwf_linop_array_w.h,v 3.1 2006-10-19 16:01:30 edwards Exp $
/*! \file
 *  \brief 4D Even Odd preconditioned domain-wall fermion linear operator
 */

#ifndef __prec_dwf_linop_array_w_h__
#define __prec_dwf_linop_array_w_h__

#include "eoprec_linop.h"
#include "actions/ferm/linop/dslash_array_w.h"
#include "actions/ferm/linop/eoprec_dwflike_linop_base_array_w.h"
#include "io/aniso_io.h"

namespace Chroma 
{ 
  //! 4D Even Odd preconditioned domain-wall Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   */
  class EvenOddPrecDWLinOpArray : public EvenOddPrecDWLikeLinOpBaseArray<LatticeFermion, 
				  multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Full constructor
    EvenOddPrecDWLinOpArray(Handle< FermState<T,P,Q> > fs,
			    const Real& WilsonMass_, const Real& m_q, int N5_,
			    const AnisoParam_t& aniso_);

    //! Destructor is automatic
    ~EvenOddPrecDWLinOpArray() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}

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


    //! Apply the Dminus operator on a lattice fermion.
    void Dminus(LatticeFermion& chi,
		const LatticeFermion& psi,
		enum PlusMinus isign,
		int s5) const;


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
    //! Partial constructor
    EvenOddPrecDWLinOpArray() {}


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
      unsigned long cbsite_flops = (4*N5+2)*Nc*Ns;
      return cbsite_flops*(Layout::sitesOnNode()/2);
    }

    //! Return flops performed by the diagonal part
    unsigned long offDiagNFlops(void) const { 
      unsigned long cbsite_flops = N5*(1320+2*Nc*Ns);
      return cbsite_flops*(Layout::sitesOnNode()/2);
    }

    //! Return flops performed by the diagonal part
    unsigned long diagInvNFlops(void) const { 
      unsigned long cbsite_flops = (10*N5-8)*Nc*Ns;
      return cbsite_flops*(Layout::sitesOnNode()/2);
    }

  private:
    Real WilsonMass;
    Real m_q;
    Real a5;
    int  N5;

    Real InvTwoKappa ;
    Real TwoKappa ;
    Real Kappa;
    Real invDfactor ;

    WilsonDslashArray  D;
  };


} // End Namespace Chroma


#endif
