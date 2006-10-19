// -*- C++ -*-
// $Id: eoprec_ovdwf_linop_array_w.h,v 3.1 2006-10-19 16:01:30 edwards Exp $
/*! \file
 *  \brief 4D Even Odd preconditioned Overlap-DWF (Borici) linear operator
 */

#ifndef __prec_ovdwf_linop_array_w_h__
#define __prec_ovdwf_linop_array_w_h__

#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/eoprec_dwflike_linop_base_array_w.h"


namespace Chroma 
{ 
  //! 4D Even Odd preconditioned Overlap-DWF (Borici) linear operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   */

  class EvenOddPrecOvDWLinOpArray : public EvenOddPrecDWLikeLinOpBaseArray<LatticeFermion, 
				    multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Full constructor
    EvenOddPrecOvDWLinOpArray(Handle< FermState<T,P,Q> > fs,
			      const Real& WilsonMass_, const Real& m_q, int N5_)
    {create(fs,WilsonMass_,m_q,N5_);}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const Real& WilsonMass_, const Real& m_q_, int N5_);

    //! Destructor is automatic
    ~EvenOddPrecOvDWLinOpArray() {}

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}

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

    //! Apply the Dminus operator on a lattice fermion. See my notes ;-)
    void Dminus(LatticeFermion& chi,
		const LatticeFermion& psi,
		enum PlusMinus isign,
		int s5) const;

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

  protected:
    //! Partial constructor
    EvenOddPrecOvDWLinOpArray() {}
    //! Hide =
    void operator=(const EvenOddPrecOvDWLinOpArray&) {}

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

} // End Namespace Chroma


#endif
