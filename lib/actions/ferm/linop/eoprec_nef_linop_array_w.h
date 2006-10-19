// -*- C++ -*-
// $Id: eoprec_nef_linop_array_w.h,v 3.1 2006-10-19 16:01:30 edwards Exp $
/*! \file
 *  \brief 4D Even Odd preconditioned NEF domain-wall fermion linear operator
 */

#ifndef __prec_nef_linop_array_w_h__
#define __prec_nef_linop_array_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_array_w.h"
#include "actions/ferm/linop/eoprec_dwflike_linop_base_array_w.h"


namespace Chroma
{
  //! 4D Even Odd preconditioned NEF domain-wall Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   */
  class EvenOddPrecNEFDWLinOpArray : public EvenOddPrecDWLikeLinOpBaseArray<LatticeFermion, 
				     multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    // ***** HACK *****
    EvenOddPrecNEFDWLinOpArray(Handle< FermState<T,P,Q> > fs,
			       const Real& WilsonMass_, const multi1d<Real>& b5_, 
			       const multi1d<Real>& c5_, const Real& m_q, int N5_)
    {QDP_error_exit("not supported");}
    // ***** HACK *****


    /*!
      Full constructor 
      Set b5 = 1.0 and c5=0.0 to get Shamir DWF with a5=1.
      Set b5 = 1.0 and c5=1.0 to get Borichi DWF.
    */
    EvenOddPrecNEFDWLinOpArray(Handle< FermState<T,P,Q> > fs,
			       const Real& WilsonMass_, const Real &b5_, 
			       const Real &c5_, const Real& m_q, int N5_)
    {create(fs,WilsonMass_,b5_,c5_,m_q,N5_);}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const Real& WilsonMass_, const Real &b5_, 
		const Real &c5_, const Real& m_q_, int N5_);

    //! set b5 and c5 given kappa and a5
    //void set_b5_c5(const Real &kappa_, const Real &a5_);

    //! Destructor is automatic
    ~EvenOddPrecNEFDWLinOpArray() {}

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
 
    //! Apply the Dminus operator on a lattice fermion. See my notes ;-)
    void Dminus(LatticeFermion& chi,
		const LatticeFermion& psi,
		enum PlusMinus isign,
		int s5) const;

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
    

    //! Apply the even-odd (odd-even) coupling piece of the NEF operator
    /*!
     * \ingroup linop
     *
     * The operator acts on the entire lattice
     *
     * \param psi 	  Pseudofermion field     	       (Read)
     * \param isign   Flag ( PLUS | MINUS )   	       (Read)
     * \param cb      checkerboard ( 0 | 1 )               (Read)
     */
    void  applyOffDiag(multi1d<LatticeFermion>& chi, 
		       const multi1d<LatticeFermion>& psi, 
		       enum PlusMinus isign,
		       const int cb) const ;

  protected:
    //! Partial constructor
    EvenOddPrecNEFDWLinOpArray() {}
    //! Partial constructor
    void operator=(const EvenOddPrecNEFDWLinOpArray&) {}

  private:
    Real WilsonMass;
    Real c5;
    Real b5;
    Real m_q;
    int  N5;

    Real c5TwoKappa ;
    Real c5InvTwoKappa ;
    Real b5TwoKappa ;
    Real b5InvTwoKappa ;

    //Real InvTwoKappa ;
    Real TwoKappa ;
    Real Kappa;
    Real invDfactor ;

    WilsonDslashArray  D;
  };

} // End Namespace Chroma



#endif
