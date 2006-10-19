
// -*- C++ -*-
// $Id: eoprec_nef_general_linop_array_w.h,v 3.1 2006-10-19 16:01:30 edwards Exp $
/*! \file
 *  \brief 4D Even Odd preconditioned NEF domain-wall fermion linear operator
 *            generalised to take array of b_5 and c_5
 */

#ifndef __prec_nef_general_linop_array_w_h__
#define __prec_nef_general_linop_array_w_h__

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
  class EvenOddPrecGenNEFDWLinOpArray : public EvenOddPrecDWLikeLinOpBaseArray<LatticeFermion, 
					multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    // ***** HACK *****
    EvenOddPrecGenNEFDWLinOpArray(Handle< FermState<T,P,Q> > fs,
				  const Real& WilsonMass_, 
				  const multi1d<Real>& b5_, 
				  const multi1d<Real>& c5_, 
				  const Real& m_q_, 
				  int N5_) 
    {create(fs, WilsonMass_, m_q_, b5_, c5_, N5_);}
    // ***** HACK *****



    /*!
      Full constructor 
      Set b5 = 1.0 and c5=0.0 to get Shamir DWF with a5=1.
      Set b5 = 1.0 and c5=1.0 to get Borichi DWF.
    */
    EvenOddPrecGenNEFDWLinOpArray(Handle< FermState<T,P,Q> > fs,
				  const Real& WilsonMass_, 
				  const Real &b5_, 
				  const Real &c5_, 
				  const Real& m_q_, 
				  int N5_)
    {
      multi1d<Real> b5_arr(N5_);
      multi1d<Real> c5_arr(N5_);
      for(int i=0; i < N5_; i++) { 
	b5_arr[i] = b5_;
	c5_arr[i] = c5_;
      }

      create(fs,WilsonMass_,m_q_, b5_arr, c5_arr, N5_);
    }

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const Real& WilsonMass_,
		const Real& m_q_, 
		const multi1d<Real>& b5_, 
		const multi1d<Real>& c5_, 
		int N5_);


    //! set b5 and c5 given kappa and a5
    //void set_b5_c5(const Real &kappa_, const Real &a5_);

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}

    //! Destructor is automatic
    ~EvenOddPrecGenNEFDWLinOpArray() {}

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


  protected:

    //! Apply the even-even (odd-odd) coupling piece of the domain-wall fermion operator
    /*!
     * \param chi     result                 (Modify)
     * \param psi     source                 (Read)
     * \param isign   Flag ( PLUS | MINUS )  (Read)
     * \param cb      checkerboard ( 0 | 1 ) (Read)
     */
    void applyDiag(multi1d<LatticeFermion>& chi, 
		   const multi1d<LatticeFermion>& psi, 
		   enum PlusMinus isign,
		   int cb) const;

    //! Apply the inverse even-even (odd-odd) coupling piece of the domain-wall fermion operator
    /*!
     * \param chi     result                    (Modify)
     * \param psi     source                    (Read)
     * \param isign   Flag ( PLUS | MINUS )     (Read)
     * \param cb      checkerboard ( 0 | 1 )    (Read)
     */
    void applyDiagInv(multi1d<LatticeFermion>& chi, 
		      const multi1d<LatticeFermion>& psi, 
		      enum PlusMinus isign,
		      int cb) const;
    

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
    void applyOffDiag(multi1d<LatticeFermion>& chi, 
		      const multi1d<LatticeFermion>& psi, 
		      enum PlusMinus isign,
		      int cb) const ;

    //! Apply the Dminus operator on a lattice fermion. See my notes ;-)
    void Dminus(LatticeFermion& chi,
		const LatticeFermion& psi,
		enum PlusMinus isign,
		int s5) const;

    

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
			   int cb) const ;

  protected:
    //! Partial constructor
    EvenOddPrecGenNEFDWLinOpArray() {}
    //! Partial constructor
    void operator=(const EvenOddPrecGenNEFDWLinOpArray&) {}

  private:
    Real WilsonMass;
    multi1d<Real> c5; // Nef Param Array
    multi1d<Real> b5; // Nef Param Array
    Real m_q;
    int  N5;

    // Derived quantities and constants 
    multi1d<Real> f_plus;   // f_{+}[i] = b_5[i]*(Nd - WilsonMass) + 1;
    multi1d<Real> f_minus;  // f_{-}[i] = c_5[i]*(Nd - WilsonMass) - 1;

    multi1d<Real> l;        // lowest row of L (left) array (N5-1 elems) 
    multi1d<Real> r;        // rightmost col of R (right) array (N5-1 elems)

    multi1d<Real> a;        // Subdiag elements of L'
                            // where L' is L in the LDU decomp of Tridiagonal 
                            // T

    multi1d<Real> b;        // Superdiag elements of U'
                            // where U' is U in LDU decomp of Tridiagonal 
                            // T

    multi1d<Real> d;        // Diag elems of D, where D is diagonal matrix
                            // in the LDU decomp of Tridiag T


    WilsonDslashArray  D;
  };

} // End Namespace Chroma



#endif
