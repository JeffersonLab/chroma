// -*- C++ -*-
/*! \file
 *  \brief Symmetric even-odd preconditioned Clover fermion linear operator
 */

#ifndef __seoprec_clover_linop_w_h__
#define __seoprec_clover_linop_w_h__

#include "state.h"
#include "fermbc.h"
#include "seoprec_logdet_linop.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/clover_term_w.h"
#include "coarsenable_op.h"


namespace Chroma 
{ 
  //! Symmetric even-odd preconditioned Clover-Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * The kernel for Clover fermions is
   *
   *      M  =  A + (d+M) - (1/2) D'
   *
   */
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  class SymEvenOddPrecCloverLinOpT : public SymEvenOddPrecLogDetLinearOperator<T_, multi1d<U_>, multi1d<U_> >,
				     public MGCoarsenableOperator<T_>
  {
  public:
    // Typedefs to save typing
    typedef T_           T;
    typedef multi1d<U_>  P;
    typedef multi1d<U_>  Q;

    //! A real number in the precision of T
    typedef typename std::conditional< std::is_same< typename WordType<T>::Type_t , REAL32 >::value ,
				       RealF , RealD >::type RealT;
    static_assert( std::is_same< typename WordType<RealT>::Type_t , typename WordType<T>::Type_t >::value ,
		   "SymEvenOddPrecCloverLinOpT: RealT must have the precision of T" );

    //! Partial constructor
    SymEvenOddPrecCloverLinOpT() {}

    //! Full constructor
    SymEvenOddPrecCloverLinOpT(Handle< FermState<T,P,Q> > fs,
			      const CloverFermActParams& param_)
    {
      create(fs,param_);
    }

    //! Destructor is automatic
    ~SymEvenOddPrecCloverLinOpT() {
    }

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const override {return D.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const CloverFermActParams& param_);

    //! Apply the the even-even block onto a source std::vector
    void unprecEvenEvenLinOp(T& chi, const T& psi,
		       enum PlusMinus isign) const override;

    //! Apply the inverse of the even-even block onto a source std::vector
    void unprecEvenEvenInvLinOp(T& chi, const T& psi,
			  enum PlusMinus isign) const override;
  
    //! Apply the the odd-odd block onto a source std::vector
     void unprecOddOddLinOp(T& chi, const T& psi,
 		     enum PlusMinus isign) const override;

     //! Apply the inverse of the odd-odd block onto a source std::vector
     void unprecOddOddInvLinOp(T& chi, const T& psi,
 			enum PlusMinus isign) const override;

    //! Apply the the even-odd block onto a source std::vector
    void unprecEvenOddLinOp(T& chi, const T& psi,
		      enum PlusMinus isign) const override;

    //! Apply the the odd-even block onto a source std::vector
    void unprecOddEvenLinOp(T& chi, const T& psi,
		      enum PlusMinus isign) const override;



    // Deriv of A_ee
    virtual void derivUnprecEvenEvenLinOp(P& ds_u, const T& chi, const T& psi,
    		enum PlusMinus isign) const override;

    // Deriv of A_oo
    virtual void derivUnprecOddOddLinOp(P& ds_u, const T& chi, const T& psi,
    		enum PlusMinus isign) const override;

    // Deriv of  D_eo
    virtual void derivUnprecEvenOddLinOp(P& ds_u, const T& chi, const T& psi,
    		enum PlusMinus isign) const override;

    // Deriv of  D_oe
    virtual void derivUnprecOddEvenLinOp(P& ds_u, const T& chi, const T& psi,
    		enum PlusMinus isign) const override;
  
    //! Apply the even-even block onto a source std::vector
    void derivLogDetEvenEvenLinOp(P& ds_u,
				  enum PlusMinus isign) const override;

    //! Apply the odd-odd block onto a source std::vector
    void derivLogDetOddOddLinOp(P& ds_u,
				enum PlusMinus isign) const override;

    // Override inherited one with a few more funkies
    //
    // Optimized operator() work on it when the category default works.
    void operator()(T& chi, const T& psi,
    		enum PlusMinus isign) const override;

    void deriv(P& ds_u, const T& chi, const T& psi,
    		enum PlusMinus isign) const override;


	//! Deriv
	void derivMultipole(P& ds_u, const multi1d<T>& chi,
			const multi1d<T>& psi,
			enum PlusMinus isign) const override;

    //! Return flops performed by the operator()
    unsigned long nFlops() const override;

    //! Get the log det of the even even part
    Double logDetEvenEvenLinOp(void) const override; 

    //! Get the log det of the odd odd part
    Double logDetOddOddLinOp(void) const override; 


    MGOperatorForm mgForm() const override
    {
      return MGOperatorForm::SYMMETRIC_EVEN_ODD_PRECONDITIONED;
    }

    const Subset& paritySubset(int p) const override { return rb[p]; }

    bool hasLocal    (int p, int q) const override { return p == q; }
    bool hasDirection(int p, int q) const override { return p != q; }

    void applyLocal(T& chi, const T& psi,
                    enum PlusMinus isign, int p, int q) const override
    {
      assert(p == q);
      chi[rb[p]] = psi;
      getFermBC().modifyF(chi, rb[p]);
    }

    void applyDirection(T& chi, const T& psi,
                        enum PlusMinus isign, int dir, int p, int q) const override
    {
      assert(q == 1 - p);
      T t;
      RealT mhalf = -0.5;
      if (isign == PLUS)
        {
          D.applyDirection(t, psi, isign, dir, p);  // output cb = p
          t[rb[p]] *= mhalf;
          invclov.apply(chi, t, isign, p);
        }
      else
        {
          invclov.apply(t, psi, isign, q);
          D.applyDirection(chi, t, isign, dir, p);
          chi[rb[p]] *= mhalf;
        }
      getFermBC().modifyF(chi, rb[p]);
    }

    void liftToEven(T& v_e, const T& v_o) const override
    {
      T t1, t2;
      D.apply(t1, v_o, PLUS, 0);                  // D_eo v_o
      invclov.apply(t2, t1, PLUS, 0);             // A_ee^{-1} D_eo v_o
      RealT half = 0.5;
      v_e[rb[0]] = half * t2;                     // -A_ee^{-1} (-1/2 D_eo) v_o
      getFermBC().modifyF(v_e, rb[0]);
    }

  private:
    CloverFermActParams param;
    DslashT D;
    CloverT clov;
    CloverT invclov;  // needed for even-even and odd-odd
    mutable double clov_apply_time;
    mutable double clov_deriv_time;
    mutable StopWatch swatch;
  };


  using SymEvenOddPrecCloverLinOp = SymEvenOddPrecCloverLinOpT<LatticeFermion, LatticeColorMatrix, WilsonDslash, CloverTerm>;
  extern template class SymEvenOddPrecCloverLinOpT<LatticeFermion, LatticeColorMatrix, WilsonDslash, CloverTerm>;

#if BASE_PRECISION == 64
  using SymEvenOddPrecCloverLinOpF = SymEvenOddPrecCloverLinOpT<LatticeFermionF, LatticeColorMatrixF, WilsonDslashF, CloverTermF>;
  extern template class SymEvenOddPrecCloverLinOpT<LatticeFermionF, LatticeColorMatrixF, WilsonDslashF, CloverTermF>;
#else
  using SymEvenOddPrecCloverLinOpF = SymEvenOddPrecCloverLinOp;
#endif



} // End Namespace Chroma


#endif
