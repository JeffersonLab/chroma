// -*- C++ -*-
/*! \file
 *  \brief Even-odd preconditioned Clover fermion linear operator
 */

#ifndef __prec_clover_linop_w_h__
#define __prec_clover_linop_w_h__

#include "state.h"
#include "fermbc.h"
#include "eoprec_logdet_linop.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/clover_term_w.h"
#include "coarsenable_op.h"


namespace Chroma 
{ 
  //! Even-odd preconditioned Clover-Dirac operator
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
  class EvenOddPrecCloverLinOpT : public EvenOddPrecLogDetLinearOperator<T_, multi1d<U_>, multi1d<U_> >,
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
		   "EvenOddPrecCloverLinOpT: RealT must have the precision of T" );

    //! Partial constructor
    EvenOddPrecCloverLinOpT() {}

    //! Full constructor
    EvenOddPrecCloverLinOpT(Handle< FermState<T,P,Q> > fs,
			    const CloverFermActParams& param_)
      {
	create(fs,param_);
      }

    //! Destructor is automatic
    ~EvenOddPrecCloverLinOpT() {
      QDPIO::cout << "CLOV_LINOP: Time spent in clov deriv (total) = " << clov_deriv_time << std::endl;
      QDPIO::cout << "CLOV_LINOP: Time spent in clov apply/invapply (total) = " << clov_apply_time << std::endl;

    }

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const override {return D.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const CloverFermActParams& param_);

    //! Apply the the even-even block onto a source std::vector
    void evenEvenLinOp(T& chi, const T& psi, 
		       enum PlusMinus isign) const override;

    //! Apply the inverse of the even-even block onto a source std::vector
    void evenEvenInvLinOp(T& chi, const T& psi, 
			  enum PlusMinus isign) const override;
  
    //! Apply the the even-odd block onto a source std::vector
    void evenOddLinOp(T& chi, const T& psi, 
		      enum PlusMinus isign) const override;

    //! Apply the the odd-even block onto a source std::vector
    void oddEvenLinOp(T& chi, const T& psi, 
		      enum PlusMinus isign) const override;

    //! Apply the the odd-odd block onto a source std::vector
    void oddOddLinOp(T& chi, const T& psi, 
		     enum PlusMinus isign) const override;

    // Override inherited one with a few more funkies
    void operator()(T& chi, const T& psi, 
		    enum PlusMinus isign) const override;

    //! Apply the even-even block onto a source std::vector
    void derivEvenEvenLinOp(P& ds_u, 
			    const T& chi, const T& psi, 
			    enum PlusMinus isign) const override;

  //! Apply the even-even block onto a source std::vector
    void derivEvenEvenLinOpMP(P& ds_u, 
			    const multi1d<T>& chi, const multi1d<T>& psi, 
			    enum PlusMinus isign) const override;

    void derivLogDetEvenEvenLinOp(P& ds_u, 
				  enum PlusMinus isign) const override;

    //! Apply the the even-odd block onto a source std::vector
    void derivEvenOddLinOp(P& ds_u, 
			   const T& chi, const T& psi, 
			   enum PlusMinus isign) const override;
 
    //! Apply the the odd-even block onto a source std::vector
    void derivOddEvenLinOp(P& ds_u, 
			   const T& chi, const T& psi, 
			   enum PlusMinus isign) const override;

    //! Apply the the odd-odd block onto a source std::vector
    void derivOddOddLinOp(P& ds_u, 
			  const T& chi, const T& psi, 
			  enum PlusMinus isign) const override;

    void derivOddOddLinOpMP(P& ds_u, 
			    const multi1d<T>& chi, const multi1d<T>& psi, 
			    enum PlusMinus isign) const override;

    //! Return flops performed by the operator()
    unsigned long nFlops() const override;

    //! Get the log det of the even even part
    Double logDetEvenEvenLinOp(void) const override; 


    MGOperatorForm mgForm() const override
    {
      return MGOperatorForm::EVEN_ODD_PRECONDITIONED;
    }

    const Subset& paritySubset(int p) const override { return rb[p]; }

    bool hasLocal    (int p, int q) const override { return p == q; }
    bool hasDirection(int p, int q) const override { return p != q; }

    void applyLocal(T& chi, const T& psi,
                    enum PlusMinus isign, int p, int q) const override
    {
      assert(p == q);
      clov.apply(chi, psi, isign, p);          // A_ee (p=0) or A_oo (p=1)
      getFermBC().modifyF(chi, rb[p]);
    }

    void applyDirection(T& chi, const T& psi,
                        enum PlusMinus isign, int dir, int p, int q) const override
    {
      assert(q == 1 - p);
      D.applyDirection(chi, psi, isign, dir, p); // H_eo (p=0) or H_oe (p=1)
      RealT mhalf = -0.5;
      chi[rb[p]] *= mhalf;
      getFermBC().modifyF(chi, rb[p]);
    }

    void liftToEven(T& v_e, const T& v_o) const override
    {
      T t1, t2;
      evenOddLinOp(t1, v_o, PLUS);             // H_eo v_o   (includes -1/2)
      evenEvenInvLinOp(t2, t1, PLUS);          // A_ee^{-1} H_eo v_o
      v_e[rb[0]] = -t2;
    }

    
  private:
    mutable T tmp1;
    mutable T tmp2;
    CloverFermActParams param;
    DslashT D;
    CloverT clov;
    CloverT invclov;  // uggh, only needed for evenEvenLinOp
    mutable double clov_apply_time;
    mutable double clov_deriv_time;
    mutable StopWatch swatch;
  };


  //! In the base precision
  using EvenOddPrecCloverLinOp = EvenOddPrecCloverLinOpT<LatticeFermion, LatticeColorMatrix, WilsonDslash, CloverTerm>;
  extern template class EvenOddPrecCloverLinOpT<LatticeFermion, LatticeColorMatrix, WilsonDslash, CloverTerm>;

#if BASE_PRECISION == 64
  //! In single precision
  using EvenOddPrecCloverLinOpF = EvenOddPrecCloverLinOpT<LatticeFermionF, LatticeColorMatrixF, WilsonDslashF, CloverTermF>;
  extern template class EvenOddPrecCloverLinOpT<LatticeFermionF, LatticeColorMatrixF, WilsonDslashF, CloverTermF>;
#else
  using EvenOddPrecCloverLinOpF = EvenOddPrecCloverLinOp;
#endif


} // End Namespace Chroma


#endif
