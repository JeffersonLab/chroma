// -*- C++ -*-
// $Id: linearop.h,v 3.0 2006-04-03 04:58:44 edwards Exp $

/*! @file
 * @brief Linear Operators
 */

#ifndef __linearop_h__
#define __linearop_h__

#include "chromabase.h"
#include "fermbc.h"

using namespace QDP::Hints;

namespace Chroma
{

  //-----------------------------------------------------------------------------------
  //! Linear Operator
  /*! @ingroup linop
   *
   * Supports creation and application for linear operators that
   * hold things like Dirac operators, etc.
   */
  template<typename T>
  class LinearOperator
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~LinearOperator() {}

    //! Apply the operator onto a source vector
    virtual void operator() (T& chi, const T& psi, enum PlusMinus isign) const = 0;

    //! Apply the operator onto a source vector to some precision
    virtual void operator() (T& chi, const T& psi, enum PlusMinus isign, 
			     Real epsilon) const
    {
      (*this)(chi,psi,isign);
    }

    //! Return the subset on which the operator acts
    virtual const OrderedSubset& subset() const = 0;
    
    //! Return the number of flops performed by operator()
    //! Base linop returns 0 and this can be overridden
    virtual unsigned long nFlops() const { return 0; }
  };


  //-----------------------------------------------------------------------------------
  //! Linear Operator to arrays
  /*! @ingroup linop
   *
   * Supports creation and application for linear operators that
   * hold things like Dirac operators, etc.
   */
  template<typename T>
  class LinearOperatorArray
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~LinearOperatorArray() {}

    //! Expected length of array index
    virtual int size() const = 0;

    //! Apply the operator onto a source vector
    virtual void operator() (multi1d<T>& chi, const multi1d<T>& psi, 
			     enum PlusMinus isign) const = 0;

    //! Apply the operator onto a source vector to some precision
    virtual void operator() (multi1d<T>& chi, const multi1d<T>& psi, 
			     enum PlusMinus isign, Real epsilon) const
    {
      (*this)(chi,psi,isign);
    }

    //! Return the subset on which the operator acts
    virtual const OrderedSubset& subset() const = 0;

    //! Return the number of flops performed by operator()
    //! Base linop returns 0 and this can be overridden
    virtual unsigned long nFlops() const { return 0; };
  };


  //-----------------------------------------------------------------------------------
  //! Differentiable Linear Operator
  /*! @ingroup linop
   *
   * Supports creation and application for linear operators that
   * hold things like Dirac operators, etc. that are differentiable
   */
  template<typename T, typename P, typename Q>
  class DiffLinearOperator : public LinearOperator<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~DiffLinearOperator() {}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Apply the derivative of the operator onto a source vector
    /*! Default implementation */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
    {
      QDPIO::cerr << "deriv: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the derivative of the operator onto a source vector to some precision
    /*! Default implementation */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign, const Real& epsilon) const
    {
      deriv(ds_u,chi,psi,isign);
    }
  };


  //-----------------------------------------------------------------------------------
  //! Differentiable Linear Operator
  /*! @ingroup linop
   *
   * Supports creation and application for linear operators that
   * hold things like Dirac operators, etc. that are differentiable
   */
  template<typename T, typename P, typename Q>
  class DiffLinearOperatorArray : public LinearOperatorArray<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~DiffLinearOperatorArray() {}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Apply the derivative of the operator onto a source vector
    /*! Default implementation */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign) const
    {
      QDPIO::cerr << "deriv: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the derivative of the operator onto a source vector to some precision
    /*! Default implementation */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign, const Real& epsilon) const
    {
      deriv(ds_u,chi,psi,isign);
    }
  };


  //----------------------------------------------------------------
  //! Unpreconditioned linear operator including derivatives
  /*! @ingroup linop
   *
   * Support for unpreconditioned linear operators with derivative
   */
  template<typename T, typename P, typename Q>
  class UnprecLinearOperator : public DiffLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~UnprecLinearOperator() {}

    //! Only defined on the entire lattice
    const OrderedSubset& subset() const {return all;}
  };


  //----------------------------------------------------------------
  //! Unpreconditioned linear operator including derivatives
  /*! @ingroup linop
   *
   * Support for unpreconditioned linear operators with derivative
   */
  template<typename T, typename P, typename Q>
  class UnprecLinearOperatorArray : public DiffLinearOperatorArray<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~UnprecLinearOperatorArray() {}

    //! Only defined on the entire lattice
    const OrderedSubset& subset() const {return all;}
  };


  //----------------------------------------------------------------
  //#include "eo_prec_linop.h"
  //#include "prec_constdet_linop.h"
  //#include "prec_logdet_linop.h"
  //---------------------------------------------------------------------
  //! Even odd Linear Operator (for staggered like things )
  /*! @ingroup linop
   *
   * Support for even-odd staggered-like linear operators
   *
   *  [   D_ee        D_eo ]
   *  [   D_oe        D_oo ]
   *
   *  Usually D_ee = D_oo = 2m
   */
  template<typename T, typename P, typename Q>
  class EvenOddLinearOperator : public DiffLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddLinearOperator() {}

    //! Only defined on the even lattice
    const OrderedSubset& subset() const {return all;}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Apply the even-even block onto a source vector
    virtual void evenEvenLinOp(T& chi, const T& psi, 
			       enum PlusMinus isign) const = 0;
  
    //! Apply the the even-odd block onto a source vector
    virtual void evenOddLinOp(T& chi, const T& psi, 
			      enum PlusMinus isign) const = 0;

    //! Apply the the odd-even block onto a source vector
    virtual void oddEvenLinOp(T& chi, const T& psi, 
			      enum PlusMinus isign) const = 0;

    //! Apply the the odd-odd block onto a source vector
    virtual void oddOddLinOp(T& chi, const T& psi, 
			     enum PlusMinus isign) const = 0;

    //! Apply the operator onto a source vector
    /*! User should make sure deriv routines do a resize  */
    virtual void operator() (T& chi, const T& psi, 
			     enum PlusMinus isign) const
    {
      T   tmp1, tmp2;  // if an array is used here, the space is not reserved
      moveToFastMemoryHint(tmp1);
      moveToFastMemoryHint(tmp2);
 
      /*  Chi   =   D    Psi   +    D    Psi   */
      /*     E       E,E    E        E,O    O  */
      evenEvenLinOp(tmp1, psi, isign);
      evenOddLinOp(tmp2, psi, isign);
      chi[rb[0]] = tmp1 + tmp2;

      /*  Chi   =  D    Psi    +  D    Psi  */
      /*     O      O,E    E       O,O    O */
      oddEvenLinOp(tmp1, psi, isign);
      oddOddLinOp(tmp2, psi, isign);
      chi[rb[1]] = tmp1 + tmp2;

      getFermBC().modifyF(chi, rb[1]);
    }

    //! Apply the even-even block onto a source vector
    virtual void derivEvenEvenLinOp(P& ds_u, 
				    const T& chi, const T& psi, 
				    enum PlusMinus isign) const 
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
 
    //! Apply the the even-odd block onto a source vector
    virtual void derivEvenOddLinOp(P& ds_u, 
				   const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the the odd-even block onto a source vector
    virtual void derivOddEvenLinOp(P& ds_u, 
				   const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the the odd-odd block onto a source vector
    virtual void derivOddOddLinOp(P& ds_u, 
				  const T& chi, const T& psi, 
				  enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the operator onto a source vector
    virtual void deriv(P& ds_u, 
		       const T& chi, const T& psi, 
		       enum PlusMinus isign) const
    {
      // Need deriv of  chi_e^dag * (D_ee * psi_e + D_eo * psi_i)
      // Need deriv of  chi_o^dag * (D_oe * psi_e + D_oo * psi_i)

      //
      // NOTE: even with even-odd decomposition, the ds_u will still have contributions
      // on all cb. So, no adding of ds_1 onto ds_u under a subset
      //
      T   tmp1, tmp2;  // if an array is used here, the space is not reserved
      moveToFastMemoryHint(tmp1);
      moveToFastMemoryHint(tmp2);

      P   ds_1;  // routines should resize

      // ds_u = chi_e ^dag * D'_ee * psi_e
      derivEvenEvenLinOp(ds_u, chi, psi, isign);

      // ds_u += chi_e ^dag * D'_eo * psi_o
      derivEvenOddLinOp(ds_1, chi, psi, isign);
      ds_u += ds_1;

      // ds_u += chi_o ^dag * D'_oe * psi_e
      derivOddEvenLinOp(ds_1, chi, psi, isign);
      ds_u += ds_1;

      // ds_u += chi_o ^dag * D'_oo * psi_o
      derivOddOddLinOp(ds_1, chi, psi, isign);
      ds_u += ds_1;

      getFermBC().zero(ds_u);
    }

    //! Return flops performed by the evenEvenLinOp
    virtual unsigned long evenEvenNFlops() const { return 0; }
    
    //! Return flops performed by the evenOddLinOp
    virtual unsigned long evenOddNFlops() const { return 0; }

    //! Return flops performed by the oddEvenLinOp
    virtual unsigned long oddEvenNFlops() const { return 0; }

    //! Return flops performed by the oddOddLinOp
    virtual unsigned long oddOddNFlops() const { return 0; }


    //! Return flops performed by the operator()
    virtual unsigned long nFlops() const { 
      return this->oddOddNFlops()
	+this->oddEvenNFlops()
	+this->evenEvenNFlops()
	+this->evenOddNFlops();
    }
  };


  //-----------------------------------------------------------------------------------
  //! Dslash-like Linear Operator
  /*! @ingroup linop
   *
   * These are concessions/optimizations for red-black checkboarding 
   */
  template<typename T, typename P, typename Q>
  class DslashLinearOperator : public DiffLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help in cleanup
    virtual ~DslashLinearOperator() {}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Apply operator on both checkerboards (entire lattice)
    virtual void operator() (T& d, const T& psi, enum PlusMinus isign) const
    {
      apply(d, psi, isign, 0);
      apply(d, psi, isign, 1);
    }

    //! Apply checkerboarded linear operator
    /*! 
     * To avoid confusion (especially of the compilers!), call the checkerboarded
     * apply instead of operator()
     */
    virtual void apply (T& chi, const T& psi, enum PlusMinus isign, int cb) const = 0;


    //! Take deriv of D
    /*!
     * \param chi     left vector                                 (Read)
     * \param psi     right vector                                (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     *
     * \return Computes   chi^dag * \dot(D} * psi  
     */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
    {
      QDPIO::cerr << "deriv: not implemented" << endl;
      QDP_abort(1);
    }

    //! Take deriv of D
    /*!
     * \param chi     left vector on cb                           (Read)
     * \param psi     right vector on 1-cb                        (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of chi vector                  (Read)
     *
     * \return Computes   chi^dag * \dot(D} * psi  
     */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign, int cb) const
    {
      QDPIO::cerr << "deriv: not implemented" << endl;
      QDP_abort(1);
    }
  };



  //-----------------------------------------------------------------------------------
  //! Dslash-like Linear Operator for arrays
  /*! @ingroup linop
   *
   * These are concessions/optimizations for red-black checkboarding 
   */
  template<typename T, typename P, typename Q>
  class DslashLinearOperatorArray : public DiffLinearOperatorArray<T,P,Q>
  {
  public:
    //! Virtual destructor to help in cleanup
    virtual ~DslashLinearOperatorArray() {}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Expected length of array index
    virtual int size() const = 0;

    //! Apply operator on both checkerboards (entire lattice)
    virtual void operator() (multi1d<T>& d, const multi1d<T>& psi, enum PlusMinus isign) const
    {
      d.resize(size());
      apply(d, psi, isign, 0);
      apply(d, psi, isign, 1);
    }

    //! Apply checkerboarded linear operator
    /*! 
     * To avoid confusion (especially of the compilers!), call the checkerboarded
     * apply instead of operator()
     */
    virtual void apply (multi1d<T>& chi, const multi1d<T>& psi, enum PlusMinus isign, int cb) const = 0;

    //! Take deriv of D
    /*!
     * \param chi     left vector                                 (Read)
     * \param psi     right vector                                (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     *
     * \return Computes   chi^dag * \dot(D} * psi  
     */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign) const
    {
      QDPIO::cerr << "deriv: not implemented" << endl;
      QDP_abort(1);
    }

    //! Take deriv of D
    /*!
     * \param chi     left vector on cb                           (Read)
     * \param psi     right vector on 1-cb                        (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of chi vector                  (Read)
     *
     * \return Computes   chi^dag * \dot(D} * psi  
     */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign, int cb) const
    {
      QDPIO::cerr << "deriv: not implemented" << endl;
      QDP_abort(1);
    }
  };

}



#endif
