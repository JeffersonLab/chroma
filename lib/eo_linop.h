// -*- C++ -*-
// $Id: eo_linop.h,v 3.2 2007-02-22 21:11:44 bjoo Exp $

/*! @file
 * @brief Linear Operators
 */

#ifndef __eo_linop_h__
#define __eo_linop_h__

#include "linearop.h"

namespace Chroma
{

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
    const Subset& subset() const {return all;}

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

}


#endif
