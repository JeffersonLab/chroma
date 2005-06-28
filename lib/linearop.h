// -*- C++ -*-
// $Id: linearop.h,v 1.15 2005-06-28 15:28:15 bjoo Exp $

/*! @file
 * @brief Linear Operators
 */

#ifndef __linearop_h__
#define __linearop_h__

#include "chromabase.h"

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
    virtual const unsigned long nFlops() { return 0; }
  };


  //! Partial specialization of Linear Operator to arrays
  /*! @ingroup linop
   *
   * Supports creation and application for linear operators that
   * hold things like Dirac operators, etc.
   */
  template<typename T>
  class LinearOperator< multi1d<T> >
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~LinearOperator() {}

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
    virtual const unsigned long nFlops() const { return 0; };
  };


  //-----------------------------------------------------------------------------------
  //! Differentiable Linear Operator
  /*! @ingroup linop
   *
   * Supports creation and application for linear operators that
   * hold things like Dirac operators, etc. that are differentiable
   */
  template<typename T, typename P>
  class DiffLinearOperator : public LinearOperator<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~DiffLinearOperator() {}

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


  //----------------------------------------------------------------
  //! Unpreconditioned linear operator including derivatives
  /*! @ingroup linop
   *
   * Support for unpreconditioned linear operators with derivative
   */
  template<typename T, typename P>
  class UnprecLinearOperator : public DiffLinearOperator<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~UnprecLinearOperator() {}

    //! Only defined on the entire lattice
    const OrderedSubset& subset() const {return all;}
  };


  //! Even-odd preconditioned linear operator
  /*! @ingroup linop
   *
   * Support for even-odd preconditioned linear operators
   * Given a matrix M written in block form:
   *
   *      [      A             D        ]
   *      [       E,E           E,O     ]
   *      [                             ]
   *      [      D             A        ]
   *      [       O,E           O,O     ]
   *
   * The preconditioning consists of using the triangular matrices
   *
   *      [      1              0        ]
   *      [       E,E            E,O     ]
   *  L = [                              ]
   *      [     D     A^(-1)    1        ]
   *      [      O,E   E,E        O,O    ]
   *
   * and
   *
   *      [      A              D       ]
   *      [       E,E            E,O    ]
   *  U = [                             ]
   *      [      0              1       ]
   *      [       O,E            O,O    ]
   *
   * The preconditioned matrix is formed from
   *
   *  M'   =  L^-1 * M * U^-1
   *
   * where
   *
   *           [      1              0        ]
   *           [       E,E            E,O     ]
   *  L^(-1) = [                              ]
   *           [   - D     A^(-1)    1        ]
   *           [      O,E   E,E        O,O    ]
   *
   * and
   *
   *           [      A^(-1)       - A^(-1) D       ]
   *           [       E,E            E,E    E,O    ]
   *  U^(-1) = [                                    ]
   *           [      0                1            ]
   *           [       O,E              O,O         ]
   *
   * Resulting in a new  M
   *
   *      [      1                    0                      ]
   *      [       E,E                  E,O                   ]
   *      [                                                  ]
   *      [      0                A     -  D    A^(-1)  D    ]
   *      [       O,E              O,O      O,E   E,E    E,O ]
   *
   *
   * This class is used to implement the resulting linear operator
   *
   *      ~
   *      M  =  A(o,o) - D(o,e) . A^-1(e,e) . D(e,o)
   *
   * By construction, the linear operator is ONLY defined on the odd subset
   *
   */

  template<typename T, typename P>
  class EvenOddPrecLinearOperator : public DiffLinearOperator<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecLinearOperator() {}

    //! Only defined on the odd lattice
    const OrderedSubset& subset() const {return rb[1];}

    //! Apply the even-even block onto a source vector
    /*! This does not need to be optimized */
    virtual void evenEvenLinOp(T& chi, const T& psi, 
			       enum PlusMinus isign) const = 0;
  
    //! Apply the inverse of the even-even block onto a source vector
    virtual void evenEvenInvLinOp(T& chi, const T& psi, 
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
    virtual void operator() (T& chi, const T& psi, 
			     enum PlusMinus isign) const
    {
      T   tmp1, tmp2; moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

      /*  Tmp1   =  D     A^(-1)     D    Psi  */
      /*      O      O,E        E,E   E,O    O */
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      oddEvenLinOp(tmp1, tmp2, isign);

      /*  Chi   =  A    Psi  -  Tmp1  */
      /*     O      O,O    O        O */
      oddOddLinOp(chi, psi, isign);
      chi[rb[1]] -= tmp1;
    }


    //! Apply the UNPRECONDITIONED operator onto a source vector
    /*! Mainly intended for debugging */
    virtual void unprecLinOp(T& chi, const T& psi, 
			     enum PlusMinus isign) const
    {
      T   tmp1, tmp2;  moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

      /*  Chi   =   A    Psi   +    D    Psi   */
      /*     E       E,E    O        E,O    O  */
      evenEvenLinOp(tmp1, psi, isign);
      evenOddLinOp(tmp2, psi, isign);
      chi[rb[0]] = tmp1 + tmp2;

      /*  Chi   =  A    Psi  -  Tmp1  */
      /*     O      O,O    O        O */
      oddEvenLinOp(tmp1, psi, isign);
      oddOddLinOp(tmp2, psi, isign);
      chi[rb[1]] = tmp1 + tmp2;
    }


    //! Apply the even-even block onto a source vector
    virtual void derivEvenEvenLinOp(P& ds_u, const T& chi, const T& psi, 
				    enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
  
    //! Apply the the even-odd block onto a source vector
    virtual void derivEvenOddLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
 
    //! Apply the the odd-even block onto a source vector
    virtual void derivOddEvenLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the the odd-odd block onto a source vector
    virtual void derivOddOddLinOp(P& ds_u, const T& chi, const T& psi, 
				  enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the derivative of the operator onto a source vector
    /*! User should make sure deriv routines do a resize  */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
    {
      // Need deriv of  (A_oo - D_oe*Ainv_ee*D_eo*psi_e)
      enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

      //
      // Make sure the deriv routines do a resize !!!
      //
      T   tmp1, tmp2, tmp3;  // if an array is used here, the space is not reserved
      moveToFastMemoryHint(tmp1);
      moveToFastMemoryHint(tmp2);
      moveToFastMemoryHint(tmp3);

      P   ds_1;  // deriv routines should resize

      //
      // NOTE: even with even-odd decomposition, the ds_u will still have contributions
      // on all cb. So, no adding of ds_1 onto ds_u under a subset
      //
      //  ds_u  =  chi^dag * A'_oo * psi
      derivOddOddLinOp(ds_u, chi, psi, isign);

      //  ds_u  -=  chi^dag * D'_oe * Ainv_ee * D_eo * psi_o
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      derivOddEvenLinOp(ds_1, chi, tmp2, isign);
      ds_u -= ds_1;

      //  ds_u  +=  chi^dag * D_oe * Ainv_ee * A'_ee * Ainv_ee * D_eo * psi_o
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      evenOddLinOp(tmp1, chi, msign);
      evenEvenInvLinOp(tmp3, tmp1, msign);
      derivEvenEvenLinOp(ds_1, tmp3, tmp2, isign);
      ds_u += ds_1;

      //  ds_u  -=  chi^dag * D_oe * Ainv_ee * D'_eo * psi_o
      evenOddLinOp(tmp1, chi, msign);
      evenEvenInvLinOp(tmp3, tmp1, msign);
      derivEvenOddLinOp(ds_1, tmp3, psi, isign);
      ds_u -= ds_1;
    }

    //! Apply the derivative of the UNPRECONDITIONED operator onto a source vector
    /*! Mainly intended for debugging */
    virtual void derivUnprecLinOp(P& ds_u, const T& chi, const T& psi, 
				  enum PlusMinus isign) const
    {
      T   tmp1, tmp2;  moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);
      P   ds_tmp;  // deriv routines should resize

      //  ds_u  =  chi^dag * A'_ee * psi
      derivEvenEvenLinOp(ds_u, chi, psi, isign);

      //  ds_u +=  chi^dag * D'_eo * psi
      derivEvenOddLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;

      //  ds_u +=  chi^dag * D'_oe * psi
      derivOddEvenLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;

      //  ds_u +=  chi^dag * A'_oo * psi
      derivOddOddLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;
    }

    //! Return flops performed by the evenEvenLinOp
    virtual const unsigned long evenEvenNFlops() const { return 0; }
    
    //! Return flops performed by the evenOddLinOp
    virtual const unsigned long evenOddNFlops() const { return 0; }

    //! Return flops performed by the oddEvenLinOp
    virtual const unsigned long oddEvenNFlops() const { return 0; }

    //! Return flops performed by the oddOddLinOp
    virtual const unsigned long oddOddNFlops() const { return 0; }

    //! Return flops performed by the evenEvenInvLinOp
    virtual const unsigned long evenEvenInvNFlops() const { return 0; }

    //! Return flops performed by the operator()
    virtual const unsigned long nFlops() const { 
      return 0;
    }

  };


  //! Partial specialization of even-odd preconditioned linear operator including derivatives
  /*! @ingroup linop
   *
   * Support for even-odd preconditioned linear operators with derivatives
   *
   * Given a matrix M written in block form:
   *
   *      [      A             D        ]
   *      [       E,E           E,O     ]
   *      [                             ]
   *      [      D             A        ]
   *      [       O,E           O,O     ]
   *
   * The preconditioning consists of using the triangular matrices
   *
   *      [      1              0        ]
   *      [       E,E            E,O     ]
   *  L = [                              ]
   *      [     D     A^(-1)    1        ]
   *      [      O,E   E,E        O,O    ]
   *
   * and
   *
   *      [      A              D       ]
   *      [       E,E            E,O    ]
   *  U = [                             ]
   *      [      0              1       ]
   *      [       O,E            O,O    ]
   *
   * The preconditioned matrix is formed from
   *
   *  M'   =  L^-1 * M * U^-1
   *
   * where
   *
   *           [      1              0        ]
   *           [       E,E            E,O     ]
   *  L^(-1) = [                              ]
   *           [   - D     A^(-1)    1        ]
   *           [      O,E   E,E        O,O    ]
   *
   * and
   *
   *           [      A^(-1)       - A^(-1) D       ]
   *           [       E,E            E,E    E,O    ]
   *  U^(-1) = [                                    ]
   *           [      0                1            ]
   *           [       O,E              O,O         ]
   *
   * Resulting in a new  M
   *
   *      [      1                    0                      ]
   *      [       E,E                  E,O                   ]
   *      [                                                  ]
   *      [      0                A     -  D    A^(-1)  D    ]
   *      [       O,E              O,O      O,E   E,E    E,O ]
   *
   *
   * This class is used to implement the resulting linear operator
   *
   *      ~
   *      M  =  A(o,o) - D(o,e) . A^-1(e,e) . D(e,o)
   *
   * By construction, the linear operator is ONLY defined on the odd subset
   *
   */

  template<typename T, typename P>
  class EvenOddPrecLinearOperator< multi1d<T>, P > : public DiffLinearOperator< multi1d<T>, P >
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecLinearOperator() {}

    //! Only defined on the odd lattice
    const OrderedSubset& subset() const {return rb[1];}

    //! Expected length of array index
    virtual int size() const = 0;

    //! Apply the even-even block onto a source vector
    /*! This does not need to be optimized */
    virtual void evenEvenLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
			       enum PlusMinus isign) const = 0;
  
    //! Apply the inverse of the even-even block onto a source vector
    virtual void evenEvenInvLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
				  enum PlusMinus isign) const = 0;
  
    //! Apply the the even-odd block onto a source vector
    virtual void evenOddLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
			      enum PlusMinus isign) const = 0;

    //! Apply the the odd-even block onto a source vector
    virtual void oddEvenLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
			      enum PlusMinus isign) const = 0;

    //! Apply the the odd-odd block onto a source vector
    virtual void oddOddLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
			     enum PlusMinus isign) const = 0;

    //! Apply the operator onto a source vector
    virtual void operator() (multi1d<T>& chi, const multi1d<T>& psi, 
			     enum PlusMinus isign) const
    {
      multi1d<T>  tmp1(size());  moveToFastMemoryHint(tmp1);
      multi1d<T>  tmp2(size());  moveToFastMemoryHint(tmp2);

      /*  Tmp1   =  D     A^(-1)     D    Psi  */
      /*      O      O,E        E,E   E,O    O */
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      oddEvenLinOp(tmp1, tmp2, isign);

      /*  Chi   =  A    Psi  -  Tmp1  */
      /*     O      O,O    O        O */
      oddOddLinOp(chi, psi, isign);
      for(int n=0; n < size(); ++n)
	chi[n][rb[1]] -= tmp1[n];
    }


    //! Apply the UNPRECONDITIONED operator onto a source vector
    /*! Mainly intended for debugging */
    virtual void unprecLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
			     enum PlusMinus isign) const
    {
      multi1d<T>  tmp1(size()); moveToFastMemoryHint(tmp1);
      multi1d<T>  tmp2(size()); moveToFastMemoryHint(tmp2);

      /*  Chi   =   A    Psi   +    D    Psi   */
      /*     E       E,E    O        E,O    O  */
      evenEvenLinOp(tmp1, psi, isign);
      evenOddLinOp(tmp2, psi, isign);
      for(int n=0; n < size(); ++n)
	chi[n][rb[0]] = tmp1[n] + tmp2[n];

      /*  Chi   =   D    Psi    +    A    Psi   */
      /*     O       O,E    E         O,O    O  */
      oddEvenLinOp(tmp1, psi, isign);
      oddOddLinOp(tmp2, psi, isign);
      for(int n=0; n < size(); ++n)
	chi[n][rb[1]] = tmp1[n] + tmp2[n];
    }


    //! Apply the even-even block onto a source vector
    virtual void derivEvenEvenLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				    enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
  
    //! Apply the the even-odd block onto a source vector
    virtual void derivEvenOddLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
 
    //! Apply the the odd-even block onto a source vector
    virtual void derivOddEvenLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the the odd-odd block onto a source vector
    virtual void derivOddOddLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				  enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the operator onto a source vector
    /*! User should make sure deriv routines do a resize  */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign) const
    {
      // Need deriv of  (A_oo - D_oe*Ainv_ee*D_eo*psi_e)
      enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

      //
      // Make sure the deriv routines do a resize !!!
      //
      multi1d<T>   tmp1(size()), tmp2(size()), tmp3(size());  // if an array is used here, the space is not reserved
      moveToFastMemoryHint(tmp1);
      moveToFastMemoryHint(tmp2);
      moveToFastMemoryHint(tmp3);

      P            ds_1;  // deriv routines should resize

      //
      // NOTE: even with even-odd decomposition, the ds_u will still have contributions
      // on all cb. So, no adding of ds_1 onto ds_u under a subset
      //
      //  ds_u  =  chi^dag * A'_oo * psi
      derivOddOddLinOp(ds_u, chi, psi, isign);

      //  ds_u  -=  chi^dag * D'_oe * Ainv_ee * D_eo * psi_o
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      derivOddEvenLinOp(ds_1, chi, tmp2, isign);
      ds_u -= ds_1;

      //  ds_u  +=  chi^dag * D_oe * Ainv_ee * A'_ee * Ainv_ee * D_eo * psi_o
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      evenOddLinOp(tmp1, chi, msign);
      evenEvenInvLinOp(tmp3, tmp1, msign);
      derivEvenEvenLinOp(ds_1, tmp3, tmp2, isign);
      ds_u += ds_1;

      //  ds_u  -=  chi^dag * D_oe * Ainv_ee * D'_eo * psi_o
      evenOddLinOp(tmp1, chi, msign);
      evenEvenInvLinOp(tmp3, tmp1, msign);
      derivEvenOddLinOp(ds_1, tmp3, psi, isign);
      ds_u -= ds_1;
    }


    //! Apply the derivative of the UNPRECONDITIONED operator onto a source vector
    /*! Mainly intended for debugging */
    virtual void derivUnprecLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				  enum PlusMinus isign) const
    {
      T   tmp1, tmp2;  // if an array is used here, the space is not reserved
      moveToFastMemoryHint(tmp1);
      moveToFastMemoryHint(tmp2);

      P   ds_tmp;  // deriv routines should resize

      //  ds_u  =  chi^dag * A'_ee * psi
      derivEvenEvenLinOp(ds_u, chi, psi, isign);

      //  ds_u +=  chi^dag * D'_eo * psi
      derivEvenOddLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;

      //  ds_u +=  chi^dag * D'_oe * psi
      derivOddEvenLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;

      //  ds_u +=  chi^dag * A'_oo * psi
      derivOddOddLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;
    }

    //! Return flops performed by the evenEvenLinOp
    virtual const unsigned long evenEvenNFlops() const { return 0; }
    
    //! Return flops performed by the evenOddLinOp
    virtual const unsigned long evenOddNFlops() const { return 0; }

    //! Return flops performed by the oddEvenLinOp
    virtual const unsigned long oddEvenNFlops() const { return 0; }

    //! Return flops performed by the oddOddLinOp
    virtual const unsigned long oddOddNFlops() const { return 0; }

    //! Return flops performed by the evenEvenInvLinOp
    virtual const unsigned long evenEvenInvNFlops() const { return 0; }

    //! Return flops performed by the operator()
    virtual const unsigned long nFlops() const { 
      return (this->oddOddNFlops()
	      +this->oddEvenNFlops()
	      +this->evenEvenInvNFlops()
	      +this->evenOddNFlops());
    }


  };


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
  template<typename T, typename P>
  class EvenOddLinearOperator : public DiffLinearOperator<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddLinearOperator() {}

    //! Only defined on the even lattice
    const OrderedSubset& subset() const {return all;}

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
    }

    //! Return flops performed by the evenEvenLinOp
    virtual const unsigned long evenEvenNFlops() const { return 0; }
    
    //! Return flops performed by the evenOddLinOp
    virtual const unsigned long evenOddNFlops() const { return 0; }

    //! Return flops performed by the oddEvenLinOp
    virtual const unsigned long oddEvenNFlops() const { return 0; }

    //! Return flops performed by the oddOddLinOp
    virtual const unsigned long oddOddNFlops() const { return 0; }


    //! Return flops performed by the operator()
    virtual const unsigned long nFlops() const { 
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
  template<typename T, typename P>
  class DslashLinearOperator : public DiffLinearOperator<T,P>
  {
  public:
    //! Virtual destructor to help in cleanup
    virtual ~DslashLinearOperator() {}

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

}



#endif
