// -*- C++ -*-
// $Id: linearop.h,v 1.7 2004-12-12 21:22:14 edwards Exp $

/*! @file
 * @brief Linear Operators
 */

#ifndef __linearop_h__
#define __linearop_h__

using namespace QDP;

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

    //! Virtual destructor to help with cleanup;
    virtual ~LinearOperator() {}
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

    //! Virtual destructor to help with cleanup;
    virtual ~LinearOperator() {}
  };


#if 0
  //-----------------------------------------------------------------------------------
  //! Approximate Linear Operator
  /*!
   * @ingroup linop
   *
   * These are linear operators that in some way require an accuracy
   * parameter in their operator() for some solvers. Examples are
   * the multi-pole approximated sign function. Where the epsilon
   * would refer to the accuracy of the multi shift solver
   */
  template<typename T>
  class ApproxLinearOperator : public LinearOperator<T>
  {
  public:
    //! Return the subset on which the operator acts
    virtual const OrderedSubset& subset() const = 0;

    //! Apply the operator onto a source vector
    virtual void operator() (T& chi, const T& psi, enum PlusMinus isign) const = 0;

    virtual ~ApproxLinearOperator() {}
  };
#endif

  //-----------------------------------------------------------------------------------
  //! Unpreconditioned linear operator
  /*! @ingroup linop
   *
   * Support for unpreconditioned linear operators
   */
  template<typename T>
  class UnprecLinearOperatorBase : public LinearOperator<T>
  {
  public:
    //! Only defined on the entire lattice
    const OrderedSubset& subset() const {return all;}

    //! Virtual destructor to help with cleanup;
    virtual ~UnprecLinearOperatorBase() {}
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

  template<typename T>
  class EvenOddPrecLinearOperatorBase : public LinearOperator<T>
  {
  public:
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
      T   tmp1, tmp2;  // if an array is used here, the space is not reserved

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
      T   tmp1, tmp2;  // if an array is used here, the space is not reserved

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

    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecLinearOperatorBase() {}
  };


  //! Partial specialization of even-odd preconditioned linear operator to arrays
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

  template<typename T>
  class EvenOddPrecLinearOperatorBase< multi1d<T> > : public LinearOperator< multi1d<T> >
  {
  public:
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
      multi1d<T>  tmp1(size());
      multi1d<T>  tmp2(size());

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
      multi1d<T>  tmp1(size());
      multi1d<T>  tmp2(size());

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

    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecLinearOperatorBase() {}
  };


#if 0
  //-----------------------------------------------------------------------------------
  //! Approximate versions of unpreconditioned linear operators
  template<typename T>
  class ApproxUnprecLinearOperatorBase : public virtual UnprecLinearOperatorBase<T>, public virtual ApproxLinearOperator<T>
  {
  public:
  };


  //! Approximate versions of even-odd preconditioned linear operators
  template<typename T>
  class ApproxEvenOddPrecLinearOperatorBase : public virtual EvenOddPrecLinearOperatorBase<T>, public virtual ApproxLinearOperator<T>
  {
  public:
  };
#endif


  //-----------------------------------------------------------------------------------
  //! Unpreconditioned linear operator including derivatives
  /*! @ingroup linop
   *
   * Support for unpreconditioned linear operators with derivative
   */
  template<typename T, typename P>
  class UnprecLinearOperator : public UnprecLinearOperatorBase<T>
  {
  public:
    //! Apply the operator onto a source vector
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
    {
      QDPIO::cerr << "deriv: not implemented" << endl;
      QDP_abort(1);
    }

    //! Virtual destructor to help with cleanup;
    virtual ~UnprecLinearOperator() {}
  };


  //! Even-odd preconditioned linear operator including derivatives
  /*! @ingroup linop
   *
   * Defines derivative of linear operator
   */
  template<typename T, typename P>
  class EvenOddPrecLinearOperator : public EvenOddPrecLinearOperatorBase<T>
  {
  public:
    //! Only defined on the odd lattice
    const OrderedSubset& subset() const {return rb[1];}

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

    //! Apply the operator onto a source vector
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecLinearOperator() {}
  };


  //! Even-odd preconditioned linear operator including derivatives
  /*! @ingroup linop
   *
   * Defines derivative of linear operator
   */
  template<typename T>
  class EvenOddPrecLinearOperator<T, multi1d<LatticeColorMatrix> > : public EvenOddPrecLinearOperatorBase<T>
  {
  public:
    //! Only defined on the odd lattice
    const OrderedSubset& subset() const {return rb[1];}

    //! Apply the even-even block onto a source vector
    virtual void derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, const T& chi, const T& psi, 
				    enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
  
    //! Apply the the even-odd block onto a source vector
    virtual void derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
 
    //! Apply the the odd-even block onto a source vector
    virtual void derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the the odd-odd block onto a source vector
    virtual void derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, const T& chi, const T& psi, 
				  enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the operator onto a source vector
    /*! User should make sure deriv routines do a resize  */
    virtual void deriv(multi1d<LatticeColorMatrix>& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
    {
      // Need deriv of  (A_oo - D_oe*Ainv_ee*D_eo*psi_e)
      enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

      //
      // Make sure the deriv routines do a resize !!!
      //
      ds_u.resize(Nd);

      T   tmp1, tmp2, tmp3;  // if an array is used here, the space is not reserved
      multi1d<LatticeColorMatrix>   ds_1;  // deriv routines should resize

      //  ds_o  =  chi^dag * A'_oo * psi
      derivOddOddLinOp(ds_u, chi, psi, isign);

      //  ds_o  -=  chi^dag * D'_oe * Ainv_ee * D_eo * psi_o
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      derivOddEvenLinOp(ds_1, chi, tmp2, isign);
      for(int mu=0; mu < ds_u.size(); ++mu)
	ds_u[mu][rb[1]] -= ds_1[mu];

      //  ds_o  +=  chi^dag * D_oe * Ainv_ee * A'_ee * Ainv_ee * D_eo * psi_o
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      evenOddLinOp(tmp1, chi, msign);
      evenEvenInvLinOp(tmp3, tmp1, msign);
      derivEvenEvenLinOp(ds_1, tmp3, tmp2, isign);
      for(int mu=0; mu < ds_u.size(); ++mu)
	ds_u[mu][rb[1]] += ds_1[mu];

      //  ds_o  -=  chi^dag * D_oe * Ainv_ee * D'_eo * psi_o
      evenOddLinOp(tmp1, chi, msign);
      evenEvenInvLinOp(tmp3, tmp1, msign);
      derivEvenOddLinOp(ds_1, tmp3, psi, isign);
      for(int mu=0; mu < ds_u.size(); ++mu)
	ds_u[mu][rb[1]] -= ds_1[mu];
    }

    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecLinearOperator() {}
  };


#if 0
  //-----------------------------------------------------------------------------------
  //! Approximate versions of unpreconditioned linear operators with derivatives
  template<typename T, typename P>
  class ApproxUnprecLinearOperator : public virtual UnprecLinearOperator<T,P>, public virtual ApproxLinearOperator<T>
  {
  public:
  };


  //! Approximate versions of even-odd preconditioned linear operators with derivatives
  template<typename T, typename P>
  class ApproxEvenOddPrecLinearOperator : public virtual EvenOddPrecLinearOperator<T,P>, public virtual ApproxLinearOperator<T>
  {
  public:
  };
#endif


  //-----------------------------------------------------------------------------------
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
  template<typename T>
  class EvenOddLinearOperatorBase : public LinearOperator<T>
  {
  public:
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

    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddLinearOperatorBase() {}
  };


  //! Even odd Linear Operator (for staggered like things ) including derivatives
  /*! @ingroup linop
   *
   * Support for even-odd staggered-like linear operators including derivatives
   *
   *  [   D_ee        D_eo ]
   *  [   D_oe        D_oo ]
   *
   *  Usually D_ee = D_oo = 2m
   */
  template<typename T, typename P>
  class EvenOddLinearOperator : public EvenOddLinearOperatorBase<T>
  {
  public:
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

    //! Apply the operator onto a source vector
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddLinearOperator() {}
  };


  //! Even odd Linear Operator (for staggered like things ) including derivatives
  /*! @ingroup linop
   *
   * Support for even-odd staggered-like linear operators including derivatives
   *
   *  [   D_ee        D_eo ]
   *  [   D_oe        D_oo ]
   *
   *  Usually D_ee = D_oo = 2m
   */
  template<typename T>
  class EvenOddLinearOperator<T, multi1d<LatticeColorMatrix> > : public EvenOddLinearOperatorBase<T>
  {
  public:
    //! Apply the even-even block onto a source vector
    virtual void derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
				    const T& chi, const T& psi, 
				    enum PlusMinus isign) const 
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
 
    //! Apply the the even-odd block onto a source vector
    virtual void derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
				   const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the the odd-even block onto a source vector
    virtual void derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
				   const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the the odd-odd block onto a source vector
    virtual void derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
				  const T& chi, const T& psi, 
				  enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the operator onto a source vector
    virtual void deriv(multi1d<LatticeColorMatrix>& ds_u, 
		       const T& chi, const T& psi, 
		       enum PlusMinus isign) const
    {
      // Need deriv of  chi_e^dag * (D_ee * psi_e + D_eo * psi_i)
      // Need deriv of  chi_o^dag * (D_oe * psi_e + D_oo * psi_i)

      ds_u.resize(Nd);
      T   tmp1, tmp2;  // if an array is used here, the space is not reserved
      multi1d<LatticeColorMatrix>   ds_1;  // routines should resize

      // ds_e = chi_e ^dag * D'_ee * psi_e
      derivEvenEvenLinOp(ds_u, chi, psi, isign);

      // ds_e += chi_e ^dag * D'_eo * psi_o
      derivEvenOddLinOp(ds_1, chi, psi, isign);
      for(int mu=0; mu < Nd; ++mu)
	ds_u[mu][rb[0]] += ds_1[mu];

      // ds_o += chi_o ^dag * D'_oe * psi_e
      derivOddEvenLinOp(ds_1, chi, psi, isign);
      for(int mu=0; mu < Nd; ++mu)
	ds_u[mu][rb[1]] += ds_1[mu];

      // ds_o += chi_o ^dag * D'_oo * psi_o
      derivOddOddLinOp(ds_1, chi, psi, isign);
      for(int mu=0; mu < Nd; ++mu)
	ds_u[mu][rb[1]] += ds_1[mu];
    }

    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddLinearOperator() {}
  };


  //-----------------------------------------------------------------------------------
  //! Dslash-like Linear Operator
  /*! @ingroup linop
   *
   * These are concessions/optimizations for red-black checkboarding 
   */
  template<typename T>
  class DslashLinearOperator : public LinearOperator<T>
  {
  public:
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

    //! Virtual destructor to help in cleanup
    virtual ~DslashLinearOperator() {}
  };

}

using namespace Chroma;


#endif
