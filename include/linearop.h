// -*- C++ -*-
// $Id: linearop.h,v 1.16 2003-11-22 21:35:20 edwards Exp $

/*! @file
 * @brief Linear Operators
 */

#ifndef __linearop_h__
#define __linearop_h__

using namespace QDP;

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

  //! Return the subset on which the operator acts
  virtual const OrderedSubset& subset() const = 0;

  //! Virtual destructor to help with cleanup;
  virtual ~LinearOperator() {}
};


//! Even-odd preconditioned linear operator
/*! @ingroup linop
 *
 * Support for even-odd preconditioned linear operators
 * Given a matrix M written in block form:
 *
 *      [      A             D        ]
 *      [        E,E           E,O    ]
 *      [                             ]
 *      [      D             A        ]
 *      [        O,E           O,O    ]

 * The preconditioning consists of left multiplying by
 *
 *      [      1      	       0       ]
 *      [        E,E            E,O    ]
 *      [                              ]
 *      [     D     A^(-1)     1       ]
 *      [      O,E   E,E        O,O    ]
 *
 *
 * Resulting in a new  M
 *
 *      [      A       	            D                      ]
 *      [        E,E                 E,O                   ]
 *      [                                                  ]
 *      [      0       	        A     -  D    A^(-1)  D    ]
 *      [        O,E             O,O      O,E   E,E    E,O ]
 *
 *
 * This class is used to implement the resulting linear opeator
 *
 *      ~
 *      M  =  A(o,o) - D(o,e) . A^-1(e,e) . D(e,o)
 *
 * By construction, the linear operator is ONLY defined on the odd subset
 *
 */

template<typename T>
class EvenOddPrecLinearOperator : public LinearOperator<T>
{
public:
  //! Only defined on the odd lattice
  const OrderedSubset& subset() const {return rb[1];}

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

  //! Virtual destructor to help with cleanup;
  virtual ~EvenOddPrecLinearOperator() {}
};



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


#endif
