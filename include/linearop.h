// -*- C++ -*-
// $Id: linearop.h,v 1.18 2003-12-02 15:41:31 edwards Exp $

/*! @file
 * @brief Linear Operators
 */

#ifndef __linearop_h__
#define __linearop_h__

using namespace QDP;

#include "state.h"


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

  //! Return the subset on which the operator acts
  virtual const OrderedSubset& subset() const = 0;

  //! Virtual destructor to help with cleanup;
  virtual ~LinearOperator() {}
};



//! Proxy for linear operator abstract class
/*! @ingroup linop
 *
 * Supports creation and application for linear operators that
 * hold things like Dirac operators, etc.
 */

template<typename T>
class LinearOperatorProxy : public LinearOperator<T>
{
public:
  //! Initialize pointer with existing pointer
  /*! Requires that the pointer p is a return value of new */
  explicit LinearOperatorProxy(const LinearOperator<T>* p=0) : linop(p) {}

  //! Copy pointer (one more owner)
  LinearOperatorProxy(const LinearOperatorProxy& a) : linop(a.linop) {}

  //! Assignment
  LinearOperatorProxy& operator=(const LinearOperatorProxy& a) 
    {linop = a.linop; return *this;}

  //! Access the value to which the pointer refers
  const LinearOperator<T>& operator*() const {return linop.operator*();}
  const LinearOperator<T>* operator->() const {return linop.operator->();}

  //! Apply the operator onto a source vector
  void operator() (T& chi, const T& psi, enum PlusMinus isign) const
    {linop->operator()(chi, psi, isign);}

  //! Return the subset on which the operator acts
  const OrderedSubset& subset() const
    {return linop->subset();}

private:
  Handle<const LinearOperator<T> >  linop;
};



//! Partial specialization of Linear Operator to arrays
/*! @ingroup linop
 *
 * Supports creation and application for linear operators that
 * hold things like Dirac operators, etc.
 */
template<typename T>
class LinearOperatorProxy< multi1d<T> > : public LinearOperator< multi1d<T> >
{
public:
  //! Initialize pointer with existing pointer
  /*! Requires that the pointer p is a return value of new */
  explicit LinearOperatorProxy(const LinearOperator< multi1d<T> >* p=0) : linop(p) {}

  //! Copy pointer (one more owner)
  LinearOperatorProxy(const LinearOperatorProxy& a) : linop(a.linop) {}

  //! Assignment
  LinearOperatorProxy& operator=(const LinearOperatorProxy& a) 
    {linop = a.linop; return *this;}

  //! Access the value to which the pointer refers
  const LinearOperator< multi1d<T> >& operator*() const {return linop.operator*();}
  const LinearOperator< multi1d<T> >* operator->() const {return linop.operator->();}

  //! Expected length of array index
  int size() const {return linop->size();}

  //! Apply the operator onto a source vector
  void operator() (multi1d<T>& chi, const multi1d<T>& psi, 
		   enum PlusMinus isign) const
    {linop->operator()(chi, psi, isign);}

  //! Return the subset on which the operator acts
  const OrderedSubset& subset() const {return linop->subset();}

private:
  Handle<const LinearOperator< multi1d<T> > >  linop;
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
 *      [      1      	      0        ]
 *      [       E,E            E,O     ]
 *  L = [                              ]
 *      [     D     A^(-1)    1        ]
 *      [      O,E   E,E        O,O    ]
 *
 * and
 *
 *      [      A      	      D       ]
 *      [       E,E            E,O    ]
 *  U = [                             ]
 *      [      0              1       ]
 *      [       O,E            O,O    ]
 *
 * The preconditioned matrix is formed from
 *
 *  M'   =  L^-1 * M * U^-1
 *
 * Resulting in a new  M
 *
 *      [      1       	            0                      ]
 *      [       E,E                  E,O                   ]
 *      [                                                  ]
 *      [      0       	        A     -  D    A^(-1)  D    ]
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
class EvenOddPrecLinearOperator : public LinearOperator<T>
{
public:
  //! Only defined on the odd lattice
  const OrderedSubset& subset() const {return rb[1];}

  //! Another way of saying only defined on the odd lattice
  virtual const int subsetCB() const {return 1;}

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
  virtual ~EvenOddPrecLinearOperator() {}
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
 *      [      1      	      0        ]
 *      [       E,E            E,O     ]
 *  L = [                              ]
 *      [     D     A^(-1)    1        ]
 *      [      O,E   E,E        O,O    ]
 *
 * and
 *
 *      [      A      	      D       ]
 *      [       E,E            E,O    ]
 *  U = [                             ]
 *      [      0              1       ]
 *      [       O,E            O,O    ]
 *
 * The preconditioned matrix is formed from
 *
 *  M'   =  L^-1 * M * U^-1
 *
 * Resulting in a new  M
 *
 *      [      1       	            0                      ]
 *      [       E,E                  E,O                   ]
 *      [                                                  ]
 *      [      0       	        A     -  D    A^(-1)  D    ]
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
class EvenOddPrecLinearOperator< multi1d<T> > : public LinearOperator< multi1d<T> >
{
public:
  //! Only defined on the odd lattice
  const OrderedSubset& subset() const {return rb[1];}

  //! Another way of saying only defined on the odd lattice
  virtual const int subsetCB() const {return 1;}

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
  virtual ~EvenOddPrecLinearOperator() {}
};



//! Proxy for even-odd preconditioned linear operator abstract class
/*! @ingroup linop
 *
 * Support for even-odd preconditioned linear operators
 */
template<typename T>
class EvenOddPrecLinearOperatorProxy : public EvenOddPrecLinearOperator<T>
{
public:
  //! Initialize pointer with existing pointer
  /*! Requires that the pointer p is a return value of new */
  explicit EvenOddPrecLinearOperatorProxy(const EvenOddPrecLinearOperator<T>* p=0) : linop(p) {}

  //! Copy pointer (one more owner)
  EvenOddPrecLinearOperatorProxy(const EvenOddPrecLinearOperatorProxy& a) : linop(a.linop) {}

  //! Assignment
  EvenOddPrecLinearOperatorProxy& operator=(const EvenOddPrecLinearOperatorProxy& a)
    {linop = a.linop; return *this;}

  //! Access the value to which the pointer refers
  const EvenOddPrecLinearOperator<T>& operator*() const {return linop.operator*();}
  const EvenOddPrecLinearOperator<T>* operator->() const {return linop.operator->();}

  //! Only defined on the odd lattice
  const OrderedSubset& subset() const {return linop->subset();}

  //! Another way of saying only defined on the odd lattice
  virtual const int subsetCB() const {return linop->subsetCB();}

  //! Apply the even-even block onto a source vector
  /*! This does not need to be optimized */
  virtual void evenEvenLinOp(T& chi, const T& psi, 
			     enum PlusMinus isign) const
    {linop->evenEvenLinOp(chi, psi, isign);}
  
  //! Apply the inverse of the even-even block onto a source vector
  virtual void evenEvenInvLinOp(T& chi, const T& psi, 
				enum PlusMinus isign) const
    {linop->evenEvenInvLinOp(chi, psi, isign);}
  
  //! Apply the the even-odd block onto a source vector
  virtual void evenOddLinOp(T& chi, const T& psi, 
			    enum PlusMinus isign) const
    {linop->evenOddLinOp(chi, psi, isign);}

  //! Apply the the odd-even block onto a source vector
  virtual void oddEvenLinOp(T& chi, const T& psi, 
			    enum PlusMinus isign) const
    {linop->oddEvenLinOp(chi, psi, isign);}

  //! Apply the the odd-odd block onto a source vector
  virtual void oddOddLinOp(T& chi, const T& psi, 
			   enum PlusMinus isign) const
    {linop->oddOddLinOp(chi, psi, isign);}

  //! Apply the operator onto a source vector
  virtual void operator() (T& chi, const T& psi, 
			   enum PlusMinus isign) const
    {linop->operator()(chi, psi, isign);}

  //! Apply the UNPRECONDITIONED operator onto a source vector
  /*! Mainly intended for debugging */
  virtual void unprecLinOp(T& chi, const T& psi, 
			   enum PlusMinus isign) const
    {linop->unprecLinOp(chi, psi, isign);}

private:
  Handle<const EvenOddPrecLinearOperator<T> >  linop;
};


//! Proxy for even-odd preconditioned linear operator abstract class
/*! @ingroup linop
 *
 * Support for even-odd preconditioned linear operators
 */
template<typename T>
class EvenOddPrecLinearOperatorProxy< multi1d<T> > : public EvenOddPrecLinearOperator< multi1d<T> >
{
public:
  //! Initialize pointer with existing pointer
  /*! Requires that the pointer p is a return value of new */
  explicit EvenOddPrecLinearOperatorProxy(const EvenOddPrecLinearOperator< multi1d<T> >* p=0) : linop(p) {}

  //! Copy pointer (one more owner)
  EvenOddPrecLinearOperatorProxy(const EvenOddPrecLinearOperatorProxy& p) : linop(p.linop) {}

  //! Assignment
  EvenOddPrecLinearOperatorProxy& operator=(const EvenOddPrecLinearOperatorProxy& a)
    {linop = a.linop; return *this;}

  //! Access the value to which the pointer refers
  const EvenOddPrecLinearOperator< multi1d<T> >& operator*() const {return linop.operator*();}
  const EvenOddPrecLinearOperator< multi1d<T> >* operator->() const {return linop.operator->();}

  //! Expected length of array index
  int size() const {return linop->size();}

  //! Only defined on the odd lattice
  const OrderedSubset& subset() const {return linop->subset();}

  //! Another way of saying only defined on the odd lattice
  virtual const int subsetCB() const {return linop->subsetCB();}

  //! Apply the even-even block onto a source vector
  /*! This does not need to be optimized */
  virtual void evenEvenLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
			     enum PlusMinus isign) const
    {linop->evenEvenLinOp(chi, psi, isign);}
  
  //! Apply the inverse of the even-even block onto a source vector
  virtual void evenEvenInvLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
				enum PlusMinus isign) const
    {linop->evenEvenInvLinOp(chi, psi, isign);}
  
  //! Apply the the even-odd block onto a source vector
  virtual void evenOddLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
			    enum PlusMinus isign) const
    {linop->evenOddLinOp(chi, psi, isign);}

  //! Apply the the odd-even block onto a source vector
  virtual void oddEvenLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
			    enum PlusMinus isign) const
    {linop->oddEvenLinOp(chi, psi, isign);}

  //! Apply the the odd-odd block onto a source vector
  virtual void oddOddLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
			   enum PlusMinus isign) const
    {linop->oddOddLinOp(chi, psi, isign);}

  //! Apply the operator onto a source vector
  virtual void operator() (multi1d<T>& chi, const multi1d<T>& psi, 
			   enum PlusMinus isign) const
    {linop->operator()(chi, psi, isign);}

  //! Apply the UNPRECONDITIONED operator onto a source vector
  /*! Mainly intended for debugging */
  virtual void unprecLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
			   enum PlusMinus isign) const
    {linop->unprecLinOp(chi, psi, isign);}

private:
  Handle<const EvenOddPrecLinearOperator< multi1d<T> > >  linop;
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
