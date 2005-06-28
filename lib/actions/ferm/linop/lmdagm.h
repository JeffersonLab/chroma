// -*- C++ -*-
// $Id: lmdagm.h,v 1.6 2005-06-28 15:28:15 bjoo Exp $

#ifndef __lmdagm_w_h__
#define __lmdagm_w_h__

#include "handle.h"
#include "linearop.h"

using namespace QDP::Hints;

namespace Chroma 
{ 
//! M^dag.M linear operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 *
 * Linear operator forming M^dag.M from an operator M
 */
template<typename T>
class lmdagm : public LinearOperator<T>
{
public:
  //! Initialize pointer with existing pointer
  /*! Requires that the pointer p is a return value of new */
  lmdagm(const LinearOperator<T>* p) : A(p) {}

  //! Copy pointer (one more owner)
  lmdagm(Handle<const LinearOperator<T> > p) : A(p) {}

  //! Destructor
  ~lmdagm() {}

  //! Subset comes from underlying operator
  inline const OrderedSubset& subset() const {return A->subset();}

  //! Apply the operator onto a source vector
  /*! For this operator, the sign is ignored */
  inline void operator() (T& chi, const T& psi, enum PlusMinus isign) const
    {
      T  tmp;  moveToFastMemoryHint(tmp);

      (*A)(tmp, psi, PLUS);
      (*A)(chi, tmp, MINUS);
    }

private:
  const Handle< const LinearOperator<T> > A;
};




//! Partial specialization of M^dag.M linear operator over arrays
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 *
 * Linear operator forming M^dag.M from an operator M
 */
template<typename T>
class lmdagm< multi1d<T> > : public LinearOperator< multi1d<T> >
{
public:
  //! Initialize pointer with existing pointer
  /*! Requires that the pointer p is a return value of new */
  lmdagm(const LinearOperator< multi1d<T> >* p) : A(p) {}

  //! Copy pointer (one more owner)
  lmdagm(Handle<const LinearOperator< multi1d<T> > > p) : A(p) {}

  //! Destructor
  ~lmdagm() {}

  //! Length of array index
  int size() const {return A->size();}

  //! Subset comes from underlying operator
  inline const OrderedSubset& subset() const {return A->subset();}

  //! Apply the operator onto a source vector
  /*! For this operator, the sign is ignored */
  inline void operator() (multi1d<T>& chi, const multi1d<T>& psi, enum PlusMinus isign) const
    {
      multi1d<T>  tmp(size()); moveToFastMemoryHint(tmp);
      (*A)(tmp, psi, PLUS);
      (*A)(chi, tmp, MINUS);
    }

private:
  const Handle< const LinearOperator< multi1d<T> > > A;
};

//! M^dag.M linear operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 *
 * Linear operator forming M^dag.M from an operator M
 */
template<typename T>
class approx_lmdagm : public LinearOperator<T>
{
public:
  //! Initialize pointer with existing pointer
  /*! Requires that the pointer p is a return value of new */
  approx_lmdagm(const LinearOperator<T>* p) : A(p) {}

  //! Copy pointer (one more owner)
  approx_lmdagm(Handle<const LinearOperator<T> > p) : A(p) {}

  //! Destructor
  ~approx_lmdagm() {}

  //! Subset comes from underlying operator
  inline const OrderedSubset& subset() const {return A->subset();}

  //! Apply the operator onto a source vector
  /*! For this operator, the sign is ignored */
  inline void operator() (T& chi, const T& psi, enum PlusMinus isign) const
    {
      T  tmp;  moveToFastMemoryHint(tmp);

      (*A)(tmp, psi, PLUS);
      (*A)(chi, tmp, MINUS);
    }

  //! Apply the operator onto a source vector
  /*! For this operator, the sign is ignored */
  inline void operator() (T& chi, const T& psi, enum PlusMinus isign, Real epsilon) const
    {
      T  tmp; moveToFastMemoryHint(tmp);

      (*A)(tmp, psi, PLUS, epsilon/Real(2));
      (*A)(chi, tmp, MINUS, epsilon/Real(2));
    }

private:
  const Handle< const LinearOperator<T> > A;
};


}; // End Namespace Chroma


#endif
