// -*- C++ -*-
// $Id: lopscl.h,v 1.1 2004-01-12 15:59:45 bjoo Exp $

#ifndef __lopscl_h__
#define __lopscl_h__

#include "handle.h"
#include "linearop.h"

using namespace QDP;

//! M^dag.M linear operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 *
 * This operator scales its input operator
 */
template<typename T, class C>
class lopscl : public LinearOperator<T>
{
public:
  //! Initialize pointer with existing pointer
  /*! Requires that the pointer p is a return value of new */
  lopscl(const LinearOperator<T>* p, const C& scale_fact_) : A(p), scale_fact(scale_fact_)  {}

  //! Copy pointer (one more owner)
  lopscl(Handle<const LinearOperator<T> > p, const C& scale_fact_) : A(p), scale_fact(scale_fact_) {}

  //! Destructor
  ~lopscl() {}

  //! Subset comes from underlying operator
  inline const OrderedSubset& subset() const {return A->subset();}

  //! Apply the operator onto a source vector
  /*! For this operator, the sign is ignored */
  inline void operator() (T& chi, const T& psi, enum PlusMinus isign) const
    {
      T  tmp;
      const OrderedSubset& sub = A->subset();
      (*A)(chi, psi, isign);
      chi[sub] *= scale_fact;
    }

private:
  const Handle< const LinearOperator<T> > A;
  const C scale_fact;
};



//! Partial specialization of scaled M  operator over arrays
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 *
 * Linear operator forming M^dag.M from an operator M
 */
template<typename T, class C>
class lopscl< multi1d<T>, C > : public LinearOperator< multi1d<T> >
{
public:
  //! Initialize pointer with existing pointer
  /*! Requires that the pointer p is a return value of new */
  lopscl(const LinearOperator< multi1d<T> >* p, const C& scale_fact_) : A(p), scale_fact(scale_fact_)  {}

  //! Copy pointer (one more owner)
  lopscl(Handle<const LinearOperator< multi1d<T> > > p, const C& scale_fact_) : A(p), scale_fact(scale_fact_)  {}

  //! Destructor
  ~lopscl() {}

  //! Length of array index
  int size() const {return A->size();}

  //! Subset comes from underlying operator
  inline const OrderedSubset& subset() const {return A->subset();}

  //! Apply the operator onto a source vector
  /*! For this operator, the sign is ignored */
  inline void operator() (multi1d<T>& chi, const multi1d<T>& psi, enum PlusMinus isign) const
    {
      multi1d<T>  tmp(size());
      const OrderedSubset& sub = A->subset();

      (*A)(chi, psi, isign);
      for(int i = 0; i < size(); i++) { 
	chi[i][sub] *= scale_fact;
      }

    }

private:
  const Handle< const LinearOperator< multi1d<T> > > A;
  const C scale_fact;
};

#endif
