// -*- C++ -*-
// $Id: llincomb.h,v 3.2 2007-02-22 21:11:46 bjoo Exp $

#ifndef __llincomb_h__
#define __llincomb_h__

#include "handle.h"
#include "linearop.h"


namespace Chroma
{
  //! Linear combination of a Linear Operator
  /*!
   * \ingroup linop
   *
   * This operator does the linear combination  (z1 + z2*B)*psi
   */
  template<typename T, class C>
  class llincomb : public LinearOperator<T>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    llincomb(const LinearOperator<T>* p, const C& add_const_, const C& scale_fact_) : 
      A(p), add_const(add_const_), scale_fact(scale_fact_)  {}

    //! Copy pointer (one more owner)
    llincomb(Handle<const LinearOperator<T> > p, const C& add_const_, const C& scale_fact_) : 
      A(p), add_const(add_const_), scale_fact(scale_fact_) {}

    //! Destructor
    ~llincomb() {}

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return A->subset();}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (T& chi, const T& psi, enum PlusMinus isign) const
      {
	const Subset& sub = A->subset();
	(*A)(chi, psi, isign);
	chi[sub] *= scale_fact;
	chi[sub] += add_const*psi;
      }

  private:
    const Handle< const LinearOperator<T> > A;
    const C add_const;
    const C scale_fact;
  };
}


#endif
