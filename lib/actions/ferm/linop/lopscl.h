// -*- C++ -*-
// $Id: lopscl.h,v 3.1 2007-02-22 21:11:46 bjoo Exp $

#ifndef __lopscl_h__
#define __lopscl_h__

#include "handle.h"
#include "linearop.h"


namespace Chroma 
{ 
  //! Scaled Linear Operator
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
    lopscl(LinearOperator<T>* p, const C& scale_fact_) : A(p), scale_fact(scale_fact_)  {}

    //! Copy pointer (one more owner)
    lopscl(Handle< LinearOperator<T> > p, const C& scale_fact_) : A(p), scale_fact(scale_fact_) {}

    //! Destructor
    ~lopscl() {}

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return A->subset();}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (T& chi, const T& psi, enum PlusMinus isign) const
      {
	const Subset& sub = A->subset();
	(*A)(chi, psi, isign);
	chi[sub] *= scale_fact;
      }

  private:
    Handle< LinearOperator<T> > A;
    const C scale_fact;
  };

  //! Scaled Linear Operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * This operator scales its input operator
   */
  template<typename T, class C>
  class approx_lopscl : public LinearOperator<T>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    approx_lopscl(LinearOperator<T>* p, const C& scale_fact_) : A(p), scale_fact(scale_fact_)  {}

    //! Copy pointer (one more owner)
    approx_lopscl(Handle< LinearOperator<T> > p, const C& scale_fact_) : A(p), scale_fact(scale_fact_) {}

    //! Destructor
    ~approx_lopscl() {}

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return A->subset();}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (T& chi, const T& psi, enum PlusMinus isign) const
      {
	const Subset& sub = A->subset();
	(*A)(chi, psi, isign);
	chi[sub] *= scale_fact;
      }

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (T& chi, const T& psi, enum PlusMinus isign, Real epsilon) const
      {
	const Subset& sub = A->subset();
	(*A)(chi, psi, isign, epsilon);
	chi[sub] *= scale_fact;
      }

  private:
    Handle< LinearOperator<T> > A;
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
  class lopsclArray : public LinearOperatorArray<T>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    lopsclArray(LinearOperatorArray<T>* p, const C& scale_fact_) : A(p), scale_fact(scale_fact_)  {}

    //! Copy pointer (one more owner)
    lopsclArray(Handle< LinearOperatorArray<T> > p, const C& scale_fact_) : A(p), scale_fact(scale_fact_)  {}

    //! Destructor
    ~lopsclArray() {}

    //! Length of array index
    int size() const {return A->size();}

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return A->subset();}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (multi1d<T>& chi, const multi1d<T>& psi, enum PlusMinus isign) const
      {
	const Subset& sub = A->subset();

	(*A)(chi, psi, isign);
	for(int i = 0; i < size(); i++) { 
	  chi[i][sub] *= scale_fact;
	}

      }

  private:
    Handle< LinearOperatorArray<T> > A;
    const C scale_fact;
  };


} // End Namespace Chroma


#endif
