// -*- C++ -*-
// $Id: ldag.h,v 1.1 2007-04-11 03:43:41 edwards Exp $
/*! \file
 * \brief M^dag composition of a linear operator
 */

#ifndef __ldag_h__
#define __ldag_h__

#include "handle.h"
#include "linearop.h"

namespace Chroma 
{ 
  //! M^dag linear operator
  /*!
   * \ingroup linop
   *
   * Linear operator forming M^dag from an operator M
   */
  template<typename T>
  class MdagLinOp : public LinearOperator<T>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    MdagLinOp(LinearOperator<T>* p) : A(p) {}

    //! Copy pointer (one more owner)
    MdagLinOp(Handle< LinearOperator<T> > p) : A(p) {}

    //! Destructor
    ~MdagLinOp() {}

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return A->subset();}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (T& chi, const T& psi, enum PlusMinus isign) const
      {
	T  tmp;  QDP::Hints::moveToFastMemoryHint(tmp);
	enum PlusMinus misign = (isign == PLUS) ? MINUS : PLUS;

	(*A)(chi, psi, misign);
      }

    unsigned long nFlops(void) const {
      unsigned long nflops=A->nFlops();
      return nflops;
    }

  private:
    Handle< LinearOperator<T> > A;
  };




  //! M^dag linear operator over arrays
  /*!
   * \ingroup linop
   *
   * Linear operator forming M^dag from an operator M
   */
  template<typename T>
  class MdagLinOpArray : public LinearOperatorArray<T>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    MdagLinOpArray(LinearOperatorArray<T>* p) : A(p) {}

    //! Copy pointer (one more owner)
    MdagLinOpArray(Handle< LinearOperatorArray<T> > p) : A(p) {}

    //! Destructor
    ~MdagLinOpArray() {}

    //! Length of array index
    int size() const {return A->size();}

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return A->subset();}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (multi1d<T>& chi, const multi1d<T>& psi, enum PlusMinus isign) const
      {
	multi1d<T>  tmp(size()); QDP::Hints::moveToFastMemoryHint(tmp);
	enum PlusMinus misign = (isign == PLUS) ? MINUS : PLUS;

	(*A)(chi, psi, misign);
      }

    unsigned long nFlops(void) const {
      unsigned long nflops=A->nFlops();
      return nflops;
    }

  private:
    Handle< LinearOperatorArray<T> > A;
  };


  //! M^dag linear operator
  /*!
   * \ingroup linop
   *
   * Linear operator forming M^dag from an operator M
   */
  template<typename T>
  class approx_lmdag : public LinearOperator<T>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    approx_lmdagm(LinearOperator<T>* p) : A(p) {}

    //! Copy pointer (one more owner)
    approx_lmdagm(Handle< LinearOperator<T> > p) : A(p) {}

    //! Destructor
    ~approx_lmdagm() {}

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return A->subset();}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (T& chi, const T& psi, enum PlusMinus isign) const
      {
	T  tmp;  QDP::Hints::moveToFastMemoryHint(tmp);
	enum PlusMinus misign = (isign == PLUS) ? MINUS : PLUS;

	(*A)(chi, psi, misign);
      }

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (T& chi, const T& psi, enum PlusMinus isign, Real epsilon) const
      {
	T  tmp; QDP::Hints::moveToFastMemoryHint(tmp);
	enum PlusMinus misign = (isign == PLUS) ? MINUS : PLUS;

	(*A)(chi, psi, misign, epsilon/Real(2));
      }

  private:
    Handle< LinearOperator<T> > A;
  };


  //! Differentiable M^dag linear operator
  /*!
   * \ingroup linop
   *
   * Diff. Linear operator forming M^dag from an operator M
   */
  template<typename T, typename P, typename Q>
  class DiffMdagLinOp : public DiffLinearOperator<T,P,Q>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    DiffMdagLinOp(DiffLinearOperator<T,P,Q>* p) : A(p) {}

    //! Copy pointer (one more owner)
    DiffMdagLinOp(Handle< DiffLinearOperator<T,P,Q> > p) : A(p) {}

    //! Destructor
    ~DiffMdagLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return A->getFermBC();}

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return A->subset();}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (T& chi, const T& psi, enum PlusMinus isign) const
      {
	T  tmp;  QDP::Hints::moveToFastMemoryHint(tmp);
	enum PlusMinus misign = (isign == PLUS) ? MINUS : PLUS;

	(*A)(chi, psi, misign);
      }

    //! Apply the derivative of the operator
    /*! 
     * Deriv of   chi^dag * A^dag * psi
     *
     * Default implementation
     */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
      {
	T   tmp; QDP::Hints::moveToFastMemoryHint(tmp);
	enum PlusMinus misign = (isign == PLUS) ? MINUS : PLUS;

	A->deriv(ds_u, chi, psi, misign);
      }


    //! Apply the derivative of the operator
    /*! 
     * Deriv of   chi^dag * A^dag * psi
     *
     * Default implementation
     */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign, const Real& epsilon) const
      {
	T   tmp;
	enum PlusMinus misign = (isign == PLUS) ? MINUS : PLUS;

	A->deriv(ds_u, chi, psi, misign, epsilon);
      }

    //! Return the number of flops performed by operator()
    unsigned long nFlops(void) const {return A->nFlops();}

  private:
    Handle< DiffLinearOperator<T,P,Q> > A;
  };




  //! M^dag linear operator over arrays
  /*!
   * \ingroup linop
   *
   * Linear operator forming M^dag from an operator M
   */
  template<typename T, typename P, typename Q>
  class DiffMdagLinOpArray : public DiffLinearOperatorArray<T,P,Q>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    DiffMdagLinOpArray(DiffLinearOperatorArray<T,P,Q>* p) : A(p) {}

    //! Copy pointer (one more owner)
    DiffMdagLinOpArray(Handle< DiffLinearOperatorArray<T,P,Q> > p) : A(p) {}

    //! Destructor
    ~DiffMdagLinOpArray() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return A->getFermBC();}

    //! Length of array index
    int size() const {return A->size();}

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return A->subset();}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (multi1d<T>& chi, const multi1d<T>& psi, enum PlusMinus isign) const
      {
	multi1d<T>  tmp(size()); QDP::Hints::moveToFastMemoryHint(tmp);
	enum PlusMinus misign = (isign == PLUS) ? MINUS : PLUS;

	(*A)(chi, psi, misign);
      }

    //! Apply the derivative of the operator
    /*! 
     * Deriv of   chi^dag * A^dag * psi
     *
     * Default implementation
     */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign) const
      {
	multi1d<T>  tmp(size()); QDP::Hints::moveToFastMemoryHint(tmp);
	enum PlusMinus misign = (isign == PLUS) ? MINUS : PLUS;

	A->deriv(ds_u, chi, psi, misign);
      }


    //! Apply the derivative of the operator
    /*! 
     * Deriv of   chi^dag * A^dag * psi
     *
     * Default implementation
     */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign, const Real& epsilon) const
      {
	multi1d<T>  tmp(size());
	enum PlusMinus misign = (isign == PLUS) ? MINUS : PLUS;

	A->deriv(ds_u, chi, psi, misign, epsilon);
      }

    //! Return the number of flops performed by operator()
    unsigned long nFlops(void) const {return A->nFlops();}

  private:
    Handle< DiffLinearOperatorArray<T,P,Q> > A;
  };


} // End Namespace Chroma


#endif
