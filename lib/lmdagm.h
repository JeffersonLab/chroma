// -*- C++ -*-
// $Id: lmdagm.h,v 3.3 2007-04-11 03:44:00 edwards Exp $
/*! \file
 * \brief M^dag*M composition of a linear operator
 */

#ifndef __lmdagm_h__
#define __lmdagm_h__

#include "handle.h"
#include "linearop.h"

namespace Chroma 
{ 
  //! M^dag.M linear operator
  /*!
   * \ingroup linop
   *
   * Linear operator forming M^dag.M from an operator M
   */
  template<typename T>
  class MdagMLinOp : public LinearOperator<T>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    MdagMLinOp(LinearOperator<T>* p) : A(p) {}

    //! Copy pointer (one more owner)
    MdagMLinOp(Handle< LinearOperator<T> > p) : A(p) {}

    //! Destructor
    ~MdagMLinOp() {}

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return A->subset();}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (T& chi, const T& psi, enum PlusMinus isign) const
      {
	T  tmp;  QDP::Hints::moveToFastMemoryHint(tmp);

	(*A)(tmp, psi, PLUS);
	(*A)(chi, tmp, MINUS);
      }

    unsigned long nFlops(void) const {
      unsigned long nflops=2*A->nFlops();
      return nflops;
    }

  private:
    Handle< LinearOperator<T> > A;
  };




  //! M^dag.M linear operator over arrays
  /*!
   * \ingroup linop
   *
   * Linear operator forming M^dag.M from an operator M
   */
  template<typename T>
  class MdagMLinOpArray : public LinearOperatorArray<T>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    MdagMLinOpArray(LinearOperatorArray<T>* p) : A(p) {}

    //! Copy pointer (one more owner)
    MdagMLinOpArray(Handle< LinearOperatorArray<T> > p) : A(p) {}

    //! Destructor
    ~MdagMLinOpArray() {}

    //! Length of array index
    int size() const {return A->size();}

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return A->subset();}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (multi1d<T>& chi, const multi1d<T>& psi, enum PlusMinus isign) const
      {
	multi1d<T>  tmp(size()); QDP::Hints::moveToFastMemoryHint(tmp);
	(*A)(tmp, psi, PLUS);
	(*A)(chi, tmp, MINUS);
      }

    unsigned long nFlops(void) const {
      unsigned long nflops=2*A->nFlops();
      return nflops;
    }

  private:
    Handle< LinearOperatorArray<T> > A;
  };


  //! M^dag.M linear operator
  /*!
   * \ingroup linop
   *
   * Linear operator forming M^dag.M from an operator M
   */
  template<typename T>
  class approx_lmdagm : public LinearOperator<T>
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

	(*A)(tmp, psi, PLUS);
	(*A)(chi, tmp, MINUS);
      }

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (T& chi, const T& psi, enum PlusMinus isign, Real epsilon) const
      {
	T  tmp; QDP::Hints::moveToFastMemoryHint(tmp);

	(*A)(tmp, psi, PLUS, epsilon/Real(2));
	(*A)(chi, tmp, MINUS, epsilon/Real(2));
      }

  private:
    Handle< LinearOperator<T> > A;
  };


  //! Differentiable M^dag.M linear operator
  /*!
   * \ingroup linop
   *
   * Diff. Linear operator forming M^dag.M from an operator M
   */
  template<typename T, typename P, typename Q>
  class DiffMdagMLinOp : public DiffLinearOperator<T,P,Q>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    DiffMdagMLinOp(DiffLinearOperator<T,P,Q>* p) : A(p) {}

    //! Copy pointer (one more owner)
    DiffMdagMLinOp(Handle< DiffLinearOperator<T,P,Q> > p) : A(p) {}

    //! Destructor
    ~DiffMdagMLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return A->getFermBC();}

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return A->subset();}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (T& chi, const T& psi, enum PlusMinus isign) const
      {
	T  tmp;  QDP::Hints::moveToFastMemoryHint(tmp);

	(*A)(tmp, psi, PLUS);
	(*A)(chi, tmp, MINUS);
      }

    //! Apply the derivative of the operator
    /*! 
     * Deriv of   chi^dag * A^dag * A * psi
     *
     * Default implementation
     */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
      {
	T   tmp; QDP::Hints::moveToFastMemoryHint(tmp);

	(*A)(tmp, psi, PLUS);
	A->deriv(ds_u, chi, tmp, MINUS);
      
	P  ds_tmp; // deriv routines should resize
	(*A)(tmp, chi, PLUS);  // note using PLUS
	A->deriv(ds_tmp, tmp, psi, PLUS);
	ds_u += ds_tmp;
      }


    //! Apply the derivative of the operator
    /*! 
     * Deriv of   chi^dag * A^dag * A * psi
     *
     * Default implementation
     */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign, const Real& epsilon) const
      {
	T   tmp;

	(*A)(tmp, psi, PLUS, epsilon);
	A->deriv(ds_u, chi, tmp, MINUS, epsilon);
      
	P  ds_tmp; // deriv routines should resize
	(*A)(tmp, chi, PLUS, epsilon);  // note using PLUS
	A->deriv(ds_tmp, tmp, psi, PLUS, epsilon);
	ds_u += ds_tmp;
      }

    //! Return the number of flops performed by operator()
    unsigned long nFlops(void) const {return 2*A->nFlops();}

  private:
    Handle< DiffLinearOperator<T,P,Q> > A;
  };




  //! M^dag.M linear operator over arrays
  /*!
   * \ingroup linop
   *
   * Linear operator forming M^dag.M from an operator M
   */
  template<typename T, typename P, typename Q>
  class DiffMdagMLinOpArray : public DiffLinearOperatorArray<T,P,Q>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    DiffMdagMLinOpArray(DiffLinearOperatorArray<T,P,Q>* p) : A(p) {}

    //! Copy pointer (one more owner)
    DiffMdagMLinOpArray(Handle< DiffLinearOperatorArray<T,P,Q> > p) : A(p) {}

    //! Destructor
    ~DiffMdagMLinOpArray() {}

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
	(*A)(tmp, psi, PLUS);
	(*A)(chi, tmp, MINUS);
      }

    //! Apply the derivative of the operator
    /*! 
     * Deriv of   chi^dag * A^dag * A * psi
     *
     * Default implementation
     */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign) const
      {
	multi1d<T>  tmp(size()); QDP::Hints::moveToFastMemoryHint(tmp);

	(*A)(tmp, psi, PLUS);
	A->deriv(ds_u, chi, tmp, MINUS);
      
	P  ds_tmp; // deriv routines should resize
	(*A)(tmp, chi, PLUS);  // note using PLUS
	A->deriv(ds_tmp, tmp, psi, PLUS);
	ds_u += ds_tmp;
      }


    //! Apply the derivative of the operator
    /*! 
     * Deriv of   chi^dag * A^dag * A * psi
     *
     * Default implementation
     */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign, const Real& epsilon) const
      {
	multi1d<T>  tmp(size());

	(*A)(tmp, psi, PLUS, epsilon);
	A->deriv(ds_u, chi, tmp, MINUS, epsilon);
      
	P  ds_tmp; // deriv routines should resize
	(*A)(tmp, chi, PLUS, epsilon);  // note using PLUS
	A->deriv(ds_tmp, tmp, psi, PLUS, epsilon);
	ds_u += ds_tmp;
      }

    //! Return the number of flops performed by operator()
    unsigned long nFlops(void) const {return 2*A->nFlops();}

  private:
    Handle< DiffLinearOperatorArray<T,P,Q> > A;
  };


} // End Namespace Chroma


#endif
