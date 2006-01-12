// -*- C++ -*-
// $Id: lmdagm.h,v 2.1 2006-01-12 05:45:16 edwards Exp $

#ifndef __lmdagm_h__
#define __lmdagm_h__

#include "handle.h"
#include "linearop.h"

using namespace QDP::Hints;

namespace Chroma 
{ 
  //! M^dag.M linear operator
  /*!
   * \ingroup linop
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

    unsigned long nFlops(void) const {
      unsigned long nflops=2*A->nFlops();
      return nflops;
    }

  private:
    const Handle< const LinearOperator<T> > A;
  };




  //! Partial specialization of M^dag.M linear operator over arrays
  /*!
   * \ingroup linop
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

    unsigned long nFlops(void) const {
      unsigned long nflops=2*A->nFlops();
      return nflops;
    }

  private:
    const Handle< const LinearOperator< multi1d<T> > > A;
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


  //! Differentiable M^dag.M linear operator
  /*!
   * \ingroup linop
   *
   * Diff. Linear operator forming M^dag.M from an operator M
   */
  template<typename T, typename P>
  class DiffMdagMLinOp : public DiffLinearOperator<T,P>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    DiffMdagMLinOp(const DiffLinearOperator<T,P>* p) : A(p) {}

    //! Copy pointer (one more owner)
    DiffMdagMLinOp(Handle<const DiffLinearOperator<T,P> > p) : A(p) {}

    //! Destructor
    ~DiffMdagMLinOp() {}

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

    //! Apply the derivative of the operator
    /*! 
     * Deriv of   chi^dag * A^dag * A * psi
     *
     * Default implementation
     */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
      {
	T   tmp; moveToFastMemoryHint(tmp);

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
    const Handle< const DiffLinearOperator<T,P> > A;
  };




  //! Partial specialization of M^dag.M linear operator over arrays
  /*!
   * \ingroup linop
   *
   * Linear operator forming M^dag.M from an operator M
   */
  template<typename T, typename P>
  class DiffMdagMLinOp< multi1d<T>, P > : public DiffLinearOperator< multi1d<T>, P >
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    DiffMdagMLinOp(const DiffLinearOperator< multi1d<T>, P >* p) : A(p) {}

    //! Copy pointer (one more owner)
    DiffMdagMLinOp(Handle<const DiffLinearOperator< multi1d<T>, P > > p) : A(p) {}

    //! Destructor
    ~DiffMdagMLinOp() {}

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

    //! Apply the derivative of the operator
    /*! 
     * Deriv of   chi^dag * A^dag * A * psi
     *
     * Default implementation
     */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign) const
      {
	multi1d<T>  tmp(size()); moveToFastMemoryHint(tmp);

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
    const Handle< const DiffLinearOperator< multi1d<T>, P > > A;
  };


} // End Namespace Chroma


#endif
