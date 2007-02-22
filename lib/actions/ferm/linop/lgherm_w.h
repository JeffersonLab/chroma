// -*- C++ -*-
// $Id: lgherm_w.h,v 3.1 2007-02-22 21:11:46 bjoo Exp $

#ifndef __lgherm_h__
#define __lgherm_h__

#include "handle.h"
#include "linearop.h"


namespace Chroma 
{ 
  //! Gamma(5) hermitian linear operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * This operator scales its input operator
   */
  template<typename T>
  class lgherm : public LinearOperator<T>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    lgherm(LinearOperator<T>* p) : A(p) {}

    //! Copy pointer (one more owner)
    lgherm(Handle< LinearOperator<T> > p): A(p) {}

    //! Destructor
    ~lgherm() {}

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return A->subset();}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (T& chi, const T& psi, enum PlusMinus isign) const
      {
	const int G5=Ns*Ns-1;

	T  tmp;
	const Subset& sub = A->subset();

	// [ Gamma(5) D ]^{dag} = Gamma(5) D
	// using D = gamma_5 D^{dag} gamma_5

	(*A)(tmp, psi, PLUS);
	chi[sub] = Gamma(G5)*tmp;
      }

  private:
    Handle< LinearOperator<T> > A;
  };



  //! Partial specialization of scaled M  operator over arrays
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * Linear operator forming M^dag.M from an operator M
   */
  template<typename T>
  class lghermArray : public LinearOperatorArray<T>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    lghermArray(LinearOperatorArray<T>* p) : A(p) {}

    //! Copy pointer (one more owner)
    lghermArray(Handle< LinearOperatorArray<T> > p) : A(p) {}

    //! Destructor
    ~lghermArray() {}

    //! Length of array index
    int size() const {return A->size();}

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return A->subset();}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (multi1d<T>& chi, const multi1d<T>& psi, enum PlusMinus isign) const
      {
	const int G5=Ns*Ns-1;

	multi1d<T>  tmp(size());
	const Subset& sub = A->subset();

	// [ Gamma_5 D ]^{dag} = Gamma_5 D 
	// using D = gamma_5 D^{dag} gamma_5
	(*A)(tmp, psi, PLUS);
	for(int i = 0; i < size(); i++) { 
	  chi[i][sub] = Gamma(G5)*tmp[i];
	}
      }

  private:
    Handle< LinearOperatorArray<T> > A;
  };


} // End Namespace Chroma


#endif
