// -*- C++ -*-
// $Id: lgherm_w.h,v 1.3 2004-12-12 21:22:16 edwards Exp $

#ifndef __lgherm_h__
#define __lgherm_h__

#include "handle.h"
#include "linearop.h"

using namespace QDP;

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
  lgherm(const LinearOperator<T>* p) : A(p) {}

  //! Copy pointer (one more owner)
  lgherm(Handle<const LinearOperator<T> > p): A(p) {}

  //! Destructor
  ~lgherm() {}

  //! Subset comes from underlying operator
  inline const OrderedSubset& subset() const {return A->subset();}

  //! Apply the operator onto a source vector
  /*! For this operator, the sign is ignored */
  inline void operator() (T& chi, const T& psi, enum PlusMinus isign) const
  {
      const int G5=Ns*Ns-1;

      T  tmp;
      const OrderedSubset& sub = A->subset();

      // [ Gamma(5) D ]^{dag} = Gamma(5) D
      // using D = gamma_5 D^{dag} gamma_5

      (*A)(tmp, psi, PLUS);
      chi[sub] = Gamma(G5)*tmp;
  }

private:
  const Handle< const LinearOperator<T> > A;
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
class lgherm< multi1d<T> > : public LinearOperator< multi1d<T> >
{
public:
  //! Initialize pointer with existing pointer
  /*! Requires that the pointer p is a return value of new */
  lgherm(const LinearOperator< multi1d<T> >* p) : A(p) {}

  //! Copy pointer (one more owner)
  lgherm(Handle<const LinearOperator< multi1d<T> > > p) : A(p) {}

  //! Destructor
  ~lgherm() {}

  //! Length of array index
  int size() const {return A->size();}

  //! Subset comes from underlying operator
  inline const OrderedSubset& subset() const {return A->subset();}

  //! Apply the operator onto a source vector
  /*! For this operator, the sign is ignored */
  inline void operator() (multi1d<T>& chi, const multi1d<T>& psi, enum PlusMinus isign) const
    {
      const int G5=Ns*Ns-1;

      multi1d<T>  tmp(size());
      const OrderedSubset& sub = A->subset();

      // [ Gamma_5 D ]^{dag} = Gamma_5 D 
      // using D = gamma_5 D^{dag} gamma_5
      (*A)(tmp, psi, PLUS);
      for(int i = 0; i < size(); i++) { 
	chi[i][sub] = Gamma(G5)*tmp[i];
      }
    }

private:
  const Handle< const LinearOperator< multi1d<T> > > A;
};


}; // End Namespace Chroma

using namespace Chroma;

#endif
