// -*- C++ -*-
#ifndef __lopishift_h__
#define __lopishift_h__

/*! \file A linear operator scaled by i sigma where sigma is a Real shift */

#include "handle.h"
#include "linearop.h"

namespace Chroma {

  template<typename T, typename C>
  class lopishift : public LinearOperator<T>
  {
  public:

    //! Initialized from pointer
    lopishift(LinearOperator<T>* p, const C& s) : M(p), shift_fact(s) {}
    
    lopishift(Handle<LinearOperator<T> > p, const C& s) : M(p), shift_fact(s) {}

    ~lopishift(){};

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return M->subset();}
    
    virtual unsigned long nFlops() const {
      return M->nFlops() + 4*Nc*Ns*(M->subset()).siteTable().size();
    }

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    inline void operator() (T& chi, const T& psi, enum PlusMinus isign) const
      {
	const Subset& sub = M->subset();
	(*M)(chi, psi, isign);
	T tmp;
	tmp[sub] = Gamma(15)*timesI(psi);
	if( isign == PLUS) {
	  chi[sub] += shift_fact*tmp; // (M + i shift) psi
	}
	else {
	  chi[sub] -= shift_fact*tmp; // (M - i shift) psi
	}
      }

  private:
    Handle< LinearOperator<T> > M;
    const C& shift_fact;
  };
}
#endif
