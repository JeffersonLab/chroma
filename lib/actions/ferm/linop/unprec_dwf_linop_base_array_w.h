// -*- C++ -*-
// $Id: unprec_dwf_linop_base_array_w.h,v 1.2 2004-09-19 02:37:06 edwards Exp $
/*! \file
 *  \brief Base class for unpreconditioned domain-wall-like fermion linear operator
 */

#ifndef __unprec_dwf_linop_base_array_w_h__
#define __unprec_dwf_linop_base_array_w_h__

#include "linearop.h"

using namespace QDP;

//namespace Chroma
//{
  //! Unpreconditioned domain-wall Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   */
  template<typename T>
  class UnprecDWLinOpBaseArray : public LinearOperator< multi1d<T> >
  {
  public:
    //! Length of DW flavor index/space
    virtual int size() const = 0;

#if 0
    //! Apply the Dplus operator on a vector in Ls.
    virtual
    void Dplus(multi1d<T>& chi,
	       const multi1d<T>& psi,
	       enum PlusMinus isign) const = 0;

    //! Apply the Dplus operator on a lattice fermion.
    virtual
    void Dplus(T& chi,
	       const T& psi,
	       enum PlusMinus isign) const = 0;
#endif

    //! Apply the Dminus operator on a vector in Ls.
    virtual
    void Dminus(multi1d<T>& chi,
		const multi1d<T>& psi,
		enum PlusMinus isign) const = 0;

    //! Apply the Dminus operator on a lattice fermion.
    virtual
    void Dminus(T& chi,
		const T& psi,
		enum PlusMinus isign) const = 0;
  };
 
//}

//using namespace Chroma;

#endif
