// -*- C++ -*-
// $Id: unprec_dwflike_linop_base_array_w.h,v 2.0 2005-09-25 21:04:30 edwards Exp $
/*! \file
 *  \brief Base class for unpreconditioned domain-wall-like fermion linear operator
 */

#ifndef __unprec_dwflike_linop_base_array_w_h__
#define __unprec_dwflike_linop_base_array_w_h__

#include "linearop.h"


namespace Chroma
{
  //! Unpreconditioned domain-wall Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   */
  template<typename T, typename P>
  class UnprecDWLikeLinOpBaseArray : public UnprecLinearOperator< multi1d<T>, P >
  {
  public:
    //! Length of DW flavor index/space
    virtual int size() const = 0;

#if 0
    //! Apply the Dplus operator on a lattice fermion.
    virtual
    void Dplus(T& chi,
	       const T& psi,
	       enum PlusMinus isign,
	       int s5) const = 0;
#endif

    //! Apply the Dminus operator on a lattice fermion.
    virtual
    void Dminus(T& chi,
		const T& psi,
		enum PlusMinus isign,
		int s5) const = 0;
  };
 
}


#endif
