// -*- C++ -*-
// $Id: unprec_dwf_linop_base_array_w.h,v 1.6 2005-01-14 20:13:06 edwards Exp $
/*! \file
 *  \brief Base class for unpreconditioned domain-wall-like fermion linear operator
 */

#ifndef __unprec_dwf_linop_base_array_w_h__
#define __unprec_dwf_linop_base_array_w_h__

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
  class UnprecDWLinOpBaseArray : public UnprecLinearOperator< multi1d<T>, P >
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
