// -*- C++ -*-
// $Id: prec_dwflike_linop_base_array_w.h,v 1.1 2005-05-28 22:37:42 edwards Exp $
/*! \file
 *  \brief Base class for even-odd preconditioned domain-wall-like linops
 */

#ifndef __prec_dwflike_linop_base_array_w_h__
#define __prec_dwflike_linop_base_array_w_h__

#include "linearop.h"


namespace Chroma
{
  //! 4D Even Odd preconditioned domain-wall Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   */
  template<typename T, typename P>
  class EvenOddPrecDWLikeLinOpBaseArray : public EvenOddPrecLinearOperator< multi1d<T>, P >
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

}; // End Namespace Chroma


#endif
