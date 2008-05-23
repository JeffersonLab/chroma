// -*- C++ -*-
// $Id: rat_approx.h,v 3.1 2008-05-23 21:31:34 edwards Exp $
/*! \file
 *  \brief Base class for rational approximations
 */

#ifndef __rat_approx_h__
#define __rat_approx_h__

#include "update/molecdyn/monomial/remez_coeff.h"

namespace Chroma
{

  //! Base class for rational approximations
  /*! @ingroup monomial
   *
   */
  class RationalApprox
  {
  public:
    //! Virtual destructor
    virtual ~RationalApprox() {}

    //! Produce the partial-fraction-expansion (PFE) and its inverse (IPFE)
    virtual void operator()(RemezCoeff_t& pfe, RemezCoeff_t& ipfe) const = 0;
  };

}  // end namespace Chroma

#endif  // Include guard



