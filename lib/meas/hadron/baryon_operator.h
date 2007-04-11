// -*- C++ -*-
// $Id: baryon_operator.h,v 1.3 2007-04-11 04:55:29 juge Exp $
/*! \file
 *  \brief Construct baryon operator
 */

#ifndef __baryon_operator_h__
#define __baryon_operator_h__

#include "chromabase.h"

namespace Chroma
{
  //! Construct baryon operators
  /*! @ingroup hadron
   *
   */
  template<typename T>
  class BaryonOperator
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~BaryonOperator() {}

    //! Construct the operator (do the contractions)
    virtual multi1d<LatticeComplex> operator()(const T& quark1, 
					       const T& quark2, 
					       const T& quark3,
					       enum PlusMinus isign) const = 0;

    virtual LatticeComplex operator()(const T& quark1, 
				      const T& quark2, 
				      const T& quark3) const = 0;

    virtual multi1d<LatticeComplex> operator()(const T& quark1, 
					       const T& quark2, 
					       const T& quark3,
								 int* qindices,
					       enum PlusMinus isign) const = 0;

    virtual LatticeComplex operator()(const T& quark1, 
					       const T& quark2, 
					       const T& quark3,
					       int* qindices) const = 0;

  };

}


#endif
