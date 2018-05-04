// -*- C++ -*-
/*! @file
 * @brief Quark displacement
 */

#ifndef __quark_displacement_h__
#define __quark_displacement_h__

#include "chromabase.h"

namespace Chroma
{
  //! Base class for quark displacement
  /*! @ingroup smear
   *
   * Supports displacement of quarks
   */
  template<typename T>
  class QuarkDisplacement
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~QuarkDisplacement() {}

    //! Displace the quark
    /*!
     * \param obj      Object to displace ( Modify )
     * \param u        Link field ( Read )
     * \param isign    PLUS is the operator and MINUS is the dagger ( Read )
     */
    virtual void operator()(T& obj, 
			    const multi1d<LatticeColorMatrix>& u, 
			    enum PlusMinus isign) const = 0;
  };

}


#endif
