// -*- C++ -*-
// $Id: quark_smearing.h,v 3.0 2006-04-03 04:59:05 edwards Exp $

/*! @file
 * @brief Quark smearing
 */

#ifndef __quark_smearing_h__
#define __quark_smearing_h__

#include "chromabase.h"

namespace Chroma
{
  //! Base class for quark smearing
  /*! @ingroup smear
   *
   * Supports smearing of quarks
   */
  template<typename T>
  class QuarkSmearing
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~QuarkSmearing() {}

    //! Smear the quark
    /*!
     * \param obj      Object to smear ( Modify )
     * \param u        Link field ( Read )
     */
    virtual void operator()(T& obj, const multi1d<LatticeColorMatrix>& u) const = 0;
  };

}


#endif
