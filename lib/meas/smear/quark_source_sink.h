// -*- C++ -*-
// $Id: quark_source_sink.h,v 2.1 2005-11-07 06:40:55 edwards Exp $

/*! @file
 * @brief Quark source or sink smearing
 */

#ifndef __quark_source_sink_h__
#define __quark_source_sink_h__

#include "chromabase.h"

namespace Chroma
{
  //! Base class for quark source and sink smearing
  /*! @ingroup smear
   *
   * Supports creation and application of smearing (with link smearing)
   * on quarks, and potentially displacements. Basically the construction
   * of a "source" or "sink" state on a pre-existing quark object
   */
  template<typename T>
  class QuarkSourceSink
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~QuarkSourceSink() {}

    //! Smear the quark
    /*!
     * \param obj      Object to source or sink smear ( Modify )
     */
    virtual void operator()(T& obj) const = 0;
  };

}


#endif
