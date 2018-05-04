// -*- C++ -*-
/*! @file
 * @brief Gauge initialization
 */

#ifndef __gauge_init_h__
#define __gauge_init_h__

#include "chromabase.h"

namespace Chroma
{
  //! Base class for gauge initialization
  /*! @ingroup gauge
   *
   * Supports initialization of gauge fields
   */
  class GaugeInit
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~GaugeInit() {}

    //! Initialize the gauge field
    virtual void operator()(XMLReader& file_xml, XMLReader& record_xml, 
			    multi1d<LatticeColorMatrix>& u) const = 0;
  };

}


#endif
