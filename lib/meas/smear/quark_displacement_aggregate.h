// -*- C++ -*-
/*! \file
 *  \brief All quark displacement constructors
 */

#ifndef __quark_displacement_aggregate_w_h__
#define __quark_displacement_aggregate_w_h__

#include "chromabase.h"
#include "io/xml_group_reader.h"

namespace Chroma
{
  //! Registration aggregator
  /*! \ingroup smear */
  namespace QuarkDisplacementEnv
  {
    bool registerAll();

    //! Returns a no-displacement group
    GroupXML_t   nullXMLGroup();
  }
}

#endif
