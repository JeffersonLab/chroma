// -*- C++ -*-
// $Id: quark_displacement_aggregate.h,v 3.1 2006-06-10 16:28:19 edwards Exp $
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
    extern const bool registered;

    //! Returns a no-displacement group
    GroupXML_t   nullXMLGroup();
  }
}

#endif
