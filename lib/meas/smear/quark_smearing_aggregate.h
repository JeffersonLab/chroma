// -*- C++ -*-
/*! \file
 *  \brief All quark smearing constructors
 */

#ifndef __quark_smearing_aggregate_w_h__
#define __quark_smearing_aggregate_w_h__

#include "chromabase.h"
#include "io/xml_group_reader.h"

namespace Chroma
{
  //! Registration aggregator
  /*! \ingroup smear */
  namespace QuarkSmearingEnv
  {
    bool registerAll();

    //! Returns a no-smearing group
    GroupXML_t   nullXMLGroup();
  }
}

#endif
