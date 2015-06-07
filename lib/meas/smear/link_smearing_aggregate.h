// -*- C++ -*-
/*! \file
 *  \brief All link smearing constructors
 */

#ifndef __link_smearing_aggregate_w_h__
#define __link_smearing_aggregate_w_h__

#include "chromabase.h"
#include "io/xml_group_reader.h"

namespace Chroma
{
  //! Registration aggregator
  /*! \ingroup smear */
  namespace LinkSmearingEnv
  {
    bool registerAll();

    //! Returns a no-linksmearing group
    GroupXML_t   nullXMLGroup();
  }
}

#endif
