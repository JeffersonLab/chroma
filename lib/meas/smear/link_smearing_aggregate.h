// -*- C++ -*-
// $Id: link_smearing_aggregate.h,v 3.1 2006-06-10 16:28:19 edwards Exp $
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
    extern const bool registered;

    //! Returns a no-linksmearing group
    GroupXML_t   nullXMLGroup();
  }
}

#endif
