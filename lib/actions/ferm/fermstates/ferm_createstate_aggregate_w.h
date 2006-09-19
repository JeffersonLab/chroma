// -*- C++ -*-
// $Id: ferm_createstate_aggregate_w.h,v 1.1 2006-09-19 17:53:36 edwards Exp $
/*! \file
 *  \brief All ferm create-state method
 */

#ifndef __ferm_createstate_aggregate_w_h__
#define __ferm_createstate_aggregate_w_h__

#include "chromabase.h"
#include "io/xml_group_reader.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup fermstates */
  namespace CreateFermStateEnv
  {
    extern const bool registered;

    //! Returns a periodic createstate group
    GroupXML_t   nullXMLGroup();
  }

}

#endif
