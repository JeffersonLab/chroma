// -*- C++ -*-
// $Id: ferm_createstate_aggregate_s.h,v 1.3 2006-11-17 02:17:31 edwards Exp $
/*! \file
 *  \brief All ferm create-state method
 */

#ifndef __ferm_createstate_aggregate_s_h__
#define __ferm_createstate_aggregate_s_h__

#include "chromabase.h"
#include "io/xml_group_reader.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup fermstates */
  namespace StaggeredCreateFermStateEnv
  {
    bool registerAll();

    //! Returns a periodic createstate group
    GroupXML_t   nullXMLGroup();
  }

}

#endif
