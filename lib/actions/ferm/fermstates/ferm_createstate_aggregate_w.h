// -*- C++ -*-
// $Id: ferm_createstate_aggregate_w.h,v 1.2 2006-09-20 20:31:41 edwards Exp $
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
    bool registerAll();

    //! Returns a periodic createstate group
    GroupXML_t   nullXMLGroup();
  }

}

#endif
