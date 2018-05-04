// -*- C++ -*-
/*! \file
 *  \brief All Wilson-type fermion actions
 */

#ifndef __fermactss_aggregate_w_h__
#define __fermactss_aggregate_w_h__

#include "chromabase.h"

namespace Chroma
{
  //! Registration aggregator
  /*! Wilson-like 4D */
  namespace WilsonTypeFermActs4DEnv
  {
    bool registerAll();
  }

  //! Registration aggregator
  /*! Wilson-like 5D actions */
  namespace WilsonTypeFermActs5DEnv
  {
    bool registerAll();
  }

  //! Registration aggregator
  /*! All Wilson-like 4D and 5D actions */
  namespace WilsonTypeFermActsEnv
  {
    bool registerAll();
  }

}

#endif
