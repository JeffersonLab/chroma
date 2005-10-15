// -*- C++ -*-
// $Id: fermacts_aggregate_w.h,v 2.1 2005-10-15 01:35:26 edwards Exp $
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
    extern const bool registered;
  }

  //! Registration aggregator
  /*! Wilson-like 5D actions */
  namespace WilsonTypeFermActs5DEnv
  {
    extern const bool registered;
  }

  //! Registration aggregator
  /*! All Wilson-like 4D and 5D actions */
  namespace WilsonTypeFermActsEnv
  {
    extern const bool registered;
  }

}

#endif
