// -*- C++ -*-
// $Id: quark_smearing_aggregate.h,v 2.1 2005-11-07 06:40:55 edwards Exp $
/*! \file
 *  \brief All quark smearing constructors
 */

#ifndef __quark_smearing_aggregate_w_h__
#define __quark_smearing_aggregate_w_h__

#include "chromabase.h"

namespace Chroma
{
  //! Registration aggregator
  /*! \ingroup smear */
  namespace PropSmearingEnv
  {
    extern const bool registered;
  }

  //! Registration aggregator
  /*! \ingroup smear */
  namespace FermSmearingEnv
  {
    extern const bool registered;
  }

  //! Registration aggregator
  /*! \ingroup smear */
  namespace ColorVecSmearingEnv
  {
    extern const bool registered;
  }
}

#endif
