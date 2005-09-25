// -*- C++ -*-
// $Id: global_metropolis_accrej.h,v 2.0 2005-09-25 21:04:41 edwards Exp $
/*! \file
 * \brief Global metropolis
 *
 * Global metropolis
 */

#ifndef GLOBAL_METROP_ACCREJ
#define GLOBAL_METROP_ACCREJ

#include "chromabase.h"


namespace Chroma 
{ 
  // This needs to move to a .cc file at some stage
  /*! @ingroup hmc */
  bool globalMetropolisAcceptReject(const Double& DeltaH);
};
#endif
