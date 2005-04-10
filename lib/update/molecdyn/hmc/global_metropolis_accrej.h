// -*- C++ -*-
// $Id: global_metropolis_accrej.h,v 1.3 2005-04-10 21:46:42 edwards Exp $
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
