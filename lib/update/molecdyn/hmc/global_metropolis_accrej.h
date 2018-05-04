// -*- C++ -*-
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
  // Metropolis accept/reject
  /*! @ingroup hmc */
  bool globalMetropolisAcceptReject(const Double& DeltaH);
}
#endif
