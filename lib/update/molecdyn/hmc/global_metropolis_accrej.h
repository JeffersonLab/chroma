// -*- C++ -*-
// $Id: global_metropolis_accrej.h,v 3.0 2006-04-03 04:59:07 edwards Exp $
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
};
#endif
