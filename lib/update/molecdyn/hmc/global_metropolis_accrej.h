// -*- C++ -*-
// $Id: global_metropolis_accrej.h,v 3.1 2006-08-26 02:08:41 edwards Exp $
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
