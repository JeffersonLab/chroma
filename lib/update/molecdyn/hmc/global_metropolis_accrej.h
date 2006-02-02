// -*- C++ -*-
// $Id: global_metropolis_accrej.h,v 2.1 2006-02-02 18:46:02 edwards Exp $
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
