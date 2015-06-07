// -*- C++ -*-

/*! @file
 * @brief Periodic ferm state and a creator
 */

#ifndef __periodic_fermstate_w_h__
#define __periodic_fermstate_w_h__

#include "chromabase.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"

namespace Chroma
{
  /*! @ingroup fermstates */
  namespace CreatePeriodicFermStateEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }

}


#endif
