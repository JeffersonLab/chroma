// -*- C++ -*-
// $Id: periodic_fermstate_s.h,v 1.2 2006-09-20 20:31:41 edwards Exp $

/*! @file
 * @brief Periodic ferm state and a creator
 */

#ifndef __periodic_fermstate_s_h__
#define __periodic_fermstate_s_h__

#include "chromabase.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"

namespace Chroma
{
  /*! @ingroup fermstates */
  namespace StaggeredCreatePeriodicFermStateEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }

}


#endif
