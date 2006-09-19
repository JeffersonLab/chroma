// -*- C++ -*-
// $Id: periodic_fermstate_s.h,v 1.1 2006-09-19 17:53:37 edwards Exp $

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
    extern const bool registered;
  }

}


#endif
