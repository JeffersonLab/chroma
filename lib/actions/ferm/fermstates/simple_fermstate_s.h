// -*- C++ -*-
// $Id: simple_fermstate_s.h,v 1.2 2006-09-20 20:31:41 edwards Exp $

/*! @file
 * @brief Simple ferm state and a creator
 */

#ifndef __simple_fermstate_s_h__
#define __simple_fermstate_s_h__

#include "chromabase.h"

namespace Chroma
{
  /*! @ingroup fermstates */
  namespace StaggeredCreateSimpleFermStateEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }

}


#endif
