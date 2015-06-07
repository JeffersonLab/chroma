// -*- C++ -*-
/*! \file
 *  \brief Simple fermionic BC
 */

#ifndef __simple_fermbc_s_h__
#define __simple_fermbc_s_h__

#include "actions/ferm/fermbcs/simple_fermbc.h"

namespace Chroma
{
  //! Name and registration
  /*! \ingroup fermbcs */
  namespace StaggeredTypeSimpleFermBCEnv
  {
    extern const std::string name;
    bool registerAll();
  }
}

#endif
