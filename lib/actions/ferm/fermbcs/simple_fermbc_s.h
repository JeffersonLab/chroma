// -*- C++ -*-
// $Id: simple_fermbc_s.h,v 1.1 2004-12-24 04:23:20 edwards Exp $
/*! \file
 *  \brief Simple fermionic BC
 */

#ifndef __simple_fermbc_s_h__
#define __simple_fermbc_s_h__

#include "actions/ferm/fermbcs/simple_fermbc.h"

namespace Chroma
{
  //! Name and registration
  namespace StaggeredTypeSimpleFermBCEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
}

#endif
