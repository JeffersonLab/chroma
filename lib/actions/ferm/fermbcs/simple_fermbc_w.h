// -*- C++ -*-
// $Id: simple_fermbc_w.h,v 1.1 2004-12-24 04:23:20 edwards Exp $
/*! \file
 *  \brief Simple fermionic BC
 */

#ifndef __simple_fermbc_w_h__
#define __simple_fermbc_w_h__

#include "actions/ferm/fermbcs/simple_fermbc.h"

namespace Chroma
{
  //! Name and registration
  namespace WilsonTypeSimpleFermBCEnv
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Name and registration
  namespace WilsonTypeSimpleFermBCArrayEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
}

#endif
