// -*- C++ -*-
// $Id: simple_fermbc_w.h,v 2.0 2005-09-25 21:04:27 edwards Exp $
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
