// -*- C++ -*-
// $Id: simple_fermbc_w.h,v 2.3 2006-03-16 03:00:13 edwards Exp $
/*! \file
 *  \brief Simple fermionic BC
 */

#ifndef __simple_fermbc_w_h__
#define __simple_fermbc_w_h__

#include "actions/ferm/fermbcs/simple_fermbc.h"

namespace Chroma
{
  //! Name and registration 
  /*! \ingroup fermbcs */
  namespace WilsonTypeSimpleFermBCEnv
  {
    extern const std::string name;
    extern const bool registered;
  }

}

#endif
