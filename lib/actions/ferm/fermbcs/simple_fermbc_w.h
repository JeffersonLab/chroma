// -*- C++ -*-
// $Id: simple_fermbc_w.h,v 2.2 2006-02-26 03:47:52 edwards Exp $
/*! \file
 *  \brief Simple fermionic BC
 */

#ifndef __simple_fermbc_w_h__
#define __simple_fermbc_w_h__

#include "actions/ferm/fermbcs/simple_fermbc.h"

namespace Chroma
{
  //! Name and registration 
  /*! \ingroup fermbc */
  namespace WilsonTypeSimpleFermBCEnv
  {
    extern const std::string name;
    extern const bool registered;
  }

}

#endif
