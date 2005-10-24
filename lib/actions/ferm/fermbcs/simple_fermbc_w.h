// -*- C++ -*-
// $Id: simple_fermbc_w.h,v 2.1 2005-10-24 05:53:55 edwards Exp $
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


  //! Name and registration
  /*! \ingroup fermbc */
  namespace WilsonTypeSimpleFermBCArrayEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
}

#endif
