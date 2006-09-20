// -*- C++ -*-
// $Id: simple_fermbc_s.h,v 3.1 2006-09-20 20:28:00 edwards Exp $
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
