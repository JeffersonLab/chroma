// -*- C++ -*-
// $Id: periodic_fermbc_w.h,v 2.2 2006-03-16 03:00:12 edwards Exp $
/*! \file
 *  \brief Periodic fermionic BC
 */

#ifndef __periodic_fermbc_w_h__
#define __periodic_fermbc_w_h__

#include "actions/ferm/fermbcs/periodic_fermbc.h"

namespace Chroma
{
  //! Name and registration 
  /*! \ingroup fermbcs */
  namespace WilsonTypePeriodicFermBCEnv
  {
    extern const std::string name;
    extern const bool registered;
  }

}

#endif
