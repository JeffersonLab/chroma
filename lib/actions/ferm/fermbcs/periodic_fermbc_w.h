// -*- C++ -*-
// $Id: periodic_fermbc_w.h,v 3.1 2006-09-20 20:28:00 edwards Exp $
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
    bool registerAll();
  }

}

#endif
