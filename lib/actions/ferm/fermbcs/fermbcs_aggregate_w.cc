// $Id: fermbcs_aggregate_w.cc,v 2.1 2006-02-26 03:47:51 edwards Exp $
/*! \file
 *  \brief All Wilson-type fermion boundary conditions
 */

#include "actions/ferm/fermbcs/fermbcs_aggregate_w.h"
#include "actions/ferm/fermbcs/simple_fermbc_w.h"
#include "actions/ferm/fermbcs/periodic_fermbc_w.h"
#include "actions/ferm/fermbcs/twisted_fermbc_w.h"
#include "actions/ferm/fermbcs/schr1link_fermbc_w.h"
#include "actions/ferm/fermbcs/schr2link_fermbc_w.h"

namespace Chroma
{

  //! Registration aggregator
  namespace WilsonTypeFermBCEnv
  {
    bool registerAll() 
    {
      bool success = true;

      success &= WilsonTypeSimpleFermBCEnv::registered;
      success &= WilsonTypePeriodicFermBCEnv::registered;
      success &= WilsonTypeTwistedFermBCEnv::registered;
//      success &= WilsonTypeSchr1LinkFermBCEnv::registered;
//      success &= WilsonTypeSchr2LinkFermBCEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
