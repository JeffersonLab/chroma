// $Id: fermbcs_aggregate_w.cc,v 3.0 2006-04-03 04:58:48 edwards Exp $
/*! \file
 *  \brief All Wilson-type fermion boundary conditions
 */

#include "actions/ferm/fermbcs/fermbcs_aggregate_w.h"
#include "actions/ferm/fermbcs/simple_fermbc_w.h"
#include "actions/ferm/fermbcs/periodic_fermbc_w.h"
#include "actions/ferm/fermbcs/twisted_fermbc_w.h"
#include "actions/ferm/fermbcs/schr_triv_fermbc_w.h"
#include "actions/ferm/fermbcs/schr_nonpert_fermbc_w.h"
#include "actions/ferm/fermbcs/schr_coupling_fermbc_w.h"
#include "actions/ferm/fermbcs/schr_chromomag_fermbc_w.h"
#include "actions/ferm/fermbcs/schr_dirich_fermbc_w.h"

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
      success &= SchrTrivialFermBCEnv::registered;
      success &= SchrNonPertFermBCEnv::registered;
      success &= SchrCouplingFermBCEnv::registered;
      success &= SchrChromoMagFermBCEnv::registered;
      success &= SchrDirichletFermBCEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
