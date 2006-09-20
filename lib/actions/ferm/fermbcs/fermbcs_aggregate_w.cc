// $Id: fermbcs_aggregate_w.cc,v 3.2 2006-09-20 20:28:00 edwards Exp $
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
    static bool registered = false;

    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= WilsonTypeSimpleFermBCEnv::registerAll();
	success &= WilsonTypePeriodicFermBCEnv::registerAll();
	success &= WilsonTypeTwistedFermBCEnv::registerAll();
	success &= SchrTrivialFermBCEnv::registerAll();
	success &= SchrNonPertFermBCEnv::registerAll();
	success &= SchrCouplingFermBCEnv::registerAll();
	success &= SchrChromoMagFermBCEnv::registerAll();
	success &= SchrDirichletFermBCEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }

}
