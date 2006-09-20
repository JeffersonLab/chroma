// $Id: fermacts_aggregate_w.cc,v 3.8 2006-09-20 20:27:58 edwards Exp $
/*! \file
 *  \brief All Wilson-type fermion actions
 */

#include "actions/ferm/fermbcs/fermbcs_aggregate_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

#include "actions/ferm/fermacts/unprec_clover_fermact_w.h"
#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/unprec_parwilson_fermact_w.h"
#include "actions/ferm/fermacts/unprec_hamberwu_fermact_w.h"
#include "actions/ferm/fermacts/unprec_dwftransf_fermact_w.h"
#include "actions/ferm/fermacts/unprec_w12_fermact_w.h"

#include "actions/ferm/fermacts/prec_clover_fermact_w.h"
#include "actions/ferm/fermacts/prec_clover_extfield_fermact_w.h"
#include "actions/ferm/fermacts/prec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/prec_parwilson_fermact_w.h"
#include "actions/ferm/fermacts/prec_slic_fermact_w.h"

#include "actions/ferm/fermacts/ovlap_partfrac4d_fermact_w.h"

#include "actions/ferm/fermacts/poly_cheb_fermact_w.h"

#include "actions/ferm/fermacts/unprec_dwf_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_ovdwf_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_ovext_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_ovlap_contfrac5d_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_zolo_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_ht_contfrac5d_fermact_array_w.h"

#include "actions/ferm/fermacts/prec_dwf_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_ovdwf_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_kno_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_zolo_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_ovlap_contfrac5d_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_ht_contfrac5d_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_ovext_fermact_array_w.h"


#include "actions/ferm/fermacts/ovext_tuning_strategy_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_aggregate.h"

#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"


namespace Chroma
{

  //! Registration aggregator
  namespace WilsonTypeFermActs4DEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// All system solvers
	success &= LinOpSysSolverEnv::registerAll();
	success &= MdagMSysSolverEnv::registerAll();
	success &= MdagMMultiSysSolverEnv::registerAll();

	// All 4D bcs
	success &= WilsonTypeFermBCEnv::registerAll();

	// All fermstates
	success &= CreateFermStateEnv::registerAll();

	// 4D actions
	success &= EvenOddPrecWilsonFermActEnv::registerAll();
	success &= UnprecWilsonFermActEnv::registerAll();
	success &= OvlapPartFrac4DFermActEnv::registerAll();
	success &= EvenOddPrecParWilsonFermActEnv::registerAll();
	success &= UnprecParWilsonFermActEnv::registerAll();

	success &= EvenOddPrecCloverFermActEnv::registerAll();
	success &= UnprecCloverFermActEnv::registerAll();
	success &= EvenOddPrecSLICFermActEnv::registerAll();

//      success &= EvenOddPrecCloverExtFieldFermActEnv::registerAll();

	success &= UnprecHamberWuFermActEnv::registerAll();
//      success &= UnprecW12FermActEnv::registerAll();

	success &= PolyChebFermActEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }


  //! Registration aggregator
  namespace WilsonTypeFermActs5DEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// All 5D system solvers
	success &= LinOpSysSolverArrayEnv::registerAll();
	success &= MdagMSysSolverArrayEnv::registerAll();
	success &= MdagMMultiSysSolverArrayEnv::registerAll();

	// All 5D bcs
	success &= WilsonTypeFermBCEnv::registerAll();

	// All fermstates
	success &= CreateFermStateEnv::registerAll();

	// 5D actions
	success &= EvenOddPrecDWFermActArrayEnv::registerAll();
	success &= UnprecDWFermActArrayEnv::registerAll();
	success &= EvenOddPrecNEFFermActArrayEnv::registerAll();
	success &= UnprecNEFFermActArrayEnv::registerAll();
	success &= UnprecOvlapContFrac5DFermActArrayEnv::registerAll();
	success &= UnprecHTContFrac5DFermActArrayEnv::registerAll();
	success &= EvenOddPrecHtContFrac5DFermActArrayEnv::registerAll();
	success &= EvenOddPrecOvlapContFrac5DFermActArrayEnv::registerAll();
	success &= UnprecOvDWFermActArrayEnv::registerAll();
	success &= EvenOddPrecOvDWFermActArrayEnv::registerAll();
	success &= UnprecOvExtFermActArrayEnv::registerAll();
	success &= EvenOddPrecOvExtFermActArrayEnv::registerAll();
	success &= UnprecZoloNEFFermActArrayEnv::registerAll();
	success &= EvenOddPrecZoloNEFFermActArrayEnv::registerAll();
	success &= EvenOddPrecKNOFermActArrayEnv::registerAll();
	success &= UnprecDWFTransfFermActEnv::registerAll();
    
	// Tuning Strategies
	success &= OvExtTuningStrategyAggregateEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }


  //! Registration aggregator
  /*! All Wilson-like 4D and 5D actions */
  namespace WilsonTypeFermActsEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// All 4D actions
	success &= WilsonTypeFermActs4DEnv::registerAll();

	// 5D actions
	success &= WilsonTypeFermActs5DEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
