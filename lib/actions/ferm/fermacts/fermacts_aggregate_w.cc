// $Id: fermacts_aggregate_w.cc,v 3.21 2009-02-10 04:22:42 edwards Exp $
/*! \file
 *  \brief All Wilson-type fermion actions
 */

#include "actions/ferm/fermbcs/fermbcs_aggregate_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

#include "actions/ferm/fermacts/unprec_clover_fermact_w.h"
#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/unprec_parwilson_fermact_w.h"
#include "actions/ferm/fermacts/unprec_graphene_fermact_w.h"
#include "actions/ferm/fermacts/unprec_hamberwu_fermact_w.h"
#include "actions/ferm/fermacts/unprec_dwftransf_fermact_w.h"
#include "actions/ferm/fermacts/unprec_w12_fermact_w.h"

#include "actions/ferm/fermacts/eoprec_clover_fermact_w.h"
#include "actions/ferm/fermacts/eoprec_clover_orbifold_fermact_w.h"
#include "actions/ferm/fermacts/eoprec_clover_extfield_fermact_w.h"
#include "actions/ferm/fermacts/eoprec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/eoprec_wilson_coarse_fine_fermact_w.h"
#include "actions/ferm/fermacts/eoprec_parwilson_fermact_w.h"
#include "actions/ferm/fermacts/eoprec_slic_fermact_w.h"
#include "actions/ferm/fermacts/eoprec_slrc_fermact_w.h"
#include "actions/ferm/fermacts/unprec_s_cprec_t_wilson_fermact_w.h"
#include "actions/ferm/fermacts/iluprec_s_cprec_t_wilson_fermact_w.h"
#include "actions/ferm/fermacts/iluprec_s_cprec_t_clover_fermact_w.h"
#include "actions/ferm/fermacts/ilu2prec_s_cprec_t_wilson_fermact_w.h"
#include "actions/ferm/fermacts/ilu2prec_s_cprec_t_clover_fermact_w.h"
#include "actions/ferm/fermacts/eo3dprec_s_cprec_t_wilson_fermact_w.h"
#include "actions/ferm/fermacts/eo3dprec_s_cprec_t_clover_fermact_w.h"

#include "actions/ferm/fermacts/ovlap_partfrac4d_fermact_w.h"

#include "actions/ferm/fermacts/poly_cheb_fermact_w.h"

#include "actions/ferm/fermacts/unprec_dwf_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_ovdwf_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_ovext_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_ovlap_contfrac5d_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_zolo_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_ht_contfrac5d_fermact_array_w.h"

#include "actions/ferm/fermacts/eoprec_dwf_fermact_array_w.h"
#include "actions/ferm/fermacts/eoprec_ovdwf_fermact_array_w.h"
#include "actions/ferm/fermacts/eoprec_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/eoprec_kno_fermact_array_w.h"
#include "actions/ferm/fermacts/eoprec_zolo_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/eoprec_ovlap_contfrac5d_fermact_array_w.h"
#include "actions/ferm/fermacts/eoprec_ht_contfrac5d_fermact_array_w.h"
#include "actions/ferm/fermacts/eoprec_ovext_fermact_array_w.h"


#include "actions/ferm/fermacts/ovext_tuning_strategy_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_accumulate_aggregate.h"
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
	success &= MdagMMultiSysSolverAccumulateEnv::registerAll();

	// All 4D bcs
	success &= WilsonTypeFermBCEnv::registerAll();

	// All fermstates
	success &= CreateFermStateEnv::registerAll();

	// 4D actions
	success &= EvenOddPrecWilsonFermActEnv::registerAll();
	success &= EvenOddPrecWilsonCoarseFineFermActEnv::registerAll();
	success &= UnprecGrapheneFermActEnv::registerAll();
	success &= UnprecWilsonFermActEnv::registerAll();
	success &= OvlapPartFrac4DFermActEnv::registerAll();
	success &= EvenOddPrecParWilsonFermActEnv::registerAll();
	success &= UnprecParWilsonFermActEnv::registerAll();

	success &= EvenOddPrecCloverFermActEnv::registerAll();
	success &= UnprecCloverFermActEnv::registerAll();
	success &= EvenOddPrecCloverOrbifoldFermActEnv::registerAll();
	success &= EvenOddPrecSLICFermActEnv::registerAll();
	success &= EvenOddPrecSLRCFermActEnv::registerAll();
	
//      success &= EvenOddPrecCloverExtFieldFermActEnv::registerAll();

	success &= UnprecHamberWuFermActEnv::registerAll();
//      success &= UnprecW12FermActEnv::registerAll();

	success &= PolyChebFermActEnv::registerAll();

#if QDP_NS == 4
#if QDP_NC == 3
#if QDP_ND == 4
	success &= UnprecSpaceCentralPrecTimeWilsonFermActEnv::registerAll();
	success &= ILUPrecSpaceCentralPrecTimeWilsonFermActEnv::registerAll();
	success &= ILUPrecSpaceCentralPrecTimeCloverFermActEnv::registerAll();
	success &= ILU2PrecSpaceCentralPrecTimeWilsonFermActEnv::registerAll();
	success &= ILU2PrecSpaceCentralPrecTimeCloverFermActEnv::registerAll();
	success &= EO3DPrecSpaceCentralPrecTimeWilsonFermActEnv::registerAll();
	success &= EO3DPrecSpaceCentralPrecTimeCloverFermActEnv::registerAll();
#endif
#endif
#endif
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
	success &= MdagMMultiSysSolverAccumulateArrayEnv::registerAll();

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
