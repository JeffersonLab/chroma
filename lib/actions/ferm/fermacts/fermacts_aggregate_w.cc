// $Id: fermacts_aggregate_w.cc,v 2.1 2005-10-02 03:08:49 bjoo Exp $
/*! \file
 *  \brief All Wilson-type fermion actions
 */

#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/unprec_parwilson_fermact_w.h"
#include "actions/ferm/fermacts/unprec_dwftransf_fermact_w.h"
#include "actions/ferm/fermacts/prec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/prec_parwilson_fermact_w.h"

#include "actions/ferm/fermacts/ovlap_partfrac4d_fermact_w.h"

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

#include "actions/ferm/fermacts/unprec_stout_fermact_w.h"

namespace Chroma
{

  //! Registration aggregator
  namespace WilsonTypeFermActs4DEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // 4D actions
      success &= EvenOddPrecWilsonFermActEnv::registered;
      success &= UnprecWilsonFermActEnv::registered;
      success &= OvlapPartFrac4DFermActEnv::registered;
      success &= EvenOddPrecParWilsonFermActEnv::registered;
      success &= UnprecParWilsonFermActEnv::registered;
      success &= UnprecStoutWilsonTypeFermActEnv::registered;
      return success;
    }

    const bool registered = registerAll();
  }


  //! Registration aggregator
  namespace WilsonTypeFermActs5DEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // 5D actions
      success &= EvenOddPrecDWFermActArrayEnv::registered;
      success &= UnprecDWFermActArrayEnv::registered;
      success &= EvenOddPrecNEFFermActArrayEnv::registered;
      success &= UnprecNEFFermActArrayEnv::registered;
      success &= UnprecOvlapContFrac5DFermActArrayEnv::registered;
      success &= UnprecHTContFrac5DFermActArrayEnv::registered;
      success &= EvenOddPrecHtContFrac5DFermActArrayEnv::registered;
      success &= EvenOddPrecOvlapContFrac5DFermActArrayEnv::registered;
      success &= UnprecOvDWFermActArrayEnv::registered;
      success &= EvenOddPrecOvDWFermActArrayEnv::registered;
      success &= UnprecOvExtFermActArrayEnv::registered;
      success &= EvenOddPrecOvExtFermActArrayEnv::registered;
      success &= UnprecZoloNEFFermActArrayEnv::registered;
      success &= EvenOddPrecZoloNEFFermActArrayEnv::registered;
      success &= EvenOddPrecKNOFermActArrayEnv::registered;
      success &= UnprecDWFTransfFermActEnv::registered;

      // Tuning Strategies
      success &= OvExtTuningStrategyAggregateEnv::registered;
      return success;
    }

    const bool registered = registerAll();
  }


  //! Registration aggregator
  /*! All Wilson-like 4D and 5D actions */
  namespace WilsonTypeFermActsEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // All 4D actions
      success &= WilsonTypeFermActs4DEnv::registered;

      // 5D actions
      success &= WilsonTypeFermActs5DEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
