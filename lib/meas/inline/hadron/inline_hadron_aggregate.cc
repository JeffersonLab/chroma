// $Id: inline_hadron_aggregate.cc,v 3.58 2009-09-16 21:23:09 colin Exp $
/*! \file
 *  \brief Inline hadron measurement aggregator
 */

#include "chroma_config.h"

#include "meas/inline/hadron/inline_hadron_aggregate.h"
#include "meas/inline/hadron/inline_apply_fermstate_w.h"
//#include "meas/inline/hadron/inline_spectrumQll.h"
#include "meas/inline/hadron/inline_create_colorvecs.h"
#include "meas/inline/hadron/inline_create_colorvecs.h"

#ifdef BUILD_LAPACK
#include "meas/inline/hadron/inline_laplace_eigs.h"
#else
#warning "Not Building Inline Laplace Eigs"
#endif

#include "meas/inline/hadron/inline_prop_3pt_w.h"
#include "meas/inline/hadron/inline_disco_w.h"
#include "meas/inline/hadron/inline_disco_eoprec_w.h"
#include "meas/inline/hadron/inline_disco_eo_eigcg_w.h"
#include "meas/inline/hadron/inline_disco_eigcg_w.h"
#include "meas/inline/hadron/inline_static_light_spec_w.h"
#include "meas/inline/hadron/inline_heavy_light_cont_w.h"
#include "meas/inline/hadron/inline_heavyhadspec_w.h"
#include "meas/inline/hadron/inline_make_source_w.h"
#include "meas/inline/hadron/inline_make_source_ferm_w.h"
#include "meas/inline/hadron/inline_propagator_w.h"
#include "meas/inline/hadron/inline_propagator_ferm_w.h"
#include "meas/inline/hadron/inline_multi_propagator_w.h"
#include "meas/inline/hadron/inline_seqsource_w.h"
#include "meas/inline/hadron/inline_seqprop_test_w.h"
#include "meas/inline/hadron/inline_hadspec_w.h"
#include "meas/inline/hadron/inline_mesonspec_w.h"
#include "meas/inline/hadron/inline_hadron_contract.h"
#include "meas/inline/hadron/inline_stag_to_wils.h"
//#include "meas/inline/hadron/inline_spectrum_w.h"
#include "meas/inline/hadron/inline_sink_smear_w.h"
#include "meas/inline/hadron/inline_diquark_w.h"
#include "meas/inline/hadron/inline_qqq_w.h"
#include "meas/inline/hadron/inline_qqq_diquark_w.h"
#include "meas/inline/hadron/inline_qqbar_w.h"
#include "meas/inline/hadron/inline_building_blocks_w.h"
#include "meas/inline/hadron/inline_noisy_building_blocks_w.h"
#include "meas/inline/hadron/inline_bar3ptfn_w.h"
//#include "meas/inline/hadron/inline_multipole_w.h"
#include "meas/inline/hadron/inline_npr_vertex_w.h"
#include "meas/inline/hadron/inline_npr_w.h"
#include "meas/inline/hadron/inline_distillution_noise.h"
#include "meas/inline/hadron/inline_prop_distillution_w.h"
#include "meas/inline/hadron/inline_prop_distillation_w.h"
#include "meas/inline/hadron/inline_prop_colorvec_w.h"
#include "meas/inline/hadron/inline_static_prop_colorvec_w.h"
#include "meas/inline/hadron/inline_annih_prop_matelem_colorvec_w.h"
#include "meas/inline/hadron/inline_prop_matelem_colorvec_w.h"
#include "meas/inline/hadron/inline_prop_matelem_pt_colorvec_w.h"
#include "meas/inline/hadron/inline_prop_and_matelem_colorvec_w.h"
#include "meas/inline/hadron/inline_prop_matelem_lm_colorvec_w.h"
#include "meas/inline/hadron/inline_baryon_matelem_colorvec_w.h"
#include "meas/inline/hadron/inline_meson_matelem_colorvec_w.h"
#include "meas/inline/hadron/inline_genprop_matelem_colorvec_w.h"
#include "meas/inline/hadron/inline_genprop_matelem_pt_colorvec_w.h"
#include "meas/inline/hadron/inline_mres_w.h"
#include "meas/inline/hadron/inline_qpropqio_w.h"
#include "meas/inline/hadron/inline_qpropadd_w.h"
#include "meas/inline/hadron/inline_qqqNucNuc_w.h"
#include "meas/inline/hadron/inline_stoch_meson_w.h"
#include "meas/inline/hadron/inline_stoch_baryon_w.h"
#include "meas/inline/hadron/inline_stoch_hadron_w.h"
#include "meas/inline/hadron/inline_stoch_group_baryon_w.h"
#include "meas/inline/hadron/inline_stoch_group_meson_w.h"
#include "meas/inline/hadron/inline_gauge_transf_obj.h"
#include "meas/inline/hadron/inline_rotate_spin_w.h"
//#include "meas/inline/hadron/inline_stoch_laph_quark_w.h"
//#include "meas/inline/hadron/inline_stoch_laph_baryon_w.h"

#include "meas/inline/hadron/inline_barspec_db_w.h"

// Grab all fermacts to make sure they are registered
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineHadronAggregateEnv
  {
    namespace
    {
      //! Local registration flag
      bool registered = false;
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Grab the fermacts
	success &= WilsonTypeFermActsEnv::registerAll();

	// Hadron stuff
	success &= InlineFermStateEnv::registerAll();

	success &= InlineMakeSourceEnv::registerAll();
	success &= InlinePropagatorEnv::registerAll();

	success &= InlineMakeSourceFermEnv::registerAll();
	success &= InlinePropagatorFermEnv::registerAll();

	success &= InlineMultiPropagatorEnv::registerAll();  // save space
	success &= InlineSeqSourceEnv::registerAll();
#ifdef BUILD_LAPACK
	success &= InlineLaplaceEigsEnv::registerAll();
#endif
	success &= InlineSeqPropTestEnv::registerAll();
	success &= InlineHadSpecEnv::registerAll();
	success &= InlineMesonSpecEnv::registerAll();
	success &= InlineHadronContractEnv::registerAll();
//	success &= InlineSpectrumEnv::registerAll();
	success &= InlineCreateColorVecsEnv::registerAll();
	success &= InlineProp3ptEnv::registerAll();
	success &= InlineDiscoEnv::registerAll();
	success &= InlineDiscoEOPrecEnv::registerAll();
	success &= InlineDiscoEoEigCGEnv::registerAll();
	success &= InlineDiscoEigCGEnv::registerAll();
	success &= InlineStagToWilsEnv::registerAll();
	success &= InlineSinkSmearEnv::registerAll();
	success &= InlineDiquarkEnv::registerAll();
	success &= InlineQQQDiquarkEnv::registerAll();
	success &= InlineQQQEnv::registerAll();
	success &= InlineQQbarEnv::registerAll();
	success &= InlineBuildingBlocksEnv::registerAll();
	success &= InlineNoisyBuildingBlocksEnv::registerAll();
	success &= InlineBar3ptfnEnv::registerAll();
//      success &= InlineMultipoleEnv::registerAll();  // not being used
	success &= InlineNprVertexEnv::registerAll();
	success &= InlineNprEnv::registerAll();
	success &= InlineDistillutionNoiseEnv::registerAll();
	success &= InlinePropDistillutionEnv::registerAll();
	success &= InlinePropDistillationEnv::registerAll();
	success &= InlinePropColorVecEnv::registerAll();
	success &= InlineStaticPropColorVecEnv::registerAll();
	success &= InlineAnnihPropMatElemColorVecEnv::registerAll();
	success &= InlinePropMatElemColorVecEnv::registerAll();
	success &= InlinePropMatElemPtColorVecEnv::registerAll();
	success &= InlinePropAndMatElemColorVecEnv::registerAll();
	success &= InlinePropMatElemLowMemoryColorVecEnv::registerAll();
	success &= InlineBaryonMatElemColorVecEnv::registerAll();
	success &= InlineMesonMatElemColorVecEnv::registerAll();
	success &= InlineGenPropMatElemColorVecEnv::registerAll();
	success &= InlineGenPropMatElemPtColorVecEnv::registerAll();
	success &= InlineMresEnv::registerAll();
	success &= InlineGaugeTransfNamedObjEnv::registerAll();
	success &= InlineRotateSpinEnv::registerAll();
	success &= InlineQpropQIOEnv::registerAll();
	success &= InlineQpropAddEnv::registerAll();
	success &= InlineQQQNucNucEnv::registerAll();
	success &= InlineBarSpecEnv::registerAll();
//	success &= InlineSpectrumQllEnv::registerAll();
	success &= InlineStaticLightSpecEnv::registerAll();
	success &= InlineHeavyLightContEnv::registerAll();
	success &= InlineHeavyHadSpecEnv::registerAll();
	success &= InlineStochMesonEnv::registerAll();
	success &= InlineStochBaryonEnv::registerAll();
	success &= InlineStochHadronEnv::registerAll();
	success &= InlineStochGroupBaryonEnv::registerAll();
	success &= InlineStochGroupMesonEnv::registerAll();
//	success &= InlineStochLaphQuarkEnv::registerAll();
//	success &= InlineStochLaphBaryonEnv::registerAll();

	registered = true;
      }
      return success;
    }

  }

}
