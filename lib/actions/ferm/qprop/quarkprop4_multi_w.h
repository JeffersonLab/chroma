// -*- C++ -*-
/*! \file
 *  \brief Full quark propagator solver
 *
 *  Given a complete propagator as a source, this does all the inversions needed
 */

#ifndef __quarkprop4_multi_w_h__
#define __quarkprop4_multi_w_h__

#include "chromabase.h"

#include "fermact.h"
#include "actions/ferm/fermacts/overlap_fermact_base_w.h"



namespace Chroma 
{

  //! Given a complete propagator as a source, this does all the inversions needed
  /*! \ingroup qprop */
  void 
  OverlapFermActBase::multiQuarkProp4(multi1d<LatticePropagator>& q_sol,
				      XMLWriter& xml_out,
				      const LatticePropagator& q_src,
				      Handle< FermState<LatticeFermion,
				      multi1d<LatticeColorMatrix>,
				      multi1d<LatticeColorMatrix> > state,
				      const multi1d<Real>& masses,
				      const Multi1InvertParam_t& invParam,
				      const int n_soln,
				      int& ncg_had);

}
#endif
