// -*- C++ -*-
// $Id: quarkprop4_s.h,v 3.5 2006-10-19 16:01:34 edwards Exp $
/*! \file
 *  \brief Full quark propagator solver
 *
 *  Given a complete propagator as a source, this does all the inversions needed
 */

#ifndef __quarkprop4_s_h__
#define __quarkprop4_s_h__

#include "stagtype_fermact_s.h"

namespace Chroma
{
  //! Given a complete propagator as a source, this does all the inversions needed
  /*! \ingroup qprop
   *
   * This routine is actually generic to all Wilson-like fermions
   *
   * \param q_sol    quark propagator ( Write )
   * \param q_src    source ( Read )
   * \param invParam inverter parameters ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */

  void quarkProp4(LatticeStaggeredPropagator& q_sol, 
		  XMLWriter& xml_out,
		  const LatticeStaggeredPropagator& q_src,
		  const StaggeredTypeFermAct<LatticeStaggeredFermion, 
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >& S_f,
		  Handle< FermState<LatticeStaggeredFermion, 
		  multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state,
		  const GroupXML_t& invParam,
		  QuarkSpinType quarkSpinType,
		  int& ncg_had);

}
#endif
