// $Id: dwf_quarkprop4_w.h,v 2.0 2005-09-25 21:04:30 edwards Exp $
/*! \file
 * \brief Full quark propagator solver for domain wall fermions
 *
 * Given a complete propagator as a source, this does all the inversions needed
 */

#ifndef __dwf_quarkprop4_w_h__
#define __dwf_quarkprop4_w_h__

#include "fermact.h"

namespace Chroma
{
  //! Given a complete propagator as a source, this does all the inversions needed
  /*! \ingroup qprop
   *
   * This routine is actually generic to Domain Wall fermions (Array) fermions
   *
   * \param q_sol    quark propagator ( Write )
   * \param q_src    source ( Read )
   * \param t_src    time slice of source ( Read )
   * \param j_decay  direction of decay ( Read )
   * \param invParam inverter parameters ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */
  void
  dwf_quarkProp4(LatticePropagator& q_sol, 
		 XMLWriter& xml_out,
		 const LatticePropagator& q_src,
		 int t_src, int j_decay,
		 Handle< const SystemSolver< multi1d<LatticeFermion> > > qpropT,
		 Handle<const ConnectState> state,
		 const Real& m_q,
		 int& ncg_had);

}

#endif
