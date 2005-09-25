// $Id: quarkprop4_w.h,v 2.0 2005-09-25 21:04:31 edwards Exp $
/*! \file
 *  \brief Full quark propagator solver
 *
 *  Given a complete propagator as a source, this does all the inversions needed
 */

#ifndef __quarkprop4_w_h__
#define __quarkprop4_w_h__

#include "fermact.h"

namespace Chroma
{
  //! Given a complete propagator as a source, this does all the inversions needed
  /*! \ingroup qprop
   *
   * This routine is actually generic to all Wilson-like fermions
   *
   * \param q_sol    quark propagator ( Write )
   * \param q_src    source ( Read )
   * \param qprop    inverter ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */

  void quarkProp4(LatticePropagator& q_sol, 
		  XMLWriter& xml_out,
		  const LatticePropagator& q_src,
		  Handle<const SystemSolver<LatticeFermion> > qprop,
		  bool nonRelProp,
		  int& ncg_had);
};
#endif
