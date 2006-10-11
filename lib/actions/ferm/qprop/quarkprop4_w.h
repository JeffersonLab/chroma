// -*- C++ -*-
// $Id: quarkprop4_w.h,v 3.3 2006-10-11 15:42:26 edwards Exp $
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
		  Handle< SystemSolver<LatticeFermion> > qprop,
		  QuarkSpinType quarkSpinType,
		  int& ncg_had);
}
#endif
