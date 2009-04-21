// -*- C++ -*-
// $Id: mespbp_w.h,v 3.1 2009-04-21 01:21:42 eneil Exp $

/*! \file
 * \brief Psibar-psi measurements
 *
 * Support for computing psibar-psi in various ways and fermion actions
 */

/*! \defgroup pbp Psibar-psi measurements
 * \ingroup meas
 *
 * Support for computing psibar-psi in various ways and fermion actions
 */

#ifndef __mespbp_w_h__
#define __mespbp_w_h__

#include "fermact.h"

namespace Chroma
{
	void MesPbp(
		Handle< SystemSolver<LatticeFermion> > qprop,
		Handle< FermState<LatticeFermion, multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> > > state,
		const multi1d<Real>& Mass,
		const int ichiral,
		XMLWriter& xml_out,
		const string& xml_group,
		const string& FermAct);

} // end namespace Chroma

#endif
