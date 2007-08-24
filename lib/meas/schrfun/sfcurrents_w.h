// -*- C++ -*-
// $Id: sfcurrents_w.h,v 3.1 2007-08-24 19:23:04 edwards Exp $
/*! \file
 *  \brief Current renormalizations within Schroedinger functional
 */

#ifndef __sfcurrents_w_h__
#define __sfcurrents_w_h__

#include "fermact.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{
  
  //! Compute the kprop used in PCAC
  /*! @ingroup schrfun */
  Propagator SFKprop(const LatticePropagator& quark_prop_f,
		     Handle< FermState<LatticeFermion, multi1d<LatticeColorMatrix>,
		     multi1d<LatticeColorMatrix> > > state,
		     const SftMom& phases);

  //! Compute Z_V
  /*! @ingroup schrfun */
  void SFCurrentZV(XMLWriter& xml_out, 
		   const string& xml_group,
		   const LatticePropagator& quark_prop_f,
		   const LatticePropagator& quark_prop_b,
		   Handle< SystemSolver<LatticeFermion> > qprop,
		   Handle< FermState<LatticeFermion, multi1d<LatticeColorMatrix>,
		   multi1d<LatticeColorMatrix> > > state,
		   const SftMom& phases);

  //! Compute Z_V
  /*! 
   * @ingroup schrfun 
   *
   * @return number of inverter iterations
   */
  int SFCurrentZA(XMLWriter& xml_out, 
		  const string& xml_group,
		  const multi1d<Real>& pseudo_prop_f,
		  const multi1d<Real>& axial_prop_f,
		  const multi1d<Real>& pseudo_prop_b,
		  const multi1d<Real>& axial_prop_b,
		  const LatticePropagator& quark_prop_f,
		  const LatticePropagator& quark_prop_b,
		  Handle< SystemSolver<LatticeFermion> > qprop,
		  Handle< FermState<LatticeFermion, multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> > > state,
		  const SftMom& phases,
		  int x0, int y0);

}

#endif
