// -*- C++ -*-
// $Id: sfcurcor_w.h,v 3.1 2006-04-10 21:18:23 edwards Exp $
/*! \file
 * \brief Schroedinger functional correlation functions
 */

#ifndef __sfcurcor_w_h__
#define __sfcurcor_w_h__

#include "chromabase.h"
#include "util/ft/sftmom.h"

namespace Chroma
{

  //! Schroedinger functional correlation functions
  /*!
   * Construct 'current correlators' and axial density used for the PCAC determination
   * in the Schroedinger Functional
   *
   * \param quark_propagator    quark propagator ( Read )
   * \param pseudo_prop         pion correlator ( Write )
   * \param axial_prop          axial-current to pion_1 correlators ( Modify )
   * \param phases        object holds list of momenta and Fourier phases ( Read )
   *
   *
   *         \
   * cc(t) =  >  < m(0, 0) c(t + L, x) > 
   *         / 
   *         ----
   *           x
   */

  void SFcurcor(const LatticePropagator& quark_propagator, 
		multi1d<Real>& pseudo_prop, 
		multi1d<Real>& axial_prop, 
		const SftMom& phases);

}
#endif
