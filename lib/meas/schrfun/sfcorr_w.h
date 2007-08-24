// -*- C++ -*-
// $Id: sfcorr_w.h,v 3.1 2007-08-24 19:23:04 edwards Exp $
/*! \file
 * \brief Schroedinger functional correlation functions
 */

#ifndef __sfcorr_w_h__
#define __sfcorr_w_h__

#include "chromabase.h"
#include "util/ft/sftmom.h"

namespace Chroma
{

  //! Schroedinger functional correlation functions
  /*!
   * @ingroup schrfun
   *
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

  void SFcorr(multi1d<Real>& pseudo_prop, 
	      multi1d<Real>& axial_prop, 
	      const LatticePropagator& quark_propagator, 
	      const SftMom& phases);

}
#endif
