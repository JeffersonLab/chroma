// -*- C++ -*-
// $IdR$
/*! \file
 *  \brief Construct all components of a meson propagator
 */

#ifndef __mescomp_h__
#define __mescomp_h__

#include "io/qprop_io.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{

  //! Convert generalized correlator object
  void convertMescomp(multi1d<Complex>& mesprop1d, const multiNd<Complex>& mesprop, 
		      const int j_decay);

  //! Construct all components of a meson propagator
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * In all mesons the colour components are contracted leaving only the
   * spin components.
   *
   * \param mesprop                  meson correlation function (in real space) ( Write )
   * \param quark_propagator_1       quark propagator ( Read )
   * \param quark_propagator_2       quark propagator ( Read )
   * \param phases                   object holds list of momenta ( Read )
   * \param t0                       coordinates of source in decay direction ( Read )
   */

  void mescomp(multiNd<Complex>& mesprop,
	       const LatticePropagator& quark_propagator_1, 
	       const LatticePropagator& quark_propagator_2,
	       const SftMom& phases,
	       int t0);

}  // end namespace Chroma

#endif
