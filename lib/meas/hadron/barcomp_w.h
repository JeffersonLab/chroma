// -*- C++ -*-
// $Id: barcomp_w.h,v 3.0 2006-04-03 04:58:59 edwards Exp $
/*! \file
 *  \brief Construct all components of a baryon propagator
 */

#ifndef __barcomp_h__
#define __barcomp_h__

#include "io/qprop_io.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{

  //! Convert generalized correlator object
  void convertBarcomp(multi1d<Complex>& barprop1d, const multiNd<Complex>& barprop, 
		      int j_decay);


  //! Construct all components of a baryon propagator
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * In all baryons the colour components are contracted with the totally
   * antisymmetric 'tensor' eps(a,b,c) = antisym_tensor(a,b,c).
   *
   * \param barprop                  baryon correlation function (in real space) ( Write )
   * \param quark_propagator_1       quark propagator ( Read )
   * \param quark_propagator_2       quark propagator ( Read )
   * \param quark_propagator_3       quark propagator ( Read )
   * \param phases                   object holds list of momenta ( Read )
   * \param t0                       coordinates of source in decay direction ( Read )
   * \param bc_spec                  boundary condition for spectroscopy ( Read )
   */

  void barcomp(multiNd<Complex>& barprop,
	       const LatticePropagator& quark_propagator_1, 
	       const LatticePropagator& quark_propagator_2,
	       const LatticePropagator& quark_propagator_3,
	       const SftMom& phases,
	       int t0, int bc_spec);

}  // end namespace Chroma

#endif
