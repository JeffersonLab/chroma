// -*- C++ -*-
// $Id: barcomp_w.h,v 1.6 2004-01-05 01:00:33 edwards Exp $
/*! \file
 *  \brief Construct all components of a baryon propagator
 */

#ifndef __barcomp_h__
#define __barcomp_h__

#include "io/qprop_io.h"

//! Construct all components of a baryon propagator
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 * In all baryons the colour components are contracted with the totally
 * antisymmetric 'tensor' eps(a,b,c) = antisym_tensor(a,b,c).
 *
 * \param quark_propagator_{1,2,3} quark propagators (read)
 * \param header_{1,2,3}           second quark propagator ( Read )
 * \param quark_propagator_3       third quark propagator  ( Read )
 * \param barprop                  baryon correlation function (in real space) ( Write )
 * \param num_mom                  number of non-zero momenta ( Read )
 * \param t0       cartesian coordinates of the source in the j_decay direction ( Read )
 * \param j_decay                  direction of the exponential decay ( Read )
 * \param bc_spec                  boundary condition for spectroscopy ( Read )
 * \param file                     the filename of the binary output file
 */
void barcomp(const LatticePropagator& quark_propagator_1, 
	     const PropHead& header_1,
	     const LatticePropagator& quark_propagator_2,
	     const PropHead& header_2,
	     const LatticePropagator& quark_propagator_3,
	     const PropHead& header_3,
	     int t0, int j_decay, int bc_spec,
	     const char file[],
	     NmlWriter& nml);

#endif
