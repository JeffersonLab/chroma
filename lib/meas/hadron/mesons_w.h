// -*- C++ -*-
// $Id: mesons_w.h,v 1.2 2003-02-26 03:19:36 edwards Exp $

#ifndef MESONS_INCLUDE
#define MESONS_INCLUDE

//! Meson 2-pt functions
/* This routine is specific to Wilson fermions!
 *
 * Construct meson propagators
 * The two propagators can be identical or different.
 *
 * \param quark_prop_1 -- first quark propagator ( Read )
 * \param quark_prop_2 -- second (anti-) quark propagator ( Read )
 * \param meson_propagator -- Ns^2 mesons ( Modify )
 * \param t_source -- cartesian coordinates of the source ( Read )
 * \param j_decay -- direction of the exponential decay ( Read )
 *
 *        ____
 *        \
 * m(t) =  >  < m(t_source, 0) m(t + t_source, x) >
 *        /
 *        ----
 *          x
 */
void mesons(const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2, 
	    multi2d<Real>& meson_propagator, 
	    const multi1d<int>& t_source, int j_decay);

#endif
