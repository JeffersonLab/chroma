// -*- C++ -*-
// $Id: mesons_w.h,v 1.3 2003-03-06 00:30:14 flemingg Exp $

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
 * \param t_source -- cartesian coordinates of the source ( Read )
 * \param sink_mom2_max -- max sink hadron mom squared ( Read )
 * \param j_decay -- direction of the exponential decay ( Read )
 * \param nml -- namelist file object ( Read )
 *
 *        ____
 *        \
 * m(t) =  >  < m(t_source, 0) m(t + t_source, x) >
 *        /
 *        ----
 *          x
 */
void mesons(const LatticePropagator& quark_prop_1,
            const LatticePropagator& quark_prop_2,
            const multi1d<int>& t_source,
            int sink_mom2_max,
            int j_decay,
            NmlWriter& nml);

#endif
