// -*- C++ -*-
// $Id: mesons_w.h,v 1.4 2003-03-14 05:14:32 flemingg Exp $

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
 * \param t0 -- timeslice coordinate of the source ( Read )
 * \param phases -- object holds list of momenta and Fourier phases ( Read )
 * \param nml -- namelist file object ( Read )
 * \param nml_group -- string used for writing nml data ( Read )
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
            SftMom& phases,
            int t0,
            NmlWriter& nml,
            char* nml_group) ;

#endif
