// -*- C++ -*-
// $Id: mesons2_w.h,v 1.1 2006-07-10 19:53:36 edwards Exp $
/*! \file
 *  \brief Meson 2-pt functions
 */

#ifndef __mesons2_h__
#define __mesons2_h__

namespace Chroma {

//! Meson 2-pt functions
/* This routine is specific to Wilson fermions!
 *
 * Construct meson propagators and writes in COMPLEX 
 *
 * The two propagators can be identical or different.
 *
 * \param quark_prop_1  first quark propagator ( Read )
 * \param quark_prop_2  second (anti-) quark propagator ( Read )
 * \param t0            timeslice coordinate of the source ( Read )
 * \param phases        object holds list of momenta and Fourier phases ( Read )
 * \param xml           xml file object ( Write )
 * \param xml_group     string used for writing xml data ( Read )
 *
 *        ____
 *        \
 * m(t) =  >  < m(t_source, 0) m(t + t_source, x) >
 *        /
 *        ----
 *          x
 */

void mesons2(const LatticePropagator& quark_prop_1,
	     const LatticePropagator& quark_prop_2,
	     const SftMom& phases,
	     int t0,
	     XMLWriter& xml,
	     const string& xml_group) ;
  
}  // end namespace Chroma

#endif
