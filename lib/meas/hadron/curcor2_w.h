// -*- C++ -*-
// $Id: curcor2_w.h,v 1.2 2003-10-10 03:46:46 edwards Exp $

#ifndef __curcor2_h__
#define __curcor2_h__

//! Construct current correlators 
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 *  The two propagators can be identical or different.

 * This includes the "rho_1--rho_2" correlators used for O(a) improvement

 * For use with "rotated" propagators we added the possibility of also
 * computing the local vector current, when no_vec_cur = 4. In this
 * case the 3 local currents come last.

 * \param u               gauge field ( Read )
 * \param quark_prop_1    first quark propagator ( Read )
 * \param quark_prop_2    second (anti-) quark propagator ( Read )
 * \param vector_current  {3,4}*(Nd-1) vector-current to rho_1 correlators ( Modify )
 * \param axial_current   2 axial-current to pion_1 correlators ( Modify )
 * \param t0              timeslice coordinate of the source ( Read )
 * \param j_decay         direction of the exponential decay ( Read )
 * \param no_vec_cur      number of vector current types, 3 or 4 ( Read )
 * \param xml             namelist file object ( Read )
 * \param xml_group       string used for writing xml data ( Read )
 *
 *         ____
 *         \
 * cc(t) =  >  < m(t_source, 0) c(t + t_source, x) >
 *         /                    
 *         ----
 *           x
 */

void curcor2(const multi1d<LatticeColorMatrix>& u, 
	     const LatticePropagator& quark_prop_1, 
	     const LatticePropagator& quark_prop_2, 
	     const SftMom& phases,
	     int t0,
	     int j_decay, int no_vec_cur,
	     XMLWriter& xml,
	     const string& xml_group);

#endif
