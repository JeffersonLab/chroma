// -*- C++ -*-
// $Id: baryon_w.h,v 3.0 2006-04-03 04:58:59 edwards Exp $
/*! \file
 *  \brief Baryon 2-pt functions
 */

#ifndef __baryon_h__
#define __baryon_h__

namespace Chroma {

//! Baryon 2-pt functions
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions! 
 *
 * Construct baryon propagators for the Proton and the Delta^+ with
 * degenerate "u" and "d" quarks, as well as the Lambda for, in
 * addition, a degenerate "s" quark. For these degenerate quarks, the
 * Lambda is degenerate with the Proton, but we keep it for compatibility
 * with the sister routine that treats non-degenerate quarks.
  multi3d<DComplex>& bardisp1;

 * The routine optionally computes time-charge reversed baryons and adds them
 * in for increased statistics.

 * \param quark_propagator   quark propagator ( Read )
 * \param t0         cartesian coordinates of the source ( Read )
 * \param bc_spec    boundary condition for spectroscopy ( Read )
 * \param time_rev   add in time reversed contribution if true ( Read )
 * \param phases     object holds list of momenta and Fourier phases ( Read )
 * \param xml        xml file object ( Read )
 * \param xml_group  group name for xml data ( Read )
 *
 */

void baryon(const LatticePropagator& quark_propagator, 
            const SftMom& phases,
            int t0, int bc_spec, bool time_rev,
            XMLWriter& xml,
	    const string& xml_group);


//! Baryon 2-pt functions
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions! 
 *
 * Construct baryon propagators for the Proton and the Delta^+ with
 * degenerate "u" and "d" quarks, as well as the Lambda for, in
 * addition, a degenerate "s" quark. For these degenerate quarks, the
 * Lambda is degenerate with the Proton, but we keep it for compatibility
 * with the sister routine that treats non-degenerate quarks.
 *
 * \param quark_propagator  quark propagator ( Read )
 * \param barprop    baryon propagator ( Modify )
 * \param phases     object holds list of momenta and Fourier phases ( Read )
 *
 *        ____
 *        \
 * b(t) =  >  < b(t_source, 0) b(t + t_source, x) >
 *        /                    
 *        ----
 *          x

 * For the Proton we take

 * |P_1, s_z=1/2> = (d C gamma_5 u) "u_up"

 * for the Lambda

 * |L_1, s_z=1/2> = 2*(u C gamma_5 d) "s_up" + (s C gamma_5 d) "u_up"
 *                  + (u C gamma_5 s) "d_up"

 * and for the Delta^+

 * |D_1, s_z=3/2> = 2*(d C gamma_- u) "u_up" + (u C gamma_- u) "d_up".

 * We have put "q_up" in quotes, since this is meant in the Dirac basis,
 * not in the 'DeGrand-Rossi' chiral basis used in the program!

 * For all baryons we compute a 'B_2' that differs from the 'B_1' above
 * by insertion of a gamma_4 between C and the gamma_{5,-}.
 * And finally, we also compute the non-relativistic baryons, 'B_3',
 * which up to a factor 1/2 are just the difference B_1 - B_2, as can
 * be seen by projecting to the "upper" components in the Dirac basis,
 * achieved by (1 + gamma_4)/2 q, for quark q.

 * The Proton_k is baryon 3*(k-1), the Lambda_k is baryon 3*(k-1)+1
 * and the Delta^+_k is baryon 3*(k-1)+2. 
 */

void baryon(const LatticePropagator& quark_propagator, 
            const SftMom& phases,
	    multi3d<DComplex>& barprop);

}  // end namespace Chroma

#endif
