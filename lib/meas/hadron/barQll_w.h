$Id: barQll_w.h,v 1.1 2004-09-08 14:00:13 kostas Exp $
/*! \file
 *  \brief Heavy Baryon (Qll)  2-pt function : Orginos and Savage
 */

#ifndef __barQll_h__
#define __barQll_h__

//! Lambdaq and SigmaQ 2-pt functions
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions! 
 *
 * Construct baryon propagators for the LambdaQ and SigmaQ(*) with
 * degenerate "u" and "d" quarks.  In the heavy quark limit the Sigma
 * and Sigma* are degenerate.  The heavy quark is inserted in the infinitely heavy quark limit
 * by a Wilson-Line without spin indices.
 * We are effectively propagating a spin-0 diquark and a spin-1 diquark.
 *
 * \param u                  gauge field (Read) 
 * \param quark_propagator   quark propagator ( Read )
 * \param src_coord          cartesian coordinates of the source ( Read )
 * \param phases             object holds list of momenta and Fourier phases ( Read )
 * \param xml                xml file object ( Read )
 * \param xml_group          group name for xml data ( Read )
 *
 */
void Qll(const multi1d<LatticeColorMatrix>& u, 
	 const LatticePropagator& quark_propagator,
	 const multi1d<int>& src_coord, 
	 const SftMom& phases,
       	 XMLWriter& xml,
	 const string& xml_group);

//! Heavy Quark Propagator
/*!
 * \ingroup hadron
 *
 * 
 * This constructs the propagator for a spinless Wilson-Line propagating from the 
 * point src_coord forward in time, and vanishing on previous time-slices.
 *
 * \param Qprop              Wilson-Line (write)
 * \param u                  Gauge Field (Read) 
 * \param src_coord          cartesian coordinates of the source ( Read )
 *
 */
void HeavyQuarkProp(LatticeColorMatrix& Qprop,
		    const multi1d<LatticeColorMatrix>& u,
		    const multi1d<int>& src_coord);
#end
