// -*- C++ -*-
// $Id: qqbar_w.h,v 3.2 2008-12-21 21:22:37 edwards Exp $
/*! \file
 *  \brief constructs 2 quark propagators contracted at the sink
 */

#ifndef __qqbar_h__
#define __qqbar_h__


namespace Chroma 
{

  //! Meson-Meson 4-pt functions
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * Construct meson-meson propagators
   * The two propagators can be identical or different.
   *
   * \param qqbar           -- the 2-quark propagator ( Write )
   * \param quark_prop_1 -- first quark propagator ( Read )
   * \param quark_prop_2 -- second (anti-) quark propagator ( Read )
   * \param t0 -- timeslice coordinate of the source ( Read )
   * \param phases -- object holds list of momenta and Fourier phases ( Read )
   *
   *          ____
   *          \                               +
   *qqbar(p,t)=> [g5 q2(t_src;t + t_src,x) g5] g5 q1(t+t_src,x;t_src) *  exp(ipx)
   *          /
   *          ----
   *            x
   */

  void compute_qqbar( multi2d<DPropagator>& qqbar,
		      const LatticePropagator& quark_prop_1,
		      const LatticePropagator& quark_prop_2, 
		      const SftMom& phases,
		      int t0) ;

  // gamma is the gamma matrix at the sink 
  void compute_qqbar( multi2d<DPropagator>& qqbar, const int gamma,
		      const LatticePropagator& quark_prop_1,
		      const LatticePropagator& quark_prop_2, 
		      const SftMom& phases,
		      int t0) ;
  
  void write_qqbar(QDPFileWriter& to,
		   multi2d<DPropagator>& qqbar, 
		   const SftMom& phases,
		   string type,
		   string sink);

} //namespace Chroma

#endif
