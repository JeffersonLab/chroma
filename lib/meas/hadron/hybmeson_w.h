// -*- C++ -*-
// $Id: hybmeson_w.h,v 3.1 2008-09-05 13:59:18 edwards Exp $
/*! \file
 *  \brief Hybrid meson 2-pt functions
 */

#ifndef __hybmeson_h__
#define __hybmeson_h__

namespace Chroma 
{

  //! Hybrid meson 2-pt functions
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * First we construct a hybrid pion and 3 hybrid rho's, followed by
   * an exotic 0^{+-}, an exotic 0^{--} and finally 2*3 exotic 1^{-+}'s.
   *
   * \param f             field strength tensor ( Read )
   * \param u_smr         the SMEARED gauge field, used in constructing the f's
   * \param quark_prop_1  first quark propagator ( Read )
   * \param quark_prop_2  second (anti-) quark propagator ( Read )
   * \param t_source      cartesian coordinates of the source ( Read )
   * \param phases        object holds list of momenta and Fourier phases ( Read )
   * \param xml           xml file object ( Read )
   * \param xml_group     string used for writing xml data ( Read )
   *
   *        ____
   *        \
   * m(t) =  >  < m(t_source, 0) m(t + t_source, x) >
   *        /                    
   *        ----
   *          x 
   */

  void hybmeson(const multi1d<LatticeColorMatrix>& f, 
		const multi1d<LatticeColorMatrix>& u_smr, 
		const LatticePropagator& quark_prop_1,
		const LatticePropagator& quark_prop_2, 
		const SftMom& phases,
		multi1d<int> t_source,
		XMLWriter& xml,
		const string& xml_group);

}  // end namespace Chroma

#endif
