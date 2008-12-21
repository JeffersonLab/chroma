// -*- C++ -*-
//$Id: mesQl_w.h,v 1.3 2008-12-21 21:22:37 edwards Exp $
/*! \file
 *  \brief Heavy light meson (Qlbar)  2-pt function : Orginos and Savage
 */

#ifndef __mesQl_h__
#define __mesQl_h__

#include "chromabase.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{

  //! Heavy-light meson 2-pt function
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions! 
   *
   * Construct propagators for a heavy-light pseudoscalar meson.
   * In the heavy quark limit the D and the D* are degenerate.  
   * The heavy quark is inserted in the infinitely heavy quark limit
   * by a Wilson-Line without spin indices.
   * We are effectively propagating a spin-1/2 light quark
   *
   * \param u                  gauge field (Read) 
   * \param quark_propagator   quark propagator ( Read )
   * \param src_coord          cartesian coordinates of the source ( Read )
   * \param phases             object holds list of momenta and Fourier phases ( Read )
   * \param xml                xml file object ( Read )
   * \param xml_group          group name for xml data ( Read )
   *
   */
  void Qlbar(const multi1d<LatticeColorMatrix>& u, 
	     const LatticePropagator& quark_propagator,
	     const multi1d<int>& src_coord, 
	     const SftMom& phases,
	     XMLWriter& xml,
	     const string& xml_group);

}
#endif
