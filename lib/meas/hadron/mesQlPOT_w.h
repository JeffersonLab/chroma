// -*- C++ -*-
//$Id: mesQlPOT_w.h,v 1.2 2008-12-21 21:22:37 edwards Exp $
/*! \file
 *  \brief Potential between 2 heavy mesons : Orginos and Savage
 */

#ifndef __mesQlPOT_h__
#define __mesQlPOT_h__

#include "chromabase.h"
#include "util/ft/sftmom.h"

namespace Chroma
{

  //! Heavy-light meson potential
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions! 
   *
   * Construct propagators for two heavy mesons in all combinations of spin up and 
   * spin down light degrees of freedom..
   * In the heavy quark limit the D and the D* are degenerate.  
   * The heavy quark is inserted in the infinitely heavy quark limit
   * by a Wilson-Line without spin indices.
   * We are effectively propagating two  spin-1/2 light degrees of freedom.
   *
   * \param u                  gauge field (Read) 
   * \param quark1             quark propagator ( Read )
   * \param quark2             quark propagator ( Read )
   * \param src1               cartesian coordinates of one source ( Read )
   * \param src2               cartesian coordinates of the other source ( Read )
   * \param phases             object holds list of momenta and Fourier phases ( Read )
   * \param xml                xml file object ( Read )
   * \param xml_group          group name for xml data ( Read )
   *
   */

  void QlQlPOT(const multi1d<LatticeColorMatrix>& u, 
	       const LatticePropagator& quark1,
	       const LatticePropagator& quark2,
	       const multi1d<int>& src1, 
	       const multi1d<int>& src2, 
	       const SftMom& phases,
	       XMLWriter& xml,
	       const string& xml_group);

  //!  Spin Transpose Function
  /*!
   * \ingroup hadron
   *
   * This is a dumb way of taking the spin transpose of a propagator,
   * while leaving all other indices untouched...suggested to 
   * us by Bob Edwards
   * \param prop                input propagator
   * \param STprop              spin transposed propagator
   *
   */
  
  void SpinTranspose(const LatticePropagator& prop, 
		     LatticePropagator& STprop);


}

#endif
