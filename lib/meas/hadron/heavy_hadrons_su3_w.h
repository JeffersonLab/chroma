// -*- C++ -*-
// $Id: heavy_hadrons_su3_w.h,v 1.2 2008-12-21 21:22:36 edwards Exp $ 
/*! \file
 *  \brief Potential between 2 heavy hadrons : Detmold
 *  Correlators checked independentely by Savage
 *

 * Includes Lambda_b's etc

 * Soon also will have the BBB "potential"

 */

#include "chromabase.h"
#include "util/ft/sftmom.h"

namespace Chroma
{

  //! Heavy hadron spectrum for SU(3) isospin limit
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions! 
   *
   
   * \param u                  gauge field (Read) 
   * \param quark1             light quark propagator ( Read )
   * \param quark2             strange quark propagator ( Read )
   * \param src               cartesian coordinates of one source "0"( Read )
   * \param phases             object holds list of momenta and Fourier phases ( Read )
   * \param xml                xml file object ( Read )
   * \param xml_group          group name for xml data ( Read )
   *
   */

  void static_light_su3(const multi1d<LatticeColorMatrix>& u, 
			const LatticePropagator& quark1,
			const LatticePropagator& quark2,
			const multi1d<int>& src,
			const SftMom& phases,
			XMLWriter& xml,
			const string& xml_group);


}

