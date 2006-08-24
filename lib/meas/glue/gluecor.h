// -*- C++ -*-
// $Id: gluecor.h,v 3.1 2006-08-24 02:33:52 edwards Exp $
/*! \file
 *  \brief Construct 0++, 2++ and 1+- glueball correlation functions from fuzzy links
 */

#ifndef __gluecor_h__
#define __gluecor_h__

#include "util/ft/sftmom.h"

namespace Chroma 
{

  //! Construct 0++, 2++ and 1+- glueball correlation functions from fuzzy links
  /*! 
   * \ingroup glue
   *
   * Construct 0++, 2++ and 1+- glueball correlation functions from
   * fuzzy links at blocking level bl_level and Write them in 
   * XML format.
   *
   * Warning: this works only for Nd = 4 !
   *
   * \param xml_out       xml file object ( Write )
   * \param xml_group     string used for writing xml data ( Read )
   * \param u             (blocked) gauge field ( Read )
   * \param bl_level      blocking level ( Read )
   * \param phases        object holds list of momenta and Fourier phases ( Read )
   */

  void gluecor(XMLWriter& xml_out, const string& xml_group,
	       const multi1d<LatticeColorMatrix>& u, 
	       const SftMom& phases,
	       int bl_level);

}  // end namespace Chroma

#endif
