// -*- C++ -*-
// $Id: polycor.h,v 3.1 2006-08-24 02:33:52 edwards Exp $
/*! \file
 *  \brief Construct Polyakov loop correlation functions from fuzzy links
 */

#ifndef __polycor_h__
#define __polycor_h__

#include "util/ft/sftmom.h"

namespace Chroma 
{

  //! Construct Polyakov loop correlation functions from fuzzy links
  /*!
   * \ingroup glue
   *
   * Construct Polyakov loop correlation functions from fuzzy links at 
   * blocking level bl_level in the directions orthogonal to j_decay
   * and Write them in (pseudo) XML format.
   *
   * \param xml_out      xml file object ( Write )
   * \param xml_group    string used for writing xml data ( Read )
   * \param u            (blocked) gauge field ( Read )
   * \param block_latt    block lattice size ( Read )
   * \param bl_level      blocking level ( Read )
   * \param phases        object holds list of momenta and Fourier phases ( Read )
   */

  void polycor(XMLWriter& xml_out, const string& xml_group,
	       const multi1d<LatticeColorMatrix>& u,
	       const SftMom& phases,
	       int bl_level);

}  // end namespace Chroma

#endif
