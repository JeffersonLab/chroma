// -*- C++ -*-
// $Id: fuzglue.h,v 3.1 2006-08-24 02:33:52 edwards Exp $
/*! \file
 *  \brief Compute 'fuzzy' (blocked) glueball correlation functions
 */

#ifndef __fuzglue_h__
#define __fuzglue_h__

namespace Chroma 
{

  //! Compute 'fuzzy' (blocked) glueball correlation functions
  /*! 
   * \ingroup glue
   *
   * Driver for computation of 'fuzzy' (blocked) glueball correlation functions
   * using Tepers 'fuzzying' method.
   *
   * Warning: this works only for Nd = 4 ! (Construction of glueball states)
   * Warning: this works only for Nc = 2 and 3 ! (Projection of blocked links)
   *
   * \param xml_out    xml file object ( Write )
   * \param xml_group  string used for writing xml data ( Read )
   * \param u          gauge field ( Read )
   * \param j_decay    direction along which the exponentrial is allowed to decay ( Read )
   * \param BlkAccu    accuracy in fuzzy link projection ( Read )
   * \param BlkMax     maximum number of iterations in fuzzy link projection ( Read ) 
   */

  void fuzglue(XMLWriter& xml_out, const string& xml_group,
	       const multi1d<LatticeColorMatrix>& u, 
	       int j_decay, const Real& BlkAccu, int BlkMax);

}  // end namespace Chroma

#endif
