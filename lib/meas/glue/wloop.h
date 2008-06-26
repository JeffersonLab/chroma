// -*- C++ -*-
// $Id: wloop.h,v 1.1 2008-06-26 14:58:35 mcneile Exp $

#ifndef __wloop_h__
#define __wloop_h__

namespace Chroma 
{


  //! Print the value of the average plaquette normalized to 1
  /*!
   * \ingroup glue
   *
   * \param xml    plaquette average (Write)
   * \param u      gauge field (Read)
   */
  void Wloop(XMLWriter& xml,
	      const string& xml_group,
	      const multi1d<LatticeColorMatrix>& u);

}  // end namespace Chroma

#endif
