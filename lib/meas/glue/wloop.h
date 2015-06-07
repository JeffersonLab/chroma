// -*- C++ -*-

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
	      const std::string& xml_group,
	      const multi1d<LatticeColorMatrix>& u);

}  // end namespace Chroma

#endif
