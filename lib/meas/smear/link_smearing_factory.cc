// -*- C++ -*-
// $Id: link_smearing_factory.cc,v 2.2 2005-11-22 19:16:04 edwards Exp $
/*! \file
 *  \brief Factory for producing link smearing objects
 */

#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "handle.h"

#include "meas/glue/mesplq.h"

namespace Chroma
{

  // Convenience function to smear link
  void linkSmear(multi1d<LatticeColorMatrix>& u, const std::string& link_xml, const std::string& link_type)
  {
    bool foo = LinkSmearingEnv::registered;  // make sure all link smearings are loaded

    try
    {
      if (link_xml != "")
      {
	std::istringstream  xml_s(link_xml);
	XMLReader  linktop(xml_s);
	const string link_path = "/LinkSmearing";
	
	Handle< LinkSmearing >
	  linkSmearing(TheLinkSmearingFactory::Instance().createObject(link_type,
								       linktop,
								       link_path));
	(*linkSmearing)(u);
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": Caught Exception smearing: " << e << endl;
      QDP_abort(1);
    }


    // Paranoia check
    if (link_xml != "")
    {
      Double w_plaq, s_plaq, t_plaq, link;
      MesPlq(u, w_plaq, s_plaq, t_plaq, link);

      QDPIO::cout << __func__ << ": w_plaq=" << w_plaq << " s_plaq=" << s_plaq
		  << " t_plaq=" << t_plaq
		  << " link=" << link
		  << endl;
    }

  }


}
