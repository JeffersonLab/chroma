// $Id: link_smearing_aggregate.cc,v 3.7 2008-01-20 03:07:51 edwards Exp $
/*! \file
 *  \brief All link smearing applicators
 */

#include "meas/smear/link_smearing_aggregate.h"

#include "meas/smear/ape_link_smearing.h"
#include "meas/smear/hyp_link_smearing.h"
#include "meas/smear/no_link_smearing.h"
#include "meas/smear/stout_link_smearing.h"
#include "meas/smear/stout_link_smearing_old.h"

namespace Chroma
{

  // Registration aggregator
  namespace LinkSmearingEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// link smearing
	success &= APELinkSmearingEnv::registerAll();
	success &= HypLinkSmearingEnv::registerAll();
	success &= NoLinkSmearingEnv::registerAll();
	success &= StoutLinkSmearingEnv::registerAll();
	success &= StoutLinkSmearingOldEnv::registerAll();

	registered = true;
      }
      return success;
    }


    // Returns a no-smearing group
    GroupXML_t   nullXMLGroup()
    {
      GroupXML_t nope;

      XMLBufferWriter xml_tmp;
      NoLinkSmearingEnv::Params  non;
      write(xml_tmp, "LinkSmearing", non);
      nope.xml = xml_tmp.printCurrentContext();
      nope.id = NoLinkSmearingEnv::name;
      nope.path = "/LinkSmearing";

      return nope;
    }

  }

}
