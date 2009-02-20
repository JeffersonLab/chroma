// $Id: link_smearing_aggregate.cc,v 3.10 2009-02-20 15:10:24 edwards Exp $
/*! \file
 *  \brief All link smearing applicators
 */

#include "meas/smear/link_smearing_aggregate.h"

#include "meas/smear/ape_link_smearing.h"
#include "meas/smear/hyp_link_smearing.h"
#include "meas/smear/no_link_smearing.h"
#include "meas/smear/stout_link_smearing.h"

namespace Chroma
{

  // Registration aggregator
  namespace LinkSmearingEnv
  {
    namespace
    {
      //! Local registration flag
      bool registered = false;
    }

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

	registered = true;
      }
      return success;
    }


    // Returns a no-smearing group
    GroupXML_t nullXMLGroup()
    {
      GroupXML_t nope;

      XMLBufferWriter xml_tmp;
      NoLinkSmearingEnv::Params  non;
      write(xml_tmp, "LinkSmearing", non);
      nope.xml = xml_tmp.printCurrentContext();
      nope.id = NoLinkSmearingEnv::getName();
      nope.path = "/LinkSmearing";

      return nope;
    }

  }

}
