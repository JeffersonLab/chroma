// $Id: link_smearing_aggregate.cc,v 3.2 2006-06-10 16:28:19 edwards Exp $
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
    bool registerAll() 
    {
      bool success = true;

      // link smearing
      success &= APELinkSmearingEnv::registered;
      success &= HypLinkSmearingEnv::registered;
      success &= NoLinkSmearingEnv::registered;
      success &= StoutLinkSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();



    // Returns a no-smearing group
    GroupXML_t   nullXMLGroup()
    {
      GroupXML_t nope;

      XMLBufferWriter xml_tmp;
      NoLinkSmearingEnv::Params  non;
      write(xml_tmp, "LinkSmearing", non);
      nope.xml = xml_tmp.str();
      nope.id = NoLinkSmearingEnv::name;

      return nope;
    }

  }

}
