// $Id: quark_smearing_aggregate.cc,v 3.2 2006-06-10 16:28:19 edwards Exp $
/*! \file
 *  \brief All quark smearing
 */

#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/no_quark_smearing.h"
#include "meas/smear/gaus_quark_smearing.h"

namespace Chroma
{

  // Registration aggregator
  namespace QuarkSmearingEnv
  {
    bool registerAll() 
    {
      bool success = true;

      success &= NoQuarkSmearingEnv::registered;
      success &= GausQuarkSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();


    // Returns a no-smearing group
    GroupXML_t   nullXMLGroup()
    {
      GroupXML_t nope;

      XMLBufferWriter xml_tmp;
      NoQuarkSmearingEnv::Params  non;
      write(xml_tmp, "SmearingParam", non);
      nope.xml = xml_tmp.str();
      nope.id = NoQuarkSmearingEnv::name;

      return nope;
    }

  }

}
