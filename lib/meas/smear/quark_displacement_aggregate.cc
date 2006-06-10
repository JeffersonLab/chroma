// $Id: quark_displacement_aggregate.cc,v 3.2 2006-06-10 16:28:19 edwards Exp $
/*! \file
 *  \brief All quark displacements
 */

#include "meas/smear/quark_displacement_aggregate.h"

#include "meas/smear/no_quark_displacement.h"
#include "meas/smear/simple_quark_displacement.h"
#include "meas/smear/deriv_quark_displacement_w.h"
#include "meas/smear/gamma_displacement_w.h"

namespace Chroma
{

  // Registration aggregator
  namespace QuarkDisplacementEnv
  {
    bool registerAll() 
    {
      bool success = true;

      success &= NoQuarkDisplacementEnv::registered;
      success &= SimpleQuarkDisplacementEnv::registered;
      success &= DerivQuarkDisplacementEnv::registered;
      success &= GammaDisplacementEnv::registered;

      return success;
    }

    const bool registered = registerAll();


    // Returns a no-displacement group
    GroupXML_t   nullXMLGroup()
    {
      GroupXML_t nope;

      XMLBufferWriter xml_tmp;
      NoQuarkDisplacementEnv::Params  non;
      write(xml_tmp, "Displacement", non);
      nope.xml = xml_tmp.str();
      nope.id = NoQuarkDisplacementEnv::name;

      return nope;
    }

  }

}
