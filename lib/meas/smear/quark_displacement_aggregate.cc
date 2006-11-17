// $Id: quark_displacement_aggregate.cc,v 3.4 2006-11-17 02:55:11 edwards Exp $
/*! \file
 *  \brief All quark displacements
 */

#include "meas/smear/quark_displacement_aggregate.h"

#include "meas/smear/no_quark_displacement.h"
#include "meas/smear/simple_quark_displacement.h"
#include "meas/smear/deriv_quark_displacement_w.h"
#include "meas/smear/gamma_displacement_w.h"

#include "meas/smear/deriv_quark_displacement_s.h"

namespace Chroma
{

  // Registration aggregator
  namespace QuarkDisplacementEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Wilson-type
	success &= NoQuarkDisplacementEnv::registerAll();
	success &= SimpleQuarkDisplacementEnv::registerAll();
	success &= DerivQuarkDisplacementEnv::registerAll();
	success &= GammaDisplacementEnv::registerAll();

	// Staggered-type
	success &= StaggeredDerivQuarkDisplacementEnv::registerAll();

	registered = true;
      }
      return success;
    }


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
