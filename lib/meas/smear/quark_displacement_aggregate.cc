// $Id: quark_displacement_aggregate.cc,v 3.11 2009-02-20 15:10:24 edwards Exp $
/*! \file
 *  \brief All quark displacements
 */

#include "meas/smear/quark_displacement_aggregate.h"

#include "meas/smear/no_quark_displacement.h"
#include "meas/smear/simple_quark_displacement.h"
#include "meas/smear/deriv_quark_displacement_w.h"
#include "meas/smear/gamma_displacement_w.h"

#include "meas/smear/deriv_quark_displacement_s.h"
#include "meas/smear/quark_flavor_s.h"

namespace Chroma
{

  // Registration aggregator
  namespace QuarkDisplacementEnv
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
	// Wilson-type
	success &= NoQuarkDisplacementEnv::registerAll();
	success &= SimpleQuarkDisplacementEnv::registerAll();
	success &= DerivQuarkDisplacementEnv::registerAll();
	success &= GammaDisplacementEnv::registerAll();

	// Staggered-type
	success &= StaggeredDerivQuarkDisplacementEnv::registerAll();
	success &= StaggeredQuarkFlavorOpEnv::registerAll();

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
      nope.xml  = xml_tmp.str();
      nope.id   = NoQuarkDisplacementEnv::getName();
      nope.path = "/Displacement";

      return nope;
    }

  }

}
