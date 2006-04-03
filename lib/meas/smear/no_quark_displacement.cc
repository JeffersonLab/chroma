// $Id: no_quark_displacement.cc,v 3.0 2006-04-03 04:59:05 edwards Exp $
/*! \file
 *  \brief No quark displacement
 */

#include "chromabase.h"

#include "meas/smear/quark_displacement_factory.h"
#include "meas/smear/no_quark_displacement.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, NoQuarkDisplacementEnv::Params& param)
  {
    NoQuarkDisplacementEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const NoQuarkDisplacementEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace NoQuarkDisplacementEnv
  {
    //! Callback function
    QuarkDisplacement<LatticePropagator>* createProp(XMLReader& xml_in,
						     const std::string& path)
    {
      return new QuarkDisplace<LatticePropagator>(Params(xml_in, path));
    }

    //! Callback function
    QuarkDisplacement<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path)
    {
      return new QuarkDisplace<LatticeFermion>(Params(xml_in, path));
    }
    
    //! Callback function
    QuarkDisplacement<LatticeColorVector>* createColorVec(XMLReader& xml_in,
							  const std::string& path)
    {
      return new QuarkDisplace<LatticeColorVector>(Params(xml_in, path));
    }
    
    //! Name to be used
    const std::string name = "NONE";

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(name, createProp);
      foo &= Chroma::TheFermDisplacementFactory::Instance().registerObject(name, createFerm);
      foo &= Chroma::TheColorVecDisplacementFactory::Instance().registerObject(name, createColorVec);
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();


    //! Parameters for running code
    Params::Params(XMLReader& xml, const string& path)
    {
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      write(xml, "DisplacementType", name);

      pop(xml);
    }

  }  // end namespace
}  // end namespace Chroma

