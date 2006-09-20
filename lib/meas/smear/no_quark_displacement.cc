// $Id: no_quark_displacement.cc,v 3.1 2006-09-20 20:28:04 edwards Exp $
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

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(name, createProp);
	success &= Chroma::TheFermDisplacementFactory::Instance().registerObject(name, createFerm);
	success &= Chroma::TheColorVecDisplacementFactory::Instance().registerObject(name, createColorVec);
	registered = true;
      }
      return success;
    }


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

