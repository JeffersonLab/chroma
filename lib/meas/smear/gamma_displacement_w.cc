// $Id: gamma_displacement_w.cc,v 3.3 2008-11-04 18:43:57 edwards Exp $
/*! \file
 *  \brief Gamma insertion/displacement
 */

#include "chromabase.h"

#include "meas/smear/quark_displacement_factory.h"
#include "meas/smear/gamma_displacement_w.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, GammaDisplacementEnv::Params& param)
  {
    GammaDisplacementEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const GammaDisplacementEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace GammaDisplacementEnv
  {
    namespace
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
    
      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name = "GAMMA_INSERTION";
    }

    //! Return the name
    std::string getName() {return name;}

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(name, createProp);
	success &= Chroma::TheFermDisplacementFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }


    //! Parameters for running code
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "gamma", gamma);
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      write(xml, "gamma", gamma);

      pop(xml);
    }

  }  // end namespace
}  // end namespace Chroma

