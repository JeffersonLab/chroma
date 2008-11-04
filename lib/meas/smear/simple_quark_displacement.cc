// $Id: simple_quark_displacement.cc,v 3.3 2008-11-04 18:43:58 edwards Exp $
/*! \file
 *  \brief Simple quark displacement
 */

#include "chromabase.h"

#include "meas/smear/quark_displacement_factory.h"
#include "meas/smear/simple_quark_displacement.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, SimpleQuarkDisplacementEnv::Params& param)
  {
    SimpleQuarkDisplacementEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const SimpleQuarkDisplacementEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace SimpleQuarkDisplacementEnv
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
      QuarkDisplacement<LatticeStaggeredPropagator>* createStagProp(XMLReader& xml_in,
								    const std::string& path)
      {
	return new QuarkDisplace<LatticeStaggeredPropagator>(Params(xml_in, path));
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
    
      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name = "SIMPLE_DISPLACEMENT";
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
	success &= Chroma::TheStagPropDisplacementFactory::Instance().registerObject(name, createStagProp);
	success &= Chroma::TheFermDisplacementFactory::Instance().registerObject(name, createFerm);
	success &= Chroma::TheColorVecDisplacementFactory::Instance().registerObject(name, createColorVec);
	registered = true;
      }
      return success;
    }


    //! Parameters for running code
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "disp_length", disp_length);
      read(paramtop, "disp_dir", disp_dir);
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      write(xml, "disp_length", disp_length);
      write(xml, "disp_dir", disp_dir);

      pop(xml);
    }

  }  // end namespace
}  // end namespace Chroma

