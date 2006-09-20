// $Id: no_link_smearing.cc,v 3.2 2006-09-20 20:28:04 edwards Exp $
/*! \file
 *  \brief No link smearing
 */

#include "chromabase.h"

#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/no_link_smearing.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const string& path, NoLinkSmearingEnv::Params& param)
  {
    NoLinkSmearingEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const NoLinkSmearingEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace NoLinkSmearingEnv
  {
    //! Callback function
    LinkSmearing* createSource(XMLReader& xml_in,
			       const std::string& path)
    {
      return new LinkSmear(Params(xml_in, path));
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
	success &= Chroma::TheLinkSmearingFactory::Instance().registerObject(name, createSource);
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

      int version = 1;
      write(xml, "version", version);
      write(xml, "LinkSmearingType", name);

      pop(xml);
    }


    //! Smear the links
    void
    LinkSmear::operator()(multi1d<LatticeColorMatrix>& u) const
    {
    }

  }
}
