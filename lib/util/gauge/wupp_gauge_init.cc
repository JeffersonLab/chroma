// $Id: milc_gauge_init.cc,v 3.3 2009-05-13 13:12:19 edwards Exp $
/*! \file
 *  \brief Read an expanded BMW config
 */

#include "util/gauge/gauge_init_factory.h"
#include "util/gauge/gauge_init_aggregate.h"

#include "util/gauge/wupp_gauge_init.h"
#include "io/readwupp.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const string& path, WUPPGaugeInitEnv::Params& param)
  {
    WUPPGaugeInitEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const WUPPGaugeInitEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace WUPPGaugeInitEnv
  {
    //! Callback function
    GaugeInit* createSource(XMLReader& xml_in,
			    const std::string& path)
    {
      return new GaugeIniter(Params(xml_in, path));
    }

    //! Name to be used
    const std::string name = "WUPP";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheGaugeInitFactory::Instance().registerObject(name, createSource);
	registered = true;
      }
      return success;
    }


    // Parameters for running code
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "cfg_file", cfg_file);
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      int version = 1;
      write(xml, "cfg_type", WUPPGaugeInitEnv::name);
      write(xml, "cfg_file", cfg_file);

      pop(xml);
    }


    // Initialize the gauge field
    void
    GaugeIniter::operator()(XMLReader& gauge_file_xml,
			    XMLReader& gauge_xml,
			    multi1d<LatticeColorMatrix>& u) const
    {
      readWupp(gauge_xml, u, params.cfg_file);
    }
  }
}
