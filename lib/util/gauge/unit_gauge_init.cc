// $Id: unit_gauge_init.cc,v 3.1 2007-02-04 22:06:42 edwards Exp $
/*! \file
 *  \brief Create a unit config
 */

#include "util/gauge/gauge_init_factory.h"
#include "util/gauge/gauge_init_aggregate.h"

#include "util/gauge/unit_gauge_init.h"

namespace Chroma
{

  //! Hooks to register the class
  namespace UnitGaugeInitEnv
  {
    //! Callback function
    GaugeInit* createSource(XMLReader& xml_in,
			    const std::string& path)
    {
      return new GaugeIniter(Params(xml_in, path));
    }

    //! Name to be used
    const std::string name = "UNIT";

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
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      int version = 1;
      write(xml, "cfg_type", UnitGaugeInitEnv::name);

      pop(xml);
    }


    // Initialize the gauge field
    void
    GaugeIniter::operator()(XMLReader& gauge_file_xml,
			    XMLReader& gauge_xml,
			    multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Starting up unit gauge (free) config" << endl;
      u.resize(Nd);
      u = 1;

      XMLBufferWriter file_xml, record_xml;
      push(file_xml, "gauge");
      write(file_xml, "id", int(0));
      pop(file_xml);
      push(record_xml, "unit");
      pop(record_xml);

      gauge_file_xml.open(file_xml);
      gauge_xml.open(record_xml);
    }
  }
}
