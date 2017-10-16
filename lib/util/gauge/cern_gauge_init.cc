/*! \file
 *  \brief Read a CERN config
 */

#include "util/gauge/gauge_init_factory.h"
#include "util/gauge/gauge_init_aggregate.h"

#include "util/gauge/cern_gauge_init.h"
#include "qdp_iogauge.h"

#include "io/cern_io.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const std::string& path, CERNGaugeInitEnv::Params& param)
  {
    CERNGaugeInitEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const std::string& path, const CERNGaugeInitEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace CERNGaugeInitEnv
  {
    //! Callback function
    GaugeInit* createSource(XMLReader& xml_in,
			    const std::string& path)
    {
      return new GaugeIniter(Params(xml_in, path));
    }

    //! Name to be used
    const std::string name = "CERN";

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
    Params::Params(XMLReader& xml, const std::string& path)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "cfg_file", cfg_file);
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const std::string& path) const
    {
      push(xml, path);
    
      int version = 1;
      write(xml, "cfg_type", CERNGaugeInitEnv::name);
      write(xml, "cfg_file", cfg_file);

      pop(xml);
    }


    // Initialize the gauge field
    void
    GaugeIniter::operator()(XMLReader& gauge_file_xml,
			    XMLReader& gauge_xml,
			    multi1d<LatticeColorMatrix>& u) const
    {
      u.resize(Nd);
      readCERN(u, params.cfg_file);

      XMLBufferWriter file_xml, record_xml;
      push(file_xml, "gauge");
      write(file_xml, "info", "A CERN Gauge Field");
      pop(file_xml);

      push(record_xml, "record");
      write(record_xml, "info", "Spooled from CERN Format File");
      pop(file_xml);

      gauge_file_xml.open(file_xml);
      gauge_xml.open(record_xml);

    }
  }
}
