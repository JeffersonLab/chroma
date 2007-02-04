// $Id: sf_gauge_init.cc,v 3.1 2007-02-04 22:06:42 edwards Exp $
/*! \file
 *  \brief Initialize a Schroedinger BC config
 */

#include "util/gauge/gauge_init_factory.h"
#include "util/gauge/gauge_init_aggregate.h"

#include "util/gauge/sf_gauge_init.h"
#include "actions/gauge/gaugebcs/schr_nonpert_gaugebc.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const string& path, SFGaugeInitEnv::Params& param)
  {
    SFGaugeInitEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const SFGaugeInitEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace SFGaugeInitEnv
  {
    //! Callback function
    GaugeInit* createSource(XMLReader& xml_in,
			    const std::string& path)
    {
      return new GaugeIniter(Params(xml_in, path));
    }

    //! Name to be used
    const std::string name = "CLASSICAL_SF";

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

      read(paramtop, "SFParams", sf);
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      int version = 1;
      write(xml, "cfg_type", SFGaugeInitEnv::name);
      write(xml, "SFParams", sf);

      pop(xml);
    }


    // Initialize the gauge field
    void
    GaugeIniter::operator()(XMLReader& gauge_file_xml,
			    XMLReader& gauge_xml,
			    multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Starting up a classical Schroedinger functional config" << endl;
      SchrNonPertGaugeBC gaugebc(params.sf);
      
      u = gaugebc.SFBndFld();

      XMLBufferWriter file_xml, record_xml;
      push(file_xml, "gauge");
      write(file_xml, "id", int(0));
      pop(file_xml);
      write(record_xml, "SF_classical", params.sf);

      gauge_file_xml.open(file_xml);
      gauge_xml.open(record_xml);
    }
  }
}
