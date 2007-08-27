// $Id: sf_gauge_init.cc,v 3.2 2007-08-27 20:06:39 uid3790 Exp $
/*! \file
 *  \brief Initialize a Schroedinger BC config
 */

#include "handle.h"
#include "util/gauge/sf_gauge_init.h"
#include "util/gauge/gauge_init_factory.h"
#include "util/gauge/gauge_init_aggregate.h"

#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"
#include "actions/gauge/gaugebcs/schroedinger_gaugebc.h"

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
    GaugeInit* createCfg(XMLReader& xml_in,
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
	success &= CreateGaugeStateEnv::registerAll();
	success &= Chroma::TheGaugeInitFactory::Instance().registerObject(name, createCfg);
	registered = true;
      }
      return success;
    }


    // Parameters for running code
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      cgs = readXMLGroup(paramtop, "GaugeState", "Name");
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      int version = 1;
      write(xml, "cfg_type", SFGaugeInitEnv::name);
      xml << cgs.xml;

      pop(xml);
    }


    // Initialize the gauge field
    void
    GaugeIniter::operator()(XMLReader& gauge_file_xml,
			    XMLReader& gauge_xml,
			    multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Starting up a classical Schroedinger functional config" << endl;

      try
      {
	//
	// Create the GaugeState object
	//
	std::istringstream  xml_s(params.cgs.xml);
	XMLReader  top(xml_s);
        QDPIO::cout << "GaugeState type = " << params.cgs.id << endl;
	
	Handle< CreateGaugeState< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	  cgs(TheCreateGaugeStateFactory::Instance().createObject(params.cgs.id,
								  top,
								  params.cgs.path));

	// Need to downcast to the appropriate BC
	const SchrGaugeBC& gaugebc = dynamic_cast<const SchrGaugeBC&>(cgs->getBC());

	// Set the fields to the appropriate background field
	u = gaugebc.SFBndFld();

	// Munge the fields with the state
	Handle<GaugeState< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	  state((*cgs)(u));

	// Pull the u fields back out from the state since they might have been
	// munged with gaugeBC's
	u = state->getLinks();
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception in creating GaugeState: " << e << endl;
	QDP_abort(1);
      }

      XMLBufferWriter file_xml, record_xml;
      push(file_xml, "gauge");
      write(file_xml, "id", int(0));
      pop(file_xml);
      params.writeXML(record_xml, "SF_classical");

      gauge_file_xml.open(file_xml);
      gauge_xml.open(record_xml);
    }
  }
}
