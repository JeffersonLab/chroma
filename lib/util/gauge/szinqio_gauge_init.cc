/*! \file
 *  \brief Read a SZINQIO config
 */

#include "util/gauge/gauge_init_factory.h"
#include "util/gauge/gauge_init_aggregate.h"

#include "util/gauge/szinqio_gauge_init.h"
#include "io/gauge_io.h"
#include "util/gauge/reunit.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const std::string& path, SZINQIOGaugeInitEnv::Params& param)
  {
    SZINQIOGaugeInitEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const std::string& path, const SZINQIOGaugeInitEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace SZINQIOGaugeInitEnv
  {
    //! Callback function
    GaugeInit* createSource(XMLReader& xml_in,
			    const std::string& path)
    {
      return new GaugeIniter(Params(xml_in, path));
    }

    //! Name to be used
    const std::string name = "SZINQIO";
    const std::string alternate_name = "SCIDAC";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheGaugeInitFactory::Instance().registerObject(name, createSource);
	success &= Chroma::TheGaugeInitFactory::Instance().registerObject(alternate_name, createSource);
	registered = true;
      }
      return success;
    }


    // Parameters for running code
    Params::Params(XMLReader& xml, const std::string& path)
    {
    	XMLReader paramtop(xml, path);

    	read(paramtop, "cfg_file", cfg_file);
    	// Default
    	cfg_pario = QDPIO_SERIAL;

    	// If the IO Node grid has more than 1 node then do parallel io
    	// as a default
    	bool pario = Layout::isIOGridDefined() && ( Layout::numIONodeGrid() > 1 );

    	// Attempt to read option with either tag
    	if ( paramtop.count("parallel_io") > 0 ) {
    		read(paramtop, "parallel_io", pario);
    	}
    	else {
    		if ( paramtop.count("ParallelIO") > 0 ) {
    			read(paramtop, "ParallelIO", pario);
    		}
    	}

	if ( paramtop.count("reunit") > 0 ) {
	     read(paramtop, "reunit", reunitP );
        }
        else { 
             reunitP = false;
        }
 
    	if ( pario )  {
    		cfg_pario = QDPIO_PARALLEL;
    	}
    	else {
    		cfg_pario = QDPIO_SERIAL;
    	}


    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const std::string& path) const
    {
      push(xml, path);
    
      int version = 1;
      write(xml, "cfg_type", SZINQIOGaugeInitEnv::name);
      write(xml, "cfg_file", cfg_file);

      bool pario = false;

      if ( cfg_pario == QDPIO_PARALLEL ) { 
    	  pario = true;
      }
      write(xml, "parallel_io", pario);
      write(xml, "reunit", reunitP);

      pop(xml);
    }


    //! Returns a link smearing group with these params
    GroupXML_t   createXMLGroup(const Params& p)
    {
      GroupXML_t foo;

      XMLBufferWriter xml_tmp;
      write(xml_tmp, "Cfg", p);
      foo.xml = xml_tmp.printCurrentContext();
      foo.id = name;
      foo.path = "/Cfg";

      return foo;
    }


    // Initialize the gauge field
    void
    GaugeIniter::operator()(XMLReader& gauge_file_xml,
			    XMLReader& gauge_xml,
			    multi1d<LatticeColorMatrix>& u) const
    {
      u.resize(Nd);
      if( params.cfg_pario == QDPIO_PARALLEL ) { 
	QDPIO::cout << "Parallel IO read" << std::endl;
      }
      readGauge(gauge_file_xml, gauge_xml, u, params.cfg_file, params.cfg_pario);
      if( params.reunitP == true ) { 
	QDPIO::cout << "Reunitarizing After read" << std::endl;
	int numbad=0;
	for(int mu=0; mu < Nd; ++mu) {
	  int numbad_mu=0;
	  reunit(u[mu], numbad_mu, REUNITARIZE_LABEL);
	  numbad += numbad_mu;
	} 
	QDPIO::cout << "Reunitarize reported " << numbad << " unitarity violations" << std::endl;
      }
    }
  }
}
