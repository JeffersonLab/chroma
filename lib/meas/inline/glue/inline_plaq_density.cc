/*! \file
 *  \brief Inline plaquette density
 */

#include "meas/inline/glue/inline_plaq_density.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"

#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma 
{ 

  namespace InlinePlaqDenEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path)
      {
	Params p(xml_in, path);
	return new InlineMeas(p);
      }

      //! Local registration flag
      bool registered = false;

      const std::string name = "PLAQ_DENSITY";
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= CreateGaugeStateEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
 
  
    namespace
    {
      //! Compute the site-level action
      void siteAction(multi1d<LatticeReal>& site_act, const multi1d<LatticeColorMatrix>& u)
      {
	START_CODE();

	// Initialize
	site_act.resize(Nd*(Nd-1)/2);

	Real one = 1.0;
	Real third = Real(1) / Real(Nc);

	// Compute the average plaquettes
	int cnt = 0;
	for(int mu=1; mu < Nd; ++mu)
	{
	  for(int nu=0; nu < mu; ++nu)
	  {
	    /* tmp_0 = u(x+mu,nu)*u_dag(x+nu,mu) */
	    /* tmp_1 = tmp_0*u_dag(x,nu)=u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
	    /* wplaq_tmp = tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)) */
	    site_act[cnt++] = one - third*real(trace(u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu])));
	  }
	}

	END_CODE();
      }
 
    }
    

    //! Parameter input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	break;

      default:
	QDPIO::cerr << "Params::Param_t: " << version 
		    << " unsupported." << std::endl;
	QDP_abort(1);
      }
    }

    //! Parameter output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& param)
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      pop(xml);
    }


    //! Parameter input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "out_file", input.out_file);
    }

    //! Parameter output
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "out_file", input.out_file);

      pop(xml);
    }


    // Params
    Params::Params()
    { 
      frequency = 0; 
    }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Params
	read(paramtop, "Param", param);

	// Ids
	read(paramtop, "NamedObject", named_obj);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << "Caught Exception reading XML: " << e << std::endl;
	QDP_abort(1);
      }
    }

    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      // Grab the object
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      push(xml_out, "PlaqDen");
      write(xml_out, "update_no", update_no);

      multi1d<LatticeReal> plaq_site;

      siteAction(plaq_site, u);

      QDPIO::cout << name << ": write plaquette density to xml file = " << params.named_obj.out_file << std::endl;
      XMLFileWriter txt(params.named_obj.out_file);

      write(txt, "plaq_site", plaq_site);

      pop(xml_out);

      END_CODE();
    } 

  }

}
