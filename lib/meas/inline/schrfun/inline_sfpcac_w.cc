// $Id: inline_sfpcac_w.cc,v 1.6 2008-02-06 18:55:18 edwards Exp $
/*! \file
 * \brief Inline Schroedinger functional measurements
 */

#include "fermact.h"
#include "meas/inline/schrfun/inline_sfpcac_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

#include "meas/schrfun/sfpcac_w.h"

namespace Chroma 
{ 
  //! SFpcac input
  void read(XMLReader& xml, const string& path, InlineSFpcacEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
  }

  //! SFpcac output
  void write(XMLWriter& xml, const string& path, const InlineSFpcacEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);

    pop(xml);
  }


  //! SFpcac input
  void read(XMLReader& xml, const string& path, InlineSFpcacEnv::Params::SFpcac_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "decay_dir", input.decay_dir);
    read(inputtop, "ZVfactP", input.ZVfactP);
    read(inputtop, "ZAfactP", input.ZAfactP);
    read(inputtop, "x0", input.x0);
    read(inputtop, "y0", input.y0);
  }

  //! SFpcac output
  void write(XMLWriter& xml, const string& path, const InlineSFpcacEnv::Params::SFpcac_t& input)
  {
    push(xml, path);

    write(xml, "decay_dir", input.decay_dir);
    write(xml, "ZVfactP", input.ZVfactP);
    write(xml, "ZAfactP", input.ZAfactP);
    write(xml, "x0", input.x0);
    write(xml, "y0", input.y0);

    pop(xml);
  }


  namespace InlineSFpcacEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "SCHROEDINGER_FUNCTIONAL_PCAC";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= WilsonTypeFermActsEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    // Param stuff
    Params::Params() { frequency = 0; }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Parameters for propagator solver
	read(paramtop, "Param", param);

	// Parameters for Schroedinger functional
	read(paramtop, "SFpcac", sfpcac);

	// Read in the output propagator/source configuration info
	read(paramtop, "NamedObject", named_obj);

	// Possible alternate XML file pattern
	if (paramtop.count("xml_file") != 0) 
	{
	  read(paramtop, "xml_file", xml_file);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }


    void
    Params::writeXML(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
    
      write(xml_out, "Param", param);
      write(xml_out, "SFpcac", sfpcac);
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }


    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "SFpcac");
	write(xml_out, "update_no", update_no);
	write(xml_out, "xml_file", xml_file);
	pop(xml_out);

	XMLFileWriter xml(xml_file);
	func(update_no, xml);
      }
      else
      {
	func(update_no, xml_out);
      }
    }


    // Real work done here
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      QDPIO::cout << name << ": Schroedinger functional calculation" << endl;

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      // Test and grab a reference to the gauge field
      multi1d<LatticeColorMatrix> u;
      XMLBufferWriter gauge_xml;
      try
      {
	u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }

      push(xml_out, "SFpcac");
      write(xml_out, "update_no", update_no);

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      params.writeXML(xml_out, "Input");

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);

      //
      // Initialize fermion action
      //
      std::istringstream  xml_s(params.param.fermact.xml);
      XMLReader  fermacttop(xml_s);
      QDPIO::cout << "FermAct = " << params.param.fermact.id << endl;

 
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, params.sfpcac.decay_dir);

      //
      // Try the factories
      //
      try
      {
	StopWatch swatch;
	swatch.reset();
	QDPIO::cout << "Try the various factories" << endl;

	// Typedefs to save typing
	typedef LatticeFermion               T;
	typedef multi1d<LatticeColorMatrix>  P;
	typedef multi1d<LatticeColorMatrix>  Q;

	// Generic Wilson-Type stuff
	Handle< WilsonTypeFermAct<T,P,Q> >
	  S_f(TheWilsonTypeFermActFactory::Instance().createObject(params.param.fermact.id,
								   fermacttop,
								   params.param.fermact.path));


	Handle< FermState<T,P,Q> > state(S_f->createState(u));

	QDPIO::cout << "Suitable factory found: do the measurements" << endl;
	Handle< SystemSolver<T> > qprop(S_f->qprop(state, params.param.invParam));

	swatch.start();

	// SF measurements
	SFpcac(qprop, state, phases, 
	       params.sfpcac.ZVfactP,
	       params.sfpcac.ZAfactP,
	       params.sfpcac.x0, params.sfpcac.y0,
	       xml_out, string("sfpcac"));

	swatch.stop();	
QDPIO::cout << "SFpcac computed: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << endl;
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const std::string& e) 
      {
	QDPIO::cout << name << ": caught exception with fermion action: " << e << endl;
      }


      pop(xml_out);  // sfpcac

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } 

  } // namespace InlineSFpcacEnv

} // namespace Chroma
