// $Id: inline_apply_fermstate_w.cc,v 3.2 2007-08-27 20:00:23 uid3790 Exp $
/*! \file
 *  \brief Inline ferm state application
 */

#include "meas/inline/hadron/inline_apply_fermstate_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/inline/io/named_objmap.h"

#include "meas/inline/io/default_gauge_field.h"

#include "actions/ferm/fermstates/ferm_createstate_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"

namespace Chroma 
{ 

  //! FermState input
  void read(XMLReader& xml, const string& path, InlineFermStateEnv::Params::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 1:
      param.cfs = readXMLGroup(paramtop, "FermState", "Name");
      break;

    default:
      QDPIO::cerr << "InlineFermStateEnv::Params::Param_t: " << version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }
  }

  //! FermState output
  void write(XMLWriter& xml, const string& path, const InlineFermStateEnv::Params::Param_t& param)
  {
    push(xml, path);

    int version = 2;
    write(xml, "version", version);
    xml << param.cfs.xml;

    pop(xml);
  }


  //! FermState input
  void read(XMLReader& xml, const string& path, InlineFermStateEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "output_id", input.output_id);
  }

  //! FermState output
  void write(XMLWriter& xml, const string& path, const InlineFermStateEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "output_id", input.output_id);

    pop(xml);
  }


  namespace InlineFermStateEnv 
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

    const std::string name = "APPLY_FERM_STATE";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
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
	QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }


    void
    Params::writeXML(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
    
      write(xml_out, "Param", param);
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }



    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();
    
      StopWatch snoop;
      snoop.reset();
      snoop.start();

      push(xml_out, "apply_ferm_state");
      write(xml_out, "update_no", update_no);
    
      // Write out the input
      params.writeXML(xml_out, "Input");

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Grab the gauge field
      multi1d<LatticeColorMatrix> u;
      XMLBufferWriter gauge_xml;

      try
      {
	u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << InlineFermStateEnv::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << InlineFermStateEnv::name << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      // Calculate some gauge invariant observables
      MesPlq(xml_out, "Observables", u);

      // Apply the fermstate
      try
      {
	// Set the construct state and modify the fields
	{
	  std::istringstream  xml_s(params.param.cfs.xml);
	  XMLReader  fermtop(xml_s);

	  Handle<CreateFermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	    cfs(TheCreateFermStateFactory::Instance().createObject(params.param.cfs.id,
								   fermtop,
								   params.param.cfs.path));

	  Handle<FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	    state((*cfs)(u));

	  // Pull the u fields back out from the state since they might have been
	  // munged with fermBC's
	  u = state->getLinks();
	}
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << InlineFermStateEnv::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << InlineFermStateEnv::name << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }

      // Again calculate some ferm invariant_observables
      MesPlq(xml_out, "Link_observables", u);

      // Now store the configuration to a memory object
      {
	XMLBufferWriter file_xml, record_xml;
	push(file_xml, "gauge");
	write(file_xml, "id", int(0));
	pop(file_xml);
	record_xml << gauge_xml;

	// Store the gauge field
	TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(params.named_obj.output_id);
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.output_id) = u;
	TheNamedObjMap::Instance().get(params.named_obj.output_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.output_id).setRecordXML(record_xml);
      }

      pop(xml_out);  // apply_ferm_state

      snoop.stop();
      QDPIO::cout << InlineFermStateEnv::name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << InlineFermStateEnv::name << ": ran successfully" << endl;

      END_CODE();
    } 

  } // namespace InlineFermStateEnv

}  // namespace Chroma
