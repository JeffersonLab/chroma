// $Id: inline_diquark_w.cc,v 1.3 2007-02-25 22:39:28 edwards Exp $
/*! \file
 * \brief Inline construction of the diquark within a QQQ
 *
 * Diquarks for QQQ calcs
 */

#include "meas/inline/hadron/inline_diquark_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/hadron/diquark_w.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/info/unique_id.h"
#include "util/ferm/diractodr.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineDiquarkEnv 
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

    const std::string name = "DIQUARK";

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


    //! Param input
    void read(XMLReader& xml, const string& path, InlineDiquarkEnv::Params::Param_t& input)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	break;

      default:
	QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "Dirac_basis", input.Dirac_basis);
    }

    //! Param output
    void write(XMLWriter& xml, const string& path, const InlineDiquarkEnv::Params::Param_t& input)
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "Dirac_basis", input.Dirac_basis);

      pop(xml);
    }

    //! Propagator input
    void read(XMLReader& xml, const string& path, InlineDiquarkEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "prop_ids", input.prop_ids);
      read(inputtop, "diquark_id", input.diquark_id);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlineDiquarkEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "prop_ids", input.prop_ids);
      write(xml, "diquark_id", input.diquark_id);

      pop(xml);
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

	// Parameters for source construction
	read(paramtop, "Param", param);

	// Read in the output propagator/source configuration info
	read(paramtop, "NamedObject", named_obj);
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
    
      // Parameters for source construction
      write(xml_out, "Param", param);

      // Write out the output propagator/source configuration info
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }


    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      // Test and grab a reference to the gauge field
      XMLBufferWriter gauge_xml;
      try
      {
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << InlineDiquarkEnv::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << InlineDiquarkEnv::name << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      push(xml_out,"diquark");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << InlineDiquarkEnv::name << ": Generalized propagator generation" << endl;
      StopWatch swatch;

      // Write out the input
      params.writeXML(xml_out, "Input");

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      params.writeXML(xml_out, "Input");

      // Write out the config info
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);

      // Check to make sure there are 2 ids
      const int Nprops = 2;
      if (params.named_obj.prop_ids.size() != Nprops)
      {
	QDPIO::cerr << "Error on input params - expecting 3 buffers" << endl;
	QDP_abort(1);
      }

      write(xml_out, "propIds", params.named_obj.prop_ids);
      write(xml_out, "diquarkId", params.named_obj.diquark_id);

      /*
       * Read the quark propagators and extract headers
       */
      multi1d<LatticePropagator> quark_propagator(Nprops);
      QQDiquark_t  diquark_header;
      diquark_header.Dirac_basis = false;
      diquark_header.forward_props.resize(Nprops);

      for(int i=0; i < Nprops; ++i)
      {
	QDPIO::cout << InlineDiquarkEnv::name << ": parse id = " << params.named_obj.prop_ids[i] << endl;
	try
	{
	  // Snarf the data into a copy
	  quark_propagator[i] =
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_ids[i]);

	  // Snarf the prop info. This is will throw if the prop_ids is not there.
	  XMLReader prop_file_xml, prop_record_xml;
	  TheNamedObjMap::Instance().get(params.named_obj.prop_ids[i]).getFileXML(prop_file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.prop_ids[i]).getRecordXML(prop_record_xml);
	
	  // Try to invert this record XML into a ChromaProp struct
	  // Also pull out the id of this source
	  {
	    read(prop_record_xml, "/SinkSmear", diquark_header.forward_props[i]);
	  }
	}
	catch( std::bad_cast ) 
	{
	  QDPIO::cerr << InlineDiquarkEnv::name << ": caught dynamic cast error" 
		      << endl;
	  QDP_abort(1);
	}
	catch (const string& e) 
	{
	  QDPIO::cerr << InlineDiquarkEnv::name << ": error message: " << e 
		      << endl;
	  QDP_abort(1);
	}
	QDPIO::cout << InlineDiquarkEnv::name << ": object successfully parsed" << endl;
      }


      // Save prop input
      write(xml_out, "Propagator_input", diquark_header);

      // Derived from input prop
      int j_decay = diquark_header.forward_props[0].source_header.j_decay;

      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, j_decay);

      // Sanity check - write out the propagator (pion) correlator in the j_decay direction
      push(xml_out, "SinkSmearedProp_correlator");
      for(int i=0; i < Nprops; ++i)
      {
	multi1d<Double> prop_corr = sumMulti(localNorm2(TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_ids[i])), 
					     phases.getSet());

	push(xml_out, "elem");
	write(xml_out, "correlator_num", i);
	write(xml_out, "sink_smeared_prop_corr", prop_corr);
	pop(xml_out);
      }
      pop(xml_out);

      /*
       * Generalized propagator calculation
       */
      // Switch to Dirac-basis if desired.
      if (params.param.Dirac_basis)
      {
	diquark_header.Dirac_basis = true;

	// The spin basis matrix
	SpinMatrix U = DiracToDRMat();

	LatticePropagator q_tmp;
	for(int i=0; i < Nprops; ++i)
	{
	  q_tmp = adj(U) * quark_propagator[i] * U;   // DeGrand-Rossi ---> Dirac
	  quark_propagator[i] = q_tmp;
	}
      }

      // Compute diquark
      QQDiquarkContract_t diquark;

      QQDiquark(diquark,
		quark_propagator[0],
		quark_propagator[1]);

      // Now save the diquark object
      try
      {
	QDPIO::cout << "Attempt to update diquark object" << endl;

	XMLBufferWriter file_xml;
	push(file_xml, "diquark");
	write(file_xml, "id", uniqueId());  // NOTE: new ID form
	pop(file_xml);

	XMLBufferWriter record_xml;
	write(record_xml, "Diquark", diquark_header);
    
	// Store the source
	TheNamedObjMap::Instance().create<QQDiquarkContract_t>(params.named_obj.diquark_id);
	TheNamedObjMap::Instance().getData<QQDiquarkContract_t>(params.named_obj.diquark_id) = 
	  diquark;
	TheNamedObjMap::Instance().get(params.named_obj.diquark_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.diquark_id).setRecordXML(record_xml);

	QDPIO::cout << "Diquark successfully update" << endl;
      }
      catch (std::bad_cast)
      {
	QDPIO::cerr << InlineDiquarkEnv::name << ": dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << InlineDiquarkEnv::name << ": error message: " << e << endl;
	QDP_abort(1);
      }

      pop(xml_out);    // diquark


      snoop.stop();
      QDPIO::cout << InlineDiquarkEnv::name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << InlineDiquarkEnv::name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

}
