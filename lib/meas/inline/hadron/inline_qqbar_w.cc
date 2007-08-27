// $Id: inline_qqbar_w.cc,v 3.5 2007-08-27 21:04:02 edwards Exp $
/*! \file
 * \brief Inline construction of qqbar
 *
 * QQbar calcs
 */

#include "meas/inline/hadron/inline_qqbar_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/hadron/mescomp_w.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/info/unique_id.h"
#include "util/ferm/diractodr.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineQQbarEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineQQbar(InlineQQbarParams(xml_in, path));
      }
      
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "QQBAR";

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
  }


  //! Param input
  void read(XMLReader& xml, const string& path, InlineQQbarParams::Param_t& input)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 4:
      break;

    default:
      /**************************************************************************/
      QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "Dirac_basis", input.Dirac_basis);
  }

  //! Param output
  void write(XMLWriter& xml, const string& path, const InlineQQbarParams::Param_t& input)
  {
    push(xml, path);

    int version = 4;
    write(xml, "version", version);
    write(xml, "Dirac_basis", input.Dirac_basis);

    pop(xml);
  }

  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineQQbarParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "prop_ids", input.prop_ids);
    read(inputtop, "qqbar_file", input.qqbar_file);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineQQbarParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "prop_ids", input.prop_ids);
    write(xml, "qqbar_file", input.qqbar_file);

    pop(xml);
  }


  // Param stuff
  InlineQQbarParams::InlineQQbarParams() { frequency = 0; }

  InlineQQbarParams::InlineQQbarParams(XMLReader& xml_in, const std::string& path) 
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

      // Ids
      read(paramtop, "NamedObject", named_obj);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineQQbarParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "Param", param);

    // Ids
    Chroma::write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }


  void 
  InlineQQbar::operator()(unsigned long update_no,
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
      QDPIO::cerr << InlineQQbarEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineQQbarEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out,"qqbar");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineQQbarEnv::name << ": generalized propagator generation" << endl;

    // Write out the input
    params.write(xml_out, "Input");

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    // Check to make sure there are 2 files
    const int Nprops = 2;
    if (params.named_obj.prop_ids.size() != Nprops)
    {
      QDPIO::cerr << InlineQQbarEnv::name << ": Error on input params - expecting 2 filenames" << endl;
      QDP_abort(1);
    }

    write(xml_out, "propIds", params.named_obj.prop_ids);
    write(xml_out, "mescompFile", params.named_obj.qqbar_file);

    /*
     * Read the quark propagators and extract headers
     */
    multi1d<LatticePropagator> quark_propagator(Nprops);
    QQbarMescomp_t  qqbar;
    qqbar.Dirac_basis = false;
    qqbar.forward_props.resize(Nprops);
    for(int i=0; i < Nprops; ++i)
    {
      try
      {
	// Snarf the data into a copy
	quark_propagator[i] =
	  TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_ids[i]);

	// Snarf the prop info. This is will throw if the prop_id is not there
	XMLReader prop_file_xml, prop_record_xml;
	TheNamedObjMap::Instance().get(params.named_obj.prop_ids[i]).getFileXML(prop_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.prop_ids[i]).getRecordXML(prop_record_xml);
   
	// Try to invert this record XML into a ChromaProp struct
	// Also pull out the id of this source
	{
	  read(prop_record_xml, "/SinkSmear/PropSink", qqbar.forward_props[i].sink_header);
	  read(prop_record_xml, "/SinkSmear/ForwardProp", qqbar.forward_props[i].prop_header);
	  read(prop_record_xml, "/SinkSmear/PropSource", qqbar.forward_props[i].source_header);
	}
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << InlineQQbarEnv::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << InlineQQbarEnv::name << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }
    }


    // Save prop input
    write(xml_out, "Propagator_input", qqbar);

    // Derived from input prop
    int t0                = qqbar.forward_props[0].source_header.t_source;
    int j_decay           = qqbar.forward_props[0].source_header.j_decay;

    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, j_decay);

    // Sanity check - write out the propagator (pion) correlator in the j_decay direction
    for(int i=0; i < Nprops; ++i)
    {
      multi1d<Double> prop_corr = sumMulti(localNorm2(quark_propagator[i]), 
					   phases.getSet());

      push(xml_out, "SinkSmearedProp_correlator");
      write(xml_out, "correlator_num", i);
      write(xml_out, "sink_smeared_prop_corr", prop_corr);
      pop(xml_out);
    }

    /*
     * Generalized propagator calculation
     */
    multiNd<Complex> mesprop;

    // Switch to Dirac-basis if desired.
    if (params.param.Dirac_basis)
    {
      qqbar.Dirac_basis = true;

      // The spin basis matrix
      SpinMatrix U = DiracToDRMat();

      LatticePropagator q_tmp;
      for(int i=0; i < Nprops; ++i)
      {
	q_tmp = adj(U) * quark_propagator[i] * U;   // DeGrand-Rossi ---> Dirac
	quark_propagator[i] = q_tmp;
      }
    }

    // Compute generation propagator
    mescomp(mesprop,
	    quark_propagator[0],
	    quark_propagator[1],
	    phases, t0);

    // Convert the data into a mult1d
    multi1d<Complex> mesprop_1d;
    convertMescomp(mesprop_1d, mesprop, j_decay);

    // Save the qqbar output
    // ONLY SciDAC output format is supported!
    {
      XMLBufferWriter file_xml;
      push(file_xml, "qqbar");
      write(file_xml, "id", uniqueId());  // NOTE: new ID form
      pop(file_xml);

      XMLBufferWriter record_xml;
      push(record_xml, "QQbar");
      write(record_xml, ".", qqbar);  // do not write the outer group
      write(record_xml, "Config_info", gauge_xml);
      pop(record_xml);  // QQbar

      // Write the scalar data
      QDPFileWriter to(file_xml, params.named_obj.qqbar_file, 
		       QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);
      write(to,record_xml,mesprop_1d);
      close(to);
    }

    pop(xml_out);    // qqbar


    snoop.stop();
    QDPIO::cout << InlineQQbarEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineQQbarEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

};
