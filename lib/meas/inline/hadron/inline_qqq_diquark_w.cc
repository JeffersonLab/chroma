// $Id: inline_qqq_diquark_w.cc,v 1.1 2007-02-25 22:39:03 edwards Exp $
/*! \file 
 * \brief Inline construction of QQQ's using a diquark
 *
 * QQQ calcs using a diquark
 */

#include "meas/inline/hadron/inline_qqq_diquark_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/hadron/barcomp_diquark_w.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/info/unique_id.h"
#include "util/ferm/diractodr.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineQQQDiquarkEnv 
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

    const std::string name = "QQQ_DIQUARK";

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
    void read(XMLReader& xml, const string& path, InlineQQQDiquarkEnv::Params::Param_t& input)
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

      read(paramtop, "sparseP", input.sparseP);
      read(paramtop, "Dirac_basis", input.Dirac_basis);
      if (input.sparseP)
      {
	read(paramtop, "SpinIndices", input.spin_indices);
      }
    }

    //! Param output
    void write(XMLWriter& xml, const string& path, const InlineQQQDiquarkEnv::Params::Param_t& input)
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "sparseP", input.sparseP);
      write(xml, "Dirac_basis", input.Dirac_basis);
      if (input.sparseP)
      {
	write(xml, "SpinIndices", input.spin_indices);
      }

      pop(xml);
    }

    //! Propagator input
    void read(XMLReader& xml, const string& path, InlineQQQDiquarkEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "diquark_id", input.diquark_id);
      read(inputtop, "prop_id", input.prop_id);
      read(inputtop, "qqq_file", input.qqq_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlineQQQDiquarkEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "diquark_id", input.diquark_id);
      write(xml, "prop_id", input.prop_id);
      write(xml, "qqq_file", input.qqq_file);

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
	QDPIO::cerr << InlineQQQDiquarkEnv::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << InlineQQQDiquarkEnv::name << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      push(xml_out,"qqq_diquark");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << InlineQQQDiquarkEnv::name << ": Generalized propagator generation" << endl;
      StopWatch swatch;

      // Write out the input
      params.writeXML(xml_out, "Input");

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      params.writeXML(xml_out, "Input");

      // Write out the config info
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 5);
      pop(xml_out);

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);

      write(xml_out, "diquarkId", params.named_obj.diquark_id);
      write(xml_out, "propId", params.named_obj.prop_id);
      write(xml_out, "barcompFile", params.named_obj.qqq_file);

      /*
       * Read the quark propagators and extract headers
       */
      const int Nprops = 3;

      LatticePropagator quark_propagator;

      QQDiquark_t   diquark_header;
      QQQBarcomp_t  qqq_header;
      qqq_header.Dirac_basis = false;
      qqq_header.forward_props.resize(Nprops);

      try
      {
	// Read diquark
	QDPIO::cout << InlineQQQDiquarkEnv::name << ": parse id = " << params.named_obj.diquark_id << endl;
	TheNamedObjMap::Instance().getData<QQDiquarkContract_t>(params.named_obj.diquark_id);
	{
	  XMLReader diquark_file_xml, diquark_record_xml;
	  TheNamedObjMap::Instance().get(params.named_obj.diquark_id).getFileXML(diquark_file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.diquark_id).getRecordXML(diquark_record_xml);
	
	  read(diquark_record_xml, "/Diquark", diquark_header);
	}
	QDPIO::cout << InlineQQQDiquarkEnv::name << ": object successfully parsed" << endl;

	// The QQQ object is composed of the diquark headers
	qqq_header.forward_props[0] = diquark_header.forward_props[0];
	qqq_header.forward_props[1] = diquark_header.forward_props[1];

	// Read 3rd prop
	QDPIO::cout << InlineQQQDiquarkEnv::name << ": parse id = " << params.named_obj.prop_id << endl;
	quark_propagator = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_id);
	{
	  XMLReader prop_file_xml, prop_record_xml;
	  TheNamedObjMap::Instance().get(params.named_obj.prop_id).getFileXML(prop_file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.prop_id).getRecordXML(prop_record_xml);
	
	  read(prop_record_xml, "/SinkSmear", qqq_header.forward_props[Nprops-1]);
	}
	QDPIO::cout << InlineQQQDiquarkEnv::name << ": object successfully parsed" << endl;
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << InlineQQQDiquarkEnv::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << InlineQQQDiquarkEnv::name << ": error message: " << e 
		    << endl;
	QDP_abort(1);
      }
      const QQDiquarkContract_t& diquark =
	TheNamedObjMap::Instance().getData<QQDiquarkContract_t>(params.named_obj.diquark_id);

      QDPIO::cout << InlineQQQDiquarkEnv::name << ": finished with gauge, diquark and prop" << endl;

      // Save prop input
      write(xml_out, "Propagator_input", qqq_header);

      // Derived from input prop
      multi1d<int> boundary = getFermActBoundary(qqq_header.forward_props[0].prop_header.fermact);
      int t0                = qqq_header.forward_props[0].source_header.t_source;
      int j_decay           = qqq_header.forward_props[0].source_header.j_decay;
      int bc_spec = boundary[j_decay];

      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, j_decay);

      // Sanity check - write out the propagator (pion) correlator in the j_decay direction
      {
	multi1d<Double> prop_corr = sumMulti(localNorm2(TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_id)), 
					     phases.getSet());

	push(xml_out, "SinkSmearedProp_correlator");
	write(xml_out, "sink_smeared_prop_corr", prop_corr);
	pop(xml_out);
      }

      /*
       * Generalized propagator calculation
       */
      multi1d<ComplexF> barprop_1d;

      // Switch to Dirac-basis if desired.
      if (params.param.Dirac_basis)
      {
	qqq_header.Dirac_basis = true;

	// The spin basis matrix
	SpinMatrix U = DiracToDRMat();

	LatticePropagator q_tmp = adj(U) * quark_propagator * U;   // DeGrand-Rossi ---> Dirac
	quark_propagator = q_tmp;
      }

      // Decide the computational mode depending on whether the data is sparse
      // or not
      if (params.param.sparseP)
      {
	qqq_header.sparseP = true;
	qqq_header.spin_indices = params.param.spin_indices;
	QQQSparse_t barprop;

	// Compute generation propagator
	barcompDiquarkSparse(barprop,
			     diquark,
			     quark_propagator,
			     params.param.spin_indices,
			     phases, t0, bc_spec);

	// Convert the data into a mult1d
	barprop_1d = barprop.serialize();
      }
      else
      {
	qqq_header.sparseP = false;
	QQQDense_t barprop;

	// Compute generation propagator
	barcompDiquarkDense(barprop,
			    diquark,
			    quark_propagator,
			    phases, t0, bc_spec);

	// Convert the data into a mult1d
	barprop_1d = barprop.serialize();
      }

      // Save the qqq output
      // ONLY SciDAC output format is supported!
      {
	XMLBufferWriter file_xml;
	push(file_xml, "qqq");
	write(file_xml, "id", uniqueId());  // NOTE: new ID form
	pop(file_xml);

	XMLBufferWriter record_xml;
	push(record_xml, "QQQ");
	write(record_xml, ".", qqq_header);  // do not write the outer group
	write(record_xml, "Config_info", gauge_xml);
	pop(record_xml);  // QQQ

	// Write the scalar data
	QDPFileWriter to(file_xml, params.named_obj.qqq_file, 
			 QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);
	write(to,record_xml,barprop_1d);
	close(to);
      }

      pop(xml_out);    // qqq


      snoop.stop();
      QDPIO::cout << InlineQQQDiquarkEnv::name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << InlineQQQDiquarkEnv::name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

}
