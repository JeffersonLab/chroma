/*! \file
 * \brief Inline construction of QpropMatMul
 *
 * QpropMatMul calculations
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_qprop_matmul_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/info/unique_id.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"
#include "util/ferm/transf.h"

namespace Chroma 
{ 
  namespace InlineQpropMatMulEnv
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineQpropMatMul(InlineQpropMatMulParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "QPROP_MAT_MUL_4D";

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
  } // end namespace


  //! QpropMatMul input
  void read(XMLReader& xml, const std::string& path, InlineQpropMatMulParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "source_id", input.source_id);
    read(inputtop, "result_id", input.result_id);
  }

  //! QpropMatMul output
  void write(XMLWriter& xml, const std::string& path, const InlineQpropMatMulParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "source_id", input.source_id);
    write(xml, "result_id", input.result_id);

    pop(xml);
  }


  // Param stuff
  InlineQpropMatMulParams::InlineQpropMatMulParams() { frequency = 0; }

  InlineQpropMatMulParams::InlineQpropMatMulParams(XMLReader& xml_in, const std::string& path)
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
    	  read(paramtop, "Frequency", frequency);
      else
    	  frequency = 1;

      // Parameters for source construction
      fermact = readXMLGroup(paramtop, "FermionAction", "FermAct");
      // Read in the output QpropMatMul/source configuration info
      read(paramtop, "NamedObject", named_obj);

      // Possible alternate XML file pattern
      if (paramtop.count("xml_file") != 0) 
      {
	read(paramtop, "xml_file", xml_file);
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
      QDP_abort(1);
    }
  }


  void
  InlineQpropMatMulParams::writeXML(XMLWriter& xml_out, const std::string& path)
  {
    push(xml_out, path);
    
    push(xml_out, "FermionAction");
    xml_out << fermact.xml;
    pop(xml_out);
    write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }


  // Function call
  void 
  InlineQpropMatMul::operator()(unsigned long update_no,
			       XMLWriter& xml_out) 
  {

      func(update_no, xml_out);
  }


  // Real work done here
  void 
  InlineQpropMatMul::func(unsigned long update_no,
			 XMLWriter& xml_out) 
  {
    START_CODE();

    QDPIO::cout << InlineQpropMatMulEnv::name << ": QpropMatMul calculation" << std::endl;

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
      QDPIO::cerr << InlineQpropMatMulEnv::name << ": caught dynamic cast error"
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << InlineQpropMatMulEnv::name << ": std::map call failed: " << e
		  << std::endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "QpropMatMul");
    write(xml_out, "update_no", update_no);

    //
    // Read in the source along with relevant information.
    // 
    XMLReader source_file_xml, source_record_xml;


    QDPIO::cout << "Snarf the source from a named buffer" << std::endl;
    try
    {
      // Try the cast to see if this is a valid source
      LatticePropagator& source_tmp =
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.source_id);

      // Snarf the source info. This is will throw if the source_id is not there
      TheNamedObjMap::Instance().get(params.named_obj.source_id).getFileXML(source_file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.source_id).getRecordXML(source_record_xml);

      // Write out the source header
      write(xml_out, "Source_file_info", source_file_xml);
      write(xml_out, "Source_record_info", source_record_xml);
    }    
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineQpropMatMulEnv::name << ": caught dynamic cast error"
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << InlineQpropMatMulEnv::name << ": error extracting source_header: " << e << std::endl;
      QDP_abort(1);
    }

    // Should be a valid cast now
    const LatticePropagator& quark_prop_source =
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.source_id);
 
    QDPIO::cout << "Source successfully read and parsed" << std::endl;


    //
    // Loop over the source color and spin, creating the source
    // and calling the relevant QpropMatMul routines. The QDP
    // terminology is that a QpropMatMul is a matrix in color
    // and spin space
    //
    try
    {
      TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.result_id);
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineQpropMatMulEnv::name << ": caught dynamic cast error"
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << InlineQpropMatMulEnv::name << ": error creating prop: " << e << std::endl;
      QDP_abort(1);
    }

    // Cast should be valid now
    LatticePropagator& quark_propagator =
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.result_id);

    //
    // Initialize fermion action
    //
    std::istringstream  xml_s(params.fermact.xml);
    XMLReader  fermacttop(xml_s);
    QDPIO::cout << "FermAct = " << params.fermact.id << std::endl;

    try {
    	StopWatch swatch;
    	swatch.reset();

    	// Typedefs to save typing
    	typedef LatticeFermion               T;
    	typedef multi1d<LatticeColorMatrix>  P;
    	typedef multi1d<LatticeColorMatrix>  Q;

    	// Generic Wilson-Type stuff


    	Handle<FermAct4D<T,P,Q> > S4(nullptr);
    	try {
    		S4=dynamic_cast<FermAct4D<T,P,Q>*>(TheFermionActionFactory::Instance().createObject(params.fermact.id,
				       fermacttop,
				       params.fermact.path));
    	}
    	catch(...) {
    		QDPIO::cout << "Cast to 4D fermact failed" << std::endl;
    		QDP_abort(1);
    	}
    	QDPIO::cout << "Creating Ferm State" << std::endl;
    	Handle< FermState<T,P,Q> > state(S4->createState(u));

    	QDPIO::cout << "Creating Linear Operator" <<std::endl;
	

    	Handle< LinearOperator<T> > M = S4->linOp(state);

    	for(int spin=0; spin < Ns; ++spin ) {
    		for(int color=0; color < Nc; ++color ) {
    			LatticeFermion src, res;
    			PropToFerm(quark_prop_source, src, color, spin);
    			(*M)(res,src,PLUS);
    			FermToProp(res,quark_propagator, color,spin);
			}
    	}

        // Write the QpropMatMul xml info
    	TheNamedObjMap::Instance().get(params.named_obj.result_id).setFileXML(source_file_xml);
    	TheNamedObjMap::Instance().get(params.named_obj.result_id).setRecordXML(source_record_xml);

    	QDPIO::cout << "QpropMatMul successfully updated" << std::endl;
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineQpropMatMulEnv::name << ": caught dynamic cast error"
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << InlineQpropMatMulEnv::name << ": error extracting prop_header: " << e << std::endl;
      QDP_abort(1);
    }

    pop(xml_out);  // QpropMatMul

    snoop.stop();
    QDPIO::cout << InlineQpropMatMulEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << std::endl;

    QDPIO::cout << InlineQpropMatMulEnv::name << ": ran successfully" << std::endl;

    END_CODE();
  } 

}

