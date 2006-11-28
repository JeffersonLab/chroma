// $Id: inline_qpropadd_w.cc,v 3.2 2006-11-28 19:30:03 edwards Exp $
/*! \file
 * \brief Inline measurement of qpropadd
 *
 * Addition of props
 */

#include "meas/inline/hadron/inline_qpropadd_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  //! Propagator parameters
  void read(XMLReader& xml, const string& path, InlineQpropAddEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "factorA", input.factorA);
    read(inputtop, "propA", input.propA);
    read(inputtop, "factorB", input.factorB);
    read(inputtop, "propB", input.propB);
    read(inputtop, "propApB", input.propApB);
  }

  //! Propagator parameters
  void write(XMLWriter& xml, const string& path, const InlineQpropAddEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "factorA", input.factorA);
    write(xml, "propA", input.propA);
    write(xml, "factorB", input.factorB);
    write(xml, "propB", input.propB);
    write(xml, "propApB", input.propApB);

    pop(xml);
  }


  namespace InlineQpropAddEnv 
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

    const std::string name = "QPROPADD";

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
    
      // Write out the output propagator/source configuration info
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }

    //--------------------------------------------------------------


    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "qpropadd");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << "QPROPADD: propagator transformation utility" << endl;

      // Write out the input
      params.writeXML(xml_out, "Input");

      //
      // Read in the source along with relevant information.
      // 
      XMLReader propA_file_xml, propA_record_xml;
    
      LatticePropagator propA ;
      LatticePropagator propB ;
      LatticePropagator propApB ;
      QDPIO::cout << "Snarf the props from a named buffer" << endl;
      try
      {
	// Try the cast to see if this is a valid source
	propA = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.propA);

	TheNamedObjMap::Instance().get(params.named_obj.propA).getFileXML(propA_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.propA).getRecordXML(propA_record_xml);

	propB = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.propB);
      }    
      catch (std::bad_cast)
      {
	QDPIO::cerr << name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": error extracting source_header: " << e << endl;
	QDP_abort(1);
      }


      // Add the props with weights
      propApB = params.named_obj.factorA*propA  +  params.named_obj.factorB*propB;


      /*
       *  Write the a source out to a named buffer
       */
      try
      {
	QDPIO::cout << "Attempt to store sequential source" << endl;


	// Store the seqsource
	TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.propApB);
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.propApB) = propApB ;
	TheNamedObjMap::Instance().get(params.named_obj.propApB).setFileXML(propA_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.propApB).setRecordXML(propA_record_xml);

	QDPIO::cout << "Propagator sum successfully stored"  << endl;
      }
      catch (std::bad_cast)
      {
	QDPIO::cerr << name << ": dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": error storing seqsource: " << e << endl;
	QDP_abort(1);
      }


      pop(xml_out);   // qpropadd
        
      QDPIO::cout << "QpropAdd ran successfully" << endl;

      END_CODE();
    }

  }

}  
