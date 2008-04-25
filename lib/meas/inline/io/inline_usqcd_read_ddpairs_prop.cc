// $Id: inline_usqcd_read_ddpairs_prop.cc,v 3.2 2008-04-25 19:27:54 bjoo Exp $
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object writing
 */

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_usqcd_read_ddpairs_prop.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineUSQCDReadDDPairsPropEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineUSQCDReadDDPairsProp(InlineUSQCDReadDDPairsPropParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "USQCD_READ_DD_PAIRS_PROP";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Inline measurement
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);

	registered = true;
      }
      return success;
    }
  }


  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineUSQCDReadDDPairsPropParams& p)
  {
    push(xml, path);
    write(xml, "Frequency", p.frequency);
    push(xml, "Param");
    write(xml, "InputFile", p.input_file_name);
    write(xml, "SourceXML", p.source_xml);
    pop(xml); // Param
    
    push(xml, "NamedObject");
    write(xml, "source_id", p.source_id); // ID Of source to create
    write(xml, "prop_id", p.prop_id); // ID Of prop to create
    pop(xml);

    if ( p.xml_file != "") { 
      write(xml, "xml_file", p.xml_file);
    }
    pop(xml);
  }

  // Param stuff
  InlineUSQCDReadDDPairsPropParams::InlineUSQCDReadDDPairsPropParams() { frequency = 0; xml_file=""; }

  InlineUSQCDReadDDPairsPropParams::InlineUSQCDReadDDPairsPropParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;


      read(paramtop, "Param/InputFile", input_file_name);
      
      XMLReader prop_header(paramtop, "Param/SourceXML");
      ostringstream dummy;
      prop_header.printCurrentContext(dummy);
      source_xml = dummy.str();
      read(paramtop, "./NamedObject/source_id", source_id);
      read(paramtop, "./NamedObject/prop_id", prop_id);

      if( paramtop.count("xml_file") == 1 ) { 
	read(paramtop, "./xml_file", xml_file);
      }
      else { 
	xml_file="";
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }




  void 
  InlineUSQCDReadDDPairsProp::operator()(unsigned long update_no,
				    XMLWriter& xml_out) 
  {
    START_CODE();
    if( params.xml_file == "" ) { 
      func(update_no, xml_out);
    }
    else { 
      XMLFileWriter separate_output(params.xml_file);
      func(update_no, separate_output);
    }

    END_CODE();
  } 


  void InlineUSQCDReadDDPairsProp::func(unsigned long update_no, XMLWriter& xml_out)
  {
    START_CODE();

    push(xml_out, "USQCDReadDDPairsProp");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineUSQCDReadDDPairsPropEnv::name << ":" << endl;
    StopWatch swatch;
    QDPIO::cout << "Attempt to read object name = " << params.input_file_name;

    XMLBufferWriter testout;
    write(testout, "InputParams", params);
   
    QDPIO::cout << "Echo Of Input Params" << endl;
    QDPIO::cout << testout.str();

    /* Do the IO here */
    /* Create the named objects */
    try { 
      TheNamedObjMap::Instance().create<LatticePropagator>(params.source_id);
      TheNamedObjMap::Instance().create<LatticePropagator>(params.prop_id);
    }
    catch(...) { 
      QDPIO::cerr << "Failed to create the source and prop" <<endl;
      QDP_abort(1);
    }

    /* Get references to the data */
    LatticePropagator& source_ref=TheNamedObjMap::Instance().getData<LatticePropagator>(params.source_id);
    LatticePropagator& prop_ref=TheNamedObjMap::Instance().getData<LatticePropagator>(params.prop_id);

    /* Open the file */
    XMLReader file_xml;
    QDPFileReader input_file(file_xml, 
			     params.input_file_name,
			     QDPIO_SERIAL);

    QDPIO::cout << "FILE XML IS: " << endl;
    
    QDPIO::cout << InlineUSQCDReadDDPairsPropEnv::name << ": ran successfully" << endl;

    pop(xml_out);  // qio_read_named_obj

    END_CODE();
}
};
