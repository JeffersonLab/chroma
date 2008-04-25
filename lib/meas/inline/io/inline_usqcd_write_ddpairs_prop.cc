// $Id: inline_usqcd_write_ddpairs_prop.cc,v 1.2 2008-04-25 19:27:54 bjoo Exp $
/*! \file
 * \brief Inline task to write an object from a named buffer
 *
 * Named object writing
 */

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_usqcd_write_ddpairs_prop.h"
#include "meas/inline/io/named_objmap.h"
#include "io/enum_io/enum_io.h"

namespace Chroma 
{ 
  namespace InlineUSQCDWriteDDPairsPropEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineUSQCDWriteDDPairsProp(InlineUSQCDWriteDDPairsPropParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "USQCD_WRITE_DD_PAIRS_PROP";

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
  void write(XMLWriter& xml, const string& path, const InlineUSQCDWriteDDPairsPropParams& p)
  {
    push(xml, path);
    write(xml, "Frequency", p.frequency);
    push(xml, "Param");
    write(xml, "OutputFile", p.output_file_name);
    write(xml, "OutputVolfmt", p.qio_volfmt);
    pop(xml); // Param
    
    push(xml, "NamedObject");
    write(xml, "prop_id", p.prop_id); // ID Of prop to create
    pop(xml);

    if ( p.xml_file != "") { 
      write(xml, "xml_file", p.xml_file);
    }
    pop(xml);
  }

  // Param stuff
  InlineUSQCDWriteDDPairsPropParams::InlineUSQCDWriteDDPairsPropParams() { frequency = 0; xml_file=""; }

  InlineUSQCDWriteDDPairsPropParams::InlineUSQCDWriteDDPairsPropParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;


      read(paramtop, "Param/OutputFile", output_file_name);
      read(paramtop, "Param/OutputVolfmt", qio_volfmt);


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
  InlineUSQCDWriteDDPairsProp::operator()(unsigned long update_no,
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


  void InlineUSQCDWriteDDPairsProp::func(unsigned long update_no, XMLWriter& xml_out)
  {
    START_CODE();

    push(xml_out, "USQCDWriteDDPairsProp");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineUSQCDWriteDDPairsPropEnv::name << ":" << endl;
    StopWatch swatch;
    QDPIO::cout << "Attempt to write object name = " << params.prop_id << endl;
    QDPIO::cout << "Output file = " << params.output_file_name << endl;

    /* Routine to write the prop here */
 



    /* End of prop writing routine */
   
    
    QDPIO::cout << InlineUSQCDWriteDDPairsPropEnv::name << ": ran successfully" << endl;

    pop(xml_out);  // qio_write_named_obj

    END_CODE();
}
};
