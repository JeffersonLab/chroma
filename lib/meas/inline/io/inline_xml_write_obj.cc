// $Id: inline_xml_write_obj.cc,v 3.1 2006-09-20 20:28:03 edwards Exp $
/*! \file
 * \brief Inline task to write an object from a named buffer
 *
 * Named object writing
 */

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_xml_write_obj.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/xml_write_obj_funcmap.h"

namespace Chroma 
{ 
  namespace InlineXMLWriteNamedObjEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineXMLWriteNamedObj(InlineXMLWriteNamedObjParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "XML_WRITE_NAMED_OBJECT";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Datatype writer
	success &= XMLWriteObjCallMapEnv::registerAll();

	// Inline measurement
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);

	registered = true;
      }
      return success;
    }
  }


  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineXMLWriteNamedObjParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "object_id", input.object_id);
    write(xml, "object_type", input.object_type);

    pop(xml);
  }

  //! File output
  void write(XMLWriter& xml, const string& path, const InlineXMLWriteNamedObjParams::File_t& input)
  {
    push(xml, path);

    write(xml, "file_name", input.file_name);

    pop(xml);
  }


  //! Object buffer
  void read(XMLReader& xml, const string& path, InlineXMLWriteNamedObjParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "object_id", input.object_id);
    read(inputtop, "object_type", input.object_type);
  }

  //! File output
  void read(XMLReader& xml, const string& path, InlineXMLWriteNamedObjParams::File_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "file_name", input.file_name);
  }


  // Param stuff
  InlineXMLWriteNamedObjParams::InlineXMLWriteNamedObjParams() { frequency = 0; }

  InlineXMLWriteNamedObjParams::InlineXMLWriteNamedObjParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Parameters for source construction
      read(paramtop, "NamedObject", named_obj);

      // Read in the destination
      read(paramtop, "File", file);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineXMLWriteNamedObjParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "NamedObject", named_obj);

    // Write out the destination
    Chroma::write(xml_out, "File", file);

    pop(xml_out);
  }


  void 
  InlineXMLWriteNamedObj::operator()(unsigned long update_no,
				     XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out, "xml_write_named_obj");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineXMLWriteNamedObjEnv::name << ": object writer" << endl;
    StopWatch swatch;

    // Write the object
    // ONLY XML output format is supported in this task
    // Other tasks could support other formats
    QDPIO::cout << "Attempt to write object name = " << params.named_obj.object_id << endl;
    write(xml_out, "object_id", params.named_obj.object_id);
    try
    {
      swatch.reset();

      // Write the object
      swatch.start();
      XMLWriteObjCallMapEnv::TheXMLWriteObjFuncMap::Instance().callFunction(params.named_obj.object_type,
									    params.named_obj.object_id,
									    params.file.file_name);
      swatch.stop();

      QDPIO::cout << "Object successfully written: time= " 
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineXMLWriteNamedObjEnv::name << ": cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineXMLWriteNamedObjEnv::name << ": error message: " << e 
		  << endl;
      QDP_abort(1);
    }
    
    QDPIO::cout << InlineXMLWriteNamedObjEnv::name << ": ran successfully" << endl;

    pop(xml_out);  // write_named_obj

    END_CODE();
  } 

};
