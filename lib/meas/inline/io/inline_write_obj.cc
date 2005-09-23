// $Id: inline_write_obj.cc,v 1.1 2005-09-23 03:43:09 edwards Exp $
/*! \file
 * \brief Inline task to write an object from a named buffer
 *
 * Named object writing
 */

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_write_obj.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/writeobj_funcmap.h"

namespace Chroma 
{ 
  namespace InlineWriteNamedObjEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineWriteNamedObj(InlineWriteNamedObjParams(xml_in, path));
    }

    const std::string name = "WRITE_NAMED_OBJECT";

    bool registerAll() 
    {
      bool success = true; 

      // Datatype writer
      success &= WriteObjCallMapEnv::registered;

      // Inline measurement
      success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);

      return success;
    }

    const bool registered = registerAll();
  };


  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineWriteNamedObjParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "object_id", input.object_id);
    write(xml, "object_type", input.object_type);

    pop(xml);
  }

  //! File output
  void write(XMLWriter& xml, const string& path, const InlineWriteNamedObjParams::File_t& input)
  {
    push(xml, path);

    write(xml, "file_name", input.file_name);
    write(xml, "file_volfmt", input.file_volfmt);

    pop(xml);
  }


  //! Object buffer
  void read(XMLReader& xml, const string& path, InlineWriteNamedObjParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "object_id", input.object_id);
    read(inputtop, "object_type", input.object_type);
  }

  //! File output
  void read(XMLReader& xml, const string& path, InlineWriteNamedObjParams::File_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "file_name", input.file_name);
    read(inputtop, "file_volfmt", input.file_volfmt);
  }


  // Param stuff
  InlineWriteNamedObjParams::InlineWriteNamedObjParams() { frequency = 0; }

  InlineWriteNamedObjParams::InlineWriteNamedObjParams(XMLReader& xml_in, const std::string& path) 
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
  InlineWriteNamedObjParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "NamedObject", named_obj);

    // Write out the destination
    Chroma::write(xml_out, "File", file);

    pop(xml_out);
  }


  void 
  InlineWriteNamedObj::operator()(const multi1d<LatticeColorMatrix>& u,
				  XMLBufferWriter& gauge_xml,
				  unsigned long update_no,
				  XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out, "write_named_obj");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineWriteNamedObjEnv::name << ": object writer" << endl;
    StopWatch swatch;

    // Write the object
    // ONLY SciDAC output format is supported in this task
    // Other tasks could support other disk formats
    QDPIO::cout << "Attempt to write object name = " << params.named_obj.object_id << endl;
    write(xml_out, "object_id", params.named_obj.object_id);
    try
    {
      swatch.reset();

      // Write the object
      swatch.start();
      TheWriteObjFuncMap::Instance().callFunction(params.named_obj.object_type,
						  params.named_obj.object_id,
						  params.file.file_name, 
						  params.file.file_volfmt, QDPIO_SERIAL);
      swatch.stop();

      QDPIO::cout << "Object successfully written: time= " 
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineWriteNamedObjEnv::name << ": cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineWriteNamedObjEnv::name << ": error message: " << e 
		  << endl;
      QDP_abort(1);
    }
    
    QDPIO::cout << InlineWriteNamedObjEnv::name << ": ran successfully" << endl;

    pop(xml_out);  // write_named_obj

    END_CODE();
  } 

};
