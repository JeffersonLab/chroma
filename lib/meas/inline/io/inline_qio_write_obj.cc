// $Id: inline_qio_write_obj.cc,v 2.0 2005-09-25 21:04:38 edwards Exp $
/*! \file
 * \brief Inline task to write an object from a named buffer
 *
 * Named object writing
 */

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_qio_write_obj.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/qio_write_obj_funcmap.h"

namespace Chroma 
{ 
  namespace InlineQIOWriteNamedObjEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineQIOWriteNamedObj(InlineQIOWriteNamedObjParams(xml_in, path));
    }

    const std::string name = "QIO_WRITE_NAMED_OBJECT";

    bool registerAll() 
    {
      bool success = true; 

      // Datatype writer
      success &= QIOWriteObjCallMapEnv::registered;

      // Inline measurement
      success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);

      return success;
    }

    const bool registered = registerAll();
  };


  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineQIOWriteNamedObjParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "object_id", input.object_id);
    write(xml, "object_type", input.object_type);

    pop(xml);
  }

  //! File output
  void write(XMLWriter& xml, const string& path, const InlineQIOWriteNamedObjParams::File_t& input)
  {
    push(xml, path);

    write(xml, "file_name", input.file_name);
    write(xml, "file_volfmt", input.file_volfmt);

    pop(xml);
  }


  //! Object buffer
  void read(XMLReader& xml, const string& path, InlineQIOWriteNamedObjParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "object_id", input.object_id);
    read(inputtop, "object_type", input.object_type);
  }

  //! File output
  void read(XMLReader& xml, const string& path, InlineQIOWriteNamedObjParams::File_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "file_name", input.file_name);
    read(inputtop, "file_volfmt", input.file_volfmt);
  }


  // Param stuff
  InlineQIOWriteNamedObjParams::InlineQIOWriteNamedObjParams() { frequency = 0; }

  InlineQIOWriteNamedObjParams::InlineQIOWriteNamedObjParams(XMLReader& xml_in, const std::string& path) 
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
  InlineQIOWriteNamedObjParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "NamedObject", named_obj);

    // Write out the destination
    Chroma::write(xml_out, "File", file);

    pop(xml_out);
  }


  void 
  InlineQIOWriteNamedObj::operator()(const multi1d<LatticeColorMatrix>& u,
				  XMLBufferWriter& gauge_xml,
				  unsigned long update_no,
				  XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out, "qio_write_named_obj");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineQIOWriteNamedObjEnv::name << ": object writer" << endl;
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
      TheQIOWriteObjFuncMap::Instance().callFunction(params.named_obj.object_type,
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
      QDPIO::cerr << InlineQIOWriteNamedObjEnv::name << ": cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineQIOWriteNamedObjEnv::name << ": error message: " << e 
		  << endl;
      QDP_abort(1);
    }
    
    QDPIO::cout << InlineQIOWriteNamedObjEnv::name << ": ran successfully" << endl;

    pop(xml_out);  // write_named_obj

    END_CODE();
  } 

};
