// $Id: inline_szin_read_obj.cc,v 2.1 2005-11-01 22:00:01 edwards Exp $
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object writing
 */

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/inline_szin_read_obj.h"
#include "meas/inline/io/szin_read_obj_funcmap.h"

namespace Chroma 
{ 
  namespace InlineSZINReadNamedObjEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineSZINReadNamedObj(InlineSZINReadNamedObjParams(xml_in, path));
    }

    const std::string name = "SZIN_READ_NAMED_OBJECT";

    bool registerAll() 
    {
      bool success = true; 

      // Datatype reader
      success &= SZINReadObjCallMapEnv::registered;

      // Inline measurement
      success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);

      return success;
    }

    const bool registered = registerAll();
  };


  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineSZINReadNamedObjParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "object_id", input.object_id);
    write(xml, "object_type", input.object_type);

    pop(xml);
  }

  //! File output
  void write(XMLWriter& xml, const string& path, const InlineSZINReadNamedObjParams::File_t& input)
  {
    push(xml, path);

    write(xml, "file_name", input.file_name);

    pop(xml);
  }


  //! Object buffer
  void read(XMLReader& xml, const string& path, InlineSZINReadNamedObjParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "object_id", input.object_id);
    read(inputtop, "object_type", input.object_type);
  }

  //! File output
  void read(XMLReader& xml, const string& path, InlineSZINReadNamedObjParams::File_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "file_name", input.file_name);
  }


  // Param stuff
  InlineSZINReadNamedObjParams::InlineSZINReadNamedObjParams() { frequency = 0; }

  InlineSZINReadNamedObjParams::InlineSZINReadNamedObjParams(XMLReader& xml_in, const std::string& path) 
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
  InlineSZINReadNamedObjParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "NamedObject", named_obj);

    // Write out the destination
    Chroma::write(xml_out, "File", file);

    pop(xml_out);
  }


  void 
  InlineSZINReadNamedObj::operator()(const multi1d<LatticeColorMatrix>& u,
				      XMLBufferWriter& gauge_xml,
				      unsigned long update_no,
				      XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out, "szin_read_named_obj");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineSZINReadNamedObjEnv::name << ": object reader" << endl;
    StopWatch swatch;

    // Read the object
    // ONLY SZIN output format is supported in this task
    // Other tasks could support other disk formats
    QDPIO::cout << "Attempt to readobject name = " << params.named_obj.object_id << endl;
    write(xml_out, "object_id", params.named_obj.object_id);
    try
    {
      swatch.reset();

      // Read the object
      swatch.start();
      TheSZINReadObjFuncMap::Instance().callFunction(params.named_obj.object_type,
						     params.named_obj.object_id,
						     params.file.file_name);
      swatch.stop();

      QDPIO::cout << "Object successfully written: time= " 
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineSZINReadNamedObjEnv::name << ": cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineSZINReadNamedObjEnv::name << ": error message: " << e 
		  << endl;
      QDP_abort(1);
    }
    
    QDPIO::cout << InlineSZINReadNamedObjEnv::name << ": ran successfully" << endl;

    pop(xml_out);  // read_named_obj

    END_CODE();
  } 

};
