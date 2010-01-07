// $Id: inline_szin_read_obj.cc,v 3.1 2006-09-20 20:28:03 edwards Exp $
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
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;

      const std::string name = "SZIN_READ_NAMED_OBJECT";
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Datatype reader
	success &= SZINReadObjCallMapEnv::registerAll();

	// Inline measurement
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);

	registered = true;
      }
      return success;
    }


    //! Object buffer
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "object_id", input.object_id);
      read(inputtop, "object_type", input.object_type);
    }

    //! File output
    void read(XMLReader& xml, const string& path, Params::File_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "file_name", input.file_name);
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
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "szin_read_named_obj");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": object reader" << endl;
      StopWatch swatch;

      // Read the object
      // ONLY SZIN output format is supported in this task
      // Other tasks could support other disk formats
      QDPIO::cout << "Attempt to read object name = " << params.named_obj.object_id << endl;
      QDPIO::cout << "Attempt to read file name = " << params.file.file_name << endl;
      write(xml_out, "object_id", params.named_obj.object_id);
      try
      {
	swatch.reset();

//      // Read the object
//      TheSZINReadObjFuncMap::Instance().dump();

	// Read the object
	swatch.start();
	SZINReadObjCallMapEnv::TheSZINReadObjFuncMap::Instance().callFunction(params.named_obj.object_type,
									      params.named_obj.object_id,
									      params.file.file_name);
	swatch.stop();

	QDPIO::cout << "SZIN Object successfully read: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << endl;
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name << ": cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": error message: " << e 
		    << endl;
	QDP_abort(1);
      }
    
      QDPIO::cout << name << ": ran successfully" << endl;

      pop(xml_out);  // read_named_obj

      END_CODE();
    } 

  }

}
