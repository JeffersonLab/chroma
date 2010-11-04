/*! \file
 * \brief Inline task to write an object from a named buffer
 *
 * Named object writing
 */

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_qio_write_obj.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/qio_write_obj_funcmap.h"
#include "io/enum_io/enum_qdpvolfmt_io.h"

namespace Chroma 
{ 
  namespace InlineQIOWriteNamedObjEnv 
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

      const std::string name = "QIO_WRITE_NAMED_OBJECT";
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Datatype writer
	success &= QIOWriteObjCallMapEnv::registerAll();

	// Inline measurement
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);

	registered = true;
      }
      return success;
    }


    //! Object buffer
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "object_id", input.object_id);
      write(xml, "object_type", input.object_type);

      pop(xml);
    }

    //! File output
    void write(XMLWriter& xml, const string& path, const Params::File_t& input)
    {
      push(xml, path);

      write(xml, "file_name", input.file_name);
      write(xml, "file_volfmt", input.file_volfmt);
      write(xml, "parallel_io", input.parallel_io);

      pop(xml);
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
      Chroma::read(inputtop, "file_volfmt", input.file_volfmt);
      if( inputtop.count("parallel_io") > 0 ) {
	read(inputtop, "parallel_io", input.parallel_io);
      }
      else { 
	input.parallel_io = false;
      }

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
    Params::writeXML(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
    
      // Parameters for source construction
      write(xml_out, "NamedObject", named_obj);

      // Write out the destination
      write(xml_out, "File", file);

      pop(xml_out);
    }


    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "qio_write_named_obj");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": object writer" << endl;
      StopWatch swatch;

      // Write the object
      // ONLY SciDAC output format is supported in this task
      // Other tasks could support other disk formats
      QDPIO::cout << "Attempt to write object name = " << params.named_obj.object_id << endl;
      write(xml_out, "object_id", params.named_obj.object_id);
      QDP_serialparallel_t parallel_io_type = QDPIO_SERIAL;
      if ( params.file.parallel_io ) { 
	QDPIO::cout << "Attempting to write with Parallel IO" << endl;
	parallel_io_type = QDPIO_PARALLEL;
      }
      else { 
	QDPIO::cout << "Attempting to write without parallel IO" << endl;
	parallel_io_type = QDPIO_SERIAL;
      }

      try
      {
	swatch.reset();

	// Write the object
	swatch.start();
	QIOWriteObjCallMapEnv::TheQIOWriteObjFuncMap::Instance().callFunction(params.named_obj.object_type,
									      params.named_obj.object_id,
									      params.file.file_name, 
									      params.file.file_volfmt, parallel_io_type);
	swatch.stop();

	QDPIO::cout << "Object successfully written: time= " 
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

      pop(xml_out);  // write_named_obj

      END_CODE();
    } 

  }

}
