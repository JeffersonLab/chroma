/*! \file
 * \brief Inline task to write and delete an object from a named buffer
 *
 * Named object writing
 */

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_qio_write_erase_obj.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineQIOWriteEraseNamedObjEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(InlineQIOWriteNamedObjEnv::Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;

      const std::string name = "QIO_WRITE_ERASE_NAMED_OBJECT";
    }

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


    // Func
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "qio_write_erase_named_obj");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": object writer" << endl;
      StopWatch swatch;

      // Write and erase the object
      // ONLY SciDAC output format is supported in this task
      // Other tasks could support other disk formats
      QDPIO::cout << "Attempt to write then delete an object name = " << params.named_obj.object_id << endl;
      write(xml_out, "object_id", params.named_obj.object_id);
      try
      {
	swatch.reset();

	// Use the qio writer
	InlineQIOWriteNamedObjEnv::InlineMeas qio_write(params);

	// Write the object
	qio_write(update_no, xml_out);

	// Now erase the object
	TheNamedObjMap::Instance().erase(params.named_obj.object_id);

	QDPIO::cout << "Object erased" << endl;
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

      pop(xml_out);  // qio_write_erase_named_obj

      END_CODE();
    } 

  }

}
