// $Id: inline_list_obj.cc,v 2.1 2005-09-26 19:37:16 edwards Exp $
/*! \file
 * \brief Inline task to list an object from a named buffer
 *
 * Named object writing
 */

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_list_obj.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineListNamedObjEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineListNamedObj(InlineListNamedObjParams(xml_in, path));
    }

    const std::string name = "LIST_NAMED_OBJECT";

    bool registerAll() 
    {
      bool success = true; 

      // Inline measurement
      success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);

      return success;
    }

    const bool registered = registerAll();
  };


  // Param stuff
  InlineListNamedObjParams::InlineListNamedObjParams() { frequency = 0; }

  InlineListNamedObjParams::InlineListNamedObjParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineListNamedObjParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    pop(xml_out);
  }


  void 
  InlineListNamedObj::operator()(const multi1d<LatticeColorMatrix>& u,
				  XMLBufferWriter& gauge_xml,
				  unsigned long update_no,
				  XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out, "list_named_obj");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineListNamedObjEnv::name << ": object list" << endl;

    // List the object
    QDPIO::cout << "Attempt to list all object names" << endl;
    try
    {
      TheNamedObjMap::Instance().dump();
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineListNamedObjEnv::name << ": cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineListNamedObjEnv::name << ": error message: " << e 
		  << endl;
      QDP_abort(1);
    }
    
    QDPIO::cout << InlineListNamedObjEnv::name << ": ran successfully" << endl;

    pop(xml_out);  // list_named_obj

    END_CODE();
  } 

};
