/*! \file
 * \brief Inline task to erase an object from a named buffer
 *
 * Named object writing
 */

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_erase_mg_space.h"
#include "meas/inline/io/named_objmap.h"

#if BASE_PRECISION == 32
#define QDP_Precision 'F'
#define QLA_Precision 'F'
#define toReal toFloat
#elif BASE_PRECISION == 64
#define QDP_Precision 'D'
#define QLA_Precision 'D'
#define toReal toDouble
#endif

extern "C" {
  // This should be placed on the include path.
#include "wilsonmg-interface.h"

}


namespace Chroma 
{ 
  namespace InlineEraseMGSpaceEnv 
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

      const std::string name = "ERASE_MG_SUBSPACE";
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    //! Object buffer
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "object_id", input.object_id);

      pop(xml);
    }


    //! Object buffer
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "object_id", input.object_id);
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

	// Ids
	read(paramtop, "NamedObject", named_obj);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << std::endl;
	QDP_abort(1);
      }
    }


    void
    Params::writeXML(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
    
      // Ids
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }


    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "erase_mg_subspace");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": object erase" << std::endl;

      // Erase the object
      QDPIO::cout << "Attempt to erase object name = " << params.named_obj.object_id << std::endl;
      write(xml_out, "object_id", params.named_obj.object_id);
      if ( TheNamedObjMap::Instance().check(params.named_obj.object_id) ) {
	QDPIO::cout << "MG Subspace: " << params.named_obj.object_id << " found. " << std::endl;

	void *mg_subspace = TheNamedObjMap::Instance().getData<void *>(params.named_obj.object_id);
	QDPIO::cout << "Calling QOP Destroy function" << std::endl;
	MGP(destroy_subspace)(mg_subspace);
     
	try { 
	  QDPIO::cout << "Attempting to delete from named object store" << std:: endl;
	  // Now erase the object
	  TheNamedObjMap::Instance().erase(params.named_obj.object_id);
	  
	  QDPIO::cout << "Object erased" << std::endl;
	}
	catch( std::bad_cast ) {
	  QDPIO::cerr << name << ": cast error" 
		      << std::endl;
	  QDP_abort(1);
	}
	catch (const std::string& e) {
	  
	  QDPIO::cerr << name << ": error message: " << e 
		      << std::endl;
	  QDP_abort(1);
	}
      }
      else { 
	QDPIO::cout << "MG Subspace: " << params.named_obj.object_id << " is not in the map. Cannot delete" << std::endl;
      }
      QDPIO::cout << name << ": ran successfully" << std::endl;

      pop(xml_out);  // erase_named_obj

      END_CODE();
    } 

  }

}
