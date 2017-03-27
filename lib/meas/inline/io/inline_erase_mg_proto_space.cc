/*! \file
 * \brief Inline task to erase an object from a named buffer
 *
 * Named object writing
 */

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_erase_mg_proto_space.h"
#include "meas/inline/io/named_objmap.h"

#include "actions/ferm/invert/mg_proto/mg_proto_helpers.h"


namespace Chroma 
{ 
  namespace InlineEraseMGProtoSpaceEnv
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

      const std::string name = "ERASE_MG_PROTO_SUBSPACE";
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

      push(xml_out, "erase_mg_proto_subspace");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": object erase" << std::endl;

      // Erase the object
      QDPIO::cout << "Attempt to erase object name = " << params.named_obj.object_id << std::endl;
      write(xml_out, "object_id", params.named_obj.object_id);
      if ( TheNamedObjMap::Instance().check(params.named_obj.object_id) ) {
    	  QDPIO::cout << "MG PROTO Subspace: " << params.named_obj.object_id << " found. " << std::endl;
    	  MGProtoHelpers::deleteMGPreconditioner(params.named_obj.object_id);
    	  QDPIO::cout << "... deleted" << std::endl;
      }
      else { 
    	  QDPIO::cout << "MG PROTO Subspace Not Found. Aborting" <<std::endl;
    	  QDP_abort(1);
      }

      END_CODE();
    } 

  }

}
