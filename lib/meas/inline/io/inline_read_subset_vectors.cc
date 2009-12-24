/*! \file
 * \brief Inline task to read an object
 */

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_read_subset_vectors.h"
#include "meas/inline/io/named_objmap.h"
#include "util/ferm/subset_vectors.h"
#include <string>

#include "util/ferm/map_obj.h"
#include "util/ferm/map_obj/map_obj_aggregate_w.h"
#include "util/ferm/map_obj/map_obj_factory_w.h"

namespace Chroma 
{ 
  namespace InlineReadSubsetVectorsEnv 
  { 
    //! Object buffer
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "object_id", input.object_id);

      // User Specified MapObject tags
      input.object_map = readXMLGroup(inputtop, "ColorVecMapObject", "MapObjType");
    }

    //! File output
    void read(XMLReader& xml, const string& path, Params::File_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "file_name", input.file_name);
    }


    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }
    
      //! Local registration flag
      bool registered = false;
      const std::string name = "READ_SUBSET_VECTORS";
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered) 
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	success &= MapObjectWilson4DEnv::registerAll();
	registered = true;
      }
      return success;
    }
  

    // Param stuff
    Params::Params() { frequency = 0; }

    Params::Params(XMLReader& reader, 
		   const std::string& path)
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

      push(xml_out, "read_subset_vectors");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": object reader" << endl;
      StopWatch swatch;

      // Read the object
      QDPIO::cout << "Attempt to read object name = " << params.named_obj.object_id << endl;

      write(xml_out, "object_id", params.named_obj.object_id);
      write(xml_out, "file_name", params.file.file_name);

      try
      {
	swatch.reset();
	swatch.start();
	
	TheNamedObjMap::Instance().create< 
	SubsetVectors<LatticeColorVector>,
	  std::string>(params.named_obj.object_id,
		       params.file.file_name);

	SubsetVectors<LatticeColorVector>& sv = TheNamedObjMap::Instance().getData< SubsetVectors<LatticeColorVector> >(params.named_obj.object_id);
	sv.openRead();

	XMLBufferWriter file_xml_buf;
	push(file_xml_buf, "FileXML");
	write(file_xml_buf,  "object_id", params.named_obj.object_id);
	write(file_xml_buf,  "file_name", params.file.file_name);
	write(file_xml_buf,  "map_size", sv.size());
	pop(file_xml_buf);

	XMLReader file_xml(file_xml_buf);

	TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML( file_xml );

	// No particularly good record XML -- this after all is not QIO. So just use file_xml again
	TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML( file_xml );

	swatch.stop();

	QDPIO::cout << "Object successfully read: time= " 
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

  } // namespace InlineReadSubsetVectorsEnv 

}
