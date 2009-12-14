// $Id: inline_nersc_read_obj.cc,v 3.2 2006-09-20 20:28:03 edwards Exp $
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object writing
 */

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_read_subset_vectors_memory.h"
#include "meas/inline/io/named_objmap.h"
#include "util/ferm/subset_vectors.h"
#include <string>

using namespace QDP;

namespace Chroma 
{ 
  namespace InlineReadSubsetVectorsMemoryEnv 
  { 

    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineReadSubsetVectorsMemory(InlineReadSubsetVectorsMemoryParams(xml_in, path));
    }
    
    //! Local registration flag
    bool registered = false;
    const std::string name = "READ_SUBSET_VECTORS_MEMORY";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered) {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  
  };

  InlineReadSubsetVectorsMemoryParams::InlineReadSubsetVectorsMemoryParams(XMLReader& reader, 
							 const std::string& path)
  {
    try 
    {
      XMLReader paramtop(reader, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Parameters for source construction
      read(paramtop, "NamedObject/object_id", named_obj.object_id);
      read(paramtop, "File/file_name", file.file_name );
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
    
  }

  void read(XMLReader& xml_in, const std::string& path, InlineReadSubsetVectorsMemoryParams& p) 
  {
    InlineReadSubsetVectorsMemoryParams tmp(xml_in, path);
    p = tmp;
  }

  void
  write(XMLWriter& xml_out, const std::string& path, const InlineReadSubsetVectorsMemoryParams& p)
  {
    push(xml_out, path);
    
    write(xml_out, "Frequency", p.frequency);

    push(xml_out, "NamedObject");
    write(xml_out, "object_id", p.named_obj.object_id);
    pop(xml_out);

    push(xml_out, "File");
    write(xml_out, "file_name", p.file.file_name);
    pop(xml_out);

    pop(xml_out);
  }

  void 
  InlineReadSubsetVectorsMemory::operator()(unsigned long update_no,
				   XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out, "read_subset_vectors");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineReadSubsetVectorsMemoryEnv::name << ": object reader" << endl;
    StopWatch swatch;

    // Read the object
    // ONLY SciDAC output format is supported in this task
    // Other tasks could support other disk formats
    QDPIO::cout << "Attempt to read object name = " << params.named_obj.object_id << endl;

    write(xml_out, "object_id", params.named_obj.object_id);
    write(xml_out, "file_name", params.file.file_name);

    try
    {
      swatch.reset();
      swatch.start();


      TheNamedObjMap::Instance().create< 
      SubsetVectors<LatticeColorVector> > (params.named_obj.object_id);

      SubsetVectors<LatticeColorVector> diskMap(params.file.file_name);
      diskMap.openRead();

      std::vector<int> keys = diskMap.dump();


      SubsetVectors<LatticeColorVector>& sv = TheNamedObjMap::Instance().getData< SubsetVectors<LatticeColorVector> >(params.named_obj.object_id);
      sv.openWrite();
      for(int i=0; i < diskMap.size(); i++) {
	EVPair<LatticeColorVector> v;
	diskMap.lookup(keys[i], v);
	sv.insert(keys[i],v);
      }
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
      QDPIO::cerr << InlineReadSubsetVectorsMemoryEnv::name << ": cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineReadSubsetVectorsMemoryEnv::name << ": error message: " << e 
		  << endl;
      QDP_abort(1);
    }
    
    QDPIO::cout << InlineReadSubsetVectorsMemoryEnv::name << ": ran successfully" << endl;

    pop(xml_out);  // read_named_obj

    END_CODE();
  } 

};
