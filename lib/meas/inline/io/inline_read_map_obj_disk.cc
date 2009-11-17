// $Id: inline_nersc_read_obj.cc,v 3.2 2006-09-20 20:28:03 edwards Exp $
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object writing
 */

#include "chromabase.h"
#include "qdp_iogauge.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_read_map_obj_disk.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineReadMapObjDiskEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineReadMapObjDisk(InlineReadMapObjDiskParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "READ_MAP_OBJECT_DISK";

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
  }


  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineReadMapObjDiskParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "object_id", input.object_id);

    pop(xml);
  }

  //! File output
  void write(XMLWriter& xml, const std::string& path, const InlineReadMapObjDiskParams& input)
  {
    input.write(xml,path);
  }


  //! Object buffer
  void read(XMLReader& xml, const string& path, InlineReadMapObjDiskParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "object_id", input.object_id);
  }



  // Param stuff
  InlineReadMapObjDiskParams::InlineReadMapObjDiskParams(XMLReader& xml_in, const std::string& path) 
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

      read(paramtop, "MapObject", map_obj_disk_p );
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineReadMapObjDiskParams::write(XMLWriter& xml_out, const std::string& path) const
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "NamedObject", named_obj);
    
    // Write out the destination
    Chroma::write(xml_out, "MapObject", map_obj_disk_p);

    pop(xml_out);
  }


  void 
  InlineReadMapObjDisk::operator()(unsigned long update_no,
				       XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out, "read_map_object_disk");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineReadMapObjDiskEnv::name << ": object reader" << endl;
    StopWatch swatch;

    // Read the object
    // ONLY SciDAC output format is supported in this task
    // Other tasks could support other disk formats
    QDPIO::cout << "Attempt to read object name = " << params.named_obj.object_id << endl;
    write(xml_out, "object_id", params.named_obj.object_id);
    try
    {
      swatch.reset();

	// Create the object as a handle. 
	// This bit will and up changing to a Factory invocation
	Handle< MapObject<KeyPropColorVec_t, LatticeFermion> >  new_map_obj_handle(
										   new MapObject
										   )

	// Create the entry
	TheNamedObjMap::Instance().create< Handle< MapObject<KeyPropColorVec_t,LatticeFermion> > >(params.named_obj.prop_id);

	// Insert
	TheNamedObjMap::Instance().getData< Handle<MapObject<KeyPropColorVec_t,LatticeFermion> > >(params.named_obj.prop_id) = new_map_obj_handle;


      TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(params.named_obj.object_id);
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.object_id) = u;


      QDPIO::cout << "Object successfully written: time= " 
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineReadMapObjDiskEnv::name << ": cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineReadMapObjDiskEnv::name << ": error message: " << e 
		  << endl;
      QDP_abort(1);
    }
    
    QDPIO::cout << InlineReadMapObjDiskEnv::name << ": ran successfully" << endl;

    pop(xml_out);  // read_named_obj

    END_CODE();
  } 

};
