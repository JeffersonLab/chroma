/*! \file
 * \brief Inline task to read an object into a named buffer
 */

#include "chromabase.h"
#include "singleton.h"
#include "funcmap.h"

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_read_map_obj_disk.h"
#include "meas/inline/io/named_objmap.h"
#include "util/ferm/key_prop_colorvec.h"
#include "util/ferm/subset_ev_pair.h"
#include "qdp_map_obj_disk.h"
#include <string>

namespace Chroma 
{ 
  namespace InlineReadMapObjDiskEnv 
  { 
    namespace ReadMapObjCallEnv
    { 
      struct DumbDisambiguator {};

      typedef SingletonHolder< 
	FunctionMap<DumbDisambiguator,
		    int,
		    std::string,
		    TYPELIST_2(const std::string&, const std::string&),
		    int (*)(const std::string&, const std::string&),
		    StringFunctionMapError> >
      TheReadMapObjFuncMap;

      namespace 
      { 
	static bool registered = false;

	template<typename K, typename V>
	int readMapObj(const std::string& object_id,
		       const std::string& file_name)
	{
	  QDP::MapObjectDisk<K,V>* obj_obj = new QDP::MapObjectDisk<K,V>();
	  obj_obj->open(file_name);

	  Handle<QDP::MapObject<K,V> > obj_handle(obj_obj);
	  TheNamedObjMap::Instance().create< Handle<QDP::MapObject<K,V> >, Handle<QDP::MapObject<K,V> > >(object_id, obj_handle);

	  return obj_handle->size();
	}

	bool registerAll(void) 
	{
	  bool success = true; 
	  if (! registered ) 
	  { 
	    success &= TheReadMapObjFuncMap::Instance().registerFunction("KeyTKeyPropColorVec_tValTLatticeFermion",
									 readMapObj<KeyPropColorVec_t, LatticeFermion>);

	    success &= TheReadMapObjFuncMap::Instance().registerFunction("KeyTintValTEVPairLatticeColorVector",
									 readMapObj<int, EVPair<LatticeColorVector> >);

	    success &= TheReadMapObjFuncMap::Instance().registerFunction("KeyTcharValTfloat",
									 readMapObj<char, float>);
	    registered = true;
	  }
	  return success;
	}
      }
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

      const std::string name = "READ_MAP_OBJECT_DISK";
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	success &= ReadMapObjCallEnv::registerAll();
	registered = true;
      }
      return success;
    }


    //! Object buffer
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "object_type", input.object_type);
      read(inputtop, "object_id", input.object_id);
    }

    //! Object buffer
    void read(XMLReader& xml, const string& path, Params::File& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "file_name", input.file_name);
    }

    Params::Params(XMLReader& reader, const std::string& path)
    {
      try 
      {
	XMLReader paramtop(reader, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	read(paramtop, "NamedObject", named_obj);
	read(paramtop, "File", file);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }


    void 
    InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "read_map_object_disk");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": object reader" << endl;
      StopWatch swatch;

      // Read the object
      // ONLY MapObject output format is supported in this task
      // Other tasks could support other disk formats
      QDPIO::cout << "Attempt to read object name = " << params.named_obj.object_id << endl;

      write(xml_out, "object_type", params.named_obj.object_type);
      write(xml_out, "object_id", params.named_obj.object_id);
      write(xml_out, "file_name", params.file.file_name);

      try
      {
	swatch.reset();
	swatch.start();

        // Read the object
        int size = ReadMapObjCallEnv::TheReadMapObjFuncMap::Instance().callFunction(params.named_obj.object_type, params.named_obj.object_id, params.file.file_name);

	XMLBufferWriter file_xml_buf;
	push(file_xml_buf, "FileXML");
	write(file_xml_buf,  "object_type", params.named_obj.object_type);
	write(file_xml_buf,  "object_id", params.named_obj.object_id);
	write(file_xml_buf,  "file_name", params.file.file_name);
	write(file_xml_buf,  "map_size", size);
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

  }
}
