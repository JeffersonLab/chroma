/*! \file
 * \brief Inline task to read an object into a named buffer
 */

#include "chromabase.h"
#include "singleton.h"
#include "funcmap.h"

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_copy_map_obj.h"
#include "meas/inline/io/named_objmap.h"

#include "util/ferm/subset_ev_pair.h"
#include "util/ferm/map_obj/map_obj_factory_w.h"
#include "util/ferm/map_obj/map_obj_aggregate_w.h"

namespace Chroma 
{ 
  namespace InlineCopyMapObjEnv 
  { 
    namespace CopyMapObjCallEnv
    { 
      struct DumbDisambiguator {};

      //! Write object function map
      /*! \ingroup inlineio */
      typedef SingletonHolder< 
	FunctionMap<DumbDisambiguator,
		    void,
		    std::string,
		    TYPELIST_1(const Params&),
		    void (*)(const Params& named_obj),
		    StringFunctionMapError> >
      TheCopyMapObjFuncMap;


      namespace 
      { 
	static bool registered = false;

	template<typename K, typename V>
	void copyMapObj(const Params& params)
	{
	  // Input object
	  Handle< QDP::MapObject<K,V> > input_obj = TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<K,V> > >(params.named_obj.input_id);
	  std::vector<K> keys; input_obj->keys(keys);

	  // Create output object
	  std::istringstream  xml_s(params.named_obj.output_obj.xml);
	  XMLReader MapObjReader(xml_s);
	
	  // Create the entry
	  Handle< QDP::MapObject<K,V> > output_obj(
	    SingletonHolder< ObjectFactory<QDP::MapObject<K,V>, 
	    std::string,
	    TYPELIST_2(XMLReader&, const std::string&),
	    QDP::MapObject<K,V>* (*)(XMLReader&, const std::string&), 
	    StringFactoryError> >::Instance().createObject(params.named_obj.output_obj.id,
							   MapObjReader,
							   params.named_obj.output_obj.path) );

	  TheNamedObjMap::Instance().create< Handle< QDP::MapObject<K,V> >, Handle< QDP::MapObject<K,V> > >(params.named_obj.output_id, output_obj);

	  // Copy the key/value-s
	  for(int i=0; i < keys.size(); i++) 
	  {
	    V v;
	    input_obj->get(keys[i],v);
	    output_obj->insert(keys[i],v);
	  }
	  output_obj->flush();
	}


	bool registerAll(void) 
	{
	  bool success = true; 
	  if (! registered ) 
	  { 
	    success &= TheCopyMapObjFuncMap::Instance().registerFunction("KeyTKeyPropColorVec_tValTLatticeFermion",
									 copyMapObj<KeyPropColorVec_t, LatticeFermion>);

	    success &= TheCopyMapObjFuncMap::Instance().registerFunction("KeyTintValTEVPairLatticeColorVector",
									 copyMapObj<int, EVPair<LatticeColorVector> >);

	    success &= TheCopyMapObjFuncMap::Instance().registerFunction("KeyTcharValTfloat",
									 copyMapObj<char, float>);
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

      const std::string name = "COPY_MAP_OBJECT";
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	success &= CopyMapObjCallEnv::registerAll();
	registered = true;
      }
      return success;
    }


    //! Object buffer
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "object_type", input.object_type);
      read(inputtop, "input_id", input.input_id);
      read(inputtop, "output_id", input.output_id);
      input.output_obj = readXMLGroup(inputtop, "MapObject", "MapObjType");
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

      registerAll();

      push(xml_out, "read_map_object_disk");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": map object copy" << endl;
      StopWatch swatch;

      // Copy the object
      // ONLY named object format is supported in this task
      // Other tasks could support other disk formats
      QDPIO::cout << "Attempt to copy input object name = " << params.named_obj.input_id << endl;

      write(xml_out, "object_type", params.named_obj.object_type);
      write(xml_out, "input_id", params.named_obj.input_id);
      write(xml_out, "output_id", params.named_obj.output_id);

      try
      {
	swatch.reset();
	swatch.start();

	// Copy the object
	CopyMapObjCallEnv::TheCopyMapObjFuncMap::Instance().callFunction(params.named_obj.object_type, params);

	// Use the xml from the first object
	XMLReader file_xml, record_xml;
	TheNamedObjMap::Instance().get(params.named_obj.input_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.input_id).getRecordXML(record_xml);

	TheNamedObjMap::Instance().get(params.named_obj.output_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.output_id).setRecordXML(record_xml);

	swatch.stop();

	QDPIO::cout << "Object successfully copied: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << endl;
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name 
		    << ": cast error for input_id= " << params.named_obj.input_id 
		    << " with type= " << params.named_obj.object_type 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": error message: " << e << endl;
	QDP_abort(1);
      }
      catch(const char* e) 
      { 
	QDPIO::cout << name << ": Caught const char * exception: " << e << endl;
	QDP_abort(1);
      }
    
      QDPIO::cout << name << ": ran successfully" << endl;

      pop(xml_out);  // read_named_obj

      END_CODE();
    } 

  }
}
