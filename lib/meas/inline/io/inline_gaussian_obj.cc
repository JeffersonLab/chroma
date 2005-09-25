// $Id: inline_gaussian_obj.cc,v 2.0 2005-09-25 21:04:38 edwards Exp $
/*! \file
 * \brief Inline task to gaussian init a named object
 *
 * Named object initialization
 */

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/inline_gaussian_obj.h"

namespace Chroma 
{ 
  //! GaussianInit init function map
  /*! \ingroup inlineio */
  typedef SingletonHolder< 
    FunctionMap<void,
		std::string,
		TYPELIST_1(const string&),
		void (*)(const string& buffer_id),
		StringFunctionMapError> >
  TheGaussianInitObjFuncMap;


  //! IO function map
  /*! \ingroup inlineio */
  namespace GaussianInitObjCallMap
  {
    //! Init a propagator
    void GaussianInitLatProp(const string& buffer_id)
    {
      TheNamedObjMap::Instance().create<LatticePropagator>(buffer_id);
      XMLBufferWriter file_xml, record_xml;

      push(file_xml,"FileXML");
      pop(file_xml);

      push(record_xml,"RecordXML");
      pop(record_xml);

      gaussian(TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id));
      TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
    }

  }  // end namespace GaussianInitCallMap


  //! IO function map environment
  /*! \ingroup inlineio */
  namespace GaussianInitObjCallMapEnv
  { 
    bool registerAll(void) 
    {
      bool success = true;
      success &= TheGaussianInitObjFuncMap::Instance().registerFunction(string("LatticePropagator"), 
									GaussianInitObjCallMap::GaussianInitLatProp);
      return success;
    }

    bool registered = registerAll();
  };


  namespace InlineGaussianInitNamedObjEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineGaussianInitNamedObj(InlineGaussianInitNamedObjParams(xml_in, path));
    }

    const std::string name = "GAUSSIAN_INIT_NAMED_OBJECT";

    bool registerAll() 
    {
      bool success = true; 

      // Gaussian init functions
      success &= GaussianInitObjCallMapEnv::registered;

      // Inline measurement
      success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);

      return success;
    }

    const bool registered = registerAll();
  };


  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineGaussianInitNamedObjParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "object_id", input.object_id);
    write(xml, "object_type", input.object_type);

    pop(xml);
  }


  //! Object buffer
  void read(XMLReader& xml, const string& path, InlineGaussianInitNamedObjParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "object_id", input.object_id);
    read(inputtop, "object_type", input.object_type);
  }


  // Param stuff
  InlineGaussianInitNamedObjParams::InlineGaussianInitNamedObjParams() { frequency = 0; }

  InlineGaussianInitNamedObjParams::InlineGaussianInitNamedObjParams(XMLReader& xml_in, const std::string& path) 
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
      QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineGaussianInitNamedObjParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Ids
    Chroma::write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }


  void 
  InlineGaussianInitNamedObj::operator()(const multi1d<LatticeColorMatrix>& u,
					 XMLBufferWriter& gauge_xml,
					 unsigned long update_no,
					 XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out, "gaussian_init_named_obj");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineGaussianInitNamedObjEnv::name << ": object list" << endl;

    // Gaussian the object
    QDPIO::cout << "Attempt to list all object names" << endl;
    try
    {
      // Gaussian init the object
      TheGaussianInitObjFuncMap::Instance().callFunction(params.named_obj.object_type,
							 params.named_obj.object_id);
    }
    catch (std::bad_cast) 
    {
      QDPIO::cerr << InlineGaussianInitNamedObjEnv::name << ": cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineGaussianInitNamedObjEnv::name << ": error message: " << e 
		  << endl;
      QDP_abort(1);
    }
    
    QDPIO::cout << InlineGaussianInitNamedObjEnv::name << ": ran successfully" << endl;

    pop(xml_out);  // gaussian_init_named_obj

    END_CODE();
  } 

};
