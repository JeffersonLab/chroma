// $Id: inline_snarf_obj.cc,v 2.2 2006-03-20 04:22:03 edwards Exp $
/*! \file
 * \brief Inline task to copy gauge arg to a named object
 *
 * Named object writing
 */

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/inline_snarf_obj.h"
#include "meas/inline/io/default_gauge_field.h"

namespace Chroma 
{ 
  namespace InlineSnarfNamedObjEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineSnarfNamedObj(InlineSnarfNamedObjParams(xml_in, path));
    }

    const std::string name = "SNARF_NAMED_OBJECT";

    bool registerAll() 
    {
      bool success = true; 

      // Inline measurement
      success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);

      return success;
    }

    const bool registered = registerAll();
  };


  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineSnarfNamedObjParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "object_id", input.object_id);

    pop(xml);
  }


  //! Object buffer
  void read(XMLReader& xml, const string& path, InlineSnarfNamedObjParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "object_id", input.object_id);
  }


  // Param stuff
  InlineSnarfNamedObjParams::InlineSnarfNamedObjParams() { frequency = 0; }

  InlineSnarfNamedObjParams::InlineSnarfNamedObjParams(XMLReader& xml_in, const std::string& path) 
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
  InlineSnarfNamedObjParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Ids
    Chroma::write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }


  void 
  InlineSnarfNamedObj::operator()(unsigned long update_no,
				  XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out, "snarf_named_obj");
    write(xml_out, "update_no", update_no);

    // Snarf the object
    QDPIO::cout << InlineSnarfNamedObjEnv::name << ": snarf default gauge field" << endl;
    try
    {
      XMLBufferWriter file_xml, record_xml;
      multi1d<LatticeColorMatrix> u;

      // Get the default field
      InlineDefaultGaugeField::get(u, file_xml, record_xml);

      TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(params.named_obj.object_id);
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.object_id) = u;
      TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
    }
    catch (std::bad_cast) 
    {
      QDPIO::cerr << InlineSnarfNamedObjEnv::name << ": cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineSnarfNamedObjEnv::name << ": error message: " << e 
		  << endl;
      QDP_abort(1);
    }
    
    QDPIO::cout << InlineSnarfNamedObjEnv::name << ": ran successfully" << endl;

    pop(xml_out);  // snarf_named_obj

    END_CODE();
  } 

};
