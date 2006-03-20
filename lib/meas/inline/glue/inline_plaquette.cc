#include "meas/inline/glue/inline_plaquette.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/default_gauge_field.h"


namespace Chroma { 

  namespace InlinePlaquetteEnv { 

    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					   const std::string& path) {

      InlinePlaquetteParams p(xml_in, path);
      return new InlinePlaquette(p);
    }

    const std::string name = "PLAQUETTE";

    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);

  };


  // Params
  InlinePlaquetteParams::InlinePlaquetteParams() 
  { 
    frequency = 0; 
    named_obj.gauge_id = InlineDefaultGaugeField::getId();
  }

  InlinePlaquetteParams::InlinePlaquetteParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Ids
      named_obj.gauge_id = InlineDefaultGaugeField::readGaugeId(paramtop, "NamedObject/gauge_id");
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void 
  InlinePlaquette::operator()(unsigned long update_no,
			      XMLWriter& xml_out) 
  {
    START_CODE();
    
    QDPIO::cout << InlinePlaquetteEnv::name << ": gauge_id = XX" 
		<< params.named_obj.gauge_id
		<< "XX" << endl;

    // Test and grab a reference to the gauge field
    XMLBufferWriter gauge_xml;
    try
    {
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlinePlaquetteEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlinePlaquetteEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "Plaquette");
    write(xml_out, "update_no", update_no);

    Double w_plaq, s_plaq, t_plaq, link;
    MesPlq(u, w_plaq, s_plaq, t_plaq, link);
    write(xml_out, "w_plaq", w_plaq);
    write(xml_out, "s_plaq", s_plaq);
    write(xml_out, "t_plaq", t_plaq);
    write(xml_out, "link", link);
    
    pop(xml_out); // pop("Plaquette");
    
    END_CODE();
  } 

};
