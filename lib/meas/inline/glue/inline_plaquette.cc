#include "meas/inline/glue/inline_plaquette.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"


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

  void 
  InlinePlaquette::operator()(const multi1d<LatticeColorMatrix>& u,
			      const unsigned long update_no,
			      XMLWriter& xml_out) 
  {
    
    push(xml_out, "Plaquette");
    write(xml_out, "update_no", update_no);

    Double w_plaq, s_plaq, t_plaq, link;
    MesPlq(u, w_plaq, s_plaq, t_plaq, link);
    write(xml_out, "w_plaq", w_plaq);
    write(xml_out, "s_plaq", s_plaq);
    write(xml_out, "t_plaq", t_plaq);
    write(xml_out, "link", link);
    
    pop(xml_out); // pop("Plaquette");
    
  } 

};
