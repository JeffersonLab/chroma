#include "chromabase.h"
#include "gaugebc.h"

#include "actions/gauge/gaugebcs/gaugebc_schroedinger.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"

namespace Chroma {

  GaugeBCSchrParams::GaugeBCSchrParams(XMLReader& xml, const std::string& path)
  {
    try {
      read(xml, "SchrFunType", SchrFun);
      read(xml, "SchrPhiMult", SchrPhiMult);
    } catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " << e << endl;
      QDP_abort(1);
    }  
  }

  void read(XMLReader& xml, const std::string& path, GaugeBCSchrParams& p)   {
    GaugeBCSchrParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, const GaugeBCSchrParams& p) {
    push(xml, path);
    write(xml, "SchrFunType", p.SchrFun);
    write(xml, "SchrPhiMult", p.SchrPhiMult);
    pop(xml);
  }

};
