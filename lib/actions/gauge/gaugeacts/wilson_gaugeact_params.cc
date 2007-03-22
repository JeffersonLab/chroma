#include "chromabase.h"
#include "actions/gauge/gaugeacts/wilson_gaugeact_params.h"

namespace Chroma {


  WilsonGaugeActParams::WilsonGaugeActParams(XMLReader& xml_in, const std::string& path) {
    XMLReader paramtop(xml_in, path);

    try {
      read(paramtop, "./beta", beta);

      //  Read optional anisoParam.
      if (paramtop.count("AnisoParam") != 0) 
	read(paramtop, "AnisoParam", aniso);
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " <<  e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, const string& path, WilsonGaugeActParams& p) {
    WilsonGaugeActParams tmp(xml, path);
    p=tmp;
  }
}
