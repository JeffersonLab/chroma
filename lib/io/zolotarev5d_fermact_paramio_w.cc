#include "chromabase.h"
#include "io/param_io.h"
#include "io/zolotarev5d_fermact_paramio.h"

using namespace QDP;
using namespace std;


Zolotarev5DFermActParams::Zolotarev5DFermActParams(XMLReader& in)
{

  // This will bomb if it fails so no need to check for NULL after
  AuxFermActHandle = read(in, "AuxFermAct");
  if( AuxFermActHandle == 0x0 ) { 
    QDPIO::cerr << "Read of AuxFermAct returned NULL Pointer" << endl;
    QDP_abort(1);
  }

  try { 

    read(in, "Mass", Mass);
    read(in, "RatPolyDeg", RatPolyDeg);
    read(in, "StateInfo", StateInfo);
  }
  catch( const string &e ) {
    QDPIO::cerr << "Caught Exception reading Zolo5D Fermact params: " << e << endl;
    QDP_abort(1);
  }
}

void write(XMLWriter& xml_out, const string& path, const Zolotarev5DFermActParams& p)
{
  if ( path != "." ) { 
    push( xml_out, path);
  }
  
  write(xml_out, "FermAct", p.getFermActType());

  switch(p.AuxFermActHandle->getFermActType()) {
  case FERM_ACT_WILSON:
    {
      const WilsonFermActParams& wilson = dynamic_cast<WilsonFermActParams&>(*(p.AuxFermActHandle));
      write(xml_out, "AuxFermAct", wilson);
    }
    break;
  default:
    QDPIO::cerr << "Unsupported AuxFermAct in write " << endl;
    QDP_abort(1);
    break;
  }

  write(xml_out, "Mass", p.Mass);
  write(xml_out, "RatPolyDeg", p.RatPolyDeg);
  write(xml_out, "StateInfo", p.StateInfo);

  if( path != "." ) { 
    pop(xml_out);
  }
}

