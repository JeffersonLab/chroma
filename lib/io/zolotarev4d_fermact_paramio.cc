#include "chromabase.h"
#include "io/param_io.h"
#include "io/zolotarev4d_fermact_paramio.h"

using namespace QDP;
using namespace std;

void read(XMLReader& xml_in, const string& path, Zolotarev4DStateInfo& info)
{
  XMLReader in(xml_in, path);

  try { 
    read(in, "NWilsVec", info.NWilsVec);
  }
  catch(const string& e) { 
    QDPIO::cerr << "Caught exception : " << e << endl;
    QDP_abort(1);
  }

  if( info.NWilsVec == 0 ) {
    // No eigenbeasties
    info.eigen_io.eigen_file="";
    info.eigen_io.eigen_volfmt = QDPIO_SINGLEFILE;

    // Try reading approx min and optional approx max
    // if approx max is not specified set it to 2*Nd
    try { 
      read(in, "ApproxMin", info.ApproxMin);
      if( in.count("ApproxMax") == 0 ) {
	info.ApproxMax = 2*Nd;
      }
      else {
	read(in, "ApproxMax", info.ApproxMax);
      }
    }
    catch( const string& e) {
      QDPIO::cerr << "Caught exception : " << e << endl;
    }
  }
  else {
    // We have eigenbeasties
    info.ApproxMin=0;
    info.ApproxMax=0;

    // Read in the eigenvector IO params
    try { 
      read(in, "Eig", info.eigen_io);
    }
    catch( const string& e ) { 
      QDPIO::cerr << "Caught exception: " << e << endl;
    }
  }
}

void write(XMLWriter& xml_out, const string& path, const Zolotarev4DStateInfo& info)
{
  if( path != "." ) { 
    push(xml_out, path);
  }

  write(xml_out, "NWilsVec", info.NWilsVec);

  if( info.NWilsVec == 0 ) {
    write(xml_out, "ApproxMin", info.ApproxMin);
    write(xml_out, "ApproxMax", info.ApproxMax);
  }
  else {
    write(xml_out, "Eig", info.eigen_io);
  }

  if( path != "." ) { 
    pop(xml_out);
  }

}

Zolotarev4DFermActParams::Zolotarev4DFermActParams(XMLReader& in)
{

  // This will bomb if it fails so no need to check for NULL after
  AuxFermActHandle = read(in, "AuxFermAct");
  
  try { 

    read(in, "Mass", Mass);
    read(in, "RatPolyDeg", RatPolyDeg);
    read(in, "InnerSolve/MaxCG", MaxCGInner);
    read(in, "InnerSolve/RsdCG", RsdCGInner);
    if( in.count("InnerSolve/ReorthFreq") == 1 ) {
      read(in, "InnerSolve/ReorthFreq", ReorthFreqInner);
    }
    else {
      ReorthFreqInner = 10; // Some default
    }

    read(in, "StateInfo", StateInfo);
  }
  catch( const string &e ) {
    QDPIO::cerr << "Caught Exception reading Zolo4D Fermact params: " << e << endl;
    QDP_abort(1);
  }
}

void write(XMLWriter& xml_out, const string& path, const Zolotarev4DFermActParams& p)
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
  push(xml_out, "InnerSolve");
  write(xml_out, "MaxCG", p.MaxCGInner);
  write(xml_out, "RsdCG", p.RsdCGInner);
  write(xml_out, "ReorthFreq", p.ReorthFreqInner);
  pop(xml_out);

  write(xml_out, "StateInfo", p.StateInfo);

  if( path != "." ) { 
    pop(xml_out);
  }
}

