#include "chromabase.h"

#include "io/param_io.h"
#include "io/fermact_paramio.h"
#include "io/zolotarev4d_fermact_paramio.h"

#include <iostream>
#include <string>

using namespace std;
using namespace QDP;

FermActParams* read(XMLReader &reader, const string& path) 
{
  // Create reader with paths staring at input path
  XMLReader top(reader, path);

  enum FermActType f_a_type;

  try { 
    read( top, "FermAct", f_a_type );
  }
  catch( const string& e ) { 
    QDPIO::cerr << "Unable to read FermAct: " << e << endl;
    QDP_abort(1);
  }

  // Now read the rest of the FermAct depending on type
  switch( f_a_type ) { 
  case FERM_ACT_WILSON: 
    return new WilsonFermActParams(top);
    break;
  case FERM_ACT_UNPRECONDITIONED_WILSON:
    return new WilsonFermActParams(top);
    break;
  case FERM_ACT_DWF:
    return new DWFFermActParams(top);
    break;
  case FERM_ACT_UNPRECONDITIONED_DWF:
    return new DWFFermActParams(top);
    break;
  case FERM_ACT_OVERLAP_DWF:
    return new DWFFermActParams(top);
    break;
  case FERM_ACT_ZOLOTAREV_4D:
    return new Zolotarev4DFermActParams(top);
    break;
  default: 
    QDPIO::cerr << "As yet unsupported FermionAction " << f_a_type << endl;
    QDP_abort(1);
  }

}

// Write out is easy. Just call virtual function
void write(XMLWriter &xml_out, const string& path, const FermActParams& p)
{
  switch( p.getFermActType() ) { 

    // Deliberate!!!
  case FERM_ACT_WILSON: // Fall through
  case FERM_ACT_UNPRECONDITIONED_WILSON: // Wilson case
    {
      const WilsonFermActParams& wilson=dynamic_cast<const WilsonFermActParams&>(p);
      write(xml_out, path, wilson);
    }
      break;
    

  case FERM_ACT_DWF:  // Fall through
  case FERM_ACT_UNPRECONDITIONED_DWF:  // Fall through
  case FERM_ACT_OVERLAP_DWF: // DWF Case
    {
      const DWFFermActParams& dwf=dynamic_cast<const DWFFermActParams&>(p);
      write(xml_out, path, dwf);
    }
    break;
    
  case FERM_ACT_ZOLOTAREV_4D:
    {
      const Zolotarev4DFermActParams& zolo4d = dynamic_cast<const Zolotarev4DFermActParams&>(p);
      write(xml_out, path, zolo4d);
    }
    break;

    defaule:
    QDPIO::cerr << "Unknown ferm action type : " << p.getFermActType() << endl;
    QDP_abort(1);
  }
}


// Now constructors:
WilsonFermActParams::WilsonFermActParams( XMLReader& xml_in ) { 

  // Re-get the FermAct
  try { 
    read(xml_in, "FermAct", my_fermact_type );
  }
  catch ( const string& e ) { 
    QDPIO::cerr << "Unable to read FermAct" << e << endl;
    QDP_abort(1);
  }

  // Check FermAct is supported
  switch( my_fermact_type ) { 
  case FERM_ACT_WILSON: 
    break;
  case FERM_ACT_UNPRECONDITIONED_WILSON:
    break;
  default: 
    QDPIO::cerr << "As yet unsupported WilsonFermAct " << my_fermact_type << endl;
    QDP_abort(1);
  }

  // Read the stuff for the action
  if (xml_in.count("Mass") != 0) {

    read(xml_in, "Mass", Mass);

    if (xml_in.count("Kappa") != 0) {
    
      QDPIO::cerr << "Error: found both a Kappa and a Mass tag" << endl;
      QDP_abort(1);

    }
  }
  else if (xml_in.count("Kappa") != 0)  {

    Real Kappa;

    read(xml_in, "Kappa", Kappa);

    Mass = kappaToMass(Kappa);    // Convert Kappa to Mass
  }
  else {
  
    QDPIO::cerr << "Error: neither Mass or Kappa found" << endl;
    QDP_abort(1);
  }    

  // There is always an aniso Param for wilson, so set it to default
  initHeader(anisoParam);

  //  Read optional anisoParam.
  if (xml_in.count("AnisoParam") != 0) {
    read(xml_in, "AnisoParam", anisoParam);
  }
 
}

// Now constructors:
DWFFermActParams::DWFFermActParams( XMLReader& xml_in ) { 

  // Re-get the FermAct
  try { 
    read(xml_in, "FermAct", my_fermact_type );
  }
  catch ( const string& e ) { 
    QDPIO::cerr << "Unable to read FermAct" << e << endl;
    QDP_abort(1);
  }

  // Check FermAct is supported
  switch( my_fermact_type ) { 
  case FERM_ACT_DWF: 
    break;
  case FERM_ACT_UNPRECONDITIONED_DWF:
    break;
  case FERM_ACT_OVERLAP_DWF:
    break;
  default: 
    QDPIO::cerr << "As yet unsupported DWFFermAct " << my_fermact_type << endl;
    QDP_abort(1);
  }

  // Read the stuff for the action
  if (xml_in.count("Mass") != 0) {

    read(xml_in, "Mass", Mass);

    if (xml_in.count("Kappa") != 0) {
    
      QDPIO::cerr << "Error: found both a Kappa and a Mass tag" << endl;
      QDP_abort(1);

    }
  }
  else if (xml_in.count("Kappa") != 0)  {

    Real Kappa;

    read(xml_in, "Kappa", Kappa);

    Mass = kappaToMass(Kappa);    // Convert Kappa to Mass
  }
  else {
  
    QDPIO::cerr << "Error: neither Mass or Kappa found" << endl;
    QDP_abort(1);
  }    

  // There is always an aniso Param for wilson, so set it to default
  initHeader(anisoParam);

  //  Read optional anisoParam.
  if (xml_in.count("AnisoParam") != 0) {
    read(xml_in, "AnisoParam", anisoParam);
  }
 
  initHeader(chiralParam);
  // Read optional chiralParam (actually not optional for DWF)
  if (xml_in.count("ChiralParam") != 0) {
    read(xml_in, "ChiralParam", chiralParam);
  }
  else {
    QDPIO::cerr << "chiralParams missing for DWF action" << endl;
    QDP_abort(1);
  }

}
    

// Now Writers
void write(XMLWriter&  xml_out, const string& path, const WilsonFermActParams& p)
{
  if( path != "." ) { 
    push(xml_out, path);
  }

  write(xml_out, "FermAct", p.getFermActType());
  write(xml_out, "Mass", p.Mass);
  write(xml_out, "AnisoParam", p.anisoParam);

  if( path != "." ) { 
    pop(xml_out);
  }

}

void write(XMLWriter& xml_out, const string& path, const DWFFermActParams& p)
{
  if( path != "." ) { 
    push(xml_out, path);
  }

  write(xml_out, "FermAct", p.getFermActType());
  write(xml_out, "Mass", p.Mass);
  write(xml_out, "AnisoParam", p.anisoParam);
  write(xml_out, "ChiralParam",p.chiralParam);
  if( path != "." ) { 
    pop(xml_out);
  }
}


      
