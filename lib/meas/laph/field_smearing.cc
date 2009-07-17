#include "field_smearing.h"
#include "xml_help.h"
#include <sstream>
using namespace std;

namespace Chroma {
  namespace LaphEnv {



   // context node of xml_in should be <smearing_scheme>

FieldSmearingInfo::FieldSmearingInfo(XMLReader& xml_in)
{
 if (!checkXMLCurrentTagName(xml_in,"smearing_scheme")){
    QDPIO::cerr << "invalid XML input to FieldSmearingInfo";
    QDP_abort(1);}
 string s1,s2;
 try{
    read(xml_in,"./link_smear_type", s1 );
    read(xml_in,"./quark_smear_type", s2 );
    read(xml_in,"./link_iterations", linkIterations );
    read(xml_in,"./link_staple_weight", linkStapleWeight );
    read(xml_in,"./number_laph_eigvecs", laphNumEigvecs );
    }
 catch(const string& err){
    QDPIO::cerr << "could not initialize FieldSmearingInfo from XML input"<<endl;
    QDP_abort(1);}
 if ((tidyString(s1)!="STOUT")||(tidyString(s2)!="LAPH")){
    QDPIO::cerr << "STOUT and LAPH smearing is required..."<<endl;
    QDP_abort(1);}
 if ((linkIterations<0)||(linkStapleWeight<0.0)||(laphNumEigvecs<1)){
    QDPIO::cerr << "invalid smearing scheme"<<endl;
    QDP_abort(1);}
}

string FieldSmearingInfo::output() const
{
 ostringstream oss;
 oss << "<smearing_scheme>"<<endl;
 oss << "  <link_smear_type> " << "STOUT" 
     << " </link_smear_type>"<<endl;
 oss << "  <link_iterations> " << linkIterations 
     << " </link_iterations>"<<endl;
 oss << "  <link_staple_weight> " << linkStapleWeight 
     << " </link_staple_weight>"<<endl;
 oss << "  <quark_smear_type> " << "LAPH" 
     << " </quark_smear_type>"<<endl;
 oss << "  <number_laph_eigvecs> " << laphNumEigvecs 
     << " </number_laph_eigvecs>"<<endl;
 oss << "</smearing_scheme>"<<endl;
 return oss.str();
}


// *************************************************************
  }
}
