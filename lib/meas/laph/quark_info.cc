#include "quark_info.h"
#include "chromabase.h"
#include <sstream>

namespace Chroma {
  namespace LaphEnv {


// ************************************************************

QuarkInfo::QuarkInfo(XMLReader& xml_in)
{ 
 if (xml_tag_count(xml_in,"QuarkInfo")==1){
    XMLReader xmlr(xml_in, "./descendant-or-self::QuarkInfo");
    set_info1(xmlr);}
 else if (xml_tag_count(xml_in,"ValenceQuarkHeader")==1){
    XMLReader xmlr(xml_in, "./descendant-or-self::ValenceQuarkHeader");
    set_info2(xmlr);}
 else{
    QDPIO::cerr << "Bad XML input to QuarkInfo"<<endl;
    QDPIO::cerr << "Expected one <QuarkInfo> or <ValenceQuarkHeader> tag"<<endl;
    QDP_abort(1);}
}

void QuarkInfo::set_info1(XMLReader& xmlr)
{
 if (xml_tag_count(xmlr,"FermionAction")!=1){
    QDPIO::cerr << "QuarkInfo initialization error: no <FermionAction> tag"<<endl;
    QDP_abort(1);}
 ostringstream strm;
 strm << "<ValenceQuarkHeader>"<<endl;
 xmlr.printCurrentContext(strm);
 strm << "</ValenceQuarkHeader>";
 quark_header = strm.str();
 setMass(xmlr,mass_name,mass);
 xmlread(xmlr, "FermAct", action_id, "QuarkInfo");
// QDPIO::cout << endl << endl <<"QuarkInfo constructor:"<<endl<<endl;
// QDPIO::cout << "action id = "<<action_id<<endl;
// QDPIO::cout << "action_header: "<<endl<< quark_header << endl << endl;
// QDPIO::cout << mass_name <<" = "<<mass<<endl;
// QDPIO::cout << "QuarkInfo Initialized" << endl<<endl;
}

void QuarkInfo::setMass(XMLReader& xmlr, string& massName, double& massValue)
{
 if (xml_tag_count(xmlr,"Mass")==1){
    massName="Mass";
    xmlread(xmlr, "Mass", massValue, "QuarkInfo");}
 else if (xml_tag_count(xmlr,"Kappa")==1){
    massName="Kappa";
    xmlread(xmlr, "Kappa", massValue, "QuarkInfo");}
 else{
    QDPIO::cerr << "could not set mass in QuarkInfo"<<endl;
    QDP_abort(1);}
}


QuarkInfo::QuarkInfo(XMLReader& xml_in, const GaugeConfigurationInfo& U)
{ 
 if (xml_tag_count(xml_in,"ValenceQuarkHeader")==1){
    XMLReader xmlr(xml_in, "./descendant-or-self::ValenceQuarkHeader");
    set_info2(xmlr);
    return;}
 if (xml_tag_count(xml_in,"QuarkInfo")!=1){
    QDPIO::cerr << "Bad XML input to QuarkInfo"<<endl;
    QDPIO::cerr << "Expected one <QuarkInfo> tag"<<endl;
    QDP_abort(1);}
 XMLReader xmlr(xml_in, "./descendant-or-self::QuarkInfo");
 int dyn=xml_tag_count(xmlr,"Dynamical");
 if (dyn>1){
    QDPIO::cerr << "Bad XML input to QuarkInfo"<<endl;
    QDPIO::cerr << "Multiple <Dynamical> tags found"<<endl;
    QDP_abort(1);}
 if (dyn==0){
    ostringstream strm;
    strm << "<ValenceQuarkHeader>"<<endl;
    xmlr.printCurrentContext(strm);
    strm << "</ValenceQuarkHeader>";
    quark_header = strm.str();
    setMass(xmlr,mass_name,mass);
    xmlread(xmlr, "FermAct", action_id, "QuarkInfo");}
 else{
    setMass(xmlr,mass_name,mass);
    quark_header.clear();
    stringstream strm;
    strm << U.getFullRecordXML();
    XMLReader xmlu(strm);
    XMLReader xmlt(xmlu,"//HMCTrj");
    string mname;
    double mval;
    int nquarks=xml_tag_count(xmlt,"FermionAction");
    for (int k=1;k<=nquarks;k++){
       ostringstream elem;
       elem << "./descendant::FermionAction[" << k <<"]";
       XMLReader xmltt(xmlt,elem.str());
       setMass(xmltt,mname,mval);
       if ((mname==mass_name)&&(abs(mval-mass)<1e-8)){
          ostringstream temp;
          temp << "<ValenceQuarkHeader>"<<endl;
          xmltt.print(temp);
          temp << "</ValenceQuarkHeader>";
          quark_header=temp.str();
          xmlread(xmltt,"FermAct", action_id, "QuarkInfo");
          break;}}
    if (quark_header.empty()){
       QDPIO::cerr << "Could not initialize QuarkInfo:"<<endl
                   << " could not find requested dynamical quark mass"<<endl
                   << " in GaugeConfigurationInfo"<<endl;
       QDP_abort(1);}
    }
// QDPIO::cout << endl << endl <<"QuarkInfo constructor:"<<endl<<endl;
// QDPIO::cout << "action id = "<<action_id<<endl;
// QDPIO::cout << "action_header: "<<endl<< quark_header << endl << endl;
// QDPIO::cout << mass_name <<" = "<<mass<<endl;
// QDPIO::cout << "QuarkInfo Initialized" << endl<<endl;
}


   // This version of the constructor assumes that header information
   // from a quark_source_sink file, for example, is passed in and
   // quark_header string can be extracted verbatim.

QuarkInfo::QuarkInfo(const string& header)
{
 extract_xml_element(header,"ValenceQuarkHeader",quark_header,
                     "QuarkInfo");
 stringstream strm;
 strm << quark_header;
 XMLReader xmlt(strm);
 XMLReader xmlr(xmlt,"/");  // due to XMLReader bug
 setMass(xmlr,mass_name,mass);
 xmlread(xmlr, "FermAct", action_id, "QuarkInfo");
}

void QuarkInfo::set_info2(XMLReader& xmlr)
{
 setMass(xmlr,mass_name,mass);
 xmlread(xmlr, "FermAct", action_id, "QuarkInfo");
 ostringstream strm;
 xmlr.print(strm);
 quark_header = strm.str();
}

  // ************************************************************

     // copy constructor

QuarkInfo::QuarkInfo(const QuarkInfo& rhs)
                : quark_header(rhs.quark_header), mass(rhs.mass), 
                  mass_name(rhs.mass_name), action_id(rhs.action_id) {}

						
QuarkInfo& QuarkInfo::operator=(const QuarkInfo& rhs)
{
 quark_header = rhs.quark_header;
 mass = rhs.mass;
 mass_name = rhs.mass_name;
 action_id = rhs.action_id;
 return *this;
}


void QuarkInfo::checkEqual(const QuarkInfo& rhs) const
{
 if (!xmlContentIsEqual(quark_header, rhs.quark_header))
    throw string("QuarkInfo xml contents do not match");
}


void QuarkInfo::matchXMLverbatim(const QuarkInfo& rhs) const
{
 if (rhs.quark_header != quark_header)
    throw string("QuarkInfo xml strings do not match");
}


bool QuarkInfo::operator==(const QuarkInfo& rhs) const
{
 return xmlContentIsEqual(quark_header,rhs.quark_header);
}



string QuarkInfo::output(int indent) const
{
 if (indent==0) return quark_header;
 else{
    string pad(3*indent,' ');
    string temp(pad);
    int pos1=0;
    int pos2=quark_header.find('\n',pos1);
    while (pos2!=string::npos){
       temp+=quark_header.substr(pos1,pos2-pos1+1)+pad;
       pos1=pos2+1;
       pos2=quark_header.find('\n',pos1);}
    temp+=quark_header.substr(pos1,quark_header.length()-pos1+1);
    return temp;}
}

void QuarkInfo::output(XMLWriter& xmlout) const
{
 stringstream oss;
 oss << quark_header;
 XMLReader xmlr(oss);
 xmlout << xmlr;
}



// ***********************************************************
  }
}
