//The definitions for the GaugeConfiguration class

#include "inverter_info.h"
#include "chromabase.h"
#include <sstream>
using namespace std;


namespace Chroma {
  namespace LaphEnv {



InverterInfo::InverterInfo(const InverterInfo& rhs) 
             : inverter_xml(rhs.inverter_xml), 
               tol(rhs.tol), id(rhs.id) {}
  
                                        
InverterInfo& InverterInfo::operator=(const InverterInfo& rhs)
{
 inverter_xml = rhs.inverter_xml;
 tol = rhs.tol;
 id = rhs.id;
 return *this;
}


InverterInfo::InverterInfo(XMLReader& xml_in)
{
 try{
    XMLReader xmlr(xml_in, "./descendant-or-self::InvertParam");
    read(xmlr, "RsdCG", tol);
    read(xmlr, "invType", id);
    ostringstream strm;
    xmlr.print(strm);
    inverter_xml = strm.str();
    QDPIO::cout << "Inverter XML:"<<endl << inverter_xml << endl;
    QDPIO::cout << "Tolerance = " << tol << endl;
    QDPIO::cout << "Inverter ID = "<<id << endl;
    QDPIO::cout << "Inverter Info Initialized" << endl;
    }
 catch(std::string& err){
    QDPIO::cerr << "ERROR: "<<err<<endl;
    QDPIO::cerr << "could not initialize InverterInfo from XML input"<<endl;
    QDP_abort(1);}
}


void InverterInfo::checkEqualTolerance(const InverterInfo& rhs) const
{
 if (rhs.tol != tol)
    throw string(__func__)+": Inverter tolerances do not match";
}

void InverterInfo::checkEqual(const InverterInfo& rhs) const
{
 if (!xmlContentIsEqual(inverter_xml,rhs.inverter_xml))
    throw string("Inverters do not match");
}

void InverterInfo::matchXMLverbatim(const InverterInfo& rhs) const
{
 if (inverter_xml!=rhs.inverter_xml)
    throw string("Inverter XML strings do not exactly match");
}

bool InverterInfo::operator==(const InverterInfo& rhs) const
{
 return xmlContentIsEqual(inverter_xml,rhs.inverter_xml);
}


string InverterInfo::output(int indent) const
{
 if (indent==0) return inverter_xml;
 else{
    string pad(3*indent,' ');
    string temp(pad);
    int pos1=0;
    int pos2=inverter_xml.find('\n',pos1);
    while (pos2!=string::npos){
       temp+=inverter_xml.substr(pos1,pos2-pos1+1)+pad;
       pos1=pos2+1;
       pos2=inverter_xml.find('\n',pos1);}
    temp+=inverter_xml.substr(pos1,inverter_xml.length()-pos1+1);
    return temp;}
}


void InverterInfo::output(XMLWriter& xmlout) const
{
 stringstream oss;
 oss << inverter_xml;
 XMLReader xmlr(oss);
 xmlout << xmlr;
}

// *******************************************************************
  }
}
