#include "quark_action_info.h"
#include "chromabase.h"
#include <sstream>

namespace Chroma {
  namespace LaphEnv {


// ************************************************************

QuarkActionInfo::QuarkActionInfo(const QuarkActionInfo& rhs)
                : ferm_act_xml(rhs.ferm_act_xml), mass(rhs.mass), 
                  mass_name(rhs.mass_name), action_id(rhs.action_id) {}
						
QuarkActionInfo& QuarkActionInfo::operator=(const QuarkActionInfo& rhs)
{
 ferm_act_xml = rhs.ferm_act_xml;
 mass = rhs.mass;
 mass_name = rhs.mass_name;
 action_id = rhs.action_id;
 return *this;
}

QuarkActionInfo::QuarkActionInfo(XMLReader& act_rdr)
{ 
 try {
    XMLReader temp(act_rdr, "//FermionAction");
    if (temp.count("Mass") != 0){
       mass_name="Mass";
       read(temp, "Mass", mass);}
    else{
       mass_name="Kappa";
       read(temp, "Kappa", mass);}
    ostringstream strm;
    temp.print(strm);
    ferm_act_xml = strm.str();
    QDPIO::cout << "Ferm Act XML = XX" << ferm_act_xml << "XX" << endl;
    read(temp, "//FermAct", action_id);
    QDPIO::cout << "Action ID = " <<  action_id << endl;
    QDPIO::cout << "Mass/Kappa = " << mass << endl;
    QDPIO::cout << "Action Info Initialized" << endl;
    }
 catch(std::string& err){
    QDPIO::cerr << "ERROR: " << err 
         << ", Fermion Action Info could not be initialized" << endl;
    QDP_abort(1);}
}


void QuarkActionInfo::checkEqual(const QuarkActionInfo& rhs) const
{
 if (!xmlContentIsEqual(ferm_act_xml, rhs.ferm_act_xml))
    throw string(__func__)+": Fermion Action xml contents do not match";
}


void QuarkActionInfo::matchXMLverbatim(const QuarkActionInfo& rhs) const
{
 if (rhs.ferm_act_xml != ferm_act_xml)
    throw string(__func__)+": Fermion Action xml does not match";
}


bool QuarkActionInfo::operator==(const QuarkActionInfo& rhs) const
{
 return xmlContentIsEqual(ferm_act_xml,rhs.ferm_act_xml);
}



string QuarkActionInfo::output(int indent) const
{
 if (indent==0) return ferm_act_xml;
 else{
    string pad(3*indent,' ');
    string temp(pad);
    int pos1=0;
    int pos2=ferm_act_xml.find('\n',pos1);
    while (pos2!=string::npos){
       temp+=ferm_act_xml.substr(pos1,pos2-pos1+1)+pad;
       pos1=pos2+1;
       pos2=ferm_act_xml.find('\n',pos1);}
    temp+=ferm_act_xml.substr(pos1,ferm_act_xml.length()-pos1+1);
    return temp;}
}

// ***********************************************************
  }
}
