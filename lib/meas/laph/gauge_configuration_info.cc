//The definitions for the GaugeConfiguration class

#include "gauge_configuration_info.h"
#include "chroma.h"
#include <sstream>

namespace Chroma {
  namespace LaphEnv {
  
 

GaugeConfigurationInfo::GaugeConfigurationInfo(XMLReader& xml_rdr)
{
 if (xml_rdr.count("//gauge_id") == 1){
 
    read(xml_rdr, "//gauge_id", gauge_id);

        //Test reference to the cfg in the named object map
    XMLBufferWriter gauge_xml_buff;
    try{   
       TheNamedObjMap::Instance().get(gauge_id).getRecordXML(gauge_xml_buff);
       }
    catch( std::bad_cast ){
       QDPIO::cerr << __func__ << ": caught dynamic cast error" << endl;
       QDP_abort(1);
       }
    catch (const string& err){
       QDPIO::cerr << __func__ << ": map call failed: " << err << endl;
       QDP_abort(1);
       }
    XMLReader xmlg(gauge_xml_buff);
    
    if (xmlg.count("//MCControl")==1){
       XMLReader xmlg1(xmlg,"//MCControl");
       XMLReader xmlg2(xmlg,"//HMCTrj");
       ostringstream oss;
       oss << "<GaugeConfigParams>";
       xmlg1.printCurrentContext(oss);
       xmlg2.printCurrentContext(oss);
       oss << "</GaugeConfigParams>";
       gauge_xml = oss.str();}
    else{
       XMLReader temp_rdr(xmlg, "/");
       ostringstream temp_strm;
       temp_strm << "<GaugeConfigParams>";
       temp_rdr.printCurrentContext(temp_strm);
       temp_strm << "</GaugeConfigParams>";
       gauge_xml = temp_strm.str();}
    }
 else{
    gauge_id="";
    XMLReader xmlg(xml_rdr,"//GaugeConfigParams");
    stringstream oss;
    xml_rdr.print(oss);
    gauge_xml=oss.str();
    }

 stringstream tmp;
 tmp << gauge_xml;
 XMLReader gauge_rdr(tmp);   
    
 if (gauge_rdr.count("/GaugeConfigParams/MCControl/StartUpdateNum") != 0){
    read(gauge_rdr, "/GaugeConfigParams/MCControl/StartUpdateNum" , traj_num);
    read(gauge_rdr, "/GaugeConfigParams/HMCTrj/MDIntegrator/t_dir", time_dir);}
 else{
    QDPIO::cerr << "WARNING: No trajectory number, time_dir found!" << endl;
    traj_num = 1000;
    time_dir = 3;}
 
 number_dir=QDP::Nd;
 time_extent=QDP::Layout::lattSize()[time_dir];

 QDPIO::cout << "gauge_xml = XX " <<  gauge_xml << "XX" << endl;
 QDPIO::cout << "TrajNum = " << traj_num << endl;
 QDPIO::cout << "Gauge Cfg Initialized" << endl;
}


    // copy constructor
    
GaugeConfigurationInfo::GaugeConfigurationInfo(const GaugeConfigurationInfo& rhs) 
         : gauge_xml(rhs.gauge_xml), gauge_id(rhs.gauge_id), 
           traj_num(rhs.traj_num), time_dir(rhs.time_dir), 
           time_extent(rhs.time_extent), number_dir(rhs.number_dir) {}


GaugeConfigurationInfo& GaugeConfigurationInfo::operator=(const GaugeConfigurationInfo& rhs)
{
 gauge_xml = rhs.gauge_xml;
 traj_num = rhs.traj_num;
 gauge_id = rhs.gauge_id;
 time_dir = rhs.time_dir;
 time_extent = rhs.time_extent;
 number_dir = rhs.number_dir;
 return *this;
}

void GaugeConfigurationInfo::checkEqual(const GaugeConfigurationInfo& rhs) const
{
 if (!xmlContentIsEqual(gauge_xml,rhs.gauge_xml))
    throw string("gauge info contents do not match");
}

void GaugeConfigurationInfo::matchXMLverbatim(const GaugeConfigurationInfo& rhs) const
{
 if (rhs.gauge_xml!=gauge_xml)
    throw string("gauge xml strings do not match");
}


bool GaugeConfigurationInfo::operator==(const GaugeConfigurationInfo& rhs) const
{return xmlContentIsEqual(gauge_xml,rhs.gauge_xml);}


string GaugeConfigurationInfo::output(int indent) const
{
 if (indent==0) return gauge_xml;
 else{
    string pad(3*indent,' ');
    string temp(pad);
    int pos1=0;
    int pos2=gauge_xml.find('\n',pos1);
    while (pos2!=string::npos){
       temp+=gauge_xml.substr(pos1,pos2-pos1+1)+pad;
       pos1=pos2+1;
       pos2=gauge_xml.find('\n',pos1);}
    temp+=gauge_xml.substr(pos1,gauge_xml.length()-pos1+1);
    return temp;}
}


// ******************************************************************
  }
}
