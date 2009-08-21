//The definitions for the GaugeConfiguration class

#include "gauge_configuration_info.h"
#include "chroma.h"
#include <sstream>

namespace Chroma {
  namespace LaphEnv {
  

// *************************************************************

   //  This constructor expects XML of the form
   //     <GaugeConfigurationInfo>
   //       <gauge_id>....</gauge_id>
   //     </GaugeConfigurationInfo>
   //  then information about the configuration is obtained
   //  from TheNamedObjMap and the gauge_header string is formed.

GaugeConfigurationInfo::GaugeConfigurationInfo(XMLReader& xml_in)
{
 if (xml_tag_count(xml_in,"GaugeConfigurationInfo")!=1){
    QDPIO::cerr << "Bad XML input to GaugeConfigurationInfo"<<endl;
    QDPIO::cerr << "Expected one <GaugeConfigurationInfo> tag"<<endl;
    QDP_abort(1);}

      //  Get "gauge_id" from xml_in, then find it in the
      //  NamedObjMap to extract the info

 XMLReader xmlr(xml_in,"./descendant-or-self::GaugeConfigurationInfo");
 xmlread(xmlr, "gauge_id", gauge_id, "GaugeConfigurationInfo");
 gauge_id=tidyString(gauge_id);
 if (gauge_id.empty()){
     QDPIO::cerr << "empty gauge_id in GaugeConfigurationInfo"<<endl;
     QDP_abort(1);}

        //Test reference to the cfg in the named object map
 XMLBufferWriter gauge_header_buff;
 try{   
    TheNamedObjMap::Instance().get(gauge_id).getRecordXML(gauge_header_buff);
    }
 catch( std::bad_cast ){
    QDPIO::cerr << "caught dynamic cast error in GaugeConfigurationInfo" << endl;
    QDP_abort(1);
    }
 catch (const string& err){
    QDPIO::cerr << "TheNamedObjMap call failed in GaugeConfigurationInfo: " << err << endl;
    QDP_abort(1);
    }
 XMLReader xmlg0(gauge_header_buff);   
 XMLReader xmlg(xmlg0,"/");         // required due to XMLReader bug

 set_info_from_gauge_header(xmlg);

 QDPIO::cout << endl << endl <<"GaugeConfigurationInfo constructor:"<<endl<<endl;
 QDPIO::cout << output()<<endl;
 QDPIO::cout << "GaugeConfigurationInfo Initialized" << endl<<endl;
}

    // *************************************************************

void GaugeConfigurationInfo::set_info_from_gauge_header(XMLReader& xmlg)
{
 if (xmlg.count(".//StartUpdateNum") != 0){
    xmlread(xmlg,"StartUpdateNum",traj_num,"GaugeConfigurationInfo");
    xmlread(xmlg,"HMCTrj/MDIntegrator/t_dir",time_dir,"GaugeConfigurationInfo");
    xmlread(xmlg,"cfg_type",config_type,"GaugeConfigurationInfo");
    xmlread(xmlg,"cfg_file",file_name,"GaugeConfigurationInfo");}
 else{
    traj_num = 1000;
    time_dir = QDP::Nd;
    config_type = "dummy";
    file_name = "none";}
 
 config_type=tidyString(config_type);
 file_name=tidyString(file_name);
 number_dir=QDP::Nd;
 time_extent=QDP::Layout::lattSize()[time_dir];
 extents=QDP::Layout::lattSize();
}

 // *************************************************************

   // This version of the constructor assumes that header information
   // from a quark_source_sink file, for example, is passed in and
   // gauge_header string can be extracted verbatim.

GaugeConfigurationInfo::GaugeConfigurationInfo(const string& header)
{
 string gauge_header;
 extract_xml_element(header,"GaugeConfigurationInfo",gauge_header,
                     "GaugeConfigurationInfo");
 set_info_from_header_string(gauge_header);
}

void GaugeConfigurationInfo::set_info_from_header_string(const string& header)
{
 stringstream oss;
 oss << header;
 XMLReader xmlg0(oss);
 XMLReader xmlg(xmlg0,"/");   // due to XMLReader bug

 xmlread(xmlg,"HMCTrajectoryNumber",traj_num,"GaugeConfigurationInfo");
 xmlread(xmlg,"TimeDir",time_dir,"GaugeConfigurationInfo");
 xmlread(xmlg,"ConfigType",config_type,"GaugeConfigurationInfo");
 xmlread(xmlg,"FileName",file_name,"GaugeConfigurationInfo");
 xmlread(xmlg,"GaugeId",gauge_id,"GaugeConfigurationInfo"); 
 xmlread(xmlg,"TimeExtent",time_extent,"GaugeConfigurationInfo");
 xmlread(xmlg,"NumberOfDir",number_dir,"GaugeConfigurationInfo");
 xmlread(xmlg,"LatticeExtents",extents,"GaugeConfigurationInfo");
}

  // ************************************************************

string GaugeConfigurationInfo::getFullRecordXML() const
{
 string gauge_xml;
 if (gauge_id.empty()) return gauge_xml;

        //Test reference to the cfg in the named object map
 XMLBufferWriter gauge_header_buff;
 try{   
    TheNamedObjMap::Instance().get(gauge_id).getRecordXML(gauge_header_buff);
    }
 catch( std::bad_cast ){
    QDPIO::cerr << "caught dynamic cast error in GaugeConfigurationInfo" << endl;
    QDP_abort(1);
    }
 catch (const string& err){
    QDPIO::cerr << "TheNamedObjMap call failed in GaugeConfigurationInfo: " << err << endl;
    QDP_abort(1);
    }
 XMLReader xmlg0(gauge_header_buff);   
 XMLReader xmlg(xmlg0,"/");         // required due to XMLReader bug

    // make the gauge_header string (the configuration could be
    // generated from Monte Carlo updating, or could be an
    // initial configuration from a cold or hot start)

 if (xml_tag_count(xmlg,"MCControl")>0){   // from Monte Carlo generated config
    try{
       XMLReader xmlg1(xmlg,".//MCControl");
       XMLReader xmlg2(xmlg,".//HMCTrj");
       ostringstream oss;
       oss << "<GaugeConfigHeader>"<<endl;
       xmlg1.print(oss);
       xmlg2.print(oss);
       oss << "</GaugeConfigHeader>";
       gauge_xml = oss.str();}
    catch(const string& err){
       QDPIO::cerr << "Could not create gauge_header in GaugeConfigurationInfo"<<endl;
       QDP_abort(1);}}
 else{
    ostringstream temp_strm;                 // from a dummy config
    temp_strm << "<GaugeConfigHeader>"<<endl;
    xmlg.printCurrentContext(temp_strm);
    temp_strm << "</GaugeConfigHeader>";
    gauge_xml = temp_strm.str();}

 return gauge_xml;
}

// *******************************************************************

    // copy constructor
    
GaugeConfigurationInfo::GaugeConfigurationInfo(
                const GaugeConfigurationInfo& rhs) 
         : file_name(rhs.file_name), gauge_id(rhs.gauge_id), 
           config_type(rhs.config_type), extents(rhs.extents),
           traj_num(rhs.traj_num), time_dir(rhs.time_dir), 
           time_extent(rhs.time_extent), number_dir(rhs.number_dir) {}


GaugeConfigurationInfo& GaugeConfigurationInfo::operator=(
               const GaugeConfigurationInfo& rhs)
{
 file_name = rhs.file_name;
 traj_num = rhs.traj_num;
 gauge_id = rhs.gauge_id;
 config_type = rhs.config_type;
 time_dir = rhs.time_dir;
 time_extent = rhs.time_extent;
 number_dir = rhs.number_dir;
 extents = rhs.extents;
 return *this;
}

void GaugeConfigurationInfo::checkEqual(const GaugeConfigurationInfo& rhs) const
{
 if  ((file_name!=rhs.file_name)
    ||(traj_num!=rhs.traj_num)
    ||(gauge_id!=rhs.gauge_id)
    ||(config_type!=rhs.config_type)
    ||(time_dir!=rhs.time_dir)
    ||(time_extent!=rhs.time_extent)
    ||(number_dir!=rhs.number_dir)
    ||(extents!=rhs.extents))
    throw string("GaugeConfigurationInfo checkEqual failed");
}


bool GaugeConfigurationInfo::operator==(const GaugeConfigurationInfo& rhs) const
{
 return ((file_name==rhs.file_name)
       &&(traj_num==rhs.traj_num)
       &&(gauge_id==rhs.gauge_id)
       &&(config_type==rhs.config_type)
       &&(time_dir==rhs.time_dir)
       &&(time_extent==rhs.time_extent)
       &&(number_dir==rhs.number_dir)
       &&(extents==rhs.extents));
}


string GaugeConfigurationInfo::output(int indent) const
{
 string pad(3*indent,' ');
 ostringstream oss;
 oss << pad << "<GaugeConfigurationInfo>"<<endl;
 oss << pad << "   <GaugeId>" << gauge_id << "</GaugeId>"<<endl;
 oss << pad << "   <FileName>" << file_name << "</FileName>"<<endl;
 oss << pad << "   <ConfigType>" << config_type << "</ConfigType>"<<endl;
 oss << pad << "   <HMCTrajectoryNumber>" << traj_num << "</HMCTrajectoryNumber>"<<endl;
 oss << pad << "   <TimeDir>" << time_dir << "</TimeDir>"<<endl;
 oss << pad << "   <TimeExtent>" << time_extent << "</TimeExtent>"<<endl;
 oss << pad << "   <NumberOfDir>" << number_dir << "</NumberOfDir>"<<endl;
 oss << pad << "   <LatticeExtents>";
 if (extents.size()>0) oss << extents[0];
 for (int k=1;k<extents.size();k++) oss << " "<<extents[k];
 oss << "</LatticeExtents>"<<endl;
 oss << pad << "</GaugeConfigurationInfo>"<<endl;
 return oss.str();
}


// ******************************************************************
  }
}
