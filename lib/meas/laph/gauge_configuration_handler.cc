#include "gauge_configuration_handler.h"
#include "chroma.h"

namespace Chroma {
  namespace LaphEnv {



GaugeConfigurationHandler::GaugeConfigurationHandler()
     : gauge_info(0), cfg(0) {}


GaugeConfigurationHandler::GaugeConfigurationHandler(XMLReader& xml_in)
{
 setInfo(xml_in);
 QDPIO::cout << "GaugeConfigurationHandler constructed"<<endl;
}


void GaugeConfigurationHandler::setInfo(XMLReader& xml_in)
{
 clear();
 try{
    gauge_info = new GaugeConfigurationInfo(xml_in);}
 catch(...){
    QDPIO::cerr << "problem allocating GaugeConfigurationInfo"<<endl;
    QDP_abort(1);}    
 QDPIO::cout << "GaugeConfig info set in GaugeConfigurationHandler"<<endl;
}

void GaugeConfigurationHandler::setData()
{
 if (!isInfoSet()){
    QDPIO::cerr << "error in GaugeConfigurationHandler:"<<endl;
    QDPIO::cerr << "  must setInfo before setData"<<endl;
    QDP_abort(1);}
 string gauge_id = gauge_info->getGaugeId();
 if (gauge_id == ""){
    QDPIO::cerr << "empty gauge_id in GaugeConfigurationHandler" << endl;
    QDP_abort(1);}

            //Assign the cfg from the named object map
 XMLBufferWriter gauge_xml_buff;
 try{
    TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(gauge_id);
    TheNamedObjMap::Instance().get(gauge_id).getRecordXML(gauge_xml_buff);
    }
 catch( std::bad_cast ){
    QDPIO::cerr << __func__ << ": caught dynamic cast error" << endl;
    QDP_abort(1);}
 catch (const string& err){
    QDPIO::cerr << __func__ << ": map call failed: " << err << endl;
    QDP_abort(1);}

 cfg = &(TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(gauge_id));

}


GaugeConfigurationHandler::~GaugeConfigurationHandler()
{
 clear();
}
    
void GaugeConfigurationHandler::clear()
{
 try {delete gauge_info;} catch(...) {QDP_abort(1);}
 gauge_info=0;
 cfg=0;
 QDPIO::cout << "GaugeConfigurationHandler cleared"<<endl;
}


string GaugeConfigurationHandler::outputInfo() const
{
 return gauge_info->output();
}


// **********************************************************************
  }
}

