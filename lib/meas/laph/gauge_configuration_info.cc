//The definitions for the GaugeConfiguration class

#include "meas/laph/gauge_configuration_info.h"
#include "chroma.h"
#include <sstream>

namespace Chroma
{
   namespace LaphEnv
   {

      GaugeConfigurationInfo::GaugeConfigurationInfo(
				const std::string& gauge_id_in) : gauge_id(gauge_id_in)
      {
            //Test reference to the cfg in the named object map
            XMLBufferWriter gauge_xml_buff;
            try
            {
               TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(gauge_id);
               TheNamedObjMap::Instance().get(gauge_id).getRecordXML(gauge_xml_buff);
            }
            catch( std::bad_cast ) 
            {
               QDPIO::cerr << __func__ << ": caught dynamic cast error" 
                  << endl;
               QDP_abort(1);
            }
            catch (const string& e) 
            {
               QDPIO::cerr << __func__ << ": map call failed: " << e 
                  << endl;
               QDP_abort(1);
            }
            
            gauge_xml = gauge_xml_buff.str();

            XMLReader gauge_rdr(gauge_xml_buff);

            if (gauge_rdr.count("/Params/MCControl/StartUpdateNum") != 0)
            {
               read(gauge_rdr, "/Params/MCControl/StartUpdateNum" , traj_num);
            }
            else
            {
               traj_num = 1000;
            }

            QDPIO::cout << "TrajNum = " << traj_num << endl;
            QDPIO::cout << "Gauge Cfg Initialized" << endl;
							
			}
	 
	 }
}

