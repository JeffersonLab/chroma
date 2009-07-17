
//    Class "GaugeConfigurationInfo" holds information about a  
//    gauge configuration.
//    It also can check that its config is the same when compared
//    against another GaugeConfigurationInfo. This is a very simple 
//    wrapper over a gauge_xml, that will also return the trajectory number.
//
//      Usage:
//
//         GaugeConfigurationInfo u(string gauge_id);   
//               -->Checks that this id in the NamedObjectMap is valid, 
//                  and extracts the trajectory number
//
//         u.check(GaugeConfigurationInfo);  
//               --> checks that another config has the same xml 
//
//         u.output();        <-- outputs the gauge xml 
//         u.getTrajNum();    <-- returns the integer number of the 
//                                RHMC trajectory for this configuration
//         u.getGaugeId();    <-- outputs the gauge_id                              


#ifndef laph_gauge_configuration_info_h
#define laph_gauge_configuration_info_h

#include "chromabase.h"

namespace Chroma
{

   namespace LaphEnv
   {

      class GaugeConfigurationInfo
      {


         public:

            GaugeConfigurationInfo(const std::string& gauge_id_in);
           
						GaugeConfigurationInfo( 
								const GaugeConfigurationInfo& rhs) : gauge_xml(rhs.output()),
								  gauge_id(rhs.getGaugeId()), traj_num(rhs.getTrajNum()) {}
						
						GaugeConfigurationInfo& operator=(const GaugeConfigurationInfo& rhs)
						{
							gauge_xml = rhs.output();
							traj_num = rhs.getTrajNum();
							gauge_id = rhs.getGaugeId();
							
							return *this;
						}

            const std::string& output() const { return gauge_xml;}

            const int getTrajNum() const {return traj_num;}
      
						const std::string& getGaugeId() const {return gauge_id; }

						void check(const GaugeConfigurationInfo& rhs) const
						{
							if (rhs.output() != gauge_xml)
							{
								QDPIO::cerr << __func__ << 
									"gauge xml does not match original XML." << endl;
								QDP_abort(1);
							}
						}

         private:
	 
            std::string gauge_xml;

						std::string gauge_id;

            int traj_num;
      };

   }

}



#endif
