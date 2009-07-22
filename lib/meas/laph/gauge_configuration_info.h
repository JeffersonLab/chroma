
//    Class "GaugeConfigurationInfo" holds information about a  
//    gauge configuration.
//    It also can check that its config is the same when compared
//    against another GaugeConfigurationInfo. This is a very simple 
//    wrapper over a gauge_xml, that will also return the trajectory number.
//
//      Usage:
//
//         GaugeConfigurationInfo u(XMLReader xml_in);   
//               -->First, checks if 'gauge_id' is a tag. If yes, uses the 
//               -->gauge_id to extract gauge_xml from the named_obj_map
//               -->If 'gauge_id' is not a tag, checks if xml_in contians a 
//               -->valid gauge_xml. Also extracts the trajectory number from
//               --> the gauge_xml
//
//         u.check(GaugeConfigurationInfo);  
//               --> checks that another config has the same xml 
//
//         u.output();        <-- outputs the gauge xml 
//         u.getTrajNum();    <-- returns the integer number of the 
//                                RHMC trajectory for this configuration
//         u.getGaugeId();    <-- outputs the gauge_id                          //         u.getTExtent()

#ifndef laph_gauge_configuration_info_h
#define laph_gauge_configuration_info_h

#include "chromabase.h"

namespace Chroma
{
  namespace LaphEnv
  {


 class GaugeConfigurationInfo
 {

   std::string gauge_xml;    // private data
   std::string gauge_id;
   int traj_num;
	
	 //Should get this from xml somehow, but for now just hard code it
	 static const decay_dir = 3;

  public:

   GaugeConfigurationInfo(XMLReader& xml_rdr);
          
   GaugeConfigurationInfo(const GaugeConfigurationInfo& rhs) 
          : gauge_xml(rhs.output()), gauge_id(rhs.getGaugeId()), 
            traj_num(rhs.getTrajNum()) {}
                  
                  GaugeConfigurationInfo& operator=(const GaugeConfigurationInfo& rhs)
                  {
                     gauge_xml = rhs.output();
                     traj_num = rhs.getTrajNum();
                     gauge_id = rhs.getGaugeId();
                     
                     return *this;
                  }

            const std::string& output() const { return gauge_xml;}

            const int getTrajNum() const {return traj_num;}
     
						const int getTimeDir() const {return decay_dir;}

						const int getTExtent() const {
							return QDP::Layout::lattsize()[decay_dir];}

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

      };

   }

}



#endif
