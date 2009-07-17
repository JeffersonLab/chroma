
//    Class "GaugeConfigurationHandler" manages the actual gauge configuration.
//    It contains a "GaugeConfigurationInfo" that holds the gauge_xml as well
//    as a reference to the actual cfg.
//
//      Usage:
//
//         GaugeConfigurationHandler u(GaugeConfigurationInfo);   
//               --> first call, get date from NamedObjectMap; 
//                   set up static parameters.
//               --> any subsequent call to the constructor causes an abort!!
//               --> basically this is a simple singleton
//
//
//         u.GetData();       <-- reference to the actual configuration
//         u.GetInfo();       <-- reference to the GaugeConfigurationInfo
//                                       


#ifndef laph_gauge_configuration_handler_h
#define laph_gauge_configuration_handler_h

#include "meas/inline/io/named_objmap.h"
#include "meas/laph/gauge_configuration_info.h"
#include "chromabase.h"



namespace Chroma
{
   namespace LaphEnv
   {
      class GaugeConfigurationHandler
      {

         public:

            //This constructor will get the cfg from the named_obj map and
            //set the static ref to the gauge_field
            GaugeConfigurationHandler(const GaugeConfigurationInfo& gauge_in);

            const multi1d<LatticeColorMatrix>& getData() const {return *(cfg);}

            const GaugeConfigurationInfo& getInfo() const {return gauge_info;}
	    
         private:
	 
            //Static reference to the gauge field 
            static multi1d<LatticeColorMatrix>* cfg;

						//Should this be static?
						const GaugeConfigurationInfo& gauge_info;

      };

   }

}



#endif
