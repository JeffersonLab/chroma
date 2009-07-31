#ifndef LAPH_GAUGE_CONFIGURATION_HANDLER_H
#define LAPH_GAUGE_CONFIGURATION_HANDLER_H

#include "meas/inline/io/named_objmap.h"
#include "gauge_configuration_info.h"
#include "chromabase.h"


namespace Chroma {
  namespace LaphEnv {


// ******************************************************************************

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


class GaugeConfigurationHandler
{

      // pointer to the info about the gauge config (internally
      // managed by this handler)
      
    const GaugeConfigurationInfo* gauge_info;

      // pointer to the gauge field (external: in NamedObjMap)
 
    const multi1d<LatticeColorMatrix>* cfg;
    
    
      // prevent copying
    GaugeConfigurationHandler(const GaugeConfigurationHandler& u);
    GaugeConfigurationHandler& operator=(const GaugeConfigurationHandler& u);


  public:

    GaugeConfigurationHandler();
    
    GaugeConfigurationHandler(XMLReader& xml_in);

    void setInfo(XMLReader& xml_in);

    ~GaugeConfigurationHandler();
    
    void clear();
    

    void setData();
    
    
      // access to the gauge configuration and its info

    const multi1d<LatticeColorMatrix>& getData() const {return *cfg;}

    const GaugeConfigurationInfo& getInfo() const {return *gauge_info;}

    string outputInfo() const;

    bool isInfoSet() const { return (gauge_info!=0);}

    bool isDataSet() const { return (cfg!=0);}

  

};

// *************************************************************************
  }
}
#endif
