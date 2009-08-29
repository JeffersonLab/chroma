#ifndef LAPH_GAUGE_CONFIGURATION_HANDLER_H
#define LAPH_GAUGE_CONFIGURATION_HANDLER_H

#include "meas/inline/io/named_objmap.h"
#include "gauge_configuration_info.h"
#include "chromabase.h"


namespace Chroma {
  namespace LaphEnv {


// *******************************************************************
// *                                                                 *
// *  Class "GaugeConfigurationHandler" manages information about    *
// *  and access to the gauge configuration in the Laph environment. *
// *  It contains a "GaugeConfigurationInfo" that holds the gauge    *
// *  header, as well as a reference to the actual cfg.              *
// *                                                                 *
// *    Basic usage:                                                 *
// *                                                                 *
// *      (a) declare a GaugeConfigurationHandler                    *
// *      (b) set the info  -- via constructor or setInfo(..)        *
// *      (c) set the data  -- via setData(..)                       *
// *      (d) use getData(..) for access to gauge configuration      *
// *                                                                 *
// *     -- getInfo(..) provides access to configuration info        *
// *     -- can be cleared and info/data reset                       *
// *                                                                 *
// *    Example:                                                     *
// *       XMLReader xml_in;                                         *
// *       GaugeConfigurationHandler uHandler;                       *
// *       uHandler.setInfo(xml_in);                                 *
// *       uHandler.setData();                                       *
// *       const multi1d<LatticeColorMatrix>& U=uHandler.getData();  *
// *    --> access to config through U                               *
// *                                                                 *
// *******************************************************************


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

    void setInfo(const std::string& header);

    ~GaugeConfigurationHandler();
    
    void clear();
    
    void setData();
    
    
      // access to the gauge configuration and its info

    const multi1d<LatticeColorMatrix>& getData();

    const GaugeConfigurationInfo& getInfo() const;

    std::string outputInfo() const;

    std::string getGaugeConfigHeader() const;

    void outputInfo(XMLWriter& xmlout) const;

    void getGaugeConfigHeader(XMLWriter& xmlout) const;

    bool isInfoSet() const { return (gauge_info!=0);}

    bool isDataSet() const { return (cfg!=0);}


  private:

    void set_info(XMLReader& xml_in);  

};

// *************************************************************************
  }
}
#endif
