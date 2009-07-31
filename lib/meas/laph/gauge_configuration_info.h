#ifndef LAPH_GAUGE_CONFIGURATION_INFO_H
#define LAPH_GAUGE_CONFIGURATION_INFO_H

#include "chromabase.h"

namespace Chroma {
  namespace LaphEnv {

// *****************************************************************
// *                                                               *
// *  Class "GaugeConfigurationInfo" holds information about a     *
// *  gauge configuration.  It also checks that its config info    *
// *  is the same when compared against another object of class    *
// *  GaugeConfigurationInfo. This is a very simple wrapper over   *
// *  a gauge_xml, that will also return the trajectory number     *
// *  as well as some other info about the gauge configuration.    *
// *                                                               *
// *  Usage:                                                       *
// *                                                               *
// *    GaugeConfigurationInfo u(XMLReader& xml_in);               *
// *                                                               *
// *       --> If "gauge_id" is found in "xml_in", then this       *
// *           info corresponds to an actual configuration.        *
// *           "gauge_id" is assigned, as well as "gauge_xml",     *
// *           and the trajectory number is extracted from the     *
// *           named object map, as well as other info.            *
// *       --> If not found, then this is metadata from            *
// *           some quantity (such as a quark source/sink)         *
// *           other than a configuration.  "gauge_id" is set      *
// *           to the empty string but the "gauge_xml" is          *
// *           assigned.                                           *
// *                                                               *
// *    GaugeConfigurationInfo u2(xml_in);                         *
// *    u.checkEqual(u2); --> checks that u2 and u have same xml   *
// *                           content; throws exception if not    *
// *    u.matchXMLverbatim(u2);  --> same as above, but uses a     *
// *                                 character by character        *
// *                               comparison of the xml strings   *
// *                                                               *
// *    string out = u.output();   // xml output                   *
// *    string out = u.output(2);  // indented xml output          *
// *    int j = u.getTrajNum();    // returns  number of RHMC      *
// *                               // trajectory for this config   *
// *    string s = u.getGaugeId(); // gauge_id string              *              
// *    int nt = u.getTimeExtent(); // time extent of lattice      *
// *    int tdir = u.getTimeDir();  // index of time               *
// *                                                               *
// *****************************************************************


class GaugeConfigurationInfo
{

  std::string gauge_xml;    // private data
  std::string gauge_id;
  int traj_num;
  int time_dir;
  int time_extent;
  int number_dir;

 public:

  GaugeConfigurationInfo(XMLReader& xml_rdr);
         
  GaugeConfigurationInfo(const GaugeConfigurationInfo& rhs);
                 
  GaugeConfigurationInfo& operator=(const GaugeConfigurationInfo& rhs);

  ~GaugeConfigurationInfo(){}
  
  void checkEqual(const GaugeConfigurationInfo& rhs) const;

  bool operator==(const GaugeConfigurationInfo& rhs) const;

  void matchXMLverbatim(const GaugeConfigurationInfo& rhs) const;



  std::string output(int indent = 0) const;

  int getTrajNum() const { return traj_num; }

  std::string getGaugeId() const { return gauge_id; }

  int getTimeExtent() const { return time_extent; }

  int getTimeDir() const { return time_dir; }
  
  int getNumberOfDirections() const { return number_dir; }

};


// **********************************************************
  }
}
#endif
