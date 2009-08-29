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
// *  GaugeConfigurationInfo. This will also return the trajectory *
// *  number as well as some other info about the gauge            *
// *  configuration.                                               *
// *                                                               *
// *  XML input:                                                   *
// *                                                               *
// *   Recall that the chroma XML input file must have the form    *
// *                                                               *
// *    <chroma>                                                   *
// *      <annotation> Job title </annotation>                     *
// *      <RNG> .... </RNG>  (optional)                            *
// *      <Param>                                                  *
// *         <nrow>X Y Z T</nrow>                                  *
// *         <InlineMeasurements>                                  *
// *            <elem>  .... </elem>   <-- each inline measurement *
// *         </InlineMeasurements>                                 *
// *      </Param>                                                 *
// *      <Cfg>                                                    *
// *         <cfg_type>CONFIG_FILE_TYPE</cfg_type>                 *
// *         <cfg_file>cfg_file_name</cfg_file>                    *
// *      </Cfg>                                                   *
// *    </chroma>                                                  *
// *                                                               *
// *   The gauge configuration specified in the <Cfg> tag gets     *
// *   the gauge_id = "default_gauge_id".                          *
// *                                                               *
// *  There are two constructors to GaugeConfigurationInfo: one    *
// *  that takes an XMLReader, and another that takes a string.    *
// *                                                               *
// *    XMLReader xml_in(...);                                     *
// *    GaugeConfigurationInfo U(xml_in);                          *
// *                                                               *
// *      --> This constructor expects XML of the form             *
// *               <GaugeConfigurationInfo>                        *
// *                  <gauge_id>....</gauge_id>                    *
// *               </GaugeConfigurationInfo>                       *
// *          then information about the configuration is obtained *
// *          from TheNamedObjMap and the gauge_header string is   *
// *          formed.  The gauge_header string is surrounded by    *
// *          the tag <GaugeConfigHeader>...</GaugeConfigHeader>.  *
// *                                                               *
// *    string headerInfo(...);                                    *
// *    GaugeConfigurationInfo U(headerInfo);                      *
// *                                                               *
// *      --> This version of the constructor assumes that header  *
// *          information from a quark_source_sink file, for       *
// *          example, is passed in and requires the string        *
// *          contains the gauge_header string, which is extracted *
// *          verbatim and has tag <GaugeConfigHeader>.            *
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

  std::string file_name;
  std::string config_type;
  std::string gauge_id;
  int traj_num;
  int time_dir;
  int time_extent;
  multi1d<int> extents;
  int number_dir;

 public:

  GaugeConfigurationInfo(XMLReader& xml_rdr);
         
  GaugeConfigurationInfo(const std::string& header);

  GaugeConfigurationInfo(const GaugeConfigurationInfo& rhs);
                 
  GaugeConfigurationInfo& operator=(const GaugeConfigurationInfo& rhs);

  ~GaugeConfigurationInfo(){}
  
  void checkEqual(const GaugeConfigurationInfo& rhs) const;

  bool operator==(const GaugeConfigurationInfo& rhs) const;




  std::string output(int indent = 0) const;

  void output(XMLWriter& xmlout) const;

  int getTrajNum() const { return traj_num; }

  std::string getGaugeId() const { return gauge_id; }

  std::string getGaugeConfigHeader() const { return output(0); }

  void getGaugeConfigHeader(XMLWriter& xmlout) const { output(xmlout); }

  std::string getFullRecordXML() const;

  int getTimeExtent() const { return time_extent; }

  int getTimeDir() const { return time_dir; }
  
  int getNumberOfDirections() const { return number_dir; }


 private:

  void set_info1(XMLReader& xmlg);
  void set_info2(XMLReader& xmlg);

};


// **********************************************************
  }
}
#endif
