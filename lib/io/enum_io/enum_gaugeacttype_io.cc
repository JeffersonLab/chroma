// -*- C++ -*-
/*! \file
 * \brief Gauge action enum
 */
#include "enum_gaugeacttype_io.h"

/* *********!!!!!!!!!!!!!!!!!!! WARNING WARNING !!!!!!!!!!!!!!!! *********** */
// These actions are only relevant to Al Hart's topological code.
// The production HMC etc does not use these enums but works off its
// own factory
/* *********!!!!!!!!!!!!!!!!!!! WARNING WARNING !!!!!!!!!!!!!!!! *********** */

#include <string>

namespace Chroma { 

  namespace GaugeActTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theGaugeActTypeMap::Instance().registerPair(std::string("WILSON"), GAUGE_ACT_WILSON );
      success = theGaugeActTypeMap::Instance().registerPair(std::string("SYMZK_1X2"), GAUGE_ACT_SYMZK_1X2 );
      success = theGaugeActTypeMap::Instance().registerPair(std::string("IWASAKI"), GAUGE_ACT_IWASAKI );
      success = theGaugeActTypeMap::Instance().registerPair(std::string("DBW2"), GAUGE_ACT_DBW2 );
      success = theGaugeActTypeMap::Instance().registerPair(std::string("3_LOOP_IMP"), GAUGE_ACT_3_LOOP_IMP );
      success = theGaugeActTypeMap::Instance().registerPair(std::string("4_LOOP_IMP"), GAUGE_ACT_4_LOOP_IMP );
      success = theGaugeActTypeMap::Instance().registerPair(std::string("5_LOOP_IMP"), GAUGE_ACT_5_LOOP_IMP );
      return success;
    }

    //!  The number of different shaped loops in largest action
    const int No_fmn = 5;

    const std::string typeIDString = "GaugeActType";
    bool registered = registerAll();
  };
  using namespace GaugeActTypeEnv;

  //! Read an GaugeActType enum
  void read(XMLReader& xml_in,  const std::string& path, GaugeActType& t) {
    theGaugeActTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an GaugeActType enum
  void write(XMLWriter& xml_out, const std::string& path, const GaugeActType& t) {
    theGaugeActTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
};
