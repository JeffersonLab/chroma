// -*- C++ -*-
// $Id: enum_gaugeacttype_io.h,v 3.0 2006-04-03 04:58:56 edwards Exp $
/*! \file
 * \brief Gauge action enum
 */
#ifndef enum_gaugeacttype_io_h
#define enum_gaugeacttype_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"

/* *********!!!!!!!!!!!!!!!!!!! WARNING WARNING !!!!!!!!!!!!!!!! *********** */
// These actions are only relevant to Al Hart's topological code.
// The production HMC etc does not use these enums but works off its
// own factory
/* *********!!!!!!!!!!!!!!!!!!! WARNING WARNING !!!!!!!!!!!!!!!! *********** */

namespace Chroma {

  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */
  //! GaugeAct type
  enum GaugeActType {
    GAUGE_ACT_WILSON = 0,
    GAUGE_ACT_SYMZK_1X2,
    GAUGE_ACT_IWASAKI,
    GAUGE_ACT_DBW2,
    GAUGE_ACT_5_LOOP_IMP,
    GAUGE_ACT_4_LOOP_IMP,
    GAUGE_ACT_3_LOOP_IMP
    
  };

  

  namespace GaugeActTypeEnv {
    //!  The number of different shaped loops in largest action
    extern const int No_fmn;  // Al's own private definition

    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<GaugeActType> > theGaugeActTypeMap;

  // Reader and writer

  //! Read a GaugeActType enum
  void read(XMLReader& r, const string& path, GaugeActType& t);

  //! Write an GaugeActType enum
  void write(XMLWriter& w, const string& path, const GaugeActType& t);

  /*! @} */   // end of group io
};
#endif
