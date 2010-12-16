// -*- C++ -*-
// $Id: enum_cfgtype_io.h,v 3.1 2006-04-16 03:06:49 edwards Exp $
/*! \file
 * \brief CfgType enum
 */

#ifndef enum_cfgtype_io_h
#define enum_cfgtype_io_h


#include "io/enum_io/enum_type_map.h"
#include "singleton.h"
#include <string>


namespace Chroma {
 
  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */
  //! Configuration type
  enum CfgType {
      CFG_TYPE_MILC = 0,
      CFG_TYPE_NERSC,
      CFG_TYPE_SCIDAC,
      CFG_TYPE_SZIN,
      CFG_TYPE_SZINQIO,
      CFG_TYPE_KYU,
      CFG_TYPE_DISORDERED,
      CFG_TYPE_UNIT,
      CFG_TYPE_CPPACS,
      CFG_TYPE_WEAK_FIELD,
      CFG_TYPE_CLASSICAL_SF,
      CFG_TYPE_WUPP,
  };


  namespace CfgTypeEnv { 
    extern bool registered; 
    extern const string typeIDString;
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<CfgType> > theCfgTypeMap;

  //! read a configuration type enum
  void read(XMLReader& xml_in, const string& path, CfgType& t);

  //! write a configuration type enum
  void write(XMLWriter& xml_out, const string& path, const CfgType& t); 

  /*! @} */   // end of group io

};

#endif
