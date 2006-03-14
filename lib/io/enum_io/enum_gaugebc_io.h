// -*- C++ -*-
// $Id: enum_gaugebc_io.h,v 2.2 2006-03-14 04:54:18 edwards Exp $
/*! \file
 * \brief Gauge bc enum
 */

#error "THIS CODE IS OBSOLETE - DO NOT USE"

#ifndef enum_gaugebc_io_h
#define enum_gaugebc_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"
#include "meas/sources/srcsnktype.h"


namespace Chroma {

  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */
  /*************** GAUGEBC ********************************/
  //! Supported Gauge BC types
  enum GaugeBCType { 
    GAUGEBC_ALL_PERIODIC, 
    GAUGEBC_SCHROEDINGER_1LINK,
    GAUGEBC_SCHROEDINGER_2LINK,
    GAUGEBC_SIMPLE 
  };

  namespace GaugeBCTypeEnv { 
    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<GaugeBCType> > theGaugeBCTypeMap;

  //! Read an GaugeBCType enum
  void read(XMLReader& r, const string& path, GaugeBCType& t);

  //! Write an GaugeBCType enum
  void write(XMLWriter& w, const string& path, const GaugeBCType& t);

  /*************  SCHROEDINGER FUNCTIONAL TYPE  *******************/

  //! Schroedinger Functional type Boundary Conditions
  enum SchrFunType {
    SF_NONE = 0, 
    SF_TRIVIAL = 1,
    SF_NONPERT = 2,
    SF_COUPLING = 3,
    SF_CHROMOMAG = 4,
    SF_DIRICHLET = 10,
};


  namespace SchrFunTypeEnv { 
    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<SchrFunType> > theSchrFunTypeMap;

  // Reader and writer

  //! Read an SchrFunType enum
  void read(XMLReader& r, const string& path, SchrFunType& t);

  //! Write an SchrFunType enum
  void write(XMLWriter& w, const string& path, const SchrFunType& t);

  /*! @} */   // end of group io
};

#endif
