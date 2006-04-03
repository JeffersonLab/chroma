// -*- C++ -*-
// $Id: enum_simplebctype_io.h,v 3.0 2006-04-03 04:58:56 edwards Exp $
/*! \file
 * \brief Simple BC type enum
 */
#ifndef enum_simplebctype_io_h
#define enum_simplebctype_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"



namespace Chroma {

  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */
  //! SimpleBC  type
  enum SimpleBCType {
    BC_TYPE_ANTIPERIODIC = -1,  // Antiperiodic = -1
    BC_TYPE_DIRICHLET,          // Dirichlet    =  0
    BC_TYPE_PERIODIC         // Periodic     =  1
  };


  namespace SimpleBCTypeEnv { 
    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<SimpleBCType> > theSimpleBCTypeMap;

  // Reader and writer

  //! Read an SimpleBCType enum
  void read(XMLReader& r, const string& path, SimpleBCType& t);

  //! Write an SimpleBCType enum
  void write(XMLWriter& w, const string& path, const SimpleBCType& t);

  /*! @} */   // end of group io

};
#endif
