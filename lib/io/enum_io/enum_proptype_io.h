// -*- C++ -*-
/*! \file
 * \brief PropType enum
 */
#ifndef enum_proptype_io_h
#define enum_proptype_io_h

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
  //! Propagator type
  enum PropType {
    PROP_TYPE_SCIDAC = 2,
    PROP_TYPE_SZIN,
    PROP_TYPE_KYU
  };


  namespace PropTypeEnv { 
    extern const std::string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef Chroma::SingletonHolder<EnumTypeMap<PropType> > thePropTypeMap;

  // Reader and writer

  //! Read a propagator type enum
  void read(XMLReader& r, const std::string& path, PropType& t);

  //! Write a propagator type enum
  void write(XMLWriter& w, const std::string& path, const PropType& t);

  /*! @} */   // end of group io
}
#endif
