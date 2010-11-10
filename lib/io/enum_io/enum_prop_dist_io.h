// -*- C++ -*-
/*! \file
 * \brief Type of lines for distillution
 */
#ifndef enum_prop_dist_io_h
#define enum_prop_dist_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"

namespace Chroma 
{

  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */
  //! Distillution line types
  enum PropDistType 
  {
    PROP_DIST_TYPE_SOURCE = 19,
    PROP_DIST_TYPE_SOLUTION = 41,
  };


  namespace PropDistTypeEnv 
  {
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder< EnumTypeMap<PropDistType> > thePropDistTypeMap;

  //! Reader
  void read(XMLReader& r, const string& path, PropDistType& t);

  //! Writer
  void write(XMLWriter& w, const string& path, const PropDistType& t);

  //! Reader
  void read(BinaryReader& r, PropDistType& t);

  //! Writer
  void write(BinaryWriter& w, const PropDistType& t);

  /*! @} */   // end of group io

};
#endif
