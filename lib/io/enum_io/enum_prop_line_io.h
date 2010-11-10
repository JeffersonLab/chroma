// -*- C++ -*-
/*! \file
 * \brief Type of lines for distillution
 */
#ifndef enum_prop_line_io_h
#define enum_prop_line_io_h

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
  enum PropLineType 
  {
    PROP_LINE_TYPE_CONN = 43,
    PROP_LINE_TYPE_DISC = 47,
  };


  namespace PropLineTypeEnv 
  {
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder< EnumTypeMap<PropLineType> > thePropLineTypeMap;

  //! Reader
  void read(XMLReader& r, const string& path, PropLineType& t);

  //! Writer
  void write(XMLWriter& w, const string& path, const PropLineType& t);

  //! Reader
  void read(BinaryReader& r, PropLineType& t);

  //! Writer
  void write(BinaryWriter& w, const PropLineType& t);

  /*! @} */   // end of group io

};
#endif
