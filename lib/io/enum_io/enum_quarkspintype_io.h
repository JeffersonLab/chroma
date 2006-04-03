// -*- C++ -*-
// $Id: enum_quarkspintype_io.h,v 3.0 2006-04-03 04:58:56 edwards Exp $

/*! \file
 * \brief Enum for what spin components of a quark prop to compute
 *
 */

#ifndef enum_quarkspintype_io_h
#define enum_quarkspintype_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"


namespace Chroma 
{
  // QuarkSpinType --------------------------------------
  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */
 
  //! Quark spin type
  /*! \ingroup io */
  enum QuarkSpinType 
  {
    QUARK_SPIN_TYPE_FULL,
    QUARK_SPIN_TYPE_UPPER,
    QUARK_SPIN_TYPE_LOWER
  };


  //! Quark spin type env
  /*! \ingroup io */
  namespace QuarkSpinTypeEnv 
  { 
    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  //! A singleton to hold the typemap
  /*! \ingroup io */
  typedef SingletonHolder<EnumTypeMap<QuarkSpinType> > theQuarkSpinTypeMap;

  // Reader and writer
  //! Read a quark spin type enum
  /*! \ingroup io */
  void read(XMLReader& r, const string& path, QuarkSpinType& t);

  //! Write a quark spin type enum
  /*! \ingroup io */
  void write(XMLWriter& w, const string& path, const QuarkSpinType& t);

  /*! @} */   // end of group io

}
#endif
