// -*- C++ -*-
// $Id: enum_plusminus_io.h,v 1.1 2006-05-23 18:11:34 edwards Exp $
/*! \file
 * \brief Plus/Minus type enum from Chroma
 */
#ifndef enum_plusminus_io_h
#define enum_plusminus_io_h

#include "chromabase.h"   // This is where the enum is defined.
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
  //! PlusMinus type
  namespace PlusMinusEnv 
  { 
    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<PlusMinus> > thePlusMinusMap;

  // Reader and writer

  //! Read an PlusMinus enum
  void read(XMLReader& r, const string& path, PlusMinus& t);

  //! Write an PlusMinus enum
  void write(XMLWriter& w, const string& path, const PlusMinus& t);

  /*! @} */   // end of group io
};
#endif
