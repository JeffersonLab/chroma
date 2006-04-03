// -*- C++ -*-
// $Id: enum_fermtype_io.h,v 3.0 2006-04-03 04:58:56 edwards Exp $
/*! \file
 * \brief FermType enum
 */

#ifndef enum_fermtype_io_h
#define enum_fermtype_io_h

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
  //! Fermion  type
  enum FermType {
    FERM_TYPE_WILSON,
    FERM_TYPE_STAGGERED
  };


  namespace FermTypeEnv { 
    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<FermType> > theFermTypeMap;

  // Reader and writer

  //! Read an Fermion Type enum
  void read(XMLReader& r, const string& path, FermType& t);

  //! Write an Fermion Type enum
  void write(XMLWriter& w, const string& path, const FermType& t);

  /*! @} */   // end of group io

};
#endif
