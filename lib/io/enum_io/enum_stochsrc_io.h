// -*- C++ -*-
// $Id: enum_stochsrc_io.h,v 3.0 2006-04-03 04:58:56 edwards Exp $
/*! \file
 * \brief Stochastic source enum
 */
#ifndef enum_stochsrc_io_h
#define enum_stochsrc_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"
#include "meas/hadron/enum_loops_s.h" // This is where the enum is defined. We are non-intrusive


namespace Chroma {

  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */

  namespace StochSrcEnv { 
    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<VolSrc> > theStochSrc ;

  // Reader and writer

  //! Read an InvType enum

  void read(XMLReader& r, const string& path, VolSrc& t);

  //! Write an InvType enum
  void write(XMLWriter& w, const string& path, const VolSrc& t);

  /*! @} */   // end of group io

};
#endif
