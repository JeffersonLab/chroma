// -*- C++ -*-
/*! \file
 * \brief Type of contractions for stochastic meson operators
 */
#ifndef enum_meson_op_io_h
#define enum_meson_op_io_h

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
  //! Meson operator contraction orderings
  enum MesonOpType 
  {
    MESON_OP_TYPE_SOURCE_SOURCE = 0,
    MESON_OP_TYPE_SOURCE_SOLUTION,
    MESON_OP_TYPE_SOLUTION_SOURCE,
    MESON_OP_TYPE_SOLUTION_SOLUTION,
  };


  namespace MesonOpTypeEnv 
  {
    extern const std::string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<MesonOpType> > theMesonOpTypeMap;

  // Reader and writer
  //! Read an approximation coefficient type enum
  void read(XMLReader& r, const std::string& path, MesonOpType& t);

  //! Write an approximation coefficient type enum
  void write(XMLWriter& w, const std::string& path, const MesonOpType& t);

  /*! @} */   // end of group io

};
#endif
