// -*- C++ -*-
/*! \file
 * \brief Coeffs enum
 */
#ifndef enum_coeffs_io_h
#define enum_coeffs_io_h

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
  //! Coeffs  type
  enum CoeffType {
      COEFF_TYPE_ZOLOTAREV = 0,
      COEFF_TYPE_TANH,
      COEFF_TYPE_TANH_UNSCALED
  };


  namespace CoeffTypeEnv { 
    extern const std::string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<CoeffType> > theCoeffTypeMap;

  // Reader and writer
  //! Read an approximation coefficient type enum
  void read(XMLReader& r, const std::string& path, CoeffType& t);

  //! Write an approximation coefficient type enum
  void write(XMLWriter& w, const std::string& path, const CoeffType& t);

  /*! @} */   // end of group io

};
#endif
