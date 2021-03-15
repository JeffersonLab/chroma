// -*- C++ -*-
/*! \file
 * \brief Eigenstd::vector type enum
 */
#ifndef enum_eigenvectype_io_h
#define enum_eigenvectype_io_h

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
  //! Eigenstd::vector type
  enum EigenVecType {
    EVEC_TYPE_SCIDAC = 2,
    EVEC_TYPE_SZIN
  };


  namespace EigenVecTypeEnv { 
    extern const std::string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef Chroma::SingletonHolder<EnumTypeMap<EigenVecType> > theEigenVecTypeMap;

  // Reader and writer

  //! Read an eigenstd::vector enum
  void read(XMLReader& r, const std::string& path, EigenVecType& t);

  //! Write an eigenstd::vector enum
  void write(XMLWriter& w, const std::string& path, const EigenVecType& t);

  /*! @} */   // end of group io
}
#endif
