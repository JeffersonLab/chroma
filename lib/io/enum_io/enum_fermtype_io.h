#ifndef enum_fermtype_io_h
#define enum_fermtype_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"


using namespace std;
using namespace Chroma;

namespace Chroma {

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
    extern const bool registered; 
    const bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<FermType> > theFermTypeMap;

  // Reader and writer

  //! Read an Fermion Type enum
  void read(XMLReader& r, const string& path, FermType& t);

  //! Write an Fermion Type enum
  void write(XMLWriter& w, const string& path, const FermType& t);

};
#endif
