#ifndef enum_innersolver_type_io_h
#define enum_innersolver_type_io_h

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
  //! OverlapInnerSolver type
  enum OverlapInnerSolverType { 
    OVERLAP_INNER_CG_SINGLE_PASS,
    OVERLAP_INNER_CG_DOUBLE_PASS
  };

  namespace OverlapInnerSolverTypeEnv { 
    extern const string typeIDString;
    extern const bool registered; 
    const bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<OverlapInnerSolverType> > theOverlapInnerSolverTypeMap;

  // Reader and writer

  //! Read a OverlapInnerSolverType enum
  void read(XMLReader& r, const string& path, OverlapInnerSolverType& t);

  //! Write an OverlapInnerSolverType enum
  void write(XMLWriter& w, const string& path, const OverlapInnerSolverType& t);

};
#endif
