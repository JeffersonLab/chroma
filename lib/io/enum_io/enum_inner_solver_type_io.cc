#include "enum_inner_solver_type_io.h"

#include <string>
using namespace std;
using namespace Chroma;

namespace Chroma { 

  namespace OverlapInnerSolverTypeEnv { 

    const bool registerAll(void) 
    {
      bool success; 
      success = theOverlapInnerSolverTypeMap::Instance().registerPair(string("SINGLE_PASS"),
								      OVERLAP_INNER_CG_SINGLE_PASS );
      success &= theOverlapInnerSolverTypeMap::Instance().registerPair(string("DOUBLE_PASS"),
								       OVERLAP_INNER_CG_DOUBLE_PASS );
      return success;
    }
    const string typeIDString = "OverlapInnerSolverType";
    const bool registered = registerAll();
  };
  using namespace OverlapInnerSolverTypeEnv;

  //! Read an OverlapInnerSolverType enum
  void read(XMLReader& xml_in,  const string& path, OverlapInnerSolverType& t) {
    theOverlapInnerSolverTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an OverlapInnerSolverType enum
  void write(XMLWriter& xml_out, const string& path, const OverlapInnerSolverType& t) {
    theOverlapInnerSolverTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
};
