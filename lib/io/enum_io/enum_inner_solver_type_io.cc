// -*- C++ -*-
/*! \file
 * \brief Inner solver enum
 */

#include "enum_inner_solver_type_io.h"

#include <string>

namespace Chroma { 

  namespace OverlapInnerSolverTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theOverlapInnerSolverTypeMap::Instance().registerPair(std::string("SINGLE_PASS"),
								      OVERLAP_INNER_CG_SINGLE_PASS );
      success &= theOverlapInnerSolverTypeMap::Instance().registerPair(std::string("DOUBLE_PASS"),
								       OVERLAP_INNER_CG_DOUBLE_PASS );
      return success;
    }
    const std::string typeIDString = "OverlapInnerSolverType";
    bool registered = registerAll();
  };
  using namespace OverlapInnerSolverTypeEnv;

  //! Read an OverlapInnerSolverType enum
  void read(XMLReader& xml_in,  const std::string& path, OverlapInnerSolverType& t) {
    theOverlapInnerSolverTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an OverlapInnerSolverType enum
  void write(XMLWriter& xml_out, const std::string& path, const OverlapInnerSolverType& t) {
    theOverlapInnerSolverTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
};
