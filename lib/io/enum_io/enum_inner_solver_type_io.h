// -*- C++ -*-
// $Id: enum_inner_solver_type_io.h,v 3.0 2006-04-03 04:58:56 edwards Exp $
/*! \file
 * \brief Inner solver enum
 */
#ifndef enum_innersolver_type_io_h
#define enum_innersolver_type_io_h

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
  //! OverlapInnerSolver type
  enum OverlapInnerSolverType { 
    OVERLAP_INNER_CG_SINGLE_PASS,
    OVERLAP_INNER_CG_DOUBLE_PASS
  };

  namespace OverlapInnerSolverTypeEnv { 
    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<OverlapInnerSolverType> > theOverlapInnerSolverTypeMap;

  // Reader and writer

  //! Read a OverlapInnerSolverType enum
  void read(XMLReader& r, const string& path, OverlapInnerSolverType& t);

  //! Write an OverlapInnerSolverType enum
  void write(XMLWriter& w, const string& path, const OverlapInnerSolverType& t);

  /*! @} */   // end of group io
};
#endif
