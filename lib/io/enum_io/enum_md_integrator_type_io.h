// -*- C++ -*-
// $Id: enum_md_integrator_type_io.h,v 3.0 2006-04-03 04:58:56 edwards Exp $
/*! \file
 * \brief MD integrator enum
 */
#ifndef enum_md_integrator_type_io_h
#define enum_md_integrator_type_io_h

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
  //! MDIntegrator  type
  enum MDIntegratorType {
    MD_PQP_LEAPFROG,
    MD_QPQ_LEAPFROG
  };


  namespace MDIntegratorTypeEnv { 
    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<MDIntegratorType> > theMDIntegratorTypeMap;

  // Reader and writer

  //! Read an MD Integrator Type enum
  void read(XMLReader& r, const string& path, MDIntegratorType& t);

  //! Write an MD Integrator Type enum
  void write(XMLWriter& w, const string& path, const MDIntegratorType& t);

  /*! @} */   // end of group io

};
#endif
