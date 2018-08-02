#ifndef enum_quda_io_h
#define enum_quda_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"
#include <quda.h>
#include <map>
namespace Chroma
{
  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */

  //! Quda Solver type
  enum QudaSolverType { 
    CG,
    BICGSTAB,
    GCR,
	CA_GCR,
    MR
  };

  
  namespace QudaSolverTypeEnv { 
    extern const std::string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<QudaSolverType> > theQudaSolverTypeMap;
  typedef SingletonHolder<std::map<QudaSolverType,QudaInverterType> > theChromaToQudaSolverTypeMap;

  // Reader and writer

  //! Read an QudaSolverType enum
  void read(XMLReader& r, const std::string& path, QudaSolverType& t);

  //! Write an QudaSolverType enum
  void write(XMLWriter& w, const std::string& path, const QudaSolverType& t);


  //! Quda Precision type
  enum QudaPrecisionType { 
    DEFAULT,
    QUARTER,
    HALF,
    SINGLE,
    DOUBLE
  };
  
  namespace QudaPrecisionTypeEnv { 
    extern const std::string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<QudaPrecisionType> > theQudaPrecisionTypeMap;
  typedef SingletonHolder<std::map<QudaPrecisionType,QudaPrecision> > theChromaToQudaPrecisionTypeMap;
  // Reader and writer

  //! Read an QudaSolverType enum
  void read(XMLReader& r, const std::string& path, QudaPrecisionType& t);

  //! Write an QudaSolverType enum
  void write(XMLWriter& w, const std::string& path, const QudaPrecisionType& t);


  //! Quda Gauge Reconstruct type
  enum QudaReconsType { 
    RECONS_NONE,
    RECONS_8,
    RECONS_12
  };
  
  namespace QudaReconsTypeEnv { 
    extern const std::string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<QudaReconsType> > theQudaReconsTypeMap;
  typedef SingletonHolder<std::map<QudaReconsType,QudaReconstructType> > theChromaToQudaReconsTypeMap;
  // Reader and writer

  //! Read an QudaReconsType enum
  void read(XMLReader& r, const std::string& path, QudaReconsType& t);

  //! Write an QudaReconsType enum
  void write(XMLWriter& w, const std::string& path, const QudaReconsType& t);


  enum QudaSchwarzMethod { 
	  INVALID_SCHWARZ,
    ADDITIVE_SCHWARZ,
    MULTIPLICATIVE_SCHWARZ,
  };
  
  namespace QudaSchwarzMethodEnv { 
    extern const std::string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<QudaSchwarzMethod> > theQudaSchwarzMethodMap;
  typedef SingletonHolder<std::map<QudaSchwarzMethod,QudaSchwarzType> > theChromaToQudaSchwarzTypeMap;
  // Reader and writer

  //! Read an QudaSchwarzMethod enum
  void read(XMLReader& r, const std::string& path, QudaSchwarzMethod& t);

  //! Write an QudaSchwarzMethod enum
  void write(XMLWriter& w, const std::string& path, const QudaSchwarzMethod& t);



}
#endif
