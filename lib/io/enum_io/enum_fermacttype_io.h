#ifndef enum_fermacttype_io_h
#define enum_fermacttype_io_h

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
  //! Fermion Action type
//! Types of fermion actions
  enum FermActType 
    {
      FERM_ACT_WILSON,
      FERM_ACT_UNPRECONDITIONED_WILSON,
      FERM_ACT_PARITY_BREAKING_WILSON,
      FERM_ACT_CLOVER,
      FERM_ACT_UNPRECONDITIONED_CLOVER,
      FERM_ACT_DWF,                           // Precond. Shamir DWF
      FERM_ACT_UNPRECONDITIONED_DWF,          // Unprec. Shamir DWF
      FERM_ACT_PROJECTED_DWF,                 // Shamir precond. DWF with E&H projection
      FERM_ACT_ZOLOTAREV_4D,                  // Overlap pole with Zolotarev coeffs
      FERM_ACT_ZOLOTAREV_5D,                  // 5D overlap op. ZOlotarev coeffs
      FERM_ACT_OVERLAP_DWF,                   // Borici
      FERM_ACT_EXTENDED_OVERLAP,              // Unprecond. N&N 5D overlap
      FERM_ACT_UNPRECONDITIONED_EXTENDED_OVERLAP,  // Precond. N&N 5D overlap
      FERM_ACT_SMEARED_LAPLACIAN_WILSON,
      FERM_ACT_PLANAR_WILSON,
      FERM_ACT_HAMBER_WU,
      FERM_ACT_STAGGERED,
      FERM_ACT_NAIK,
      FERM_ACT_ASQTAD
    };

  


  namespace FermActTypeEnv { 
    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<FermActType> > theFermActTypeMap;

  // Reader and writer

  //! Read an Fermion Action Type enum
  void read(XMLReader& r, const string& path, FermActType& t);

  //! Write an Fermion Action Type enum
  void write(XMLWriter& w, const string& path, const FermActType& t);

};
#endif
