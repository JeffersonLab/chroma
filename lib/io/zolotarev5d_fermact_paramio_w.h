#ifndef ZOLOTAREV5D_FERMACT_PARAMIO_H
#define ZOLOTAREV5D_FERMACT_PARAMIO_H

#include <string>
#include "io/param_io.h"
#include "io/fermact_paramio.h"
#include "io/eigen_io.h"
#include "io/overlap_state_info.h"

using namespace std;
using namespace QDP;

class Zolotarev5DFermActParams : public FermActParams {
 public:
  Zolotarev5DFermActParams(XMLReader& in);
  
  // Satisfy virtual functions
  const enum FermActType getFermActType(void) const { 
    return FERM_ACT_ZOLOTAREV_5D;
  }

  // Public members
  FermActParams* AuxFermActHandle;

  Real Mass;
  int RatPolyDeg;
  ZolotarevStateInfo StateInfo;
  
  // Destructor
  ~Zolotarev5DFermActParams() {
    if ( AuxFermActHandle == 0x0 ) { 

      // This is OK, as AuxFermActHandle is of type FermActParams
      // which has virtual destructor
      delete AuxFermActHandle;
    }
  }

  // Copy
  Zolotarev5DFermActParams(const Zolotarev5DFermActParams& p) : AuxFermActHandle(p.AuxFermActHandle->clone()), Mass(p.Mass), RatPolyDeg(p.RatPolyDeg), StateInfo(p.StateInfo)  {}

  // Virtual constructor 
  Zolotarev5DFermActParams* clone(void) const { 
    return new Zolotarev5DFermActParams( *this );
  }

  Real& getMass(void)  { return Mass;  }
  void setMass(const Real& m) { Mass = m ; }

};


void write(XMLWriter& xml_out, const string& path, const Zolotarev5DFermActParams& p);

#endif
