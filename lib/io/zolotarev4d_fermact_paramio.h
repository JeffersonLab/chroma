#ifndef ZOLOTAREV4D_FERMACT_PARAMIO_H
#define ZOLOTAREV4D_FERMACT_PARAMIO_H

#include <string>
#include "io/param_io.h"
#include "io/fermact_paramio.h"
#include "io/eigen_io.h"

using namespace std;
using namespace QDP;



struct Zolotarev4DStateInfo {
  Real ApproxMin;
  Real ApproxMax;
  int  NWilsVec;
  EigenIO_t eigen_io;
};

class Zolotarev4DFermActParams : public FermActParams {
 public:
  Zolotarev4DFermActParams(XMLReader& in);
  
  // Satisfy virtual functions
  const enum FermActType getFermActType(void) const { 
    return FERM_ACT_ZOLOTAREV_4D;
  }

  // Public members
  FermActParams* AuxFermActHandle;

  Real Mass;
  int RatPolyDeg;
  Real RsdCGInner;
  int   MaxCGInner;
  int   ReorthFreqInner;
  Zolotarev4DStateInfo StateInfo;
  
  // Destructor
  ~Zolotarev4DFermActParams() {
    if ( AuxFermActHandle == 0x0 ) { 

      // This is OK, as AuxFermActHandle is of type FermActParams
      // which has virtual destructor
      delete AuxFermActHandle;
    }
  }

  // Copy
  Zolotarev4DFermActParams(const Zolotarev4DFermActParams& p) : Mass(p.Mass),
    RatPolyDeg(p.RatPolyDeg), RsdCGInner(p.RsdCGInner), MaxCGInner(p.MaxCGInner), ReorthFreqInner(p.ReorthFreqInner), StateInfo(p.StateInfo), AuxFermActHandle(p.AuxFermActHandle->clone()) {}

  // Virtual constructor 
  Zolotarev4DFermActParams* clone(void) const { 
    return new Zolotarev4DFermActParams( *this );
  }

  Real& getMass(void)  { return Mass;  }
  void setMass(const Real& m) { Mass = m ; }

};


void read(XMLReader& xml_in, const string& path, Zolotarev4DStateInfo& info);
void write(XMLWriter& xml_out, const string& path, const Zolotarev4DStateInfo& info);

void write(XMLWriter& xml_out, const string& path, const Zolotarev4DFermActParams& p);

#endif
