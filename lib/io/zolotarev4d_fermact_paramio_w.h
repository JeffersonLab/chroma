#ifndef ZOLOTAREV4D_FERMACT_PARAMIO_H
#define ZOLOTAREV4D_FERMACT_PARAMIO_H

#include <string>
#include "io/param_io.h"
#include "io/eigen_io.h"
#include "io/overlap_state_info.h"

#include "actions/ferm/fermacts/overlap_fermact_base_w.h"

using namespace std;
using namespace QDP;



class Zolotarev4DFermActParams : public FermActParams {
 public:
  Zolotarev4DFermActParams(XMLReader& in);
  
  // Satisfy virtual functions
  enum FermActType getFermActType(void) const { 
    return FERM_ACT_ZOLOTAREV_4D;
  }

  // Public members
  FermActParams* AuxFermActHandle;

  Real Mass;
  int RatPolyDeg;
  int RatPolyDegPrecond;
  Real RsdCGInner;
  int   MaxCGInner;
  int   ReorthFreqInner;

  ZolotarevStateInfo StateInfo;
  OverlapInnerSolverType InnerSolverType;

  // Destructor
  ~Zolotarev4DFermActParams() {
    if ( AuxFermActHandle == 0x0 ) { 

      // This is OK, as AuxFermActHandle is of type FermActParams
      // which has virtual destructor
      delete AuxFermActHandle;
    }
  }

  // Copy
  Zolotarev4DFermActParams(const Zolotarev4DFermActParams& p) : AuxFermActHandle(p.AuxFermActHandle->clone()), Mass(p.Mass), RatPolyDeg(p.RatPolyDeg), RatPolyDegPrecond(p.RatPolyDegPrecond), RsdCGInner(p.RsdCGInner), MaxCGInner(p.MaxCGInner), ReorthFreqInner(p.ReorthFreqInner), StateInfo(p.StateInfo), InnerSolverType(p.InnerSolverType) {}

  // Virtual constructor 
  Zolotarev4DFermActParams* clone(void) const { 
    return new Zolotarev4DFermActParams( *this );
  }

  Real& getMass(void)  { return Mass;  }
  void setMass(const Real& m) { Mass = m ; }

};


void read(XMLReader& xml_in, const string& path, OverlapInnerSolverType& p);
void write(XMLWriter& xml_out, const string& path, const OverlapInnerSolverType& p);
void write(XMLWriter& xml_out, const string& path, const Zolotarev4DFermActParams& p);

#endif
