#ifndef FERMACT_PARAMIO_H
#define FERMACT_PARAMIO_H

#include "qdp.h"
#include "io/param_io.h"

using namespace QDP;

#include <string>

using namespace std;

// Base class for fermact params
class FermActParams { 
 public:
  virtual const enum FermActType getFermActType(void) const = 0;

  // Virtual destructor
  virtual ~FermActParams() {}

  // Virtual copy constructor like thingie
  virtual FermActParams* clone(void) const = 0;

  virtual Real& getMass(void) = 0;
  virtual void setMass(const Real& m) = 0;

};


// Top level reader. Constructs appropriate from reader and
// returns as base class
FermActParams *read(XMLReader& reader, const string& path);

void write(XMLWriter &xml_out, const string& path, const FermActParams& p);


// Params for wilson ferm acts
class WilsonFermActParams : public FermActParams {
 private:
  enum FermActType my_fermact_type;
 public:
  // Struct members
  Real Mass;
  AnisoParam_t anisoParam;

  // Constructor: Read from XML Reader
  WilsonFermActParams( XMLReader& xml_in );
  
  // Copy
  WilsonFermActParams( const WilsonFermActParams& p) : my_fermact_type(p.my_fermact_type), Mass(p.Mass), anisoParam(p.anisoParam) {}
   
  // Satisfy virtual functions
  const enum FermActType getFermActType(void) const { return my_fermact_type; }
  Real& getMass(void) { return Mass;  }
  void setMass(const Real& m) { Mass = m ; }

  WilsonFermActParams* clone(void) const {
    // Use my copy constructor 
    return  new WilsonFermActParams(*this);
  }
};

void write(XMLWriter &xml_out, const string& path, const WilsonFermActParams& p);

class DWFFermActParams : public FermActParams { 
 private:
  enum FermActType my_fermact_type;
 public:
  Real Mass;
  AnisoParam_t anisoParam;
  ChiralParam_t chiralParam;

  // Constructor: Read from XML Reader
  DWFFermActParams( XMLReader& xml_in);

  // copy
  DWFFermActParams( const DWFFermActParams& p ) : my_fermact_type(p.my_fermact_type), Mass(p.Mass), anisoParam(p.anisoParam), chiralParam(p.chiralParam) {}

  Real& getMass(void) { return Mass;  }
  void setMass(const Real& m) { Mass = m ; }

  DWFFermActParams* clone(void) const {
    return new DWFFermActParams(*this);
  }

  // Satisfy virtual functions
  const enum FermActType getFermActType(void) const { return my_fermact_type; }


};

void write(XMLWriter &xml_out, const string& path, const DWFFermActParams& p);




#endif
