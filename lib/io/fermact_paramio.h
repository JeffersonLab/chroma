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

};

// Params for wilson ferm acts
class WilsonFermActParams : public FermActParams {
 private:
  enum FermActType my_fermact_type;
 public:
  // Constructor: Read from XML Reader
  WilsonFermActParams( XMLReader& xml_in );

  // Satisfy virtual functions
  const enum FermActType getFermActType(void) const { return my_fermact_type; }


  // Struct members
  Real Mass;
  AnisoParam_t anisoParam;
};

class DWFFermActParams : public FermActParams { 
 private:
  enum FermActType my_fermact_type;
 public:

  // Constructor: Read from XML Reader
  DWFFermActParams( XMLReader& xml_in);
    
  // Satisfy virtual functions
  const enum FermActType getFermActType(void) const { return my_fermact_type; }

  // Struct Members
  Real Mass;
  AnisoParam_t anisoParam;
  ChiralParam_t chiralParam;
};

// Top level reader. Constructs appropriate from reader and
// returns as base class
FermActParams *read(XMLReader& reader, const string& path);

void write(XMLWriter &xml_out, const string& path, const FermActParams& p);
void write(XMLWriter &xml_out, const string& path, const WilsonFermActParams& p);
void write(XMLWriter &xml_out, const string& path, const DWFFermActParams& p);

#endif
