#ifndef GAUGE_BC_IO_H
#define GAUGE_BC_IO_H

#include "chromabase.h"

#include <string>

using namespace QDP;
using namespace std;

// Supported Gauge BC types
enum GaugeBCType_t { 
  GAUGEBC_ALL_PERIODIC, 
  GAUGEBC_SCHROEDINGER_1LINK,
  GAUGEBC_SCHROEDINGER_2LINK,
  GAUGEBC_SIMPLE 
};

// Read and write the GaugeBCType
void read(XMLReader& xml, const string& path, GaugeBCType_t& t);
void write(XMLWriter& xml, const string& path, const GaugeBCType_t& t);

// Read and write teh SchrFunType_t
//! Schroedinger Functional type Boundary Conditions
enum SchrFunType_t {
  SF_NONE = 0, 
  SF_TRIVIAL = 1,
  SF_NONPERT = 2,
  SF_COUPLING = 3,
  SF_CHROMOMAG = 4,
  SF_DIRICHLET = 10,
};

void read(XMLReader& xml, const string& path, SchrFunType_t& t);
void write(XMLWriter& xml, const string& path, const SchrFunType_t& t);


// Base class
class GaugeBCParamsBase {
 public:
  // virtual destructor
  virtual ~GaugeBCParamsBase() {}

  // virtual copy function
  virtual GaugeBCParamsBase* clone(void) const = 0;

  // get type
  virtual const GaugeBCType_t getType(void) const = 0;

};

// read and write the base class
GaugeBCParamsBase* readGaugeBCParams(XMLReader& xml, const string& path);

void write(XMLWriter& xml, const string& path, const GaugeBCParamsBase& p);

class GaugeBCPeriodicParams : public GaugeBCParamsBase {
 public:
  GaugeBCPeriodicParams(XMLReader& xml);

  // Copy constructor
  GaugeBCPeriodicParams(const GaugeBCPeriodicParams& p) {}

  // clone 
  GaugeBCPeriodicParams* clone(void) const {
    return new GaugeBCPeriodicParams(*this);
  }
  // get type
  const GaugeBCType_t getType(void) const { 
    return GAUGEBC_ALL_PERIODIC;
  }
};

void write(XMLWriter& xml, const string& path, const GaugeBCPeriodicParams& p);

// Simple Gauge BC -- set boundary to complex constant 
class GaugeBCSimpleParams : public GaugeBCParamsBase {
 public: 
  // Initialise from XML
  GaugeBCSimpleParams(XMLReader& xml_in);

  // Copy 
  GaugeBCSimpleParams(const GaugeBCSimpleParams& p) : boundary(p.boundary) {}

  // Clone
  GaugeBCSimpleParams* clone(void) const { 
    return new GaugeBCSimpleParams(*this);
  }

  // get type
  const GaugeBCType_t getType(void) const {
    return GAUGEBC_SIMPLE;
  }

  // Get the boundary
  const multi1d<Complex>& getBoundary(void) const {
    return boundary;
  }

 private:
  multi1d<Complex> boundary;
};

void write(XMLWriter& xml, const string& path, const GaugeBCSimpleParams& p);

class GaugeBCSchrParams : public GaugeBCParamsBase {
public:

  GaugeBCSchrParams(XMLReader& xml);
  GaugeBCSchrParams(const GaugeBCSchrParams& p) : SchrFun(p.SchrFun), 
						  SchrPhiMult(p.SchrPhiMult),
						  type_t(p.type_t) {};

  GaugeBCSchrParams* clone(void) const { 
    return new GaugeBCSchrParams(*this);
  }
  
  const GaugeBCType_t getType(void) const { 
    return type_t;
  }

  const SchrFunType_t getSchrFun(void) const { return SchrFun; }
  const Real getSchrPhiMult(void) const { return SchrPhiMult; }
private:
  SchrFunType_t SchrFun;
  Real SchrPhiMult;
  GaugeBCType_t type_t;

};

void write(XMLWriter& xml, const string& path, const GaugeBCSchrParams& p);

#endif
