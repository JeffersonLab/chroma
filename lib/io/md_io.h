#ifndef MD_IO_H
#define MD_IO_H

#include "chromabase.h"
#include "enum_io/enum_io.h"


#include <string>

using namespace QDP;
using namespace std;
using namespace Chroma;

// Base class for MD Integrator parameters
class MDIntegratorParamsBase { 
public: 
  // Virtual destructor
  virtual ~MDIntegratorParamsBase(void) {}

  // Clone function to allow copying through the base class
  virtual MDIntegratorParamsBase* clone(void) const = 0;

  // Function to get at the TypeID
  virtual MDIntegratorType getType(void) const = 0;
};

// The class for parameterising a simple Leapfrog class
class LeapfrogParams : public MDIntegratorParamsBase {
public:
  // No dynamic data -- destructor does nothing
  ~LeapfrogParams(void) {};

  // Read from XML
  LeapfrogParams(XMLReader& top);

  // Copy constructor
  LeapfrogParams(const LeapfrogParams& p) : dt(p.dt), tau(p.tau) {}

  // virtual copy constructor
  LeapfrogParams* clone(void) const {
    return new LeapfrogParams(*this);
  }

  MDIntegratorType getType(void) const {
    return my_type;
  }

  Real getStepSize(void) const { return dt; }
  Real getTrajLength(void) const { return tau; }

private:
  Real dt;
  Real tau;
  MDIntegratorType my_type;
};

void write(XMLWriter& xml, const string& path, const LeapfrogParams& p);

MDIntegratorParamsBase *readMDIntegratorParams(XMLReader& xml, const string& path);

void write(XMLWriter& xml, const string& path, const MDIntegratorParamsBase& p);
#endif
