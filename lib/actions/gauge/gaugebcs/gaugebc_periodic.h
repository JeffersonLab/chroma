#ifndef GAUGEBC_PERIODIC_H
#define GAUGEBC_PERIODIC_H

#include "chromabase.h"
#include "gaugebc.h"
#include "actions/gauge/gaugebc_factory.h"

using namespace QDP;
using namespace std;
using namespace Chroma;

namespace Chroma {

  namespace PeriodicGaugeBCEnv { 
    extern const std::string name;
    extern const bool registered;
  };

  class PeriodicGaugeBC : public GaugeBC
  {
  public:
    //! Only empty constructor
    PeriodicGaugeBC() {}

    //! Copy constructor
    PeriodicGaugeBC(const PeriodicGaugeBC& a) {}

    //! Destructor is automatic
    ~PeriodicGaugeBC() {}

    //! Assignment
    PeriodicGaugeBC& operator=(const PeriodicGaugeBC&) {return *this;}

    //! Modify U fields in place
    /*! NOP */
    void modify(multi1d<LatticeColorMatrix>& u) const {}

    //! Zero the U fields in place on the masked links
    void zero(multi1d<LatticeColorMatrix>& u) const {}

    //! Says if there are non-trivial BC links
    bool nontrivialP() const {return false;}

  private:

  };
}; // End namespace Chroma

#endif
