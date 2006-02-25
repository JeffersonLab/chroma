// $Id: gaugebc_periodic.h,v 2.1 2006-02-25 19:47:46 edwards Exp $
/*! \file
 *  \brief Periodic gauge boundary conditions
 */

#ifndef __gaugebc_periodic_h__
#define __gaugebc_periodic_h__

#include "chromabase.h"
#include "gaugebc.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"

namespace Chroma {

  /*! @ingroup gaugebcs */
  namespace PeriodicGaugeBCEnv { 
    extern const std::string name;
    extern const bool registered;
  };

  //! Periodic gauge
  /*! @ingroup gaugebcs */
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
} // End namespace Chroma

#endif
