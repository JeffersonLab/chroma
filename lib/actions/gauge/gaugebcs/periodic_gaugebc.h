// -*- C++ -*-
/*! \file
 *  \brief Periodic gauge boundary conditions
 */

#ifndef __periodic_gaugebc_h__
#define __periodic_gaugebc_h__

#include "gaugebc.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"

namespace Chroma {

  /*! @ingroup gaugebcs */
  namespace PeriodicGaugeBCEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }

  //! Periodic gauge
  /*! @ingroup gaugebcs */
  template<typename P, typename Q>
  class PeriodicGaugeBC : public GaugeBC<P,Q>
  {
  public:

    //! Only empty constructor
    PeriodicGaugeBC() {}

    //! Destructor is automatic
    ~PeriodicGaugeBC() {}

    //! Modify U fields in place
    /*! NOP */
    void modify(Q& u) const {}

    //! Zero the U fields in place on the masked links
    void zero(P& u) const {}

    //! Says if there are non-trivial BC links
    bool nontrivialP() const {return false;}

  private:
    //! Hide assignment
    void operator=(const PeriodicGaugeBC<P,Q>&) {}

  };

} // End namespace Chroma

#endif
