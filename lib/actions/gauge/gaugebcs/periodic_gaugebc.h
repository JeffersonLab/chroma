// $Id: periodic_gaugebc.h,v 3.0 2006-04-03 04:58:54 edwards Exp $
/*! \file
 *  \brief Periodic gauge boundary conditions
 */

#ifndef __periodic_gaugebc_h__
#define __periodic_gaugebc_h__

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
  class PeriodicGaugeBC : public GaugeBC< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

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
    void operator=(const PeriodicGaugeBC&) {}

  };

} // End namespace Chroma

#endif
