// -*- C++ -*-
// $Id: periodic_gaugestate.h,v 1.2 2006-09-20 20:28:01 edwards Exp $

/*! @file
 * @brief Periodic gauge state and a creator
 */

#ifndef __periodic_gaugestate_h__
#define __periodic_gaugestate_h__

#include "state.h"
#include "create_state.h"
#include "handle.h"
#include "actions/gauge/gaugebcs/periodic_gaugebc.h"

namespace Chroma
{

  /*! @ingroup gaugestates */
  namespace CreatePeriodicGaugeStateEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }


  //! Periodic version of GaugeState 
  /*! @ingroup gaugestates
   *
   * Only needs to hold a gauge field and gauge bc
   */
  class PeriodicGaugeState : public GaugeState< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Full constructor
    PeriodicGaugeState(const Q& q_) : 
      gbc(new PeriodicGaugeBC()), q(q_) {}

    //! Destructor
    ~PeriodicGaugeState() {}

    //! Return the link fields needed in constructing linear operators
    const Q& getLinks() const {return q;}

    //! Return the gauge BC object for this state
    const GaugeBC<P,Q>& getBC() const {return *gbc;}

    //! Return the gauge BC object for this state
    Handle< GaugeBC<P,Q> > getGaugeBC() const {return gbc;}

  private:
    PeriodicGaugeState() {}  // hide default constructur
    void operator=(const PeriodicGaugeState&) {} // hide =

  private:
    Handle< GaugeBC<P,Q> >  gbc;
    Q q;
  };



  //! Create a periodic gauge connection state
  /*! @ingroup gaugestates
   *
   * This is a factory class for producing a connection state
   */
  class CreatePeriodicGaugeState : 
    public CreateGaugeState< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Full constructor
    CreatePeriodicGaugeState() : gbc(new PeriodicGaugeBC()) {}

    //! Destructor
    ~CreatePeriodicGaugeState() {}
   
    //! Construct a ConnectState
    PeriodicGaugeState* operator()(const Q& q) const
      {
	return new PeriodicGaugeState(q);
      }

    //! Return the gauge BC object for this state
    const GaugeBC<P,Q>& getBC() const {return *gbc;}

    //! Return the gauge BC object for this state
    Handle< GaugeBC<P,Q> > getGaugeBC() const {return gbc;}

  private:
    void operator=(const CreatePeriodicGaugeState&) {} // hide =

  private:
    Handle< GaugeBC<P,Q> >  gbc;
  };

}


#endif
