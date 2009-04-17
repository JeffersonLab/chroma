// -*- C++ -*-
// $Id: periodic_gaugestate.h,v 1.3 2009-04-17 02:05:36 bjoo Exp $

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
  template<typename P, typename Q>
  class PeriodicGaugeState : public GaugeState<P,Q>
  {
  public:
    //! Full constructor
    PeriodicGaugeState(const Q& q_) : 
      gbc(new PeriodicGaugeBC<P,Q>()), q(q_) {}

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
    void operator=(const PeriodicGaugeState<P,Q>&) {} // hide =

  private:
    Handle< GaugeBC<P,Q> >  gbc;
    Q q;
  };



  //! Create a periodic gauge connection state
  /*! @ingroup gaugestates
   *
   * This is a factory class for producing a connection state
   */
  template<typename P, typename Q>
  class CreatePeriodicGaugeState : 
    public CreateGaugeState< P, Q>
  {
  public:
    //! Full constructor
    CreatePeriodicGaugeState() : gbc(new PeriodicGaugeBC<P,Q>()) {}

    //! Destructor
    ~CreatePeriodicGaugeState() {}
   
    //! Construct a ConnectState
    PeriodicGaugeState<P,Q>* operator()(const Q& q) const
      {
	return new PeriodicGaugeState<P,Q>(q);
      }

    //! Return the gauge BC object for this state
    const GaugeBC<P,Q>& getBC() const {return *gbc;}

    //! Return the gauge BC object for this state
    Handle< GaugeBC<P,Q> > getGaugeBC() const {return gbc;}

  private:
    void operator=(const CreatePeriodicGaugeState<P,Q>&) {} // hide =

  private:
    Handle< GaugeBC<P,Q> >  gbc;
  };

}


#endif
