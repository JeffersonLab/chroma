// -*- C++ -*-
// $Id: simple_gaugestate.h,v 3.0 2006-04-03 04:58:54 edwards Exp $

/*! @file
 * @brief Simple gauge state and a creator
 */

#ifndef __simple_gaugestate_h__
#define __simple_gaugestate_h__

#include "state.h"
#include "create_state.h"
#include "handle.h"

namespace Chroma
{


  /*! @ingroup gaugeacts */
  namespace CreateSimpleGaugeStateEnv 
  { 
    extern const std::string name;
    extern const bool registered;
  }


  //! Simple version of GaugeState 
  /*! @ingroup gaugeacts
   *
   * Only needs to hold a gauge field and gauge bc
   */
  template<typename P, typename Q>
  class SimpleGaugeState : public GaugeState<P,Q>
  {
  public:
    //! Full constructor
    SimpleGaugeState(Handle< GaugeBC<P,Q> > gbc_, const Q& q_) : 
      gbc(gbc_), q(q_) 
      {
	// Apply the BC
	gbc->modify(q);
      }

    //! Destructor
    ~SimpleGaugeState() {}

    //! Return the link fields needed in constructing linear operators
    const Q& getLinks() const {return q;}

    //! Return the gauge BC object for this state
    const GaugeBC<P,Q>& getBC() const {return *gbc;}

  private:
    SimpleGaugeState() {}  // hide default constructur
    void operator=(const SimpleGaugeState&) {} // hide =

  private:
    Handle< GaugeBC<P,Q> >  gbc;
    Q q;
  };



  //! Create a simple gauge connection state
  /*! @ingroup gaugeacts
   *
   * This is a factory class for producing a connection state
   */
  template<typename P, typename Q>
  class CreateSimpleGaugeState : public CreateGaugeState<P,Q>
  {
  public:
    //! Full constructor
    CreateSimpleGaugeState(Handle< GaugeBC<P,Q> > gbc_) : gbc(gbc_) {}

    //! Destructor
    ~CreateSimpleGaugeState() {}
   
    //! Construct a ConnectState
    SimpleGaugeState<P,Q>* operator()(const Q& q) const
      {
	return new SimpleGaugeState<P,Q>(gbc, q);
      }

    //! Return the gauge BC object for this state
    const GaugeBC<P,Q>& getBC() const {return *gbc;}

  private:
    CreateSimpleGaugeState() {}  // hide default constructur
    void operator=(const CreateSimpleGaugeState&) {} // hide =

  private:
    Handle< GaugeBC<P,Q> >  gbc;
  };



}


#endif
