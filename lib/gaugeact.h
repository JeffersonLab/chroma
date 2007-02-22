// -*- C++ -*-
// $Id: gaugeact.h,v 3.1 2007-02-22 21:11:45 bjoo Exp $

/*! @file
 * @brief Class structure for gauge actions
 */

#ifndef __gaugeactt_h__
#define __gaugeactt_h__


#include "state.h"
#include "gaugebc.h"
#include "create_state.h"

namespace Chroma
{

  //! Abstract base class for gauge actions
  /*! @ingroup actions
   *
   * Supports creation and application for gauge actions
   */
  template<typename P, typename Q>
  class GaugeAction
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~GaugeAction() {}

    //! Return the factory object that produces a state
    /*! 
     * The user will supply the FermState in a derived class 
     *
     * NOTE: this function is public since we have nested gaugeacts
     * that forward this call to the member gaugeact
     */
    virtual const CreateGaugeState<P,Q>& getCreateState() const = 0;

    //! Given links, create the state
    virtual GaugeState<P,Q>* createState(const Q& q) const
    {
      return getCreateState()(q);
    }

    //! Return the gauge BC object for this action
    virtual const GaugeBC<P,Q>& getGaugeBC() const
    {
      return getCreateState().getBC();
    }

    //! Return the set on which the gauge action is defined
    virtual const Set& getSet(void) const = 0;

    //! Compute dS/dU
    /*! Default version. Derived class should override this if needed. */
    virtual void deriv(P& result, const Handle< GaugeState<P,Q> >& state) const 
    {
      QDPIO::cerr << "GaugeAction::deriv not implemented" << endl;
      QDP_abort(1);
    }

    //! Compute the action on a gauge configuration
    virtual Double S(const Handle< GaugeState<P,Q> >& state) const = 0;

  };



  //! Base class for gauge actions with links appearing linearly in the action
  /*! @ingroup actions
   *
   * Here, we are assuming that the "Q" must be a multi1d<LatticeColorMatrix>
   * to be able to use heatbath. Namely, each gauge link must appear in the 
   * Lagrangian density only once.
   *
   * To fully support heatbath, this code needs more work so that 2-link actions
   * can also be supported.
   */
  class LinearGaugeAction : public GaugeAction< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Virtual destructor to help with cleanup;
    virtual ~LinearGaugeAction() {}

    //! Compute staple
    /*! 
     * Default version. Derived class should override this if needed. 
     */
    virtual void staple(LatticeColorMatrix& result,
			const Handle< GaugeState<P,Q> >& state,
			int mu, int cb) const = 0;
  };

}


#endif
