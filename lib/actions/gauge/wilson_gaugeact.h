// -*- C++ -*-
// $Id: wilson_gaugeact.h,v 1.3 2004-07-23 12:37:12 bjoo Exp $
/*! \file
 *  \brief Wilson gauge action
 */

#ifndef __wilson_gaugeact_h__
#define __wilson_gaugeact_h__

#include "gaugeact.h"

using namespace QDP;

//! Wilson gauge action
/*! \ingroup gaugeact
 *
 * The standard Wilson gauge action
 */

class WilsonGaugeAct : public GaugeAction
{
public:
  //! General GaugeBC
  WilsonGaugeAct(Handle< GaugeBC > gbc_, 
		 const Real& beta_) : 
    gbc(gbc_), beta(beta_) {}

  //! Constructor with different MD beta
  WilsonGaugeAct(Handle< GaugeBC > gbc_,
		 const Real& beta_,
		 const Real& betaMD_) :
    gbc(gbc_), beta(beta_) {}

  //! Copy constructor
  WilsonGaugeAct(const WilsonGaugeAct& a) : 
    gbc(a.gbc), beta(a.beta) {}


  //! Assignment
  WilsonGaugeAct& operator=(const WilsonGaugeAct& a)
    {gbc=a.gbc; beta=a.beta; return *this;}

  //! Clone function for virtual copy
  WilsonGaugeAct* clone(void) const {
    return new WilsonGaugeAct(*this);
  }

  //! Is anisotropy used?
  bool anisoP() const {return false;}

  //! Anisotropy factor
  const Real anisoFactor() const {return Real(1);}

  //! Anisotropic direction
  int tDir() const {return Nd-1;}

  //! Return the set on which the gauge action is defined
  /*! Defined on the even-off (red/black) set */
  const OrderedSet& getSet() const {return rb;}

  //! Produce a gauge boundary condition object
  const GaugeBC& getGaugeBC() const {return *gbc;}

  //! Compute staple
  /*! Default version. Derived class should override this if needed. */
  void staple(LatticeColorMatrix& result,
	      Handle<const ConnectState> state,
	      int mu, int cb) const;

  //! Compute dS/dU
  void dsdu(multi1d<LatticeColorMatrix>& result,
	    const multi1d<LatticeColorMatrix>& u) const;

  //! Compute the actions
  Double S(const multi1d<LatticeColorMatrix>& u) const;

  //! Destructor is automatic
  ~WilsonGaugeAct() {}

  // Accessors -- non mutable members.
  const Real getBeta(void) const { return beta; }

private:
  Handle< GaugeBC >  gbc;  // Gauge Boundary Condition
  Real beta;               // The coupling Beta

};

#endif
