// -*- C++ -*-
// $Id: wilson_gaugeact.h,v 1.2 2004-03-29 21:33:03 edwards Exp $
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
  WilsonGaugeAct(Handle< GaugeBC<LatticeFermion> > gbc_, 
		 const Real& beta_) : 
    gbc(gbc_), beta(beta_) {}

  //! Copy constructor
  WilsonGaugeAct(const WilsonGaugeAct& a) : 
    gbc(a.gbc), beta(a.beta) {}

  //! Assignment
  WilsonGaugeAct& operator=(const WilsonGaugeAct& a)
    {gbc=a.gbc; beta=a.beta; return *this;}

  //! Is anisotropy used?
  bool anisoP() const {return false;}

  //! Anisotropy factor
  const Real& anisoFactor() const {return 1.0;}

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
  /*! Default version. Derived class should override this if needed. */
  void dsdu(multi1d<LatticeColorMatrix>& result,
	    Handle<const ConnectState> state) const;

  //! Destructor is automatic
  ~WilsonGaugeAct() {}

private:
  Handle< GaugeBC >  gbc;
  Real beta;
};

#endif
