// -*- C++ -*-
// $Id: rect_gaugeact.h,v 3.0 2006-04-03 04:58:54 edwards Exp $
/*! \file
 *  \brief Rectangle gauge action
 */

#ifndef __rect_gaugeact_h__
#define __rect_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace RectGaugeActEnv { 
    extern const string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct RectGaugeActParams {
    // Base Constructor
    RectGaugeActParams();
    
    // Read params from some root path
    RectGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real coeff;  
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, RectGaugeActParams& param);
  

  //! Rect gauge action
  /*! \ingroup gaugeacts
   *
   * The standard Rect gauge action
   */

  class RectGaugeAct : public LinearGaugeAction
  {
  public:
    //! General CreateGaugeState<P,Q>
    RectGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		 const Real& coeff_) : 
      cgs(cgs_), coeff(coeff_) {}

    //! Read rectangle coefficient from a param struct
    RectGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		 const RectGaugeActParams& p) :
      cgs(cgs_), coeff(p.coeff) {}

    //! Is anisotropy used?
    bool anisoP() const {return false;}

    //! Anisotropy factor
    const Real anisoFactor() const {return Real(1);}

    //! Anisotropic direction
    int tDir() const {return Nd-1;}

    //! Return the set on which the gauge action is defined
    /*! Defined on the even-off (red/black) set */
    const OrderedSet& getSet() const {return rb;}

    //! Compute staple
    /*! Default version. Derived class should override this if needed. */
    void staple(LatticeColorMatrix& result,
		const Handle< GaugeState<P,Q> >& state,
		int mu, int cb) const;

    //! Compute dS/dU
    void deriv(multi1d<LatticeColorMatrix>& result,
	       const Handle< GaugeState<P,Q> >& state) const;

    //! Compute the actions
    Double S(const Handle< GaugeState<P,Q> >& state) const;

    //! Produce a gauge create state object
    const CreateGaugeState<P,Q>& getCreateState() const {return *cgs;}

    //! Destructor is automatic
    ~RectGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getCoeff(void) const { return coeff; }

  protected:
    //! Partial construcor
    RectGaugeAct() {}
    //! Hide assignment
    void operator=(const RectGaugeAct& a) {}

  private:
    Handle< CreateGaugeState<P,Q> >  cgs;  // Create gauge state
    Real coeff;              // The coupling coefficient
  };

};


#endif
