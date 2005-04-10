// -*- C++ -*-
// $Id: rbc_gaugeact.h,v 1.2 2005-04-10 22:32:38 edwards Exp $
/*! \file
 *  \brief RG style plaquette + rectangle gauge action following RBC conventions
 */

#ifndef __rbc_gaugeact_h__
#define __rbc_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace RBCGaugeActEnv { 
    extern const string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct RBCGaugeActParams {
    // Base Constructor
    RBCGaugeActParams();
    
    // Read params from some root path
    RBCGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real beta;  
    Real c1;  
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, RBCGaugeActParams& param);
  

  //! RG gauge action
  /*! \ingroup gaugeacts
   *
   * A RG (plaquette + rectangle) gauge action following RBC conventions
   *
   * S = -(beta/Nc)*[(1-8*c_1)*Plaq + c_1 * Rect]
   */

  class RBCGaugeAct : public GaugeAction
  {
  public:
    //! General GaugeBC
    RBCGaugeAct(Handle< GaugeBC > gbc_, 
	       const Real& beta_, const Real& c1_) : 
      beta(beta_), c1(c1_) {init(gbc_);}

    //! Read beta from a param struct
    RBCGaugeAct(Handle< GaugeBC > gbc_, 
	       const RBCGaugeActParams& p) :
      beta(p.beta), c1(p.c1) {init(gbc_);}

    //! Copy constructor
    RBCGaugeAct(const RBCGaugeAct& a) : 
      beta(a.beta), c1(a.c1), plaq(a.plaq), rect(a.rect) {}

    //! Assignment
    RBCGaugeAct& operator=(const RBCGaugeAct& a)
    {beta=a.beta; c1=a.c1; plaq=a.plaq; rect=a.rect; return *this;}

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
    const GaugeBC& getGaugeBC() const {return plaq->getGaugeBC();}

    //! Compute staple
    /*! Default version. Derived class should override this if needed. */
    void staple(LatticeColorMatrix& result,
		Handle<const ConnectState> state,
		int mu, int cb) const
    {
      plaq->staple(result,state,mu,cb);

      LatticeColorMatrix tmp;
      rect->staple(tmp,state,mu,cb);
      result += tmp;
    }

    //! Compute dS/dU
    void dsdu(multi1d<LatticeColorMatrix>& result,
	      const Handle<const ConnectState> state) const
    {
      plaq->dsdu(result,state);

      multi1d<LatticeColorMatrix> tmp;
      rect->dsdu(tmp,state);
      result += tmp;
    }

    //! Compute the actions
    Double S(const Handle<const ConnectState> state) const
    {
      return plaq->S(state) + rect->S(state);
    }

    //! Destructor is automatic
    ~RBCGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getBeta(void) const { return beta; }

  protected:
    //! Private initializer
    void init(Handle< GaugeBC > gbc);

  private:
    Real   beta;               // The coupling Beta
    Real   c1;                 // The relative rectangle factor
    Handle<PlaqGaugeAct> plaq; // Hold a plaquette gaugeact
    Handle<RectGaugeAct> rect; // Hold a rectangle gaugeact

  };

};


#endif
