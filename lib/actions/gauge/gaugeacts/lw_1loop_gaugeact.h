// -*- C++ -*-
// $Id: lw_1loop_gaugeact.h,v 1.2 2005-01-14 20:13:06 edwards Exp $
/*! \file
 *  \brief 1-loop tadpole-improved Luscher-Weisz gauge action
 */

#ifndef __lw_1loop_gaugeact_h__
#define __lw_1loop_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "actions/gauge/gaugeacts/pg_gaugeact.h"

namespace Chroma
{

  namespace LW1LoopGaugeActEnv { 
    extern const string name;
    extern const bool registered;
  }

  // Parameter structure
  struct LW1LoopGaugeActParams {
    // Base Constructor
    LW1LoopGaugeActParams();
    
    // Read params from some root path
    LW1LoopGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real beta;  
    Real u0;  
  };
  
  void read(XMLReader& xml, const string& path, LW1LoopGaugeActParams& param);
  

  //! LW1Loop gauge action
  /*! \ingroup gaugeact
   *
   * The standard LW1Loop gauge action
   */

  class LW1LoopGaugeAct : public GaugeAction
  {
  public:
    //! General GaugeBC
    LW1LoopGaugeAct(Handle< GaugeBC > gbc_, 
		   const Real& beta_, const Real& u0_) : 
      beta(beta_), u0(u0_) {init(gbc_);}

    //! Read beta from a param struct
    LW1LoopGaugeAct(Handle< GaugeBC > gbc_, 
		   const LW1LoopGaugeActParams& p) :
      beta(p.beta), u0(p.u0) {init(gbc_);}

    //! Copy constructor
    LW1LoopGaugeAct(const LW1LoopGaugeAct& a) : 
      beta(a.beta), u0(a.u0), plaq(a.plaq), rect(a.rect), pg(a.pg) {}

    //! Assignment
    LW1LoopGaugeAct& operator=(const LW1LoopGaugeAct& a)
    {beta=a.beta; u0=a.u0; plaq=a.plaq; rect=a.rect; pg=a.pg; return *this;}

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

      pg->staple(tmp,state,mu,cb);
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

      pg->dsdu(tmp,state);
      result += tmp;
    }

    //! Compute the actions
    Double S(const Handle<const ConnectState> state) const
    {
      return plaq->S(state) + rect->S(state) + pg->S(state);
    }

    //! Destructor is automatic
    ~LW1LoopGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getBeta(void) const { return beta; }

  protected:
    //! Private initializer
    void init(Handle< GaugeBC > gbc);

  private:
    Real   beta;               // The coupling Beta
    Real   u0;                 // Tadpole factor
    Handle<PlaqGaugeAct> plaq; // Hold a plaquette gaugeact
    Handle<RectGaugeAct> rect; // Hold a rectangle gaugeact
    Handle<PgGaugeAct>   pg;   // Hold a parallelogram gaugeact

  };

};


#endif
