// -*- C++ -*-
// $Id: rect_gaugeact.h,v 1.2 2005-01-14 20:13:06 edwards Exp $
/*! \file
 *  \brief Rectangle gauge action
 */

#ifndef __rect_gaugeact_h__
#define __rect_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"

namespace Chroma
{

  namespace RectGaugeActEnv { 
    extern const string name;
    extern const bool registered;
  }

  // Parameter structure
  struct RectGaugeActParams {
    // Base Constructor
    RectGaugeActParams();
    
    // Read params from some root path
    RectGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real coeff;  
  };
  
  void read(XMLReader& xml, const string& path, RectGaugeActParams& param);
  

  //! Rect gauge action
  /*! \ingroup gaugeact
   *
   * The standard Rect gauge action
   */

  class RectGaugeAct : public GaugeAction
  {
  public:
    //! General GaugeBC
    RectGaugeAct(Handle< GaugeBC > gbc_, 
		 const Real& coeff_) : 
      gbc(gbc_), coeff(coeff_) {}

    //! Read rectangle coefficient from a param struct
    RectGaugeAct(Handle< GaugeBC > gbc_, 
		 const RectGaugeActParams& p) :
      gbc(gbc_), coeff(p.coeff) {}

    //! Copy constructor
    RectGaugeAct(const RectGaugeAct& a) : 
      gbc(a.gbc), coeff(a.coeff) {}

    //! Assignment
    RectGaugeAct& operator=(const RectGaugeAct& a)
    {gbc=a.gbc; coeff=a.coeff; return *this;}

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
	      const Handle<const ConnectState> state) const;

    //! Compute the actions
    Double S(const Handle<const ConnectState> state) const;

    //! Destructor is automatic
    ~RectGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getCoeff(void) const { return coeff; }

  private:
    Handle< GaugeBC >  gbc;  // Gauge Boundary Condition
    Real coeff;              // The coupling coefficient
  };

};


#endif
