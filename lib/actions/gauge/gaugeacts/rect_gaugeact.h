// -*- C++ -*-
// $Id: rect_gaugeact.h,v 3.6 2007-05-21 22:21:54 bjoo Exp $
/*! \file
 *  \brief Rectangle gauge action
 */

#ifndef __rect_gaugeact_h__
#define __rect_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "io/aniso_io.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace RectGaugeActEnv 
  { 
    extern const string name;
    bool registerAll();
    
    extern double getTime();
  };

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct RectGaugeActParams {
    // Base Constructor
    RectGaugeActParams() {
      no_temporal_2link=false;
    }
    
    // Read params from some root path
    RectGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real coeff_s;
    Real coeff_t1;
    Real coeff_t2;
    bool no_temporal_2link;
    AnisoParam_t aniso;
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
    //! Backward compatibility:
    RectGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		 const Real& coeff_); 


    //! General CreateGaugeState<P,Q>
    RectGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		 const Real& coeff_s_, 
		 const Real& coeff_t1_,
		 const Real& coeff_t2_,
		 const bool no_temporal_2link_,
		 const AnisoParam_t& aniso_);  

    //! Read rectangle coefficient from a param struct
    RectGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		 const RectGaugeActParams& p);

    //! Return the set on which the gauge action is defined
    /*! Defined on the even-off (red/black) set */
    const Set& getSet() const {return rb;}

    //! Compute staple
    /*! Default version. Derived class should override this if needed. */
    void staple(LatticeColorMatrix& result,
		const Handle< GaugeState<P,Q> >& state,
		int mu, int cb) const;

    //! Compute dS/dU
    void deriv(multi1d<LatticeColorMatrix>& result,
	       const Handle< GaugeState<P,Q> >& state) const;

    //! compute spatial dS/dU given a time direction
    void derivSpatial(multi1d<LatticeColorMatrix>& result,
		      const Handle< GaugeState<P,Q> >& state) const;


		      

    //! compute spatial dS/dU given a time direction
    // No t_dir param is needed as it is hidden in the params
    void derivTemporal(multi1d<LatticeColorMatrix>& result,
		       const Handle< GaugeState<P,Q> >& state) const;




    //! Compute the actions
    Double S(const Handle< GaugeState<P,Q> >& state) const;

    //! Compute the spatial part of the action given a time direction
    Double spatialS(const Handle< GaugeState<P,Q> >& state) const;

    //! Compute the temporal part of the action given a time direction
    Double temporalS(const Handle< GaugeState<P,Q> >& state) const;

    //! Produce a gauge create state object
    const CreateGaugeState<P,Q>& getCreateState() const {return *cgs;}

    //! Destructor is automatic
    ~RectGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getCoeffS(void) const { return params.coeff_s; }
    const Real getCoeffT1(void) const { return params.coeff_t1; }
    const Real getCoeffT2(void) const { return params.coeff_t2; }

    //! Is anisotropy used?
    bool anisoP() const {return params.aniso.anisoP ;}

    //! Anisotropy factor
    const Real anisoFactor() const { return params.aniso.xi_0; }

    //! Anisotropic direction
    int tDir() const {return params.aniso.t_dir;}
    const bool noTemporal21LoopsP(void) const {return params.no_temporal_2link;}

  protected:
    //! Partial construcor
    RectGaugeAct() {}
    //! Hide assignment
    void operator=(const RectGaugeAct& a) {}

  private:
    Handle< CreateGaugeState<P,Q> >  cgs;  // Create gauge state
    RectGaugeActParams params; // THe parameter struct

    // A function for computing the  contribution from one rectangle 
    // in the mu/nu plane specified.
    void deriv_part(int mu, int nu, Real c_munu,
		    multi1d<LatticeColorMatrix>& ds_u, 
		    const multi1d<LatticeColorMatrix>& u) const;

    // A function for computing the contribution to the action from 
    // one rectangle int he mu,nu plane specified
    void S_part(int mu, int nu, Real c, LatticeReal& lgimp, 
		const multi1d<LatticeColorMatrix>& u) const;

  };

};


#endif
