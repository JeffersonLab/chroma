// -*- C++ -*-
// $Id: plaq_gaugeact.h,v 3.7 2008-01-20 19:47:24 edwards Exp $
/*! \file
 *  \brief Plaquette gauge action
 */

#ifndef __plaq_gaugeact_h__
#define __plaq_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "io/aniso_io.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace PlaqGaugeActEnv 
  { 
    extern const string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct PlaqGaugeActParams 
  {
    //! Base Constructor
    PlaqGaugeActParams() {}
    
    //! Read params from some root path
    PlaqGaugeActParams(XMLReader& xml_in, const std::string& path);

    multi2d<Real>  coeffs;       /*!< Array of coefficients for aniso */
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, PlaqGaugeActParams& param);
  
//  /*! @ingroup gaugeacts */
//  void write(XMLWriter& xml, const string& path, const PlaqGaugeActParams& param);
  

  //! Plaquette gauge action
  /*! \ingroup gaugeacts
   *
   * The standard Plaquette gauge action
   */

  class PlaqGaugeAct : public LinearGaugeAction
  {
  public:
    //! General CreateGaugeState<P,Q>
    /*!< Supplied for callers with simple params */
    PlaqGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		 const Real& coeff,
		 const AnisoParam_t& aniso) : cgs(cgs_)
      {init(coeff,aniso);}

    //! General CreateGaugeState<P,Q>
    /*!< Supplied for callers with simple params */
    PlaqGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		 const Real& coeff_s,
		 const Real& coeff_t,
		 const AnisoParam_t& aniso) : cgs(cgs_)
      {init(coeff_s,coeff_t,aniso);}

    //! Read coeff from a param struct
    PlaqGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		 const PlaqGaugeActParams& p) :
      cgs(cgs_), param(p) {}

    //! Return the set on which the gauge action is defined
    /*! Defined on the even-off (red/black) set */
    const Set& getSet() const {return rb;}

    //! Compute staple
    void staple(LatticeColorMatrix& result,
		const Handle< GaugeState<P,Q> >& state,
		int mu, int cb) const;

    //! Compute staple
    void stapleSpatial(LatticeColorMatrix& result,
		       const Handle< GaugeState<P,Q> >& state,
		       int mu, int cb, int t_dir) const;

    //! Compute staple
    void stapleTemporal(LatticeColorMatrix& result,
			const Handle< GaugeState<P,Q> >& state,
			int mu, int cb, int t_dir) const;

    //! Compute dS/dU
    void deriv(multi1d<LatticeColorMatrix>& result,
	       const Handle< GaugeState<P,Q> >& state) const;

    //! compute spatial dS/dU given a time direction
    void derivSpatial(multi1d<LatticeColorMatrix>& result,
		      const Handle< GaugeState<P,Q> >& state,
		      int t_dir) const;

    //! compute spatial dS/dU given a time direction
    void derivTemporal(multi1d<LatticeColorMatrix>& result,
		       const Handle< GaugeState<P,Q> >& state,
		       int t_dir) const;

    //! Produce a gauge create state object
    const CreateGaugeState<P,Q>& getCreateState() const {return *cgs;}

    //! Compute the actions
    Double S(const Handle< GaugeState<P,Q> >& state) const;

    //! Compute the spatial part of the action given a time direction
    Double spatialS(const Handle< GaugeState<P,Q> >& state, int t_dir) const;

    //! Compute the temporal part of the action given a time direction
    Double temporalS(const Handle< GaugeState<P,Q> >& state, int t_dir) const;

    //! Destructor is automatic
    ~PlaqGaugeAct() {}

  protected:
    PlaqGaugeAct() {}
    void operator=(const PlaqGaugeAct& a) {}       //! Hide assignment

    //! Internal initializer for non-general input
    void init(const Real& coeff, 
	      const AnisoParam_t& aniso);

    //! Internal initializer for non-general input
    void init(const Real& coeff_s, 
	      const Real& coeff_t, 
	      const AnisoParam_t& aniso);

  private:
    Handle< CreateGaugeState<P,Q> >  cgs;  /*!< Create Gauge State */
    PlaqGaugeActParams  param;             /*!< The parameters */
  };

};


#endif
