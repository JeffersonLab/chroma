// -*- C++ -*-
// $Id: plaq_coarse_fine_gaugeact.h,v 3.1 2008-01-16 19:01:13 edwards Exp $
/*! \file
 *  \brief Plaquette gauge action that supports 2x2 fine/coarse style anisotropy
 */

#ifndef __plaq_coarse_fine_gaugeact_h__
#define __plaq_coarse_fine_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace PlaqCoarseFineGaugeActEnv 
  { 
    extern const string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct PlaqCoarseFineGaugeActParams 
  {
    // Base Constructor
    PlaqCoarseFineGaugeActParams() {}
    
    // Read params from some root path
    PlaqCoarseFineGaugeActParams(XMLReader& xml_in, const std::string& path);

    multi1d<bool> coarseP;
    Real coeff_ff;
    Real coeff_cc;
    Real coeff_cf;
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, PlaqCoarseFineGaugeActParams& param);
  

  //! Plaquette gauge action that supports 2x2 fine/coarse style anisotropy
  /*! \ingroup gaugeacts
   *
   * Some anisotropy is used in this class. 
   * If it is not used, then why do you need this class?
   * However, the coefficients could be the same for both coarse
   * and fine directions, so the isotropic limit is supported, but
   * not necessarily optimized for that limit.
   */

  class PlaqCoarseFineGaugeAct : public LinearGaugeAction
  {
  public:
    //! Read coeff from a param struct
    PlaqCoarseFineGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
			   const PlaqCoarseFineGaugeActParams& p) :
      cgs(cgs_), param(p) {init();}

    //! Is anisotropy used?
    bool anisoP() const {return true;}

    //! Coarse directions
    const multi1d<bool>& coarseDirs() const {return param.coarseP;}

    //! Return the set on which the gauge action is defined
    /*! Defined on the even-off (red/black) set */
    const Set& getSet() const {return rb;}

    //! Compute staple
    void staple(LatticeColorMatrix& result,
		const Handle< GaugeState<P,Q> >& state,
		int mu, int cb) const;

    //! Compute dS/dU
    void deriv(multi1d<LatticeColorMatrix>& result,
	       const Handle< GaugeState<P,Q> >& state) const;

    //! Produce a gauge create state object
    const CreateGaugeState<P,Q>& getCreateState() const {return *cgs;}

    //! Compute the actions
    Double S(const Handle< GaugeState<P,Q> >& state) const;

    //! Destructor is automatic
    ~PlaqCoarseFineGaugeAct() {}

    // Accessors -- non mutable members.
    const Real& getCoeffCC() const {return param.coeff_cf;}
    const Real& getCoeffFF() const {return param.coeff_ff;}
    const Real& getCoeffCF() const {return param.coeff_cf;}

  protected:
    PlaqCoarseFineGaugeAct() {}                         /*!< Hide default constructor */
    void operator=(const PlaqCoarseFineGaugeAct& a) {}  /*!< Hide assignment */

    //! Internal initializer
    void init();

  private:
    Handle< CreateGaugeState<P,Q> >  cgs;     /*!< Create Gauge State */
    PlaqCoarseFineGaugeActParams     param;   /*!< The parameters */
    multi2d<Real>                    coeffs;  /*!< Array of coefficients for aniso */
  };

}


#endif
