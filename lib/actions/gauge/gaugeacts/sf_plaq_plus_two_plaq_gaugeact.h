// -*- C++ -*-
/*! \file
 *  \brief Plaquette plus two-plaquette (plaquette squared) gauge action
 */

#ifndef __sf_plaq_plus_two_plaq_gaugeact_h__
#define __sf_plaq_plus_two_plaq_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace SFPlaqPlusTwoPlaqGaugeActEnv 
  { 
    bool registerAll();

    //! Parameter structure
    /*! @ingroup gaugeacts */
    struct Params 
    {
      //! Base Constructor
      Params() {}
    
      //! Read params from some root path
      Params(XMLReader& xml_in, const std::string& path);

      Real  beta_F;     /*! Coupling for fundamental plaquette */
      Real  beta_A;     /*! Coupling for adjoint plaquette */
      int   decay_dir;  /*! Decay direction */
    };
  

    //! PlaqPlusTwoPlaquette gauge action
    /*! \ingroup gaugeacts
     *
     * The standard PlaqPlusTwoPlaquette gauge action
     */
    class GaugeAct : public LinearGaugeAction
    {
    public:
      // Typedefs to save typing
      typedef multi1d<LatticeColorMatrix>  P;
      typedef multi1d<LatticeColorMatrix>  Q;

      //! General CreateGaugeState<P,Q>
      //! Read coeff from a param struct
      GaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, const Params& p);

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

      //! Compute the actions
      Double S(const Handle< GaugeState<P,Q> >& state) const;

      //! Produce a gauge create state object
      const CreateGaugeState<P,Q>& getCreateState() const {return *cgs;}

      //! Destructor is automatic
      ~GaugeAct() {}

    protected:
      //! Hide assignment
      void operator=(const GaugeAct& a) {}
      
      //! Compute the site-level action
      void siteAction(multi2d<LatticeReal>& site_act, const Handle< GaugeState<P,Q> >& state) const;

      //! Compute dS/dU
      void derivPlaqFun(multi1d<LatticeColorMatrix>& ds_u,
			const Handle< GaugeState<P,Q> >& state) const;
      
      //! Compute dS/dU
      void derivPlaqTwo(multi1d<LatticeColorMatrix>& ds_u,
			const Handle< GaugeState<P,Q> >& state) const;
      
    private:
      Handle< CreateGaugeState<P,Q> >  cgs;          /*!< Create Gauge State */
      Params                           param;        /*!< The parameters */
      multi2d<LatticeReal>             plaq_weight;  /*! SF weights for the plaquettes */
    };

  }

}


#endif
