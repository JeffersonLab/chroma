/*! \file
 *  \brief Plaquette plus two-plaquette (plaquette squared) gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/plaq_plus_two_plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"
#include "io/aniso_io.h"

namespace Chroma
{
 
  namespace PlaqPlusTwoPlaqGaugeActEnv 
  { 
    namespace
    {
      GaugeAction< multi1d<LatticeColorMatrix>, 
		   multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
								 const std::string& path) 
      {
	return new GaugeAct(CreateGaugeStateEnv::reader(xml, path), 
			    Params(xml, path));
      }
      
      const std::string name = "PLAQ_PLUS_TWO_PLAQ_GAUGEACT";

      //! Local registration flag
      static bool registered = false;
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheGaugeActFactory::Instance().registerObject(name, createGaugeAct);
	registered = true;
      }
      return success;
    }


    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      XMLReader paramtop(xml_in, path);

      try 
      {	
	read(paramtop, "beta_F", beta_F);
	read(paramtop, "beta_A", beta_A);
      }
      catch( const std::string& e ) { 
	QDPIO::cerr << "Error reading XML: " <<  e << endl;
	QDP_abort(1);
      }
    }


    //! Compute the action
    Double GaugeAct::S(const Handle< GaugeState<P,Q> >& state) const
    {
      // Action at the site level
      multi2d<LatticeReal> plaq_site;
      this->siteAction(plaq_site, state);

      // Total action
      // Fundamental part
      Double act_F = zero;

      for(int mu=1; mu < Nd; ++mu)
      {
	for(int nu=0; nu < mu; ++nu)
	{
	  // Sum over plaquettes
	  act_F += sum(plaq_site[mu][nu]);
	}
      }

      // TwoPlaq part
      Double act_A = zero;

      for(int mu=1; mu < Nd; ++mu)
      {
	for(int nu=0; nu < mu; ++nu)
	{
	  // Sum over plaquettes
	  act_A += sum(plaq_site[mu][nu] * plaq_site[mu][nu]);
	}
      }

      // Normalize
      Real act = -param.beta_F * act_F - Real(0.5)*param.beta_A * act_A;

      return act;
    }
 

    //! Compute the site-level action
    void GaugeAct::siteAction(multi2d<LatticeReal>& site_act, const Handle< GaugeState<P,Q> >& state) const
    {
      START_CODE();

      // Initialize
      site_act.resize(Nd,Nd);
      site_act = zero;

      // Handle< const GaugeState<P,Q> > u_bc(createState(u));
      // Apply boundaries
      const multi1d<LatticeColorMatrix>& u = state->getLinks();

      // Compute the average plaquettes
      for(int mu=1; mu < Nd; ++mu)
      {
	for(int nu=0; nu < mu; ++nu)
	{
	  /* tmp_0 = u(x+mu,nu)*u_dag(x+nu,mu) */
	  /* tmp_1 = tmp_0*u_dag(x,nu)=u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
	  /* wplaq_tmp = tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)) */
	  site_act[mu][nu] += real(trace(u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu])));
	  
	  // Account for normalization
	  site_act[mu][nu] *= Real(1) / Real(Nc);

	  // Keep a copy
	  site_act[nu][mu] = site_act[mu][nu];
	}
      }


      END_CODE();
    }
 

    //! Compute staple
    /*! Default version. Derived class should override this if needed. */
    void GaugeAct::staple(LatticeColorMatrix& result,
			  const Handle< GaugeState<P,Q> >& state,
			  int mu, int cb) const
    {
      QDPIO::cerr << __func__ << ": staple not possible\n";
      QDP_abort(1);
    }


    //! Compute dS/dU
    void GaugeAct::deriv(multi1d<LatticeColorMatrix>& ds_u,
			 const Handle< GaugeState<P,Q> >& state) const
    {
      START_CODE();

      // Derivate at the site level
      multi1d<LatticeColorMatrix> deriv_fun;
      this->derivPlaqFun(deriv_fun, state);

      // Derivate at the site level
      multi1d<LatticeColorMatrix> deriv_two;
      this->derivPlaqTwo(deriv_two, state);

      // Total derivative
      // Fold in normalization. The (1/2) in the two_plaq is removed in the derivative via the chain rule
      ds_u.resize(Nd);

      for(int mu = 0; mu < Nd; mu++)
      {
	ds_u[mu]  = (-param.beta_F) * deriv_fun[mu];
	ds_u[mu] += (-param.beta_A) * deriv_two[mu];
      }

      END_CODE();
    }


    //! Compute dS/dU
    void GaugeAct::derivPlaqFun(multi1d<LatticeColorMatrix>& ds_u,
				const Handle< GaugeState<P,Q> >& state) const
    {
      START_CODE();

      LatticeColorMatrix tmp_0;
      LatticeColorMatrix tmp_1;
      LatticeColorMatrix tmp_2;

      const multi1d<LatticeColorMatrix>& u = state->getLinks();

      ds_u.resize(Nd);
      ds_u =zero;

      for(int mu=0; mu < Nd; mu++) 
      {
	LatticeColorMatrix G;
	G = zero;
      
	for(int nu = 0; nu < Nd; nu++) 
	{ 
	  if (mu == nu) continue;

	  LatticeColorMatrix tmp_1 = shift(u[nu], FORWARD, mu);
	  LatticeColorMatrix tmp_2 = shift(u[mu], FORWARD, nu);

	  LatticeColorMatrix up_staple   = tmp_1*adj(tmp_2)*adj(u[nu]);
	  LatticeColorMatrix down_staple = adj(tmp_1)*adj(u[mu])*u[nu];

	  G += up_staple + shift(down_staple, BACKWARD, nu);
	}
	
	ds_u[mu] = u[mu]*G;
      }

      // It is 1/(4Nc) to account for normalisation relevant to fermions
      // in the taproj, which is a factor of 2 different from the 
      // one used here.
      for(int mu=0; mu < Nd; mu++)
      {
	ds_u[mu] *= Real(1)/(Real(2*Nc));
      }

      // Zero the force on any fixed boundaries
      getGaugeBC().zero(ds_u);

      END_CODE();
    }



    //! Compute dS/dU
    void GaugeAct::derivPlaqTwo(multi1d<LatticeColorMatrix>& ds_u,
				const Handle< GaugeState<P,Q> >& state) const
    {
      START_CODE();

      // Action at the site level
      multi2d<LatticeReal> plaq_site;
      this->siteAction(plaq_site, state);

      // Derivative part
      LatticeColorMatrix tmp_0;
      LatticeColorMatrix tmp_1;
      LatticeColorMatrix tmp_2;

      const multi1d<LatticeColorMatrix>& u = state->getLinks();

      ds_u.resize(Nd);
      ds_u =zero;

      for(int mu=0; mu < Nd; mu++) 
      {
	LatticeColorMatrix G;
	G = zero;
      
	for(int nu = 0; nu < Nd; nu++) 
	{ 
	  if (mu == nu) continue;

	  LatticeColorMatrix tmp_1 = shift(u[nu], FORWARD, mu);
	  LatticeColorMatrix tmp_2 = shift(u[mu], FORWARD, nu);

	  LatticeColorMatrix up_staple   = tmp_1*adj(tmp_2)*adj(u[nu]);
	  LatticeColorMatrix down_staple = adj(tmp_1)*adj(u[mu])*u[nu];

	  G += up_staple * plaq_site[mu][nu];
	  G += shift(down_staple * plaq_site[mu][nu], BACKWARD, nu);
	}
	
	ds_u[mu] = u[mu]*G;
      }

      // It is 1/(4Nc) to account for normalisation relevant to fermions
      // in the taproj, which is a factor of 2 different from the 
      // one used here.
      for(int mu=0; mu < Nd; mu++)
      {
	ds_u[mu] *= Real(1)/(Real(2*Nc));
      }

      // Zero the force on any fixed boundaries
      getGaugeBC().zero(ds_u);

      END_CODE();
    }

  }

}
