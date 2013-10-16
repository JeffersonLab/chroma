/*! \file
 *  \brief Plaquette gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/plaq_plus_adjoint_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"
#include "io/aniso_io.h"

namespace Chroma
{
 
  namespace PlaqPlusAdjointGaugeActEnv 
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
      
      const std::string name = "PLAQ_PLUS_ADJOINT_GAUGEACT";

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
      multi2d<LatticeComplex> plaq_site;
      this->siteAction(plaq_site, state);

      // Total action
      // Fundamental part
      Double act_F = zero;

      for(int mu=1; mu < Nd; ++mu)
      {
	for(int nu=0; nu < mu; ++nu)
	{
	  // Sum over plaquettes
	  act_F += sum(real(plaq_site[mu][nu]));
	}
      }

      // Adjoint part
      Double act_A = zero;
      Double one = 1.0;

      for(int mu=1; mu < Nd; ++mu)
      {
	for(int nu=0; nu < mu; ++nu)
	{
	  // Sum over plaquettes
	  // NOTE: do the subtraction per site to mitigate loss of precision in subtracting big numbers
	  act_A += sum(localNorm2(plaq_site[mu][nu]) - one);
	}
      }

      // Normalize
      Real act = -(param.beta_F / Real(Nc)) * act_F - (param.beta_A / Real(Nc*Nc-1)) * act_A;

      return act;
    }
 

    //! Compute the site-level action
    void GaugeAct::siteAction(multi2d<LatticeComplex>& site_act, const Handle< GaugeState<P,Q> >& state) const
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
	  site_act[mu][nu] += trace(u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu]));
	  
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
      multi1d<LatticeColorMatrix> deriv_adj;
      this->derivPlaqAdj(deriv_adj, state);

      // Total derivative
      // Fold in normalization. The (1/2) in the adjoint is removed in the derivative via the chain rule
      ds_u.resize(Nd);

      for(int mu = 0; mu < Nd; mu++)
      {
	ds_u[mu]  = (-param.beta_F) * deriv_fun[mu];
	ds_u[mu] += (-param.beta_A) * deriv_adj[mu];
      }

      END_CODE();
    }


    //! Compute dS/dU
    void GaugeAct::derivPlaqFun(multi1d<LatticeColorMatrix>& ds_u,
				const Handle< GaugeState<P,Q> >& state) const
    {
      START_CODE();

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
    void GaugeAct::derivPlaqAdj(multi1d<LatticeColorMatrix>& ds_u,
				const Handle< GaugeState<P,Q> >& state) const
    {
      START_CODE();

      // Action at the site level
      multi2d<LatticeComplex> plaq_site;
      this->siteAction(plaq_site, state);

      // Derivative part
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

	  G += up_staple * adj(plaq_site[mu][nu]);
	  G += shift(down_staple * plaq_site[mu][nu], BACKWARD, nu);
	}

	ds_u[mu] = u[mu]*G;
      }

      // Normalize with 1/(Nc^2-1) . Would normally have an additional 1/2.
      // However, only took deriv. of the tr(U_p) and not the adj(tr(U_p)),
      // so need to multiply by 2.
      for(int mu=0; mu < Nd; mu++)
      {
	ds_u[mu] *= Real(1)/ Real(1*(Nc*Nc-1));
      }

      // Zero the force on any fixed boundaries
      getGaugeBC().zero(ds_u);

      END_CODE();
    }

  }

}
