/*! \file
 *  \brief Plaquette plus a power of a adjoint gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/plaq_plus_adj_power_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"
#include "io/aniso_io.h"

namespace Chroma
{
 
  namespace PlaqPlusAdjPowerGaugeActEnv 
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
      
      const std::string name = "PLAQ_PLUS_ADJ_POWER_GAUGEACT";

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
	read(paramtop, "beta", beta);
	read(paramtop, "gamma", gamma);
	read(paramtop, "q", q);
      }
      catch( const std::string& e ) { 
	QDPIO::cerr << "Error reading XML: " <<  e << std::endl;
	QDP_abort(1);
      }
    }


    //! Compute the action
    Double GaugeAct::S(const Handle< GaugeState<P,Q> >& state) const
    {
      // Action at the site level
      multi2d<LatticeComplex> plq;
      this->siteAction(plq, state);

      // Total action
      Double act_F = zero;
      Double act_A = zero;
      Real   one   = Real(1);
      int       q2 = param.q >> 1;

      for(int mu=1; mu < Nd; ++mu)
      {
	for(int nu=0; nu < mu; ++nu)
	{
	  // Sum over plaquettes
	  act_F += sum(one - real(plq[mu][nu]));
	  act_A += sum(pow(one - localNorm2(plq[mu][nu]), q2));
	}
      }

      // Normalize
      Double act = param.beta*act_F + param.gamma*act_A;

      return act;
    }
 

    //! Compute the plaquette
    void GaugeAct::siteAction(multi2d<LatticeComplex>& plq, const Handle< GaugeState<P,Q> >& state) const
    {
      START_CODE();

      // Initialize
      plq.resize(Nd,Nd);
      plq = zero;

      // Handle< const GaugeState<P,Q> > u_bc(createState(u));
      // Apply boundaries
      const multi1d<LatticeColorMatrix>& u = state->getLinks();
      Real   third = Real(1) / Real(Nc);

      // Compute the average plaquettes
      for(int mu=1; mu < Nd; ++mu)
      {
	for(int nu=0; nu < mu; ++nu)
	{
	  /* tmp_0 = u(x+mu,nu)*u_dag(x+nu,mu) */
	  /* tmp_1 = tmp_0*u_dag(x,nu)=u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
	  /* wplaq_tmp = tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)) */
	  plq[mu][nu] += third*trace(u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu])); 
	  
	  // Keep a copy
	  plq[nu][mu] = plq[mu][nu];
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

      // Action at the site level
      multi2d<LatticeComplex> plq;
      this->siteAction(plq, state);

      const multi1d<LatticeColorMatrix>& u = state->getLinks();
      multi1d<LatticeColorMatrix> deriv_fun(Nd); 
      multi1d<LatticeColorMatrix> deriv_adj(Nd);
      ds_u.resize(Nd);
      
      Real one   = Real(1);
      int     q2 = param.q >> 1;

      for(int mu=0; mu < Nd; mu++) 
      {
	deriv_fun[mu] = zero;
	deriv_adj[mu] = zero;

	for(int nu = 0; nu < Nd; nu++) 
	{ 
	  if (mu == nu) continue;

	  const LatticeReal& plaq = pow(one - localNorm2(plq[mu][nu]), q2-1);

	  LatticeColorMatrix tmp_1 = shift(u[nu], FORWARD, mu);
	  LatticeColorMatrix tmp_2 = shift(u[mu], FORWARD, nu);
	  LatticeComplex     downp = shift(plaq, BACKWARD,nu);

	  LatticeColorMatrix up_plq   = u[mu]*tmp_1*adj(tmp_2)*adj(u[nu]);
	  LatticeColorMatrix down_plq = u[mu]*shift(adj(tmp_1)*adj(u[mu])*u[nu],BACKWARD,nu);

	  deriv_fun[mu] += up_plq + down_plq;

	  deriv_adj[mu] += plaq * up_plq*conj(trace(up_plq));
	  deriv_adj[mu] += downp * down_plq*conj(trace(down_plq));
	
	}// nu

	// Fold in the normalization from the action
	ds_u[mu]  = (-param.beta/Real(2*Nc)             ) * deriv_fun[mu];
	ds_u[mu] += (-Real(q2) * param.gamma/Real(Nc*Nc)) * deriv_adj[mu];
      }// mu

      // Zero the force on any fixed boundaries
      getGaugeBC().zero(ds_u);

      END_CODE();
    }


  }//PlaqPlusAdjPowerGaugeActEnv

} // Chroma
