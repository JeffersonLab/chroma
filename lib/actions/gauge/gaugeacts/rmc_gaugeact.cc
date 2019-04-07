/*! \file
 *  \brief Plaquette gauge action as sum of characters
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/rmc_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"
#include "io/aniso_io.h"

namespace Chroma
{
 
  namespace RMCGaugeActEnv 
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
      
      const std::string name = "RMC_GAUGEACT";

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
	read(paramtop, "alpha", alpha);
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
      multi2d<LatticeColorMatrix> plq;
      this->siteAction(plq, state);

      // Total action
      Double act_F = zero;
      
      Double three = Nc;       

      for(int mu=1; mu < Nd; ++mu)
      {
	for(int nu=0; nu < mu; ++nu)
	{
	  LatticeComplex plaq=trace(plq[mu][nu]);
	  // Sum over plaquettes
	  act_F += sum(real(plaq)-three);	
	  
	}
      }

      // Normalize
      Real act = -(param.beta / Real(Nc)) * act_F;

      return act;
    }
 

    //! Compute the plaquette
    void GaugeAct::siteAction(multi2d<LatticeColorMatrix>& plq, const Handle< GaugeState<P,Q> >& state) const
    {
      START_CODE();

      // Initialize
      plq.resize(Nd,Nd);
      plq = zero;

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
	  plq[mu][nu] += (u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu])) ; 
	  
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
     START_CODE();

    const multi1d<LatticeColorMatrix>& u = state->getLinks();
    
    result = zero;
    LatticeColorMatrix tmp1, tmp2;
    LatticeColorMatrix u_nu_mu;

    for(int nu=0; nu < Nd; ++nu) 
    {
      if( nu == mu ) continue;
      
      u_nu_mu = shift(u[nu],FORWARD,mu);

      // +forward staple
      tmp1[rb[cb]] = u_nu_mu * adj(shift(u[mu],FORWARD,nu));
      tmp2[rb[cb]] = tmp1 * adj(u[nu]);

      result[rb[cb]] += param.beta * tmp2;

      // +backward staple
      tmp1[rb[cb]] = adj(shift(u_nu_mu,BACKWARD,nu)) * adj(shift(u[mu],BACKWARD,nu));
      tmp2[rb[cb]] = tmp1 * shift(u[nu],BACKWARD,nu);

      result[rb[cb]] += param.beta * tmp2;
    }

    // NOTE: a heatbath code should be responsible for resetting links on
    // a boundary. The staple is not really the correct place.

    END_CODE();
     
	    
	/*    
      QDPIO::cerr << __func__ << ": staple not possible\n";
      QDP_abort(1);*/
    }


    //! Compute dS/dU
    void GaugeAct::deriv(multi1d<LatticeColorMatrix>& ds_u,
			 const Handle< GaugeState<P,Q> >& state) const
    {
      START_CODE();

      LatticeColorMatrix tmp_1;
      LatticeColorMatrix tmp_2;

      const multi1d<LatticeColorMatrix>& u = state->getLinks();
      multi1d<LatticeColorMatrix> deriv(Nd); 
      
     
      ds_u.resize(Nd);
      for(int mu=0; mu < Nd; mu++) 
      {
	deriv[mu] = zero ;

	
	for(int nu = 0; nu < Nd; nu++) 
	{ 
	  if (mu == nu) continue;

	  LatticeColorMatrix tmp_1 = shift(u[nu], FORWARD, mu);
	  LatticeColorMatrix tmp_2 = shift(u[mu], FORWARD, nu);

	  LatticeColorMatrix up_plq   = u[mu]*tmp_1*adj(tmp_2)*adj(u[nu]);
	  LatticeColorMatrix down_plq = u[mu]*shift(adj(tmp_1)*adj(u[mu])*u[nu],BACKWARD,nu);

	  deriv[mu] += up_plq + down_plq ;	                   
	                   
	
	}// nu
	// Fold in the normalization from the action
	ds_u[mu]  = (-param.beta/Real(2*Nc)        ) * deriv[mu];
		
      }// mu

      // Zero the force on any fixed boundaries
      getGaugeBC().zero(ds_u);

      END_CODE();
    }


     //!  Update coupling coefficient Beta using RMC
    void GaugeAct::updateBeta(multi2d<LatticeReal>& RMCBeta, const Handle< GaugeState<P,Q> >& state) const
    {
	START_CODE();

	multi2d<LatticeColorMatrix> plq;
        this->siteAction(plq, state);

	const multi1d<LatticeColorMatrix>& u = state->getLinks();
	Double M = (param.beta - param.alpha) * 1.5;

	Double three = Nc; 

	for (int mu=0; mu<Nd; mu++) {
		for (int nu=mu+1; nu<Nd; nu++) {
			LatticeReal plaq=three-real(trace(plq[mu][nu]));
			LatticeReal rnd_num;
			random(rnd_num);
			LatticeReal rnd_prob = exp((param.beta-param.alpha)*plaq-M);
			RMCBeta[mu][nu] = where(rnd_num < rnd_prob,param.alpha,param.beta);
		}
	}

	END_CODE();
     }

     //! Restore/initialize coupling coefficient Beta to original values
     void GaugeAct::restoreBeta(multi2d<LatticeReal>& RMCBeta) const
     {
	START_CODE();
	for (int mu=0; mu<Nd; mu++) {
		for(int nu=mu+1;nu<Nd;nu++) {
			RMCBeta[mu][nu] = param.beta;
		}
	}
	END_CODE();

     }

     //! Compute the RMC Action
     Double GaugeAct::RMC_S(multi2d<LatticeReal>& RMCBeta, const Handle< GaugeState<P,Q> >& state) const
     {
	START_CODE();
	
	multi2d<LatticeColorMatrix> plq;
        this->siteAction(plq, state);

	Double s_pg = zero;

	const multi1d<LatticeColorMatrix>& u = state->getLinks();

	Double M = (param.beta - param.alpha) * 1.5;

	Double three = Nc; 

 	for (int mu=0; mu<Nd; mu++) {
		for (int nu=mu+1; nu<Nd; nu++) {
			LatticeReal plaq=three-real(trace(plq[mu][nu]));
			s_pg += sum(RMCBeta[mu][nu]*plaq - (RMCBeta[mu][nu]-param.alpha)/(param.beta-param.alpha)*
				log(Double(1)-exp((param.beta-param.alpha)*plaq-M)));
		}
	}

	s_pg *= Double(-1)/Double(Nc);

	END_CODE();

	return s_pg;
      }

  //! Compute dS_{RMC}/dU
  void GaugeAct::RMC_deriv(multi2d<LatticeReal>& RMCBeta, multi1d<LatticeColorMatrix> & ds_u, const Handle< GaugeState<P,Q> >& state) const
      {
	START_CODE();

	multi2d<LatticeColorMatrix> plq;
        this->siteAction(plq, state);


    	ds_u.resize(Nd);

    	ds_u = zero;

        const multi1d<LatticeColorMatrix>& u = state->getLinks();
      
	Double M = (param.beta - param.alpha) * 1.5;

	Double three = Nc; 

	for (int mu=0; mu<Nd; mu++) {
		for (int nu=mu+1; nu<Nd; nu++) {
			LatticeColorMatrix ds_orig = Real(1.0)-plq[mu][nu];
			
			LatticeReal plaq=three-real(trace(plq[mu][nu]));
			
			LatticeColorMatrix tmp_1 = shift(u[nu], FORWARD, mu);
	  		LatticeColorMatrix tmp_2 = shift(u[mu], FORWARD, nu);

	  		LatticeColorMatrix up_plq   = u[mu]*tmp_1*adj(tmp_2)*adj(u[nu]);
	 		LatticeColorMatrix down_plq = u[mu]*shift(adj(tmp_1)*adj(u[mu])*u[nu],BACKWARD,nu);

	  		LatticeColorMatrix ds_A = (up_plq + down_plq)/Double(2);				
			
			ds_u[mu] += RMCBeta[mu][nu]*ds_orig - (RMCBeta[mu][nu]-param.alpha)/(param.beta-param.alpha)*
				(param.alpha-param.beta)/(exp((param.beta-param.alpha)*plaq-M)-Double(1))*ds_A;
		}
		ds_u[mu] *= Double(-1)/Double(Nc);
	}
	// Zero the force on any fixed boundaries
   	getGaugeBC().zero(ds_u);	
	
	END_CODE();

      }

  }//RMCGaugeActEnv

} // Chroma
