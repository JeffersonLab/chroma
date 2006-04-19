// $Id: plaq_gaugeact.cc,v 3.1 2006-04-19 02:29:45 edwards Exp $
/*! \file
 *  \brief Plaquette gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugeacts/gauge_createstate_factory.h"
#include "actions/gauge/gaugeacts/gauge_createstate_aggregate.h"
#include "meas/glue/mesplq.h"


namespace Chroma
{
 
  namespace PlaqGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, 
		 multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
							       const std::string& path) 
    {
      return new PlaqGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
			      PlaqGaugeActParams(xml, path));
    }

    const std::string name = "PLAQ_GAUGEACT";
    const bool registered = TheGaugeActFactory::Instance().registerObject(name, 
									  createGaugeAct);
  };


  PlaqGaugeActParams::PlaqGaugeActParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);

    try {
      read(paramtop, "./coeff", coeff);

      //  Read optional anisoParam.
      if (paramtop.count("AnisoParam") != 0) 
	read(paramtop, "AnisoParam", aniso);
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " <<  e << endl;
      QDP_abort(1);
    }
  }


  void read(XMLReader& xml, const string& path, PlaqGaugeActParams& p) 
  {
    PlaqGaugeActParams tmp(xml, path);
    p=tmp;
  }


  // Internal initializer
  void
  PlaqGaugeAct::init()
  {
    coeffs.resize(Nd,Nd);
    coeffs = zero;

    for(int mu = 0; mu < Nd; ++mu)
    {
      for(int nu = mu+1; nu < Nd; ++nu) 
      { 
	coeffs[mu][nu] = param.coeff;

	if( anisoP() && (mu == tDir() || nu == tDir()) )
	{
	  coeffs[mu][nu] *= param.aniso.xi_0;
	}
	else
	{
	  coeffs[mu][nu] /= param.aniso.xi_0;
	}

	coeffs[nu][mu] = coeffs[mu][nu];
      }
    }
  }


  //! Compute staple
  /*!
   * \param u_mu_staple   result      ( Write )
   * \param state         gauge field ( Read )
   * \param mu            direction for staple ( Read )
   * \param cb            subset on which to compute ( Read )
   */
  void
  PlaqGaugeAct::staple(LatticeColorMatrix& u_mu_staple,
		       const Handle< GaugeState<P,Q> >& state,
		       int mu, int cb) const
  {
    // This bit of code was taken from  chroma/lib/update/heatbath/u_staple.cc
    // Supposedly it works.

    const multi1d<LatticeColorMatrix>& u = state->getLinks();
				 
    int t_dir   = tDir();
    bool AnisoP = anisoP();
    Double xi02 = anisoFactor() * anisoFactor();

    u_mu_staple = zero;
    for(int nu=0; nu<Nd; nu++) 
    {
      if( nu == mu ) continue;

      Real xi02_tmp;
      if( AnisoP && (mu == t_dir || nu == t_dir) )
      {
	xi02_tmp = xi02;
      }
      else
      {
	xi02_tmp = 1.0;
      }

      // anisotropic lattice, time-like staple
      // +forward staple
      u_mu_staple[rb[cb]]+=
	shift(u[nu],FORWARD,mu)*
	shift(adj(u[mu]),FORWARD,nu)*adj(u[nu])*
	xi02_tmp;
      // +backward staple
      u_mu_staple[rb[cb]]+=
	shift(shift(adj(u[nu]),FORWARD,mu),BACKWARD,nu)*
	shift(adj(u[mu]),BACKWARD,nu)*
	shift(u[nu],BACKWARD,nu)*
	xi02_tmp;
    }

    END_CODE();
  }



  //! Computes the derivative of the fermionic action respect to the link field
  /*!
   *         |  dS
   * ds_u -- | ----   ( Write )
   *         |  dU  
   *
   * \param ds_u       result      ( Write )
   * \param state      gauge field ( Read )
   */
  void
  PlaqGaugeAct::deriv(multi1d<LatticeColorMatrix>& ds_u,
		      const Handle< GaugeState<P,Q> >& state) const
  {
    START_CODE();

    ds_u.resize(Nd);

    LatticeColorMatrix tmp_0;
    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;

    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    ds_u = zero;

    for(int mu = 0; mu < Nd; mu++)
    {
      for(int nu=mu+1; nu < Nd; nu++) 
      {
	for(int cb=0; cb < 2; cb++) 
	{ 
	  tmp_0[rb[cb]] = shift(u[mu], FORWARD, nu)*shift(adj(u[nu]), FORWARD, mu);
	  tmp_0[rb[cb]] *= coeffs[mu][nu];   // c[mu][nu] = c[nu][mu]
	  tmp_1[rb[cb]] = tmp_0*adj(u[mu]);
	  tmp_2[rb[cb]] = u[nu]*tmp_1;
	  ds_u[nu][rb[cb]] += tmp_2;
	  ds_u[mu][rb[cb]] += adj(tmp_2);
	  ds_u[mu][rb[1-cb]] += shift(tmp_1, BACKWARD, nu)*shift(u[nu], BACKWARD, nu);
	  tmp_1[rb[cb]] = adj(u[nu])*u[mu];
	  ds_u[nu][rb[1-cb]] += shift(adj(tmp_0),BACKWARD,mu)*shift(tmp_1, BACKWARD, mu);
	}
      }
      
      // It is 1/(4Nc) to account for normalisation relevant to fermions
      // in the taproj, which is a factor of 2 different from the 
      // one used here.
      ds_u[mu] *= Real(-1)/(Real(2*Nc));
    }


#if 0
    ds_u.resize(Nd);
    ds_u =zero;

    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    for(int mu=0; mu < Nd; mu++) 
    {
      LatticeColorMatrix G;
      G = zero;
      
      for(int nu = mu+1; nu < Nd; nu++) 
      { 
	LatticeColorMatrix up_staple;
	LatticeColorMatrix down_staple;
	
	LatticeColorMatrix tmp_1;
	LatticeColorMatrix tmp_2;
	
	tmp_1 = shift( u[nu], FORWARD, mu);
	tmp_2 = shift( u[mu], FORWARD, nu);

	up_staple = tmp_1*adj(tmp_2)*adj(u[nu]);
	down_staple = adj(tmp_1)*adj(u[mu])*u[nu];

	G += up_staple + shift(down_staple, BACKWARD, nu);
      }

      ds_u[mu] = u[mu]*G;
    }

    // Pure Gauge factor (-coeff/Nc and a factor of 2 because of the forward
    // and backward staple of the force)
    for(int mu=0; mu < Nd; mu++) { 
      ds_u[mu] *= Real(-param.coeff)/(Real(2*Nc));
    }
#endif

    // Zero the force on any fixed boundaries
    getGaugeBC().zero(ds_u);

    END_CODE();
  }

  // Get the gauge action
  //
  // S = -(coeff/(Nc) Sum Re Tr Plaq
  //
  // w_plaq is defined in MesPlq as
  //
  // w_plaq =( 2/(V*Nd*(Nd-1)*Nc)) * Sum Re Tr Plaq
  //
  // so 
  // S = -coeff * (V*Nd*(Nd-1)/2) w_plaq 
  //   = -coeff * (V*Nd*(Nd-1)/2)*(2/(V*Nd*(Nd-1)*Nc))* Sum Re Tr Plaq
  //   = -coeff * (1/(Nc)) * Sum Re Tr Plaq

  Double
  PlaqGaugeAct::S(const Handle< GaugeState<P,Q> >& state) const
  {
    Double S_pg = zero;

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
	Double tmp = 
	  sum(real(trace(u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu]))));

	S_pg += tmp * Double(coeffs[mu][nu]);
      }
    }

    // Normalize
    S_pg *= Double(-1)/Double(Nc);

    return S_pg;
  } 

}

