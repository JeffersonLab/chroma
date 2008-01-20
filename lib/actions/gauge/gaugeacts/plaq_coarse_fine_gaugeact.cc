// $Id: plaq_coarse_fine_gaugeact.cc,v 3.3 2008-01-20 16:11:43 edwards Exp $
/*! \file
 *  \brief Plaquette gauge action that supports 2x2 fine/coarse style anisotropy
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/plaq_coarse_fine_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma
{
 
  namespace PlaqCoarseFineGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, 
		 multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
							       const std::string& path) 
    {
      return new PlaqCoarseFineGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
					PlaqCoarseFineGaugeActParams(xml, path));
    }

    const std::string name = "PLAQ_COARSE_FINE_GAUGEACT";

    //! Local registration flag
    static bool registered = false;

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
  }


  PlaqCoarseFineGaugeActParams::PlaqCoarseFineGaugeActParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);

    try 
    {
      read(paramtop, "coarseP", coarseP);
      read(paramtop, "coeff_cc", coeff_cc);
      read(paramtop, "coeff_ff", coeff_ff);
      read(paramtop, "coeff_cf", coeff_cf);
    }
    catch( const std::string& e ) 
    {
      QDPIO::cerr << __func__ << ": Error reading XML: " <<  e << endl;
      QDP_abort(1);
    }

    if (coarseP.size() != Nd)
    {
      QDPIO::cerr << __func__ << ": invalid size of coarseP" << endl;
      QDP_abort(1);
    }
  }


  void read(XMLReader& xml, const string& path, PlaqCoarseFineGaugeActParams& p) 
  {
    PlaqCoarseFineGaugeActParams tmp(xml, path);
    p=tmp;
  }


  // Internal initializer
  void
  PlaqCoarseFineGaugeAct::init()
  {
    START_CODE();

    coeffs.resize(Nd,Nd);
    coeffs = zero;

    for(int mu = 0; mu < Nd; ++mu)
    {
      for(int nu = mu+1; nu < Nd; ++nu) 
      { 
	if( coarseDirs()[mu] && coarseDirs()[nu] ) 
	{
	  // Coarse-coarse
	  coeffs[mu][nu] = param.coeff_cc;
	}
	else if( (! coarseDirs()[mu]) && (! coarseDirs()[nu]) )
	{
	  // Fine-fine
	  coeffs[mu][nu] = param.coeff_ff;
	}
	else
	{
	  // Coarse-fine or Fine-coarse
	  coeffs[mu][nu] = param.coeff_cf;
	}

	coeffs[nu][mu] = coeffs[mu][nu];
      }
    }
    
    END_CODE();
  }


  //! Compute staple
  /*!
   * \param u_mu_staple   result      ( Write )
   * \param state         gauge field ( Read )
   * \param mu            direction for staple ( Read )
   * \param cb            subset on which to compute ( Read )
   */
  void
  PlaqCoarseFineGaugeAct::staple(LatticeColorMatrix& u_mu_staple,
				 const Handle< GaugeState<P,Q> >& state,
				 int mu, int cb) const
  {
    START_CODE();

    // This bit of code was taken from  chroma/lib/update/heatbath/u_staple.cc
    // Supposedly it works.

    const multi1d<LatticeColorMatrix>& u = state->getLinks();
    
    u_mu_staple = zero;
    LatticeColorMatrix tmp1, tmp2;
    LatticeColorMatrix u_nu_mu;

    for(int nu=0; nu<Nd; nu++) 
    {
      if( nu == mu ) continue;
      
      u_nu_mu = shift(u[nu],FORWARD,mu);

      // +forward staple
      tmp1[rb[cb]] = u_nu_mu * adj(shift(u[mu],FORWARD,nu));
      tmp2[rb[cb]] = tmp1 * adj(u[nu]);

      u_mu_staple[rb[cb]] += coeffs[mu][nu] * tmp1;

      // +backward staple
      tmp1[rb[1-cb]] = adj(shift(u_nu_mu,BACKWARD,nu)) * adj(shift(u[mu],BACKWARD,nu));
      tmp2[rb[cb]]   = tmp1 * shift(u[nu],BACKWARD,nu);

      u_mu_staple[rb[cb]] += coeffs[mu][nu] * tmp2;
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
  PlaqCoarseFineGaugeAct::deriv(multi1d<LatticeColorMatrix>& ds_u,
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
  PlaqCoarseFineGaugeAct::S(const Handle< GaugeState<P,Q> >& state) const
  {
    START_CODE();

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
    
    END_CODE();

    return S_pg;
  } 
  
}

