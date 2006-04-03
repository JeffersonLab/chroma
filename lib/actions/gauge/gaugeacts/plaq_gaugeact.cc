// $Id: plaq_gaugeact.cc,v 3.0 2006-04-03 04:58:54 edwards Exp $
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
    GaugeAction< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
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


  //! Compute staple
  /*!
   * \param u_staple   result      ( Write )
   * \param state      gauge field ( Read )
   * \param mu         direction for staple ( Read )
   * \param cb         subset on which to compute ( Read )
   */
  void
  PlaqGaugeAct::staple(LatticeColorMatrix& u_staple,
		       const Handle< GaugeState<P,Q> >& state,
		       int mu, int cb) const
  {
    QDPIO::cout << "PlaqGaugeAct::staple() --- Is this tested ? BJ" << endl;
    QDPIO::cout << "Use at own risk" << endl;

    const multi1d<LatticeColorMatrix>& u = state->getLinks();
				 
    LatticeColorMatrix tmp_1;

    // Get the set
    const OrderedSet& actionSet = getSet();

    // Need to have Even/Odd checkerboarding of 2 subsets
    if (actionSet.numSubsets() != 2)
    {
      QDPIO::cerr << "PlaqGaugeAct::staple  implemented only for even/odd" << endl;
      QDP_abort(1);
    }

    // Aniso^2
    const Real xi02 = anisoFactor() * anisoFactor();  
    const int t_dir = tDir();
  
    // Initialise the staple to zero
    u_staple = zero;

    for(int nu = 0; nu < Nd; ++nu)
    {
      if (nu == mu) continue;

      // Forward staple  
      // tmp_1(x) = u(x+mu,nu)*u_dag(x+nu,mu)  
      tmp_1[actionSet[cb]] = shift(u[nu], FORWARD, mu) * shift(adj(u[mu]), FORWARD, nu);

      if( anisoP() )  	
	if( mu == t_dir || nu == t_dir )
	  tmp_1[actionSet[cb]] *= xi02;

      // u_staple(x) +=  tmp_1 * u_dag(x,nu)
      //   += u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)  
      u_staple[actionSet[cb]] += tmp_1 * adj(u[nu]);
          

      // Backward staple  
      // tmp_1(x) = u(x,mu)*u(x+mu,nu)  
      tmp_1[actionSet[1-cb]] = u[mu] * shift(u[nu], FORWARD, mu);

      if( anisoP() )  	
	if( mu == t_dir || nu == t_dir )
	  tmp_1[actionSet[1-cb]] = xi02 * tmp_1;

      // u_staple(x) += tmp_1_dag(x-nu) * u(x-nu,nu)
      //  += u_dag(x+mu-nu,nu)*u_dag(x-nu,mu)*u(x-nu,nu)  
      u_staple[actionSet[cb]] += shift(adj(tmp_1), BACKWARD, nu) * shift(u[nu], BACKWARD, nu);
    }  // closes nu loop */

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

    for(int mu = 0; mu < Nd; ++mu)
    {
      for(int nu=mu+1; nu<Nd; nu++) 
      {
	for(int cb=0; cb < 2; cb++) 
	{ 
	  tmp_0[rb[cb]] = shift(u[mu], FORWARD, nu)*shift(adj(u[nu]), FORWARD, mu);
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
      ds_u[mu] *= Real(-1)*Real(param.coeff)/(Real(2*Nc));
    }


#if 0
    ds_u.resize(Nd);
    ds_u =zero;

    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    for(int mu=0; mu < Nd; mu++) { 
      LatticeColorMatrix G;
      G = zero;
      
      for(int nu = mu+1; nu < Nd; nu++) { 
	
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

	S_pg += tmp;
      }
    }

    // Normalize
    S_pg *= Double(-1)*Double(param.coeff)/Double(Nc);

    return S_pg;
  } 

}

