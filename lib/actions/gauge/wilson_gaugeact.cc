// $Id: wilson_gaugeact.cc,v 1.3 2004-07-19 16:02:59 bjoo Exp $
/*! \file
 *  \brief Wilson gauge action
 */

#include "chromabase.h"
#include "actions/gauge/wilson_gaugeact.h"
#include "meas/glue/mesplq.h"

//! Compute staple
/*!
 * \param u_staple   result      ( Write )
 * \param state      gauge field ( Read )
 * \param mu         direction for staple ( Read )
 * \param cb         subset on which to compute ( Read )
 */

void
WilsonGaugeAct::staple(LatticeColorMatrix& u_staple,
		       Handle<const ConnectState> state,
		       int mu, int cb) const
{
  START_CODE("WilsonGaugeAct::staple");

  const multi1d<LatticeColorMatrix>& u = state->getLinks();
				 
  LatticeColorMatrix tmp_1;

  // Get the set
  const OrderedSet& actionSet = getSet();

  // Need to have Even/Odd checkerboarding of 2 subsets
  if (actionSet.numSubsets() != 2)
  {
    QDPIO::cerr << "WilsonGaugeAct::staple  implemented only for even/odd" << endl;
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

  END_CODE("WilsonGaugeAct::staple");
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
WilsonGaugeAct::dsdu(multi1d<LatticeColorMatrix>& ds_u,
		     Handle<const ConnectState> state) const
{
  START_CODE("WilsonGaugeAct::dsdu");

  const multi1d<LatticeColorMatrix>& u = state->getLinks();
				 
  LatticeColorMatrix tmp_0;
  LatticeColorMatrix tmp_1;
  LatticeColorMatrix tmp_2;

  for(int mu = 0; mu < Nd; mu++) {
    ds_u[mu] = zero;
  }

  for(int mu = 0; mu < Nd; ++mu)
  {
    for(int nu=mu+1; nu<Nd; nu++)
    {
      tmp_0 = shift(u[mu], FORWARD, nu) * adj(shift(u[nu], FORWARD, mu));
      tmp_1 = tmp_0 * adj(u[mu]);
      tmp_2 = u[nu] * tmp_1;
      
      ds_u[nu] += tmp_2;
      ds_u[mu] += adj(tmp_2);

      ds_u[mu] += shift(tmp_1, BACKWARD, nu)*shift(u[nu], BACKWARD, nu);

      tmp_1 = adj(u[nu])*u[mu];

      ds_u[nu] += shift(adj(tmp_0), BACKWARD, mu)*shift(tmp_1, BACKWARD, mu);
    }      
  }


  // Pure Gauge factor
  /*for(int mu=0; mu < Nd; ++mu) { 
    ds_u[mu] *= beta/(Double(2*Nc));
    }*/

  END_CODE("WilsonGaugeAct::dsdu");
}

// Get the gauge action
//
// S = (beta/(2Nc) Sum Re Tr Plaq
//
// w_plaq is defined in MesPlq as
//
// w_plaq =( 2/(V*Nd*(Nd-1)*Nc)) * Sum Re Tr Plaq
//
// so 
// S = beta/(2) * (V*Nd*(Nd-1)/2) w_plaq 
//   = beta/(2) * (V*Nd*(Nd-1)/2)*(2/(V*Nd*(Nd-1)*Nc))* Sum Re Tr Plaq
//   = beta * (1/(Nc)) * Sum Re Tr Plaq

Double
WilsonGaugeAct::S(const multi1d<LatticeColorMatrix>& u) const
{
  Double S_pg;

  Handle< const ConnectState> u_bc(createState(u));
  // Apply boundaries
  //Handle<const ConnectState> u_bc(createState(u));
  Double w_plaq, s_plaq, t_plaq, link;

  // Measure the plaquette
  MesPlq(u_bc->getLinks(), w_plaq, s_plaq, t_plaq, link);

  // Undo Mes Plaq Normalisation
  S_pg = Double( Layout::vol()*Nd*(Nd-1) ) / Double(2);
  S_pg *= beta*w_plaq;
    
  return S_pg;
} 
