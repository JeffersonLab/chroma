// $Id: wilson_gaugeact.cc,v 1.1 2004-03-03 01:50:12 edwards Exp $
/*! \file
 *  \brief Wilson gauge action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/wilson_gaugeact.h"


//! Compute staple
/*!
 * \param u_staple   result      ( Write )
 * \param state      gauge field ( Read )
 * \param mu         direction for staple ( Read )
 * \param cb         subset on which to compute ( Read )
 */

void
WilsonGaugeAct::staple(multi1d<LatticeColorMatrix>& u_staple,
		       Handle<const ConnectState> state,
		       int mu, int cb)
{
  START_CODE("WilsonGaugeAct::staple");

  const multi1d<LatticeColorMatrix>& u = state->getLinks();
				 
  LatticeColorMatrix tmp_1;

  if (getSet() != 2)
  {
    QDPIO::cerr << "WilsonGaugeAct::dsdu  implemented only for even/odd" << endl;
    QDP_abort(1);
  }

  const Real xi02 = anisoFactor() * anisoFactor();  
  const int t_dir = tDir();
  
  u_staple = 0;

  for(int nu = 0; nu < Nd; ++nu)
  {
    if (nu == mu) continue;

    /* Forward staple */
    /* tmp_1(x) = u(x+mu,nu)*u_dag(x+nu,mu) */
    tmp_1[rb[cb]] = shift(u[nu], FORWARD, mu) * shift(adj(u[mu]), FORWARD, nu);

    if( anisoP() )  	
      if( mu == t_dir || nu == t_dir )
	tmp_1[rb[cb]] *= xi02;

    /* u_staple(x) +=  tmp_1 * u_dag(x,nu)
       += u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
    u_staple[rb[cb]] += tmp_1 * adj(u[nu]);
          

    /* Backward staple */
    /* tmp_1(x) = u(x,mu)*u(x+mu,nu) */
    tmp_1[rb[1-cb]] = u[mu] * shift(u[nu], FORWARD, mu);

    if( anisoP() )  	
      if( mu == t_dir || nu == t_dir )
	tmp_1[rb[1-cb]] = xi02 * tmp_1;

    /* u_staple(x) += tmp_1_dag(x-nu) * u(x-nu,nu)
       += u_dag(x+mu-nu,nu)*u_dag(x-nu,mu)*u(x-nu,nu) */
    u_staple[rb[cb]] += shift(adj(tmp_1), BACKWARD, nu) * shift(u[nu], BACKWARD, nu);
  }  /* closes nu loop */

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
		     Handle<const ConnectState> state)
{
  START_CODE("WilsonGaugeAct::dsdu");

  const multi1d<LatticeColorMatrix>& u = state->getLinks();
				 
  LatticeColorMatrix tmp_0;
  LatticeColorMatrix tmp_1;
  LatticeColorMatrix tmp_2;

  ds_u = 0;
  
  for(int mu = 0; mu < Nd; ++mu)
  {
    for(int nu=mu+1; nu<Nd; nu++)
    {
      tmp_0 = shift(u[mu], FORWARD, nu) * adj(shift(u[nu], FORWARD, mu));
      tmp_1 = tmp_0 * adj(u[mu]);
      tmp_2 = u[nu] * tmp_1;
      ds_u[nu] += tmp_2 + shift(adj(tmp_0)*adj(u[nu])*u[mu], BACKWARD, mu);
      ds_u[mu] += adj(tmp_2) + shift(tmp_1*u[nu], BACKWARD, nu);
    }      
  }
  
  END_CODE("WilsonGaugeAct::dsdu");
}
