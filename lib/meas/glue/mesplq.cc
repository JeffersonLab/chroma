// $Id: mesplq.cc,v 1.5 2003-02-16 04:14:37 edwards Exp $
/*! \file
 *  \brief Plaquette measurement
 */

#include "chromabase.h"
#include "meas/glue/mesplq.h"


using namespace QDP;

// Primitive way to indicate the time direction
static int tDir() {return Nd-1;}

//! Return the value of the average plaquette normalized to 1
/*!
 * \param u -- gauge field (Read)
 * \param w_plaq -- plaquette average (Write)
 * \param s_plaq -- space-like plaquette average (Write)
 * \param t_plaq -- time-like plaquette average (Write)
 * \param link   -- space-time average link (Write)
 */

void MesPlq(const multi1d<LatticeColorMatrix>& u, Double& w_plaq, Double& s_plaq, 
	    Double& t_plaq, Double& link)
{
  s_plaq = t_plaq = w_plaq = link = 0.0;

  // Compute the average plaquettes
  for(int mu=1; mu < Nd; ++mu)
  {
    for(int nu=0; nu < mu; ++nu)
    {
#if 0
      /* tmp_0 = u(x+mu,nu)*u_dag(x+nu,mu) */
      LatticeColorMatrix tmp_0 = shift(u[nu],FORWARD,mu) * adj(shift(u[mu],FORWARD,nu));

      /* tmp_1 = tmp_0*u_dag(x,nu)=u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
      LatticeColorMatrix tmp_1 = tmp_0 * adj(u[nu]);

      /* tmp = sum(tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu))) */
      Double tmp = sum(real(trace(u[mu]*tmp_1)));

#else
      /* tmp_0 = u(x+mu,nu)*u_dag(x+nu,mu) */
      /* tmp_1 = tmp_0*u_dag(x,nu)=u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
      /* wplaq_tmp = tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)) */
      Double tmp = 
	sum(real(trace(u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu]))));
#endif
      w_plaq += tmp;

      if (mu == tDir() || nu == tDir())
	t_plaq += tmp;
      else 
	s_plaq += tmp;
    }
  }
  
  // Normalize
  w_plaq *= 2.0 / double(Layout::vol()*Nd*(Nd-1)*Nc);
  
  if (Nd > 2) 
    s_plaq *= 2.0 / double(Layout::vol()*(Nd-1)*(Nd-2)*Nc);
  
  t_plaq /= double(Layout::vol()*(Nd-1)*Nc);
  

  // Compute the average link
  for(int mu=0; mu < Nd; ++mu)
    link += sum(real(trace(u[mu])));

  link /= double(Layout::vol()*Nd*Nc);
}
