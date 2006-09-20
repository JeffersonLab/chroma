// $Id: pg_gaugeact.cc,v 3.3 2006-09-20 20:28:00 edwards Exp $
/*! \file
 *  \brief Parallelogram gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/pg_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma
{
 
  namespace PgGaugeActEnv 
  {
    //! Callback
    GaugeAction< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
											    const std::string& path) 
    {
      return new PgGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
			    PgGaugeActParams(xml, path));
    }

    const std::string name = "PG_GAUGEACT";

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


  // Param constructor/reader
  PgGaugeActParams::PgGaugeActParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);

    try {
      read(paramtop, "./coeff", coeff);
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " <<  e << endl;
      QDP_abort(1);
    }
  }

  // Read params
  void read(XMLReader& xml, const string& path, PgGaugeActParams& p) 
  {
    PgGaugeActParams tmp(xml, path);
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
  PgGaugeAct::staple(LatticeColorMatrix& u_staple,
		     const Handle< GaugeState<P,Q> >& state,
		     int mu, int cb) const
  {
    QDPIO::cerr << "PgGaugeAct::staple() - not converted from szin" << endl;
    QDP_abort(1);
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
  PgGaugeAct::deriv(multi1d<LatticeColorMatrix>& ds_u,
		    const Handle< GaugeState<P,Q> >& state) const
  {
    START_CODE();

    ds_u.resize(Nd);

    LatticeColorMatrix tmp_0;
    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;
    LatticeColorMatrix tmp_3;
    LatticeColorMatrix tmp_4;
    LatticeColorMatrix tmp_5;
    LatticeColorMatrix tmp_tot;
    multi1d<int> fdir(Nd);  

    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    ds_u = zero;

         
    // It is 1/(4Nc) to account for normalisation relevant to fermions
    // in the taproj, which is a factor of 2 different from the 
    // one used here.

    Real coeff_tmp = Real(-1)*Real(coeff)/(Real(2*Nc));

    for(int mu = 0; mu < Nd; ++mu)
    {
      tmp_tot = zero;
    
      /* generalized parallelogram */
      for(int nu=mu+1; nu < Nd; nu++)
      {
	for(int j=0, k=0; k < Nd; k++)
	  if((k!=mu) && (k!=nu))
	    fdir[j++] = k;

	tmp_3 = zero;
	
	for(int cb=0; cb < 2; cb++)
	{
	  /* |_  part fitting piece is ~| */
	  tmp_0[rb[1-cb]] = shift(u[mu], FORWARD, nu) * shift(adj(u[nu]), FORWARD, mu);
	  tmp_5[rb[1-cb]] = tmp_0 * coeff_tmp;

	  for(int k=0; k < Nd-2; k++)
	  {
	    int eta = fdir[k];
	    tmp_2[rb[cb]] = shift(u[eta], FORWARD, nu) * shift(tmp_5, FORWARD, eta);

	    if(k==0)
	    {
	      tmp_4[rb[cb]] = tmp_2 * shift(adj(u[eta]), FORWARD, mu);
	    }
	    else
	    {
              tmp_4[rb[cb]] += tmp_2 * shift(adj(u[eta]), FORWARD, mu);
	    }
	  }
	  tmp_3[rb[cb]] += tmp_4 * adj(u[mu]);
	  tmp_tot[rb[cb]] += adj(tmp_4) * adj(u[nu]);

	  /* _| piece */
	  tmp_0[rb[cb]] = u[nu] * shift(u[mu], FORWARD, nu);
	  tmp_5[rb[cb]] = tmp_0 * coeff_tmp;
	  for(int k=0; k < Nd-2; k++)
	  {
	    int eta = fdir[k];
	    tmp_2[rb[1-cb]] = shift(adj(u[eta])*tmp_5, BACKWARD, mu);
	    tmp_1[rb[1-cb]] = shift(u[eta], FORWARD, nu);
	    if(k==0)
	    {
	      tmp_4[rb[cb]] = shift(tmp_2*tmp_1, BACKWARD, eta);
	    }
	    else
	    {
	      tmp_4[rb[cb]] += shift(tmp_2*tmp_1, BACKWARD, eta);
	    }
	  }
	  tmp_3[rb[cb]] += adj(tmp_4) * shift(u[mu], BACKWARD, mu);
	  tmp_tot[rb[1-cb]] += shift(u[nu]*adj(tmp_4), FORWARD, mu);

	  /* ~| part */
	  tmp_0[rb[1-cb]] = adj(u[nu]) * u[mu];
	  tmp_5[rb[cb]] = shift(tmp_0, BACKWARD, nu) * coeff_tmp;
	  for(int k=0; k < Nd-2; k++)
	  {
	    int eta = fdir[k];
	    tmp_2[rb[1-cb]] = u[eta] * shift(tmp_5, FORWARD, eta);
	    if(k==0)
	    {
	      tmp_4[rb[cb]] = shift(tmp_2, BACKWARD, mu) * adj(shift(u[eta], BACKWARD, nu));
	    }
	    else
	    {
	      tmp_4[rb[cb]] += shift(tmp_2, BACKWARD, mu) * adj(shift(u[eta], BACKWARD, nu));
	    }
	  }

	  tmp_1[rb[cb]] = shift(u[nu], BACKWARD, nu);
	  tmp_tot[rb[1-cb]] += shift(adj(tmp_1)*adj(tmp_4), FORWARD, mu);
	  tmp_1[rb[cb]] = shift(u[mu], BACKWARD, mu);
	  tmp_3[rb[1-cb]] += shift(adj(tmp_1)*tmp_4, FORWARD, nu);
	  
	  /* Finally the |~ part */
	  tmp_0[rb[cb]] = u[mu] * shift(u[nu], FORWARD, mu);
	  tmp_5[rb[cb]] = tmp_0 * coeff_tmp;
	  for(int k=0; k < Nd-2; k++)
	  {
	    int eta = fdir[k];
	    tmp_2[rb[1-cb]] = shift(adj(u[eta])*tmp_5, BACKWARD, nu);
	    tmp_1[rb[1-cb]] = shift(u[eta], FORWARD, mu);
	    if(k==0)
	    {
	      tmp_4[rb[cb]] = shift(tmp_2*tmp_1, BACKWARD, eta);
	    }
	    else
	    {
	      tmp_4[rb[cb]] += shift(tmp_2*tmp_1, BACKWARD, eta);
	    }
	  }
	  tmp_tot[rb[cb]] += adj(tmp_4) * shift(u[nu], BACKWARD, nu);
	  tmp_3[rb[1-cb]] += shift(u[mu]*adj(tmp_4), FORWARD, nu);

	} /* end loop over cb */	  

	/* mult tmp_3 by u(x,nu) matrix) */
	ds_u[nu] += u[nu] * tmp_3;
      } /* end loop over nu */
        
      /* ds_u =  u(x,mu) * tmp_tot */
      ds_u[mu] += u[mu] * tmp_tot;
    }
  
    // Zero the force on any fixed boundaries
    getGaugeBC().zero(ds_u);

    END_CODE();
  }



  // Get the gauge action
  //
  // S = -(coeff/(Nc) Sum Re Tr Pg
  //
  Double
  PgGaugeAct::S(const Handle< GaugeState<P,Q> >& state) const
  {
    START_CODE();

    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    LatticeColorMatrix tmp_0;
    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;
    LatticeReal lgimp = zero;

    // Parallelogram term
    if (Nd != 4)
    {
      QDPIO::cerr << "Parallelogram implemented only for Nd==4" << endl;
      QDP_abort(1);
    }

    for(int dir = 0; dir < Nd; ++dir) /* orthogonal direction to the 3-cube */
    {
      /*+ */
      /* There are 3 orientations for  */
      /* (mu,nu,rho,-mu,-nu,-rho) */
      /*- */
      for(int mu0 = 0; mu0 < 2; ++mu0) /* 3 orientations in a 3-cube */
      {
	int mu = mu0;
	if (mu >= dir) 
	  ++mu;

	for(int nu0 = 0; nu0 <= 2; ++nu0)
	{
	  if (nu0 == mu0) continue;

	  int nu = nu0;
	  if (nu >= dir) 
	    ++nu;

	  for(int rho0 = mu0+1; rho0 <= 2; ++rho0)
	  {
	    if (rho0 == nu0) continue;
	    
	    int rho = rho0;
	    if (rho >= dir) 
	      ++rho;

	    for(int cb = 0; cb < 2; ++cb)
	    {
	      /* 6-link parallelogram. */

	      /* tmp_1(x) = u(x+nu,rho) */
	      tmp_1[rb[1-cb]] = shift(u[rho], FORWARD, nu);

	      /* tmp_2 = u(x,nu) * tmp_1 = u(x,nu)*u(x+nu,rho) */
	      tmp_2[rb[1-cb]] = u[nu] * tmp_1;

	      /* tmp_1(x) = tmp_2(x+mu) */
	      tmp_1[rb[cb]] = shift(tmp_2, FORWARD, mu);

	      /* tmp_2 = u(x,mu) * tmp_1 = u(x,mu)*u(x+mu,nu)*u(x+mu+nu,rho) */
	      tmp_2[rb[cb]] = u[mu] * tmp_1;

	      /* tmp_0 = u_dag(x,rho) * tmp_2 */
	      /*       = u_dag(x,rho)*u(x,mu)*u(x+mu,nu)*u(x+mu+nu,rho) */
	      tmp_0[rb[cb]] = adj(u[rho]) * tmp_2;

	      /* tmp_1(x) = u(x+nu,mu) */
	      tmp_1[rb[1-cb]] = shift(u[mu], FORWARD, nu);

	      /* tmp_2 = u(x,nu) * tmp_1 = u(x,nu)*u(x+nu,mu) */
	      tmp_2[rb[1-cb]] = u[nu] * tmp_1;

	      /* tmp_1(x) = tmp_2(x+rho) */
	      tmp_1[rb[cb]] = shift(tmp_2, FORWARD, rho);

	      /* lgimp += trace(tmp_1_dag * tmp_0) */
	      /*       += trace(u_dag(x+nu+rho,mu)*u_dag(x+rho,nu)*u_dag(x,rho) */
	      /*         *u_(x,mu)*u(x+mu,nu)*u(x+mu+nu,rho)) */
	      lgimp[rb[cb]] += real(trace(adj(tmp_1) * tmp_0));
	    }
	  }
	}
      }

      /*+ */
      /* 4th orientation in a 3-cube */
      /* An arbitrary orientation is taken for the ordering */
      /* (mu,-nu,rho,-mu,nu,-rho) */
      /*- */
      int mu  = 0;
      int nu  = 1;
      int rho = 2;
      if (mu  >= dir) mu  = mu + 1;
      if (nu  >= dir) nu  = nu + 1;
      if (rho >= dir) rho = rho + 1;

      for(int cb = 0; cb < 2; ++cb)
      {
	/* 6-link parallelogram. */

	/* tmp_1 = u_dag(x,nu) * u(x,rho) */
	tmp_1[rb[cb]] = adj(u[nu]) * u[rho];

	/* tmp_2(x) = tmp_1(x-rho) */
	tmp_2[rb[1-cb]] = shift(tmp_1, BACKWARD, rho);

	/* tmp_1(x) = tmp_2(x+mu) */
	tmp_1[rb[cb]] = shift(tmp_2, FORWARD, mu);

	/* tmp_2 = tmp_1 * u_dag(x,mu) */
	/*       = u_dag(x+mu-rho,nu)*u(x+mu-rho,rho)*u_dag(x,mu) */
	tmp_2[rb[cb]] = tmp_1 * adj(u[mu]);

	/* tmp_1 = tmp_2 * u(x,nu) */
	/*       = u_dag(x+mu-rho,nu)*u(x+mu-rho,rho)*u_dag(x,mu)*u(x,nu) */
	tmp_1[rb[cb]] = tmp_2 * u[nu];

	/* tmp_2(x) = tmp_1(x-nu) */
	tmp_2[rb[1-cb]] = shift(tmp_1, BACKWARD, nu);

	/* tmp_1(x) = tmp_2(x+rho) */
	tmp_1[rb[cb]] = shift(tmp_2, FORWARD, rho);

	/* tmp_0 = u(x,mu) * tmp_1 */
	/*       = u(x,mu)*u_dag(x+mu-nu,nu)*u(x+mu-nu,rho) */
	/*        *u_dag(x-nu+rho,mu)*u(x-nu+rho,nu) */
	tmp_0[rb[cb]] = u[mu] * tmp_1;

	/* lgimp += trace(tmp_0 * u_dag(x,rho)) */
	/*       += trace(u(x,mu)*u_dag(x+mu-nu,nu)*u(x+mu-nu,rho) */
	/*         *u_dag(x-nu+rho,mu)*u(x-nu+rho,nu)*u_dag(x,rho)) */
//	lgimp[rb[cb]] += real(trace(tmp_0 * adj(u[rho])));
	lgimp[rb[cb]] += real(trace(adj(u[rho]) * tmp_0));
      }
    }

    Double S_pg = sum(lgimp);
    S_pg *= -coeff / Real(Nc);      // note sign
  
    END_CODE();

    return S_pg;
  } 

}

