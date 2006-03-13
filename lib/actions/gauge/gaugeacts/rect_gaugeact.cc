// $Id: rect_gaugeact.cc,v 2.1 2006-03-13 05:13:43 edwards Exp $
/*! \file
 *  \brief Rectangle gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugebcs/gaugebc_aggregate.h"

namespace Chroma
{
 
  namespace RectGaugeActEnv 
  {
    //! Callback
    GaugeAction* createGaugeAct(XMLReader& xml, const std::string& path) 
    {
      return new RectGaugeAct(GaugeTypeGaugeBCEnv::reader(xml, path), 
			      RectGaugeActParams(xml, path));
    }

    const std::string name = "RECT_GAUGEACT";
    const bool registered = TheGaugeActFactory::Instance().registerObject(name, 
									  createGaugeAct);
  };


  // Param constructor/reader
  RectGaugeActParams::RectGaugeActParams(XMLReader& xml_in, const std::string& path) {
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
  void read(XMLReader& xml, const string& path, RectGaugeActParams& p) 
  {
    RectGaugeActParams tmp(xml, path);
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
  RectGaugeAct::staple(LatticeColorMatrix& u_staple,
		       Handle<const ConnectState> state,
		       int mu, int cb) const
  {
    QDPIO::cerr << "RectGaugeAct::staple() - not converted from szin" << endl;
    QDP_abort(1);

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
  RectGaugeAct::dsdu(multi1d<LatticeColorMatrix>& ds_u,
		     const Handle< const ConnectState> state) const
  {
    START_CODE();

    ds_u.resize(Nd);

    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;
    LatticeColorMatrix tmp_3;
    LatticeColorMatrix tmp_4;
    LatticeColorMatrix tmp_5;
    LatticeColorMatrix tmp_6;
    LatticeColorMatrix tmp_tot;

    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    ds_u = zero;

    // It is 1/(4Nc) to account for normalisation relevant to fermions
    // in the taproj, which is a factor of 2 different from the 
    // one used here.

    Real coeff_tmp = Real(-1)*Real(coeff)/Real(2*Nc);

    for(int mu = 0; mu < Nd; ++mu)
    {
      tmp_tot = zero;
    
      /*  2*1 rectangles */

      /* calculate double_mu links */
      tmp_3 = u[mu] * shift(u[mu], FORWARD, mu);
    
      for(int cb=0; cb < 2; cb++)
      {
	tmp_6[rb[1-cb]] = tmp_3 * coeff_tmp;

	for(int j=0, nu=0; nu < Nd; nu++)
	{
	  if(nu == mu)
	    continue;

	  /* forward plaquette */
	  tmp_1[rb[cb]] = u[nu] * shift(tmp_6, FORWARD, nu);
	  tmp_5[rb[1-cb]] = shift(adj(tmp_1)*tmp_3, BACKWARD, mu);
	  ds_u[nu][rb[cb]] += u[nu] * shift(tmp_5, BACKWARD, mu);
	  tmp_5[rb[1-cb]] = shift(u[nu], FORWARD, mu);

	  /* at this point we add the nu contribution directly to ds_u */
	  /* we could make tmp_tot carry a direction index in order to avoid that */
	  if(j++ == 0)
	  {
 	    tmp_4[rb[cb]] = tmp_1 * adj(shift(tmp_5, FORWARD, mu));
	    ds_u[nu][rb[cb]] += tmp_4 * adj(tmp_3);
	  }
	  else
	  {
	    tmp_2 = tmp_1 * adj(shift(tmp_5, FORWARD, mu));
	    tmp_4 += tmp_2; /* sum to the staple */
	    ds_u[nu][rb[cb]] += tmp_2 * adj(tmp_3);
	  }
	  
	  /* backward plaquette */
	  tmp_1[rb[1-cb]] = adj(u[nu]) * tmp_6;
	  tmp_2[rb[cb]] = shift(u[nu], FORWARD, mu);
	  tmp_5[rb[1-cb]] = shift(tmp_2, FORWARD, mu);
	  tmp_4[rb[cb]] += shift(tmp_1*tmp_5, BACKWARD, nu);
	}

	tmp_tot[rb[cb]] += shift(u[mu], FORWARD, mu) * adj(tmp_4);
	tmp_tot[rb[1-cb]] += shift(adj(tmp_4)*u[mu], BACKWARD, mu);
      }
    
      /* ds_u =  u(x,mu) * tmp_tot */
      ds_u[mu] += u[mu] * tmp_tot;
    }
  
    // Zero the force on any fixed boundaries
    getGaugeBC().zero(ds_u);

    END_CODE();
  }


  // Get the gauge action
  //
  // S = -(coeff/(Nc) Sum Re Tr Rect
  //
  // w_rect is defined in MesPlq as
  //
  // w_rect =( 2/(V*Nd*(Nd-1)*Nc)) * Sum Re Tr Rect
  //
  // so 
  // S = -coeff * (V*Nd*(Nd-1)/2) w_rect 
  //   = -coeff * (V*Nd*(Nd-1)/2)*(2/(V*Nd*(Nd-1)*Nc))* Sum Re Tr Rect
  //   = -coeff * (1/(Nc)) * Sum Re Tr Rect

  Double
  RectGaugeAct::S(const Handle<const ConnectState> state) const
  {
    START_CODE();

    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    LatticeColorMatrix tmp_0;
    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;
    LatticeColorMatrix tmp_tot;
    LatticeReal lgimp = zero;

    // 2x1 rectangle piece

    for(int mu = 1; mu < Nd; ++mu)
    {
      for(int nu = 0; nu < mu; ++nu)
      {
	for(int cb = 0; cb < 2; ++cb)
	{
	  /* 1x2 (vertical) rectangle */

	  /* tmp_1(x) = u(x-nu,nu) */
	  tmp_1[rb[cb]] = shift(u[nu], BACKWARD, nu);

	  /* tmp_0 = tmp_1 * u(x,nu) = u(x,nu)*u(x+nu,nu) */
	  tmp_0[rb[cb]] = tmp_1 * u[nu];

	  /* tmp_2(x) = tmp_0(x+mu) */
	  tmp_2[rb[1-cb]] = shift(tmp_0, FORWARD, mu);

	  /* tmp_1(x) = u(x+nu,mu) */
	  tmp_1[rb[1-cb]] = shift(u[mu], FORWARD, nu);

	  /* tmp_0 = tmp_2 * tmp_1_dag = u(x+mu-nu,nu)*u(x+mu,nu)*u_dag(x+nu,mu) */
	  tmp_0[rb[1-cb]] = tmp_2 * adj(tmp_1);

	  /* tmp_1 = tmp_0 * u_dag(x,nu) */
	  /*       = u(x+mu-nu,nu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
	  tmp_1[rb[1-cb]] = tmp_0 * adj(u[nu]);

	  /* tmp_2(x) = tmp_1(x+nu) */
	  tmp_2[rb[cb]] = shift(tmp_1, FORWARD, nu);

	  /* tmp_tot = tmp_2 * u_dag(x,nu) */
	  /*         = u(x+mu,nu)*u(x+mu+nu,nu)*u_dag(x+2*nu,mu) */
	  /*          *u_dag(x+nu,nu)*u_dag(x,nu) */
	  tmp_tot[rb[cb]] = tmp_2 * adj(u[nu]);


	  /* 2x1 (horizontal) rectangle */

	  /* tmp_1(x) = u(x+mu,nu) */
	  tmp_1[rb[1-cb]] = shift(u[nu], FORWARD, mu);

	  /* tmp_0 = u(x,nu) * tmp_1 = u(x,mu)*u(x+mu,nu) */
	  tmp_0[rb[1-cb]] = u[mu] * tmp_1;

	  /* tmp_2(x) = tmp_0(x-nu) */
	  tmp_2[rb[cb]] = shift(tmp_0, BACKWARD, nu);

	  /* tmp_0 = tmp_2 * u_dag(x,mu) = u(x-nu,mu)*u(x+mu-nu,nu)*u_dag(x,mu) */
	  tmp_0[rb[cb]] = tmp_2 * adj(u[mu]);

	  /* tmp_1(x) = tmp_0(x+mu) */
	  tmp_1[rb[1-cb]] = shift(tmp_0, FORWARD, mu);

	  /* tmp_0 = tmp_1 * u_dag(x,mu) */
	  /*       = u(x+mu-nu,mu)*u(x+2*mu-nu,nu)*u_dag(x+mu,mu)*u_dag(x,mu) */
	  tmp_0[rb[1-cb]] = tmp_1 * adj(u[mu]);

	  /* tmp_2(x) = tmp_0(x+nu) */
	  tmp_2[rb[cb]] = shift(tmp_0, FORWARD, nu);

	  /* tmp_tot += tmp_2 * u_dag(x,nu) */
	  /*         += u(x+mu,mu)*u(x+2*mu,nu)*u_dag(x+mu+nu,mu)*u_dag(x+nu,mu) */
	  /*           *u_dag(x,nu) */
	  tmp_tot[rb[cb]] += tmp_2 * adj(u[nu]);

	  /* lgimp += trace(u(x,mu) * tmp_tot) */
	  lgimp[rb[cb]] += real(trace(u[mu] * tmp_tot));
	}
      }
    }
  
    Double S_rect = sum(lgimp);
    S_rect *= -coeff / Real(Nc);   // note sign here
  
    END_CODE();
    
    return S_rect;
  } 

}

