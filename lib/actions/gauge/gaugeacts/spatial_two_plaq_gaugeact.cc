// $Id: spatial_two_plaq_gaugeact.cc,v 1.1 2006-07-19 21:42:06 bjoo Exp $
/*! \file
 *  \brief Plaquette gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/spatial_two_plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugeacts/gauge_createstate_factory.h"
#include "actions/gauge/gaugeacts/gauge_createstate_aggregate.h"
#include "meas/glue/mesplq.h"


namespace Chroma
{
 
  namespace SpatialTwoPlaqGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, 
		 multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
							       const std::string& path) 
    {
      return new SpatialTwoPlaqGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
			      SpatialTwoPlaqGaugeActParams(xml, path));
    }

    const std::string name = "SPATIAL_TWO_PLAQ_GAUGEACT";
    const bool registered = TheGaugeActFactory::Instance().registerObject(name, 
									  createGaugeAct);
  };


  SpatialTwoPlaqGaugeActParams::SpatialTwoPlaqGaugeActParams(XMLReader& xml_in, const std::string& path) 
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


  void read(XMLReader& xml, const string& path, SpatialTwoPlaqGaugeActParams& p) 
  {
    SpatialTwoPlaqGaugeActParams tmp(xml, path);
    p=tmp;
  }


  // Internal initializer
  void
  SpatialTwoPlaqGaugeAct::init()
  {
 
    // THis term is only spatial. Really I just need to divide in the 
    // aniso factors
    if ( anisoP() ) { 
      param.coeff /= param.aniso.xi_0;
    }
    QDPIO::cout << "aniso.t_dir" << param.aniso.t_dir << endl;
  }


  //! Compute staple
  /*!
   * \param u_mu_staple   result      ( Write )
   * \param state         gauge field ( Read )
   * \param mu            direction for staple ( Read )
   * \param cb            subset on which to compute ( Read )
   */
  void
  SpatialTwoPlaqGaugeAct::staple(LatticeColorMatrix& u_mu_staple,
		       const Handle< GaugeState<P,Q> >& state,
		       int mu, int cb) const
  {
    QDPIO::cerr << "Not implemented " << endl;
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
  SpatialTwoPlaqGaugeAct::deriv(multi1d<LatticeColorMatrix>& ds_u,
		      const Handle< GaugeState<P,Q> >& state) const
  {
    START_CODE();

    ds_u.resize(Nd);


    const multi1d<LatticeColorMatrix>& u = state->getLinks();



    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    LatticeColorMatrix tmp, tmp2;

    ds_u = zero;
    ds_tmp = zero;

    for(int mu = 0; mu < Nd; mu++) {
      if ( mu == tDir() ) continue;
      
      	for(int nu=mu+1; nu < Nd; nu++) {
	  if( nu == tDir() ||   mu == nu ) continue;

	  // Plaquette Forces -> 
	  // For F_mu
	  //
	  //        (1)               (2)
	  //     < ----- ^      x  ^      |             
	  //     |       |         |      |
	  //     |       |    +    |      |
	  //  x  V       |         <------V

	  // For F_nu
	  //       (3)                (4)
	  //     ------->         <------
	  //            |    +    | 
	  //            |         |
	  //   x <----- V         V----->

	  // First  lets do 
	  //    <-----
	  //   |
	  //   |           (we'll use this for (1) and (4))
	  //   V
	  tmp = shift( adj(u[mu]), FORWARD, nu)*adj(u[nu]);

	  // Now we do 
	  //
	  //           |
	  //           |  (we'll use this for (2) and (3)
	  //           |
	  //   <-------v
	  tmp2 = shift( adj(u[nu]), FORWARD, mu)*adj(u[mu]);


	  // Make munu plaquette which is just adj(tmp2)*tmp1

	  LatticeColorMatrix P_munu = adj(tmp2)*tmp;
	  
	  // Take the real trace
	  LatticeReal re_tr_P_munu = real(trace(P_munu));

	  // Shift to get the t+1 slice
	  LatticeReal s_re_tr_P_munu = shift(re_tr_P_munu, FORWARD, tDir());

	  
	  // Now Term (1)
	  ds_tmp[mu] = u[mu]*shift(u[nu], FORWARD, mu)*tmp;

	  // Now term (2)
	  ds_tmp[mu] +=u[mu]*shift(tmp2*u[nu], BACKWARD, nu);

	  // Now term (3)
	  ds_tmp[nu] = u[nu]*shift( u[mu], FORWARD, nu)*tmp2;
	  
	  // Now term (4)
	  ds_tmp[nu] += u[nu]*shift(tmp*u[mu], BACKWARD, mu);

	  // Zero the force on any fixed boundaries
	  // getGaugeBC().zero(ds_tmp);

	  // Chain rule:
	  //
	  //    d/dt ( P_{ij}(t) P_{ij}(t+1) )
	  //     = [ d/dt P_{ij}(t) ] *  P_{ij}(t+1)
	  //      +  P_{ij}(t) [ d/dt P_{ij}(t+1) ]

	  // Since P_{ij} is a lattice scalar I can pull it to the front:
	  // so:
	  //
	  //  P_{ij}(t+1) [ d/dt P_{ij}(t) ] +  P_{ij}(t) [ d/dt P_{ij}(t+1) ]
	  //
	  //  =>  P_{ij}(t+1) F +  P_{ij}(t) F(t+1).
	  //
	  //  P_ij contributes to F_i and F_j
	  // so 
	  //      P_{ij}(t+1) F_i +  P_{ij}(t) F_i(t+1)
	  //    + P_{ij}(t+1) F_j +  P_{ij}(t) F_j(t+1)

	  //   P_ij(t+1) F_{j} 
	  ds_u[nu]  += s_re_tr_P_munu*ds_tmp[nu];
	  ds_u[nu]  += re_tr_P_munu*shift(ds_tmp[nu], FORWARD, tDir());

	  //            P_ij(t+1) F_{i}
	  ds_u[mu]  += s_re_tr_P_munu*ds_tmp[mu];
	  ds_u[mu]  += re_tr_P_munu*shift(ds_tmp[mu], FORWARD, tDir());
	  
	}
    }


    for(int mu=0; mu < Nd; mu++) {
      if (mu == tDir()) {
	// ds_u[mu] = zero;
      }
      else { 
	ds_u[mu] *= Real(-1)/(Real(8));
      }
    }

 
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
  SpatialTwoPlaqGaugeAct::S(const Handle< GaugeState<P,Q> >& state) const
  {
    Double S_pg = zero;

    // Handle< const GaugeState<P,Q> > u_bc(createState(u));
    // Apply boundaries
    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    // Compute the average plaquettes
    for(int mu=0; mu < Nd; ++mu)
    {
      if( mu == tDir()) continue;

      for(int nu=mu+1; nu < Nd; ++nu)
      {
	if( nu == tDir()) continue;
	/* tmp_0 = u(x+mu,nu)*u_dag(x+nu,mu) */
	/* tmp_1 = tmp_0*u_dag(x,nu)=u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
	/* wplaq_tmp = tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)) */

	LatticeReal P_munu = real(trace(u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu])));
	LatticeReal tmp = shift(P_munu, FORWARD, tDir());
#if 1
	S_pg += sum(P_munu*tmp);
#else
	S_pg += sum(P_munu);
#endif
      }
    }

    // Normalize -> 1 factor of Nc from each P_munu
    S_pg *= Double(-1)/Double(2);

    return S_pg;
  } 

}

