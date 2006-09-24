// $Id: plaq_plus_spatial_two_plaq_gaugeact.cc,v 3.7 2006-09-24 20:58:30 edwards Exp $
/*! \file
 *  \brief Plaquette gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/plaq_plus_spatial_two_plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"
#include "meas/glue/mesplq.h"
#include "util/gauge/taproj.h"

namespace Chroma
{
 
  namespace PlaqPlusSpatialTwoPlaqGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, 
		 multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
							       const std::string& path) 
    {
      return new PlaqPlusSpatialTwoPlaqGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
			      PlaqPlusSpatialTwoPlaqGaugeActParams(xml, path));
    }

    const std::string name = "PLAQ_PLUS_SPATIAL_TWO_PLAQ_GAUGEACT";

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
    
    static double time_spent = 0;
    double getTime() { return time_spent; }
  };


  PlaqPlusSpatialTwoPlaqGaugeActParams::PlaqPlusSpatialTwoPlaqGaugeActParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);

    try {
      //  Read optional anisoParam.
      if (paramtop.count("AnisoParam") != 0) 
	read(paramtop, "AnisoParam", aniso);
      
      read(paramtop, "./coeff_plaq_s", coeff_plaq_s);
      read(paramtop, "./coeff_plaq_t", coeff_plaq_t);
      read(paramtop, "./coeff_two_plaq", coeff_two_plaq);


    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " <<  e << endl;
      QDP_abort(1);
    }
  }


  void read(XMLReader& xml, const string& path, PlaqPlusSpatialTwoPlaqGaugeActParams& p) 
  {
    PlaqPlusSpatialTwoPlaqGaugeActParams tmp(xml, path);
    p=tmp;
  }


  // Internal initializer
  void
  PlaqPlusSpatialTwoPlaqGaugeAct::init()
  {
    START_CODE();
 
    if ( anisoP() ) { 
      // SPatial guys divided by xi_0
      param.coeff_plaq_s /= param.aniso.xi_0;
      param.coeff_two_plaq /= param.aniso.xi_0;

      // Temporal guy mutiliplied by xi_0
      param.coeff_plaq_t *= param.aniso.xi_0;

    }
    QDPIO::cout << "aniso.t_dir" << param.aniso.t_dir << endl;
    
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
  PlaqPlusSpatialTwoPlaqGaugeAct::staple(LatticeColorMatrix& u_mu_staple,
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
   *
   * This is a funny action. It is  P(x) P(x+t) 
   * Where P(x) is the real trace of the plaquette
   *
   * The chain rule means the derivative is:
   *
   *    \dot{ P(x) } P(x+t) + P(x) \dot{ P(x + t) }
   *
   * I can resum this, in a different way to collect all the terms
   * on a link U(x). This resumming involves x + t -> x
   * 
   * so the derivative then becomes
   *
   *  \dot{ P(x) } P(x+t) + P(x-t) \dot{ P(x) }
   *
   * = [ P(x+t) + P(x-t) ] \dot{ P(x) } = P_sum \dot{P(x)}
   * 
   * Now \dot{P(x)} contains the usual force contribugtions from the plaquette
   * and P_sum is an array of lattice reals.
   *
   * The only remaining complication is that I generate force contributions
   * in the routine for U_(x,i), U_(x+j,i), U_(x,j) and U_(x+i,j)
   *
   * and this is done by shifting staples from x to x+j and x+i
   *
   * during this shifting I have also to ensure the same shifting of P_sum
   *  
   *
   * \param ds_u       result      ( Write )
   * \param state      gauge field ( Read )
   */
  void
  PlaqPlusSpatialTwoPlaqGaugeAct::deriv(multi1d<LatticeColorMatrix>& ds_u,
		      const Handle< GaugeState<P,Q> >& state) const
  {
    START_CODE();

    QDP::StopWatch swatch;
    swatch.reset();
    swatch.start();

    ds_u.resize(Nd);

    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    LatticeColorMatrix tmp, tmp2;

    int t_dir = tDir();


    ds_tmp = zero;

    // Do the spatial plaquettes and 2 plaquettes first
    {    
      Real factor = param.coeff_two_plaq / Real(2*Nc);
      LatticeReal l_plaq_coeff_s = Real(param.coeff_plaq_s);

      for(int mu = 0; mu < Nd; mu++) {
	if ( mu == t_dir ) continue;
	
      	for(int nu=mu+1 ; nu < Nd; nu++) {
	  if( nu == t_dir ) continue;
	  
	  // Plaquette Forces -> 
	  // For F_i
	  //
	  //    U(x,i) term        U(x+j,i) term
	  //        (1)               (2)
	  //     < ----- ^      x  ^      |             
	  //     |       |         |      |
	  //     |       |    +    |      |
	  //  x  V       |         <------V

	  // For F_nu
	  //
	  //    U(x, j) term       U(x+i, j) term
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
	  //
	  //
	  LatticeColorMatrix u_mu_plus_nu = shift(u[mu], FORWARD, nu);
	  LatticeColorMatrix u_nu_plus_mu = shift(u[nu], FORWARD, mu);

	  tmp = adj(u_mu_plus_nu)*adj(u[nu]);

	  // Now we do 
	  //
	  //           |
	  //           |  (we'll use this for (2) and (3)
	  //           |
	  //   <-------v
	  tmp2 = adj(u_nu_plus_mu)*adj(u[mu]);


	  // Make munu plaquette which is just adj(tmp2)*tmp1

	 
	  // Take the real trace
	  LatticeReal P_munu = real(trace(adj(tmp2)*tmp));

	  // NOw I have c_plaq_s + (c_plaq_t/2Nc)(P_munu(x-t)+P_munu(x+t));

	  LatticeReal P_sum = l_plaq_coeff_s +
	    factor*(shift(P_munu, FORWARD, t_dir)+
	            shift(P_munu, BACKWARD, t_dir));

	  // It is tmp and tmp2 that will
	  // need to move when I generate the force contrib
	  // for U(x+i,j) and U(x+j,i) (in terms 2 and 4)
	  // It is at these times I also ned to shift the P_sum
	  // in exactly the same way. Terms (1) and (2) 
	  // have tmp and tmp2 unshifted and require the local P_sum
	  //
	  // Since the P_sum is made of lattice reals  multiplication
	  // is commutative, so it is safe to multiply them into 
	  // the tmp and tmp2 right here.
	  
	  tmp *= P_sum;
	  tmp2 *= P_sum;

	  // Now Term (1)
	  ds_tmp[mu] += u_nu_plus_mu*tmp;

	  // Now term (2)
	  ds_tmp[mu] += shift(tmp2*u[nu], BACKWARD, nu);

	  // Now term (3)
	  ds_tmp[nu] += u_mu_plus_nu*tmp2;
	  
	  // Now term (4)
	  ds_tmp[nu] += shift(tmp*u[mu], BACKWARD, mu);
	  
	}
      }
    } // l_plaq_coeff_s goes away here.

    
    // Now just do the temporal forces
    for(int nu=0; nu < Nd; nu++) {
      if( nu == t_dir ) continue;

      // Plaquette Forces -> 
      // For F_i
      //
      //    U(x,i) term        U(x+j,i) term
      //        (1)               (2)
      //     < ----- ^      x  ^      |             
      //     |       |         |      |
      //     |       |    +    |      |
      //  x  V       |         <------V
      
      // For F_nu
      //
      //    U(x, j) term       U(x+i, j) term
      //       (3)                (4)
      //     ------->         <------
      //            |    +    | 
      //            |         |
      //   x <----- V         V----->
      LatticeColorMatrix u_t_plus_nu = shift(u[t_dir], FORWARD, nu);
      LatticeColorMatrix u_nu_plus_t = shift(u[nu], FORWARD, t_dir);
      
      // First  lets do 
      //    <-----
      //   |
      //   |           (we'll use this for (1) and (4))
      //   V
      tmp = adj(u_t_plus_nu)*adj(u[nu]);
      
      // Now we do 
      //
      //           |
      //           |  (we'll use this for (2) and (3)
      //           |
      //   <-------v
      
      tmp2 = adj(u_nu_plus_t)*adj(u[t_dir]);
      

      // Now Term (1)
      // I can write this directly to ds_tmp[t_dir]
      // which has not been used before
      ds_tmp[t_dir] += u_nu_plus_t*tmp;
      
      // Now term (2)
      ds_tmp[t_dir] += shift(tmp2*u[nu], BACKWARD, nu);
      
      // Terms 3 and 4: I use ds_u[nu] as a dummy 
      // to collect both terms. Then I multiply with the temporal 
      // coefficient and add to ds_tmp[nu]

      // Now term (3)
      ds_u[nu] = u_t_plus_nu*tmp2;
      
      // Now term (4)
      ds_u[nu] += shift(tmp*u[t_dir], BACKWARD, t_dir);
      
      // Multiply in the coeff for the nu guys
      ds_tmp[nu] += param.coeff_plaq_t*ds_u[nu];
    }
    
    // Now multiply in the ds_tmp[t_dir] all at once
    ds_tmp[t_dir] *= param.coeff_plaq_t;

    // Finally multiply by u[mu] in the spatial dirs 
    for(int mu=0; mu < Nd; mu++) {
	ds_u[mu] = u[mu]*ds_tmp[mu];
	ds_u[mu] *= Real(-1)/Real(2*Nc);
    }

 
    // Zero the force on any fixed boundaries
    getGaugeBC().zero(ds_u);
    
    swatch.stop();
    PlaqPlusSpatialTwoPlaqGaugeActEnv::time_spent += swatch.getTimeInSeconds();    

    END_CODE();
  }

  // Get the gauge action
  //
  // S = -coeff_plaq_s Sum P_{ij}-(coeff_two_plaq/2) Sum P_{ij}(x) P_{ij}(x+t)
  //     -coeff_plaq_t Sum P_{t}
  //
  // where P_{ij}(x) is the plaquette on lattice point x and i,j are 
  // only spatial directions.
  //
  // P has normalisatiom (1/Nc) and so the overall normalisation
  // for the product is (1/Nc^2). Giving the final multiplier as
  //
  //  -coeff/(2 Nc^2 )
  Double
  PlaqPlusSpatialTwoPlaqGaugeAct::S(const Handle< GaugeState<P,Q> >& state) const
  {
    START_CODE();

    Double S_pg = zero;
    

    Real factor = param.coeff_two_plaq /(Real(2)*Real(Nc));

    // Handle< const GaugeState<P,Q> > u_bc(createState(u));
    // Apply boundaries
    const multi1d<LatticeColorMatrix>& u = state->getLinks();
    int t_dir = tDir();

    
    // Compute the spatial terms
    for(int mu=0; mu < Nd; ++mu)
    {
      if( mu == t_dir) continue;

      for(int nu=mu+1; nu < Nd; ++nu)
      {
	if( nu == t_dir) continue;
	/* tmp_0 = u(x+mu,nu)*u_dag(x+nu,mu) */
	/* tmp_1 = tmp_0*u_dag(x,nu)=u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
	/* wplaq_tmp = tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)) */

	LatticeReal P_munu = real(trace(u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu])));
	LatticeReal other_term=param.coeff_plaq_s+factor*shift(P_munu, FORWARD, t_dir);
	
	S_pg += sum( other_term*P_munu );

      }
    }


    Double S_pg_t = zero;
    
    for(int nu = 0; nu < Nd; nu++) { 
      if ( nu == t_dir ) continue;
         
      S_pg_t += sum ( real(trace(u[t_dir]*shift(u[nu],FORWARD,t_dir)*adj(shift(u[t_dir],FORWARD,nu))*adj(u[nu]))) );

    }

    // Normalize -> 1 factor of Nc from each P_munu
    S_pg *= Double(-1)/Double(Nc);
    S_pg += Double(-param.coeff_plaq_t)/Double(Nc)*S_pg_t;

    END_CODE();

    return S_pg;
  } 

}

