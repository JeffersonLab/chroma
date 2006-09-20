// $Id: spatial_two_plaq_gaugeact.cc,v 1.6 2006-09-20 20:28:00 edwards Exp $
/*! \file
 *  \brief Plaquette gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/spatial_two_plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"
#include "meas/glue/mesplq.h"
#include "util/gauge/taproj.h"

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
    START_CODE();
 
    // THis term is only spatial. Really I just need to divide in the 
    // aniso factors
    if ( anisoP() ) { 
      param.coeff /= param.aniso.xi_0;
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
  SpatialTwoPlaqGaugeAct::deriv(multi1d<LatticeColorMatrix>& ds_u,
		      const Handle< GaugeState<P,Q> >& state) const
  {
    START_CODE();

    ds_u.resize(Nd);


    const multi1d<LatticeColorMatrix>& u = state->getLinks();


    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    LatticeColorMatrix tmp, tmp2;

    int t_dir = tDir();


    ds_tmp = zero;

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
	  LatticeColorMatrix u_mu_plus_nu = shift(u[mu],FORWARD, nu);
	  LatticeColorMatrix u_nu_plus_mu = shift(u[nu],FORWARD, mu);
	  // First  lets do 
	  //    <-----
	  //   |
	  //   |           (we'll use this for (1) and (4))
	  //   V
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

	  LatticeReal P_sum = shift(P_munu, FORWARD, t_dir)+
	    shift(P_munu, BACKWARD, t_dir);

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


    // Finally multiply by u[mu] in the spatial dirs 
    // and zero out junk in the time dir
    for(int mu=0; mu < Nd; mu++) {
      if (mu == t_dir) {
	ds_u[mu] = zero; // This component is otherwise untouched
	                 // So we should zero junk out of it
      }
      else { 
	ds_u[mu] = u[mu]*ds_tmp[mu];
	ds_u[mu] *= Real(-param.coeff)/Real(4*Nc*Nc);
      }
    }
 
    // Zero the force on any fixed boundaries
    getGaugeBC().zero(ds_u);
    
    END_CODE();
  }

  // Get the gauge action
  //
  // S = -(coeff/2) Sum P_{ij}(x) P_{ij}(x+t)
  //
  // where P_{ij}(x) is the plaquette on lattice point x and i,j are 
  // only spatial directions.
  //
  // P has normalisatiom (1/Nc) and so the overall normalisation
  // for the product is (1/Nc^2). Giving the final multiplier as
  //
  //  -coeff/(2 Nc^2 )
  Double
  SpatialTwoPlaqGaugeAct::S(const Handle< GaugeState<P,Q> >& state) const
  {
    START_CODE();

    Double S_pg = zero;

    // Handle< const GaugeState<P,Q> > u_bc(createState(u));
    // Apply boundaries
    const multi1d<LatticeColorMatrix>& u = state->getLinks();
    int t_dir = tDir();

    
    // Compute the average plaquettes
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

	S_pg += sum(P_munu*shift(P_munu, FORWARD, t_dir));

      }
    }

    // Normalize -> 1 factor of Nc from each P_munu
    S_pg *= Double(-param.coeff)/Double(2*Nc*Nc);

    END_CODE();

    return S_pg;
  } 

}

