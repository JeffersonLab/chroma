// $Id: plaq_gaugeact.cc,v 3.10 2008-01-20 19:47:24 edwards Exp $
/*! \file
 *  \brief Plaquette gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"
#include "meas/glue/mesplq.h"


namespace Chroma
{
 
  namespace PlaqGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, 
		 multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
							       const std::string& path) 
    {
      return new PlaqGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
			      PlaqGaugeActParams(xml, path));
    }

    const std::string name = "PLAQ_GAUGEACT";

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


  PlaqGaugeActParams::PlaqGaugeActParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);

    try 
    {
      // Want a multi2d read here, but this is not supported currently in chroma. 
      // I'll go ahead and hack by reading Nd rows or columns. Yuk, because I
      // doubt anybody will use this part of the code.
      coeffs.resize(Nd,Nd);

      XMLReader coefftop(paramtop, "coeffs");
      for(int i=0; i < Nd; ++i)
      {
	std::ostringstream os;
	os << "elem[" << i+1 << "]";

	multi1d<Real> cc;
	read(coefftop, os.str(), cc);
	coeffs[i] = cc;
      }
      
//      read(paramtop, "coeffs", coeffs);
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


// A multi2d xml write doesn't exist at this moment.
//  void write(XMLWriter& xml, const string& path, const PlaqGaugeActParams& p) 
//  {
//    push(xml, path);
//    write(xml, "coeffs", p.coeffs);
//    pop(xml);
//  }


  // Internal initializer
  void PlaqGaugeAct::init(const Real& coeff, const AnisoParam_t& aniso)
  {
    START_CODE();

    param.coeffs.resize(Nd,Nd);
    param.coeffs = zero;

    int t_dir;
    Real coeff_t = coeff;
    Real coeff_s = coeff;

    if (aniso.anisoP)
    {
      t_dir = aniso.t_dir;
      coeff_t *= aniso.xi_0;
      coeff_s /= aniso.xi_0;
    }
    else
    {
      t_dir = -1;
    }

    for(int mu = 0; mu < Nd; ++mu)
    {
      for(int nu = mu+1; nu < Nd; ++nu) 
      { 
	if( mu == t_dir || nu == t_dir ) 
	{
	  // Temporal Plaquette in either mu or nu direction
	  param.coeffs[mu][nu] = coeff_t;
	}
	else {
	  // Spatial Plaquette
	  param.coeffs[mu][nu] = coeff_s;
	}

	param.coeffs[nu][mu] = param.coeffs[mu][nu];
      }
    }
    
    END_CODE();
  }


  // Internal initializer
  void PlaqGaugeAct::init(const Real& coeff_s, 
			  const Real& coeff_t, 
			  const AnisoParam_t& aniso)
  {
    START_CODE();

    param.coeffs.resize(Nd,Nd);
    param.coeffs = zero;

    int t_dir = (aniso.anisoP) ? aniso.t_dir : -1;

    for(int mu = 0; mu < Nd; ++mu)
    {
      for(int nu = mu+1; nu < Nd; ++nu) 
      { 
	if( mu == t_dir || nu == t_dir ) 
	{
	  // Temporal Plaquette in either mu or nu direction
	  param.coeffs[mu][nu] = coeff_t;
	  if( aniso.anisoP ) {
	    param.coeffs[mu][nu] *= aniso.xi_0;
	  }

	}
	else {
	  // Spatial Plaquette
	  param.coeffs[mu][nu] = coeff_s;
	  if( aniso.anisoP ) { 
	    param.coeffs[mu][nu] /= aniso.xi_0;
	  }
	}

	param.coeffs[nu][mu] = param.coeffs[mu][nu];
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
  PlaqGaugeAct::staple(LatticeColorMatrix& u_mu_staple,
		       const Handle< GaugeState<P,Q> >& state,
		       int mu, int cb) const
  {
    START_CODE();

    const multi1d<LatticeColorMatrix>& u = state->getLinks();
    
    u_mu_staple = zero;
    LatticeColorMatrix tmp1, tmp2;
    LatticeColorMatrix u_nu_mu;

    for(int nu=0; nu < Nd; ++nu) 
    {
      if( nu == mu ) continue;
      
      u_nu_mu = shift(u[nu],FORWARD,mu);

      // +forward staple
      tmp1[rb[cb]] = u_nu_mu * adj(shift(u[mu],FORWARD,nu));
      tmp2[rb[cb]] = tmp1 * adj(u[nu]);

      u_mu_staple[rb[cb]] += param.coeffs[mu][nu] * tmp2;

      // +backward staple
      tmp1[rb[cb]] = adj(shift(u_nu_mu,BACKWARD,nu)) * adj(shift(u[mu],BACKWARD,nu));
      tmp2[rb[cb]] = tmp1 * shift(u[nu],BACKWARD,nu);

      u_mu_staple[rb[cb]] += param.coeffs[mu][nu] * tmp2;
    }

    // NOTE: a heatbath code should be responsible for resetting links on
    // a boundary. The staple is not really the correct place.

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
  PlaqGaugeAct::stapleSpatial(LatticeColorMatrix& u_mu_staple,
			      const Handle< GaugeState<P,Q> >& state,
			      int mu, int cb, int t_dir) const
  {
    START_CODE();

    u_mu_staple = zero;
    if( mu == t_dir ) return; // Short circuit. 

    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    for(int nu=0; nu<Nd; nu++) 
    {
      if( (nu != mu) && (nu != t_dir) ) 
      {
	// nu is spatial
	// anisotropic lattice, time-like staple
	// +forward staple
	u_mu_staple[rb[cb]]+=
	  shift(u[nu],FORWARD,mu)*
	  shift(adj(u[mu]),FORWARD,nu)*adj(u[nu])*
	  param.coeffs[mu][nu];
	// +backward staple
	u_mu_staple[rb[cb]]+=
	  shift(shift(adj(u[nu]),FORWARD,mu),BACKWARD,nu)*
	  shift(adj(u[mu]),BACKWARD,nu)*
	  shift(u[nu],BACKWARD,nu)*
	  param.coeffs[mu][nu];
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
  PlaqGaugeAct::stapleTemporal(LatticeColorMatrix& u_mu_staple,
			       const Handle< GaugeState<P,Q> >& state,
			       int mu, int cb, int t_dir) const
  {
    START_CODE();

    // This bit of code was taken from  chroma/lib/update/heatbath/u_staple.cc
    // Supposedly it works.

    const multi1d<LatticeColorMatrix>& u = state->getLinks();
    u_mu_staple = zero;
    
    for(int nu=0; nu<Nd; nu++) 
    {
      bool doit_mu_is_time;
      bool doit_nu_is_time;

      doit_mu_is_time = ( mu == t_dir ) && (nu != t_dir); // true for ts plaquettes
      doit_nu_is_time = ( nu == t_dir ) && (mu != t_dir); // true for st plaquettes

      if( doit_mu_is_time || doit_nu_is_time ) 
      {
	// anisotropic lattice, time-like staple
	// +forward staple
	u_mu_staple[rb[cb]]+=
	  shift(u[nu],FORWARD,mu)*
	  shift(adj(u[mu]),FORWARD,nu)*adj(u[nu])*
	  param.coeffs[mu][nu];
	// +backward staple
	u_mu_staple[rb[cb]]+=
	  shift(shift(adj(u[nu]),FORWARD,mu),BACKWARD,nu)*
	  shift(adj(u[mu]),BACKWARD,nu)*
	  shift(u[nu],BACKWARD,nu)*
	  param.coeffs[mu][nu];
      }
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

    for(int mu = 0; mu < Nd; mu++)
    {
      for(int nu=mu+1; nu < Nd; nu++) 
      {
	for(int cb=0; cb < 2; cb++) 
	{ 
	  tmp_0[rb[cb]] = shift(u[mu], FORWARD, nu)*shift(adj(u[nu]), FORWARD, mu);
	  tmp_0[rb[cb]] *= param.coeffs[mu][nu];   // c[mu][nu] = c[nu][mu]
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


  //! compute spatial dS/dU given a time direction
  void 
  PlaqGaugeAct::derivSpatial(multi1d<LatticeColorMatrix>& ds_u,
			     const Handle< GaugeState<P,Q> >& state,
			     int t_dir) const
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
	// Pick out spatial plaquettes only
	if( (mu != t_dir) && (nu != t_dir) ) 
	{
	  for(int cb=0; cb < 2; cb++) 
	  { 
	    tmp_0[rb[cb]] = shift(u[mu], FORWARD, nu)*shift(adj(u[nu]), FORWARD, mu);
	    tmp_0[rb[cb]] *= param.coeffs[mu][nu];   // c[mu][nu] = c[nu][mu]
	    tmp_1[rb[cb]] = tmp_0*adj(u[mu]);
	    tmp_2[rb[cb]] = u[nu]*tmp_1;
	    ds_u[nu][rb[cb]] += tmp_2;
	    ds_u[mu][rb[cb]] += adj(tmp_2);
	    ds_u[mu][rb[1-cb]] += shift(tmp_1, BACKWARD, nu)*shift(u[nu], BACKWARD, nu);
	    tmp_1[rb[cb]] = adj(u[nu])*u[mu];
	    ds_u[nu][rb[1-cb]] += shift(adj(tmp_0),BACKWARD,mu)*shift(tmp_1, BACKWARD, mu);
	  }
	}
      }
      
      // It is 1/(4Nc) to account for normalisation relevant to fermions
      // in the taproj, which is a factor of 2 different from the 
      // one used here.
      ds_u[mu] *= Real(-1)/(Real(2*Nc));

    }

    // Zero the force on any fixed boundaries
    getGaugeBC().zero(ds_u);

    END_CODE();

  }
  
  //! compute spatial dS/dU given a time direction
  void 
  PlaqGaugeAct::derivTemporal(multi1d<LatticeColorMatrix>& ds_u,
			      const Handle< GaugeState<P,Q> >& state,
			      int t_dir) const
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
	// Pick out temporal plaquettes only
	if( (mu == t_dir) || (nu == t_dir) ) 
	{
	  for(int cb=0; cb < 2; cb++) 
	  {
	    tmp_0[rb[cb]] = shift(u[mu], FORWARD, nu)*shift(adj(u[nu]), FORWARD, mu);
	    tmp_0[rb[cb]] *= param.coeffs[mu][nu];   // c[mu][nu] = c[nu][mu]
	    tmp_1[rb[cb]] = tmp_0*adj(u[mu]);
	    tmp_2[rb[cb]] = u[nu]*tmp_1;
	    ds_u[nu][rb[cb]] += tmp_2;
	    ds_u[mu][rb[cb]] += adj(tmp_2);
	    ds_u[mu][rb[1-cb]] += shift(tmp_1, BACKWARD, nu)*shift(u[nu], BACKWARD, nu);
	    tmp_1[rb[cb]] = adj(u[nu])*u[mu];
	    ds_u[nu][rb[1-cb]] += shift(adj(tmp_0),BACKWARD,mu)*shift(tmp_1, BACKWARD, mu);
	  }
	}
      }
      
      // It is 1/(4Nc) to account for normalisation relevant to fermions
      // in the taproj, which is a factor of 2 different from the 
      // one used here.
      ds_u[mu] *= Real(-1)/(Real(2*Nc));
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
  PlaqGaugeAct::S(const Handle< GaugeState<P,Q> >& state) const
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

	S_pg += tmp * Double(param.coeffs[mu][nu]);
      }
    }

    // Normalize
    S_pg *= Double(-1)/Double(Nc);
    
    END_CODE();

    return S_pg;
  } 


  //! Compute the spatial part of the action given a time direction
  Double PlaqGaugeAct::spatialS(const Handle< GaugeState<P,Q> >& state, int t_dir) const
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
	if( (nu != t_dir) && (mu != t_dir) ) 
	{
	  Double tmp = 
	    sum(real(trace(u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu]))));
	  
	  S_pg += tmp * Double(param.coeffs[mu][nu]);
	}
      }
    }

    // Normalize
    S_pg *= Double(-1)/Double(Nc);
    
    END_CODE();

    return S_pg;

  }

    //! Compute the temporal part of the action given a time direction
  Double PlaqGaugeAct::temporalS(const Handle< GaugeState<P,Q> >& state, int t_dir) const
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
	if( (nu == t_dir) || (mu == t_dir) ) 
	{
	  Double tmp = 
	    sum(real(trace(u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu]))));
	  
	  S_pg += tmp * Double(param.coeffs[mu][nu]);
	}
      }
    }

    // Normalize
    S_pg *= Double(-1)/Double(Nc);
    
    END_CODE();

    return S_pg;

  }
  
}

