// $Id: wilson_gaugeact.cc,v 1.10 2005-01-04 14:57:25 bjoo Exp $
/*! \file
 *  \brief Wilson gauge action
 */

#include "chromabase.h"
#include "gaugebc.h"

#include "gaugeact_factory.h"
#include "gaugebc_factory.h"

#include "actions/gauge/wilson_gaugeact.h"
#include "meas/glue/mesplq.h"


using namespace std;
using namespace QDP;
using namespace Chroma;

namespace Chroma
{
 
  namespace WilsonGaugeActEnv { 
    GaugeAction* createGaugeAct(XMLReader& xml, const std::string& path) {
      return new WilsonGaugeAct(WilsonGaugeActParams(xml, path));
    }

    const std::string name = "WILSON_GAUGEACT";
    const bool registered = TheGaugeActFactory::Instance().registerObject(name, 
									  createGaugeAct);
  };


  WilsonGaugeActParams::WilsonGaugeActParams(XMLReader& xml_in, const std::string& path) {
    XMLReader paramtop(xml_in, path);

    try {
      read(paramtop, "./beta", beta);

      
      XMLReader bc_xml(paramtop, "./GaugeBC");
      read(bc_xml, "./Name", bc_xml_name);

      std::ostringstream os;
      bc_xml.print(os);
      bc_xml_string = os.str();
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " <<  e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, const string& path, WilsonGaugeActParams& p) {
    WilsonGaugeActParams tmp(xml, path);
    p=tmp;
  }

  WilsonGaugeAct::WilsonGaugeAct(const WilsonGaugeActParams& p) {
    beta = p.beta;
    std::string tmpstring = p.bc_xml_string;
    std::istringstream bc_is(tmpstring);
    XMLReader bc_xml_tmp(bc_is);
    
    gbc = TheGaugeBCFactory::Instance().createObject(p.bc_xml_name, bc_xml_tmp , "/GaugeBC");
  }

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
    QDPIO::cout << "WilsonGaugeAct::staple() --- Is this tested ? BJ" << endl;
    QDPIO::cout << "Use at own risk" << endl;

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
  WilsonGaugeAct::dsdu(multi1d<LatticeColorMatrix>& ds_u,
		       const Handle< const ConnectState> state) const
  {
    START_CODE();

    ds_u.resize(Nd);

    LatticeColorMatrix tmp_0;
    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;

    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    for(int mu = 0; mu < Nd; mu++) {
      ds_u[mu] = zero;
    }

    for(int mu = 0; mu < Nd; ++mu)
    {
      for(int nu=mu+1; nu<Nd; nu++) {
	for(int cb=0; cb < 2; cb++) { 
	  
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
      
      ds_u[mu] *= Real(1)*Real(beta)/(Real(4*Nc));
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

    // Pure Gauge factor (-beta/Nc and a factor of 2 because of the forward
    // and backward staple of the force)
    for(int mu=0; mu < Nd; mu++) { 
      ds_u[mu] *= Real(beta)/(Real(2*Nc));
    }
#endif

    END_CODE();
  }

  // Get the gauge action
  //
  // S = -(beta/(Nc) Sum Re Tr Plaq
  //
  // w_plaq is defined in MesPlq as
  //
  // w_plaq =( 2/(V*Nd*(Nd-1)*Nc)) * Sum Re Tr Plaq
  //
  // so 
  // S = -beta * (V*Nd*(Nd-1)/2) w_plaq 
  //   = -beta * (V*Nd*(Nd-1)/2)*(2/(V*Nd*(Nd-1)*Nc))* Sum Re Tr Plaq
  //   = -beta * (1/(Nc)) * Sum Re Tr Plaq

  Double
  WilsonGaugeAct::S(const Handle<const ConnectState> state) const
  {
    Double S_pg;

    // Handle< const ConnectState> u_bc(createState(u));
    // Apply boundaries
    Double w_plaq, s_plaq, t_plaq, link;

    // Measure the plaquette
    MesPlq(state->getLinks(), w_plaq, s_plaq, t_plaq, link);

    // Undo Mes Plaq Normalisation
    S_pg = Double(Layout::vol()*Nd*(Nd-1)*Nc)/Double(2);

    S_pg *= Double(-1)*Double(beta)/Double(2*Nc);

    // Took out minus sign -- may need to put back in...
    S_pg *= w_plaq;
    
    return S_pg;
  } 

}

