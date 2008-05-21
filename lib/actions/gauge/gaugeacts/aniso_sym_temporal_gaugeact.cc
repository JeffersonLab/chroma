// $Id: aniso_sym_temporal_gaugeact.cc,v 3.3 2008-05-21 17:07:50 bjoo Exp $
/*! \file
 *  \brief  Temporal part of Anisotropic Tree leve LW gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/aniso_sym_temporal_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

#include  "actions/gauge/gaugeacts/aniso_sym_shared_functions.h"
namespace Chroma
{
 
  namespace AnisoSymTemporalGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, 
		 multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
							       const std::string& path) 
    {
      return new AnisoSymTemporalGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
				       AnisoSymGaugeActParams(xml, path));
    }

    const std::string name = "ANISO_SYM_TEMPORAL_GAUGEACT";

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


  //! Compute action due to temporal \mu x 2\nu rectangles and plaquettes
  Double AnisoSymTemporalGaugeAct::S(const Handle< GaugeState<P,Q> >& state) const
  {
    START_CODE();

    const multi1d<LatticeColorMatrix>& u_bc = state->getLinks();

    LatticeReal lgimp = zero;
    
    int mu; 
    int nu;

    // length 1 in the t_dir (nu = t_dir) 
    nu = param.aniso.t_dir;
    for(mu=0; mu < Nd; mu++) { 
      if( mu != nu ) { 
	AnisoSym::S_part(mu, 
			 nu, 
			 param.aniso.t_dir,
			 plaq_c_t, 
			 rect_c_t_2, 
			 true,
			 lgimp,
			 u_bc);
      }
    }

    // length 2 in the t_dir (mu = t_dir) 
    // No rectangle contribution but still some plaquette
    // It is OK to conitnue to accumulate into lgimp
    // At this point the unwanted rectangles will be skipped
    mu = param.aniso.t_dir;
    for(int nu=0; nu < Nd; nu++) { 
      if( mu != nu ) { 
	AnisoSym::S_part(mu, 
			 nu, 
			 param.aniso.t_dir,
			 plaq_c_t, 
			 rect_c_t_2, 
			 true,
			 lgimp,
			 u_bc);
      }
    }


    // If user gave a zero point energy, subtract it off
    if( param.use_subtraction  ) { 
      LatticeReal ff;
      ff = param.sub_zero;
      lgimp -= ff;
    }

    Double ret_val = sum(lgimp);
    
    // Multiply in normalisation
    ret_val *= -Double(1)/Double(Nc);
    
    END_CODE();

    return ret_val;
  }

  //! Compute dS/dU
  void AnisoSymTemporalGaugeAct::deriv(multi1d<LatticeColorMatrix>& result,
				       const Handle< GaugeState<P,Q> >& state) const
  {
    result.resize(Nd);
    int mu;
    int nu;

    multi1d<LatticeColorMatrix> ds_tmp(Nd);


    const multi1d<LatticeColorMatrix>& u_bc = state->getLinks();
    
    for(mu=0; mu < Nd; mu++) { 
      result[mu] =zero;
      ds_tmp[mu]= zero;
    }

    // length 1 in the t_dir (nu = t_dir) 
    // Accumulate into ds_tmp
    nu = param.aniso.t_dir;
    for(mu=0; mu < Nd; mu++) { 
      if( mu != nu ) { 
	AnisoSym::deriv_part(mu, 
			     nu, 
			     param.aniso.t_dir,
			     plaq_c_t, 
			     rect_c_t_2, 
			     true,
			     ds_tmp,
			     u_bc);

      }
    }


    // length 2 in the t_dir (mu = t_dir)
    // Accumulate into ds_tmp

    mu = param.aniso.t_dir;
    for(int nu=0; nu < Nd; nu++) { 
      if( mu != nu ) { 
	AnisoSym::deriv_part(mu, 
			     nu, 
			     param.aniso.t_dir,
			     plaq_c_t, 
			     rect_c_t_2, 
			     true,
			     ds_tmp,
			     u_bc);

      }
    }     

    // Close up the loops
    for(int mu=0; mu < Nd; mu++) { 
      result[mu] = u_bc[mu]*ds_tmp[mu];
    }

    // Apply BCs
    getGaugeBC().zero(result);
      
    END_CODE();

  }

  // Private initializer
  void
  AnisoSymTemporalGaugeAct::init(void)
  {
    START_CODE();

    // Do the plaquette first. Spatial and temporal coeffs
    // anisotropy multiplied in in the terms constructor

    // Various tadpole things
    // spatial powers
    Real u_s_2 = param.u_s * param.u_s;
    Real u_s_4 = u_s_2 * u_s_2;

    // temporal powers
    Real u_t_2 = param.u_t * param.u_t;

    // Parameter needed for gauge act, but should never be used
    plaq_c_t = param.beta * Real(4) / ( Real(3) * u_s_2 * u_t_2 );
   
    // Loops that are short in the time direction
    // Param needed for rect_gaugeact driver, but never used
    rect_c_t_2 = - param.beta / ( Real(12)*u_s_4*u_t_2);

    // Fold in aniso factors
    if( param.aniso.anisoP ) { 
      plaq_c_t *= param.aniso.xi_0;
      rect_c_t_2 *= param.aniso.xi_0;
    }

    QDPIO::cout << "Real(Nc)*(u_s_2*u_t_2+u_s_4*u_t_2)="<<Real(Nc)*( Real(3)*u_s_2*u_t_2/Real(4)-Real(12)*u_s_4*u_t_2)/param.beta << endl;
    END_CODE();
  } 

}

