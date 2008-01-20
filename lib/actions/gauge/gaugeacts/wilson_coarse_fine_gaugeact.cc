// $Id: wilson_coarse_fine_gaugeact.cc,v 3.1 2008-01-20 19:46:53 edwards Exp $
/*! \file
 *  \brief Wilson gauge action supporting 2+2 anisotropy
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/wilson_coarse_fine_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma
{
 
  namespace WilsonCoarseFineGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, 
		 multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
							       const std::string& path) 
    {
      return new WilsonCoarseFineGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
					  WilsonCoarseFineGaugeActParams(xml, path));
    }

    const std::string name = "WILSON_COARSE_FINE_GAUGEACT";

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


  WilsonCoarseFineGaugeActParams::WilsonCoarseFineGaugeActParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);

    try {
      read(paramtop, "beta", beta);
      read(paramtop, "coarse_dirs", coarse_dirs);
      read(paramtop, "xi", xi);
      read(paramtop, "coeff_ff", coeff_ff);
      read(paramtop, "coeff_cc", coeff_cc);
    }
    catch( const std::string& e ) 
    { 
      QDPIO::cerr << __func__ << ": Error reading XML: " <<  e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, const string& path, WilsonCoarseFineGaugeActParams& p) 
  {
    WilsonCoarseFineGaugeActParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const string& path, const WilsonCoarseFineGaugeActParams& p) 
  {
    push(xml, path);

    int version = 1;

    write(xml, "version", version);
    write(xml, "beta", p.beta);
    write(xml, "coarse_dirs", p.coarse_dirs);
    write(xml, "xi", p.xi);
    write(xml, "coeff_ff", p.coeff_ff);
    write(xml, "coeff_cc", p.coeff_cc);

    pop(xml);
  }

  

  //! Read beta from a param struct
  WilsonCoarseFineGaugeAct::WilsonCoarseFineGaugeAct(Handle< CreateGaugeState<P,Q> > cgs,
						     const WilsonCoarseFineGaugeActParams& p) 
    : param(p) 
  {
    init(cgs);
  }


  // Private initializer
  void
  WilsonCoarseFineGaugeAct::init(Handle< CreateGaugeState<P,Q> > cgs)
  {
    START_CODE();

    // Sanity check. Insist there are 2 coarse and 2 fine directions
    if (param.coarse_dirs.size() != 4 || Nd != 4)
    {
      QDPIO::cerr << WilsonCoarseFineGaugeActEnv::name << ": coarse_dirs and Nd must be size=4" << endl;
      QDP_abort(1);
    }

    int cnt = 0;
    for(int i=0; i < param.coarse_dirs.size(); ++i)
    {
      if (param.coarse_dirs[i]) ++cnt;
    }

    if (cnt != 2)
    {
      QDPIO::cerr << WilsonCoarseFineGaugeActEnv::name << ": not 2 coarse dirs" << endl;
      QDP_abort(1);
    }


    // Fold in normalizations
    multi2d<Real> coeffs(Nd,Nd);
    coeffs = zero;

    Real g2   = 2.0*Nc / param.beta;
    Real xi2  = param.xi * param.xi;
    Real c_ff = param.beta * (1.0 + param.coeff_ff*g2) * xi2;
    Real c_cc = param.beta * (1.0 + param.coeff_cc*g2) / xi2;
    Real c_cf = param.beta;

    for(int mu = 0; mu < Nd; ++mu)
    {
      for(int nu = mu+1; nu < Nd; ++nu) 
      { 
	if( (! param.coarse_dirs[mu]) && (! param.coarse_dirs[nu]) )
	{
	  // Fine-fine Plaquette in either mu or nu direction
	  coeffs[mu][nu] = c_ff;
	}
	else if( param.coarse_dirs[mu] && param.coarse_dirs[nu] ) 
	{
	  // Coarse-coarse Plaquette in either mu or nu direction
	  coeffs[mu][nu] = c_cc;
	}
	else 
	{
	  // Coarse-fine or fine-coarse Plaquette
	  coeffs[mu][nu] = c_cf;
	}

	coeffs[nu][mu] = coeffs[mu][nu];
      }
    }
    
    // Create gauge action
    PlaqGaugeActParams plaq_param;
    plaq_param.coeffs = coeffs;
    plaq = new PlaqGaugeAct(cgs,plaq_param);
    
    END_CODE();
  } 

}

