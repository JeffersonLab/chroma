// $Id: ape_link_smearing.cc,v 1.1 2005-10-28 21:31:04 edwards Exp $
/*! \file
 *  \brief APE link smearing
 */

#include "chromabase.h"

#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/ape_link_smearing.h"
#include "meas/smear/ape_smear.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace APELinkSmearingEnv
  {
    //! Callback function
    LinkSmearing* createSource(XMLReader& xml_in,
			       const std::string& path)
    {
      return new APELinkSmearing(APELinkSmearingParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "APE_SMEAR";

    //! Register all the factories
    bool registerAll()
    {
      return Chroma::TheLinkSmearingFactory::Instance().registerObject(name, createSource);
    }

    //! Register the source construction
    const bool registered = registerAll();
  }


  //! Parameters for running code
  APELinkSmearingParams::APELinkSmearingParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 1:
      break;

    default :
      QDPIO::cerr << "APELinkSmearingParams input version " << version << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "link_smear_num", param.link_smear_num);
    if( param.link_smear_num < 0 )
    {
      QDPIO::cerr << "apesmear: invalid number of ape smearing iterations, link_smear_num = " 
		  << param.link_smear_num << endl;
      QDP_abort(1);
    }

    read(paramtop, "link_smear_fact", param.link_smear_fact);
    read(paramtop, "no_smear_dir", param.no_smear_dir);
  }

  // Read parameters
  void read(XMLReader& xml, const string& path, APELinkSmearingParams& param)
  {
    APELinkSmearingParams tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const APELinkSmearingParams::Param_t& param)
  {
    push(xml, path);
    
    int version = 1;
    write(xml, "version", version);

    write(xml, "link_smear_num", param.link_smear_num);
    write(xml, "link_smear_fact", param.link_smear_fact);
    write(xml, "no_smear_dir", param.no_smear_dir);

    pop(xml);
  }


  //! Smear the links
  multi1d<LatticeColorMatrix>
  APELinkSmearing::operator()(const multi1d<LatticeColorMatrix>& u) const
  {
    // Now ape smear
    multi1d<LatticeColorMatrix> u_ape(Nd);
    u_ape = u;

    if (params.param.link_smear_num > 0)
    {
      QDPIO::cout << "APE Smear gauge field" << endl;

      int BlkMax = 100;
      Real BlkAccu = 1.0e-5;

      for(int i=0; i < params.param.link_smear_num; ++i)
      {
	multi1d<LatticeColorMatrix> u_tmp(Nd);

	for(int mu = 0; mu < Nd; ++mu)
	  if ( mu != params.param.no_smear_dir )
	    APE_Smear(u_ape, u_tmp[mu], mu, 0,
		      params.param.link_smear_fact, BlkAccu, BlkMax,
		      params.param.no_smear_dir);
	  else
	    u_tmp[mu] = u_ape[mu];

	u_ape = u_tmp;
      }
      QDPIO::cout << "Gauge field APE-smeared!" << endl;
    }

    return u_ape;
  }

}
