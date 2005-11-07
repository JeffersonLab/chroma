// $Id: stout_link_smearing.cc,v 1.1 2005-11-07 18:05:42 edwards Exp $
/*! \file
 *  \brief Stout link smearing
 */

#include "chromabase.h"

#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/stout_link_smearing.h"
#include "meas/smear/stout_smear.h"

namespace Chroma
{
  //! Hooks to register the class
  namespace StoutLinkSmearingEnv
  {
    //! Callback function
    LinkSmearing* createSource(XMLReader& xml_in,
			       const std::string& path)
    {
      return new StoutLinkSmearing(StoutLinkSmearingParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "STOUT_SMEAR";

    //! Register all the factories
    bool registerAll()
    {
      return Chroma::TheLinkSmearingFactory::Instance().registerObject(name, createSource);
    }

    //! Register the source construction
    const bool registered = registerAll();
  }


  //! Parameters for running code
  StoutLinkSmearingParams::StoutLinkSmearingParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 2:
      break;

    default :

      QDPIO::cerr << "Input version " << version << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "link_smear_num", param.link_smear_num);
    read(paramtop, "link_smear_fact", param.link_smear_fact);
    read(paramtop, "no_smear_dir", param.no_smear_dir);
  }

  // Read parameters
  void read(XMLReader& xml, const string& path, StoutLinkSmearingParams& param)
  {
    StoutLinkSmearingParams tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const StoutLinkSmearingParams::Param_t& param)
  {
    push(xml, path);
    
    int version = 2;
    write(xml, "version", version);
    write(xml, "link_smear_num", param.link_smear_num);
    write(xml, "link_smear_fact", param.link_smear_fact);
    write(xml, "no_smear_dir", param.no_smear_dir);

    pop(xml);
  }


  //! Smear the links
  void
  StoutLinkSmearing::operator()(multi1d<LatticeColorMatrix>& u) const
  {
    // Now stout smear
    multi1d<LatticeColorMatrix> u_stout(Nd);
    u_stout = u;

    if (params.param.link_smear_num > 0)
    {
      QDPIO::cout << "Stout Smear gauge field" << endl;

      int BlkMax = 100;
      Real BlkAccu = 1.0e-5;

      for(int i=0; i < params.param.link_smear_num; ++i)
      {
	multi1d<LatticeColorMatrix> u_tmp(Nd);

	for(int mu = 0; mu < Nd; ++mu)
	  if ( mu != params.param.no_smear_dir )
	    stout_smear(u_tmp[mu], u_stout, mu,
			params.param.link_smear_fact,
			params.param.no_smear_dir);
	  else
	    u_tmp[mu] = u_stout[mu];

	u_stout = u_tmp;
      }
      QDPIO::cout << "Gauge field Stout-smeared!" << endl;
    }
  }

}
