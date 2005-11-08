// $Id: hyp_link_smearing.cc,v 1.2 2005-11-08 05:33:19 edwards Exp $
/*! \file
 *  \brief Hyp link smearing
 */

#include "chromabase.h"

#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/hyp_link_smearing.h"
#include "meas/smear/hyp_smear.h"
#include "meas/smear/hyp_smear3d.h"

namespace Chroma
{
  //! Hooks to register the class
  namespace HypLinkSmearingEnv
  {
    //! Callback function
    LinkSmearing* createSource(XMLReader& xml_in,
			       const std::string& path)
    {
      return new HypLinkSmearing(HypLinkSmearingParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "HYP_SMEAR";

    //! Register all the factories
    bool registerAll()
    {
      return Chroma::TheLinkSmearingFactory::Instance().registerObject(name, createSource);
    }

    //! Register the source construction
    const bool registered = registerAll();
  }


  //! Parameters for running code
  HypLinkSmearingParams::HypLinkSmearingParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);
    param.num_smear = 1;
    param.no_smear_dir = -1;

    switch (version) 
    {
    case 2:
      break;

    case 3:
      read(paramtop, "num_smear", param.num_smear);
      break;

    case 4:
      read(paramtop, "num_smear", param.num_smear);
      read(paramtop, "no_smear_dir", param.no_smear_dir);
      break;

    default :
      QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "alpha1", param.alpha1);
    read(paramtop, "alpha2", param.alpha2);
    read(paramtop, "alpha3", param.alpha3);
  }

  // Read parameters
  void read(XMLReader& xml, const string& path, HypLinkSmearingParams& param)
  {
    HypLinkSmearingParams tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const HypLinkSmearingParams::Param_t& param)
  {
    push(xml, path);

    int version = 4;
    write(xml, "version", version);
    write(xml, "LinkSmearingType", HypLinkSmearingEnv::name);

    /* this version allows a variable num_smear */
    write(xml, "num_smear", param.num_smear);
    write(xml, "alpha1", param.alpha1);
    write(xml, "alpha2", param.alpha2);
    write(xml, "alpha3", param.alpha3);
    write(xml, "no_smear_dir", param.num_smear);

    pop(xml);
  }


  //! Smear the links
  void
  HypLinkSmearing::operator()(multi1d<LatticeColorMatrix>& u) const
  {
    // Now hyp smear
    if (params.param.num_smear > 0)
    {
      QDPIO::cout << "Hyp Smear gauge field" << endl;

      int BlkMax = 100;
      Real BlkAccu = 1.0e-5;

      for (int n = 0; n < params.param.num_smear; n++)
      {
	multi1d<LatticeColorMatrix> u_hyp(Nd);

	if (params.param.no_smear_dir < 0 || params.param.no_smear_dir >= Nd)
	  Hyp_Smear(u, u_hyp, 
		    params.param.alpha1, params.param.alpha2, params.param.alpha3, 
		    BlkAccu, BlkMax);
	else
	  Hyp_Smear3d(u, u_hyp, 
		      params.param.alpha1, params.param.alpha2, params.param.alpha3, 
		      BlkAccu, BlkMax, params.param.no_smear_dir);

	u = u_hyp;
      }

      QDPIO::cout << "Gauge field Hyp-smeared!" << endl;
    }
  }

}
