// $Id: hyp_link_smearing.cc,v 1.3 2005-11-16 02:34:58 edwards Exp $
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

  // Read parameters
  void read(XMLReader& xml, const string& path, HypLinkSmearingEnv::Params& param)
  {
    HypLinkSmearingEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const HypLinkSmearingEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace HypLinkSmearingEnv
  {
    //! Callback function
    LinkSmearing* createSource(XMLReader& xml_in,
			       const std::string& path)
    {
      return new LinkSmear(Params(xml_in, path));
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


    //! Parameters for running code
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);
      num_smear = 1;
      no_smear_dir = -1;

      switch (version) 
      {
      case 2:
	break;

      case 3:
	read(paramtop, "num_smear", num_smear);
	break;

      case 4:
	read(paramtop, "num_smear", num_smear);
	read(paramtop, "no_smear_dir", no_smear_dir);
	break;

      default :
	QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "alpha1", alpha1);
      read(paramtop, "alpha2", alpha2);
      read(paramtop, "alpha3", alpha3);
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 4;
      write(xml, "version", version);
      write(xml, "LinkSmearingType", name);

      /* this version allows a variable num_smear */
      write(xml, "num_smear", num_smear);
      write(xml, "alpha1", alpha1);
      write(xml, "alpha2", alpha2);
      write(xml, "alpha3", alpha3);
      write(xml, "no_smear_dir", num_smear);

      pop(xml);
    }


    //! Smear the links
    void
    LinkSmear::operator()(multi1d<LatticeColorMatrix>& u) const
    {
      // Now hyp smear
      if (params.num_smear > 0)
      {
	QDPIO::cout << "Hyp Smear gauge field" << endl;

	int BlkMax = 100;
	Real BlkAccu = 1.0e-5;

	for (int n = 0; n < params.num_smear; n++)
	{
	  multi1d<LatticeColorMatrix> u_hyp(Nd);

	  if (params.no_smear_dir < 0 || params.no_smear_dir >= Nd)
	    Hyp_Smear(u, u_hyp, 
		      params.alpha1, params.alpha2, params.alpha3, 
		      BlkAccu, BlkMax);
	  else
	    Hyp_Smear3d(u, u_hyp, 
			params.alpha1, params.alpha2, params.alpha3, 
			BlkAccu, BlkMax, params.no_smear_dir);

	  u = u_hyp;
	}

	QDPIO::cout << "Gauge field Hyp-smeared!" << endl;
      }
    }

  }
}
