// $Id: hyp_link_smearing.cc,v 3.3 2008-11-04 18:43:57 edwards Exp $
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
    namespace
    {
      //! Callback function
      LinkSmearing* createSource(XMLReader& xml_in,
				 const std::string& path)
      {
	return new LinkSmear(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name = "HYP_SMEAR";
    }

    //! Return the name
    std::string getName() {return name;}

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheLinkSmearingFactory::Instance().registerObject(name, createSource);
	registered = true;
      }
      return success;
    }


    //! Parameters for running code
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);
      num_smear = 1;
      no_smear_dir = -1;
      BlkMax = 100;
      BlkAccu = 1.0e-5;

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

      case 5:
	read(paramtop, "num_smear", num_smear);
	read(paramtop, "no_smear_dir", no_smear_dir);
	read(paramtop, "BlkMax", BlkMax);
	read(paramtop, "BlkAccu", BlkAccu);
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

      int version = 5;
      write(xml, "version", version);
      write(xml, "LinkSmearingType", name);

      /* this version allows a variable num_smear */
      write(xml, "alpha1", alpha1);
      write(xml, "alpha2", alpha2);
      write(xml, "alpha3", alpha3);
      write(xml, "num_smear", num_smear);
      write(xml, "no_smear_dir", num_smear);
      write(xml, "BlkMax", BlkMax);
      write(xml, "BlkAccu", BlkAccu);

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

	for (int n = 0; n < params.num_smear; n++)
	{
	  multi1d<LatticeColorMatrix> u_hyp(Nd);

	  if (params.no_smear_dir < 0 || params.no_smear_dir >= Nd)
	    Hyp_Smear(u, u_hyp, 
		      params.alpha1, params.alpha2, params.alpha3, 
		      params.BlkAccu, params.BlkMax);
	  else
	    Hyp_Smear3d(u, u_hyp, 
			params.alpha1, params.alpha2, params.alpha3, 
			params.BlkAccu, params.BlkMax, params.no_smear_dir);

	  u = u_hyp;
	}

	QDPIO::cout << "Gauge field Hyp-smeared!" << endl;
      }
    }

  }
}
