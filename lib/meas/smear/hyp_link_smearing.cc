/*! \file
 *  \brief Hyp link smearing
 */

#include "chromabase.h"

#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/hyp_link_smearing.h"
#include "util/gauge/hyp_utils.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const std::string& path, HypLinkSmearingEnv::Params& param)
  {
    HypLinkSmearingEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const std::string& path, const HypLinkSmearingEnv::Params& param)
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
    Params::Params(XMLReader& xml, const std::string& path)
    {
      XMLReader paramtop(xml, path);

      int version = 1;
      if (paramtop.count("version") > 0)
	read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	int no_smear_dir;
	read(paramtop, "no_smear_dir", no_smear_dir);
	smear_dirs.resize(Nd);
	smear_dirs = true;
	smear_dirs[no_smear_dir] = false;
	break;

      default :
	QDPIO::cerr << "Input parameter version " << version << " unsupported." << std::endl;
	QDP_abort(1);
      }

      read(paramtop, "link_smear_num", link_smear_num);
      read(paramtop, "link_smear_fact", link_smear_fact);

      read(paramtop, "alpha1", alpha1);
      read(paramtop, "alpha2", alpha2);
      read(paramtop, "alpha3", alpha3);
      read(paramtop, "BlkMax", BlkMax);
      read(paramtop, "BlkAccu", BlkAccu);
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const std::string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "LinkSmearingType", HypLinkSmearingEnv::name);

      // Structure
      write(xml, "link_smear_num", link_smear_num);
      write(xml, "link_smear_fact", link_smear_fact);
      write(xml, "smear_dirs", smear_dirs);

      // Parameters
      write(xml, "alpha1", alpha1);
      write(xml, "alpha2", alpha2);
      write(xml, "alpha3", alpha3);
      write(xml, "BlkMax", BlkMax);
      write(xml, "BlkAccu", BlkAccu);
      
      pop(xml);
    }


    //! Smear the links
    void
    LinkSmear::operator()(multi1d<LatticeColorMatrix>& u) const
    {
#if 0
      // Now hyp smear
      multi1d<LatticeColorMatrix> u_hyp = u;
      multi1d<LatticeColorMatrix> u_tmp(Nd);

      if (params.link_smear_num > 0)
      {
	QDPIO::cout << "Hyp Smear gauge field" << std::endl;

	for(int i=0; i < params.link_smear_num; ++i)
	{
          Hyping::smear_links(u_hyp, u_tmp,
                              params.smear_dirs,
                              params.alpha1, params.alpha2, params.alpha3, 
                              params.BlkMax, params.BlkAccu);
	  u_hyp = u_tmp;
	}
	QDPIO::cout << "Gauge field Hyp-smeared!" << std::endl;
      }
      
      u = u_hyp;
      
      // Now hyp smear
      if (params.num_smear > 0)
      {
	QDPIO::cout << "Hyp Smear gauge field" << std::endl;

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

	QDPIO::cout << "Gauge field Hyp-smeared!" << std::endl;
      }
#endif

    }
  }
}
