// $Id: ape_link_smearing.cc,v 3.3 2008-11-04 18:43:57 edwards Exp $
/*! \file
 *  \brief APE link smearing
 */

#include "chromabase.h"

#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/ape_link_smearing.h"
#include "meas/smear/ape_smear.h"

namespace Chroma
{
  // Read parameters
  void read(XMLReader& xml, const string& path, APELinkSmearingEnv::Params& param)
  {
    APELinkSmearingEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const APELinkSmearingEnv::Params& param)
  {
    param.writeXML(xml, path);
  }



  //! Hooks to register the class
  namespace APELinkSmearingEnv
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
      const std::string name = "APE_SMEAR";
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

      int version = 1;
      if (paramtop.count("version") != 0)
	read(paramtop, "version", version);

      BlkMax = 100;
      BlkAccu = 1.0e-5;

      switch (version) 
      {
      case 1:
	break;

      case 2:
	read(paramtop, "BlkMax", BlkMax);
	read(paramtop, "BlkAccu", BlkAccu);
	break;

      default :
	QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "link_smear_num", link_smear_num);
      read(paramtop, "link_smear_fact", link_smear_fact);
      read(paramtop, "no_smear_dir", no_smear_dir);
    }

    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      int version = 2;
      write(xml, "version", version);
      write(xml, "LinkSmearingType", name);
      write(xml, "link_smear_num", link_smear_num);
      write(xml, "link_smear_fact", link_smear_fact);
      write(xml, "no_smear_dir", no_smear_dir);
      write(xml, "BlkMax", BlkMax);
      write(xml, "BlkAccu", BlkAccu);

      pop(xml);
    }


    //! Smear the links
    void
    LinkSmear::operator()(multi1d<LatticeColorMatrix>& u) const
    {
      // Now ape smear
      multi1d<LatticeColorMatrix> u_ape =  u;

      if (params.link_smear_num > 0)
      {
	QDPIO::cout << "APE Smear gauge field" << endl;

	for(int i=0; i < params.link_smear_num; ++i)
	{
	  multi1d<LatticeColorMatrix> u_tmp(Nd);

	  for(int mu = 0; mu < Nd; ++mu)
	    if ( mu != params.no_smear_dir )
	      APE_Smear(u_ape, u_tmp[mu], mu, 0,
			params.link_smear_fact, params.BlkAccu, params.BlkMax,
			params.no_smear_dir);
	    else
	      u_tmp[mu] = u_ape[mu];

	  u_ape = u_tmp;
	}
	QDPIO::cout << "Gauge field APE-smeared!" << endl;
      }

      u = u_ape;
    }

  }
}
