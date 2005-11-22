// $Id: stout_link_smearing.cc,v 1.6 2005-11-22 18:55:35 edwards Exp $
/*! \file
 *  \brief Stout link smearing
 */

#include "chromabase.h"

#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/stout_link_smearing.h"
#include "meas/smear/stout_smear.h"

#include "meas/glue/mesplq.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const string& path, StoutLinkSmearingEnv::Params& param)
  {
    StoutLinkSmearingEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const StoutLinkSmearingEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace StoutLinkSmearingEnv
  {
    //! Callback function
    LinkSmearing* createSource(XMLReader& xml_in,
			       const std::string& path)
    {
      return new LinkSmear(Params(xml_in, path));
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


    //! Parameters for running code
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

//      int version;
//      read(paramtop, "version", version);
//      switch (version) 
//      {
//      case 2:
//	break;
//      default :
//	QDPIO::cerr << "Input version " << version << " unsupported." << endl;
//	QDP_abort(1);
//      }

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
      write(xml, "LinkSmearingType", StoutLinkSmearingEnv::name);
      write(xml, "link_smear_num", link_smear_num);
      write(xml, "link_smear_fact", link_smear_fact);
      write(xml, "no_smear_dir", no_smear_dir);

      pop(xml);
    }


    //! Smear the links
    void
    LinkSmear::operator()(multi1d<LatticeColorMatrix>& u) const
    {
      // Now stout smear
      multi1d<LatticeColorMatrix> u_stout(Nd);
      u_stout = u;

      if (params.link_smear_num > 0)
      {
	QDPIO::cout << "Stout Smear gauge field" << endl;

	int BlkMax = 100;
	Real BlkAccu = 1.0e-5;

	for(int i=0; i < params.link_smear_num; ++i)
	{
	  multi1d<LatticeColorMatrix> u_tmp(Nd);

	  for(int mu = 0; mu < Nd; ++mu)
	    if ( mu != params.no_smear_dir )
	      stout_smear(u_tmp[mu], u_stout, mu,
			  params.link_smear_fact,
			  params.no_smear_dir);
	    else
	      u_tmp[mu] = u_stout[mu];

	  u_stout = u_tmp;
	}
	QDPIO::cout << "Gauge field Stout-smeared!" << endl;
      }

      u = u_stout;
    }
  }
}
