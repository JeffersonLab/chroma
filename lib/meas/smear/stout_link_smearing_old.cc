// $Id: stout_link_smearing_old.cc,v 3.1 2008-01-20 03:07:51 edwards Exp $
/*! \file
 *  \brief Stout link smearing using the old (non-gauge covariant) method
 */

#include "chromabase.h"

#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/stout_link_smearing_old.h"
#include "meas/smear/stout_smear_old.h"

#include "meas/glue/mesplq.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const string& path, StoutLinkSmearingOldEnv::Params& param)
  {
    StoutLinkSmearingOldEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const StoutLinkSmearingOldEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace StoutLinkSmearingOldEnv
  {
    //! Callback function
    LinkSmearing* createSource(XMLReader& xml_in,
			       const std::string& path)
    {
      return new LinkSmear(Params(xml_in, path));
    }

    //! Name to be used
    const std::string name = "STOUT_SMEAR_OLD";

    //! Local registration flag
    static bool registered = false;

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

      int version = 2;
      if (paramtop.count("version") > 0)
	read(paramtop, "version", version);

      switch (version) 
      {
      case 2:
      {
	int no_smear_dir;
	read(paramtop, "no_smear_dir", no_smear_dir);
	smear_dirs.resize(Nd);
	smear_dirs = true;
	smear_dirs[no_smear_dir] = false;
      }
      break;

      case 3:
	read(paramtop, "smear_dirs", smear_dirs);
	break;

      default:
	QDPIO::cerr << StoutLinkSmearingOldEnv::name << ": Input version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "link_smear_num", link_smear_num);
      read(paramtop, "link_smear_fact", link_smear_fact);
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      int version = 3;
      write(xml, "version", version);
      write(xml, "LinkSmearingType", StoutLinkSmearingOldEnv::name);
      write(xml, "link_smear_num", link_smear_num);
      write(xml, "link_smear_fact", link_smear_fact);
      write(xml, "smear_dirs", smear_dirs);

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
	QDPIO::cout << "Stout Smear gauge field using old (non gauge covarient) method" << endl;

	for(int i=0; i < params.link_smear_num; ++i)
	{
	  multi1d<LatticeColorMatrix> u_tmp(Nd);

	  for(int mu = 0; mu < Nd; ++mu)
	    if ( params.smear_dirs[mu] )
	      stout_smear_old(u_tmp[mu], u_stout, mu,
			      params.link_smear_fact,
			      params.smear_dirs);
	    else
	      u_tmp[mu] = u_stout[mu];

	  u_stout = u_tmp;
	}
	QDPIO::cout << "Gauge field Stout-smeared using old (non gauge covarient) method!" << endl;
      }

      u = u_stout;
    }
  }
}
