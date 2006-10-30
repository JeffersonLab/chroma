// $Id: mom_source_const.cc,v 3.1 2006-10-30 21:59:09 edwards Exp $
/*! \file
 *  \brief Momentum (wall) source construction
 */

#include "chromabase.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/mom_source_const.h"
#include "meas/sources/walfil_w.h"
#include "util/ft/sftmom.h"
#include "util/ferm/transf.h"

namespace Chroma
{
  // Read parameters
  void read(XMLReader& xml, const string& path, MomWallQuarkSourceConstEnv::Params& param)
  {
    MomWallQuarkSourceConstEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const MomWallQuarkSourceConstEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace MomWallQuarkSourceConstEnv
  {
    //! Callback function
    QuarkSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							   const std::string& path)
    {
      return new SourceConst<LatticePropagator>(Params(xml_in, path));
    }

    //! Name to be used
    const std::string name("MOMENTUM_VOLUME_SOURCE");

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(name, createProp);
	registered = true;
      }
      return success;
    }


    //! Initialize
    Params::Params()
    {
      j_decay = -1;
      t_source = -1;
    }


    //! Read parameters
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "t_dir", t_dir);
      read(paramtop, "mom", mom);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      write(xml, "mom", mom);
      write(xml, "t_dir", t_dir);
      pop(xml);
    }


    //! Construct the source
    template<>
    LatticePropagator
    SourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "MomWall source" << endl;

      // Initialize the slow Fourier transform phases
      multi1d<int> mom3(Nd-1);
      for(int mu=0,j=0; mu < Nd; ++mu)
      {
	if (mu != params.param.t_dir)
	  mom3[j++] = params.param.mom[mu];
      }
      int mom2_max = norm2(mom3);
      SftMom phases(mom2_max, true, params.param.t_dir);
      mom3 = phases.canonicalOrder(mom3);

      // Create the quark source
      LatticePropagator quark_source;

      Real fact = twopi * Real(params.param.mom[t_dir]) / Real(Layout::lattSize()[t_dir]);
      LatticeReal p_dot_t = cos(QDP::Layout::latticeCoordinate(params.param.t_dir) * fact);

      for(int color_source = 0; color_source < Nc; ++color_source)
      {
	for(int spin_source = 0; spin_source < Ns; ++spin_source)
	{
	  // MomWall fill a fermion source. Insert it into the propagator source
	  LatticeFermion chi;
	  
	  walfil(chi, 
		 -1,                  // outside space-time
		 -1,                  // outside space-time
		 color_source, spin_source);

	  chi *= p_dot_t;

	  for (int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) 
	  {
	    multi1d<int> mom = phases.canonicalOrder(phases.numToMom(sink_mom_num));

	    if (mom == mom3)
	      chi *= phases[sink_mom_num];
	  }

	  FermToProp(chi, quark_source, color_source, spin_source);
	}
      }

      return quark_source;
    }

  }

}
