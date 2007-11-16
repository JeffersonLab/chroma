// $Id: mom_source_const.cc,v 3.5 2007-11-16 22:27:33 kostas Exp $
/*! \file
 *  \brief Momentum (wall) source construction
 */

#include "chromabase.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/mom_source_const.h"
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


  //! Fill a specific color and spin index with 1.0 within a volume
  /*! \ingroup sources */
  void boxfil(LatticeFermion& a, int color_index, int spin_index)
  {
    START_CODE();

    if (color_index >= Nc || color_index < 0)
      QDP_error_exit("invalid color index", color_index);

    if (spin_index >= Ns || spin_index < 0)
      QDP_error_exit("invalid spin index", spin_index);

    // Write ONE to all field
    Real one = 1;
    Complex sitecomp = cmplx(one,0);
    ColorVector sitecolor = zero;
    Fermion sitefield = zero;

    pokeSpin(sitefield,
	     pokeColor(sitecolor,sitecomp,color_index),
	     spin_index);

    // Broadcast to all sites
    a = sitefield;  // QDP (not installed version) now supports   construct OLattice = OScalar
      
    END_CODE();
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
      t_dir = -1;
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
      read(paramtop, "t_srce", t_srce);

      read(paramtop, "mom", mom);

      if (mom.size() != Nd)
      {
	QDPIO::cerr << name << ": wrong size of mom array: expected length=" << Nd << endl;
	QDP_abort(1);
      }
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      write(xml, "mom", mom);
      write(xml, "t_dir", t_dir);
      write(xml, "t_srce", t_srce);

      pop(xml);
    }


    //! Construct the source
    template<>
    LatticePropagator
    SourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Volume Momentum Source" << endl;

      // Initialize the slow Fourier transform phases
      multi1d<int> mom3(Nd-1);
      for(int mu=0,j=0; mu < Nd; ++mu)
      {
	if (mu != params.t_dir)
	  mom3[j++] = params.mom[mu];
      }
      int mom2_max = norm2(mom3);
      SftMom phases(mom2_max, params.t_srce, true, params.t_dir);
      mom3 = phases.canonicalOrder(mom3);

      // Create the quark source
      LatticePropagator quark_source;

      Real fact = twopi * Real(params.mom[params.t_dir]) / Real(Layout::lattSize()[params.t_dir]);
      LatticeReal p_dot_t = cos(QDP::Layout::latticeCoordinate(params.t_dir) * fact);

      for(int color_source = 0; color_source < Nc; ++color_source)
      {
	for(int spin_source = 0; spin_source < Ns; ++spin_source)
	{
	  // MomWall fill a fermion source. Insert it into the propagator source
	  LatticeFermion chi;
	  boxfil(chi, color_source, spin_source);

	  // Multiply in the time direction phases (not handled in sftmom)
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
