// $Id: partwall_source_const.cc,v 3.2 2008-11-04 18:43:58 edwards Exp $
/*! \file
 *  \brief Partial wall source construction
 */

#include "chromabase.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/partwall_source_const.h"
#include "meas/sources/walfil_w.h"
#include "util/ferm/transf.h"


namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const string& path, PartialWallQuarkSourceConstEnv::Params& param)
  {
    PartialWallQuarkSourceConstEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const PartialWallQuarkSourceConstEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  // Read parameters
  void read(XMLReader& xml, const string& path, FixedDir_t& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "dir", param.dir);
    read(paramtop, "coord", param.coord);
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const FixedDir_t& param)
  {
    write(xml, "dir", param.dir);
    write(xml, "coord", param.coord);
  }


  //! Hooks to register the class
  namespace PartialWallQuarkSourceConstEnv
  {
    namespace
    {
      //! Callback function
      QuarkSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							     const std::string& path)
      {
	return new SourceConst<LatticePropagator>(Params(xml_in, path));
      }
      
      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name("PARTIAL_WALL_SOURCE");
    }

    //! Return the name
    std::string getName() {return name;}

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
//      j_decay = -1;
//      t_source = -1;
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

      read(paramtop, "FixedDirs", fixed_dirs);
//      read(paramtop, "j_decay", j_decay);
//      read(paramtop, "t_source", t_source);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      write(xml, "FixedDirs", fixed_dirs);
//      write(xml, "j_decay", j_decay);
//      write(xml, "t_source", t_source);
      pop(xml);
    }


    //! Construct the source
    template<>
    LatticePropagator
    SourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Partial Wall source" << endl;

      // Create a mask of which sites are included
      LatticeBoolean mask = true;
      for(int i=0; i < params.fixed_dirs.size(); ++i)
      {
	int dir   = params.fixed_dirs[i].dir;
	int coord = params.fixed_dirs[i].coord;

	if (dir < 0 || dir >= Nd)
	{
	  QDPIO::cerr << name << ": invalid direction, dir=" << dir << endl;
	  QDP_abort(1);
	}
	
	if (coord < 0 || coord >= QDP::Layout::lattSize()[dir])
	{
	  QDPIO::cerr << name << ": invalid coordinate, coord=" << coord << endl;
	  QDP_abort(1);
	}
	
	mask &= where(Layout::latticeCoordinate(dir) == coord, Boolean(true), Boolean(false));
      }


      // Create the quark source
      LatticePropagator quark_source = zero;

      for(int color_source = 0; color_source < Nc; ++color_source)
      {
	for(int spin_source = 0; spin_source < Ns; ++spin_source)
	{
	  // Write ONE to all field
	  Real one = 1;
	  Complex sitecomp = cmplx(one,0);
	  ColorVector sitecolor = zero;
	  Fermion sitefield = zero;

	  pokeSpin(sitefield,
		   pokeColor(sitecolor,sitecomp,color_source),
		   spin_source);

	  // Narrow the context to the desired slice.
	  LatticeFermion tmp;
	  tmp = sitefield;

	  // Partial Wall fill a fermion source. Insert it into the propagator source
	  LatticeFermion chi = where(mask, tmp, LatticeFermion(zero));

	  FermToProp(chi, quark_source, color_source, spin_source);
	}
      }

      return quark_source;
    }

  }

}
