// $Id: rndz2wall_source_const.cc,v 3.3 2008-11-04 18:43:58 edwards Exp $
/*! \file
 *  \brief Random Z2 wall source construction
 */

#include "chromabase.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/rndz2wall_source_const.h"
#include "util/ferm/transf.h"

namespace Chroma
{
  // Read parameters
  void read(XMLReader& xml, const string& path, RandZ2WallQuarkSourceConstEnv::Params& param)
  {
    RandZ2WallQuarkSourceConstEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const RandZ2WallQuarkSourceConstEnv::Params& param)
  {
    param.writeXML(xml, path);
  }



  //! Hooks to register the class
  namespace RandZ2WallQuarkSourceConstEnv
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
      const std::string name("RAND_Z2_WALL_SOURCE");
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

      read(paramtop, "ran_seed", ran_seed);

      /**
#if 0
#warning "CHECK IF SETTING SEED IS DESIRED BEHAVIOR"
      //
      // WARNING!!! Not sure this is the desired behavior
      //
      read(paramtop, "ran_seed", ran_seed);
      QDP::RNG::setrn(ran_seed);
#else
#warning "CHECK IF GRABBING SEED IS DESIRED BEHAVIOR"
      QDP::RNG::savern(ran_seed);
#endif
      **/

      read(paramtop, "j_decay", j_decay);
      read(paramtop, "t_source", t_source);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "ran_seed", ran_seed);
      write(xml, "j_decay", j_decay);
      write(xml, "t_source", t_source);
      pop(xml);
    }


    //! Construct the source
    template<>
    LatticePropagator
    SourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Rand Z2 Wall source" << endl;

      
      // Save current seed
      Seed ran_seed;
      QDP::RNG::savern(ran_seed);

      // Set the seed to desired value
      QDP::RNG::setrn(params.ran_seed);
      
      

      // Create the quark source
      LatticePropagator quark_source;

      multi1d<LatticeColorVector> tmp_color_vec(Nc);

      LatticeReal rnd;
      LatticeReal ar, ai;
      LatticeComplex z;

      random(rnd);
      ar = where( rnd>0.5, LatticeReal(sqrt(0.5)), LatticeReal(-sqrt(0.5)) );
      random(rnd);
      ai = where( rnd>0.5, LatticeReal(sqrt(0.5)), LatticeReal(-sqrt(0.5)) );
      z = cmplx(ar, ai);

      for(int i=0; i<Nc; i++) {
	tmp_color_vec[i] = zero;
	pokeColor(tmp_color_vec[i], z, i);
      }

      for(int color_source = 0; color_source < Nc; ++color_source)
      {
	QDPIO::cout << "color = " << color_source << endl; 

	LatticeColorVector src_color_vec = zero;

	int mu = params.j_decay;
	int slice = params.t_source;
	src_color_vec = where( Layout::latticeCoordinate(mu) == slice,
			       tmp_color_vec[color_source],
			       LatticeColorVector(zero));

	for(int spin_source = 0; spin_source < Ns; ++spin_source)
	{
	  QDPIO::cout << "spin = " << spin_source << endl; 

	  // Insert a ColorVector into spin index spin_source
	  // This only overwrites sections, so need to initialize first
	  LatticeFermion chi = zero;

	  CvToFerm(src_color_vec, chi, spin_source);
      

	  /*
	   *  Move the source to the appropriate components
	   *  of quark source.
	   */
	  FermToProp(chi, quark_source, color_source, spin_source);
	}
      }

      
      // Reset the seed
      QDP::RNG::setrn(ran_seed);

      return quark_source;
    }

  }

}
