// $Id: dilutez2_source_const.cc,v 2.1 2005-11-12 06:31:57 edwards Exp $
/*! \file
 *  \brief Random Z2 wall source construction
 */

#include "chromabase.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/dilutez2_source_const.h"
//#include "util/ferm/transf.h"
#include "meas/sources/z2_src.h"

namespace Chroma
{
  //! Initialize
  DiluteZ2QuarkSourceConstParams::DiluteZ2QuarkSourceConstParams()
  {
    j_decay = -1;
    t_source = -1;
  }


  //! Read parameters
  DiluteZ2QuarkSourceConstParams::DiluteZ2QuarkSourceConstParams(XMLReader& xml, const string& path)
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
    read(paramtop, "j_decay", j_decay);
    read(paramtop, "t_source", t_source);

    read(paramtop, "spatial_mask_size", spatial_mask_size);
    read(paramtop, "spatial_mask", spatial_mask);
    read(paramtop, "color_mask", color_mask);
    read(paramtop, "spin_mask", spin_mask);
  }

  // Read parameters
  void read(XMLReader& xml, const string& path, DiluteZ2QuarkSourceConstParams& param)
  {
    DiluteZ2QuarkSourceConstParams tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const DiluteZ2QuarkSourceConstParams& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "ran_seed", param.ran_seed);
    write(xml, "j_decay", param.j_decay);
    write(xml, "t_source", param.t_source);

    write(xml, "spatial_mask_size", param.spatial_mask_size);
    write(xml, "spatial_mask", param.spatial_mask);
    write(xml, "color_mask", param.color_mask);
    write(xml, "spin_mask", param.spin_mask);

    pop(xml);
  }



  //! Hooks to register the class
  namespace DiluteZ2QuarkSourceConstEnv
  {
    //! Callback function
    QuarkSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							   const std::string& path)
    {
      return new DiluteZ2QuarkSourceConst<LatticePropagator>(DiluteZ2QuarkSourceConstParams(xml_in, path));
    }

    //! Callback function
    QuarkSourceConstruction<LatticeFermion>* createFerm(XMLReader& xml_in,
							const std::string& path)
    {
      return new DiluteZ2QuarkSourceConst<LatticeFermion>(DiluteZ2QuarkSourceConstParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("RAND_DILUTE_Z2_SOURCE");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
//      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(name, createProp);
      foo &= Chroma::TheFermSourceConstructionFactory::Instance().registerObject(name, createFerm);
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();
  }


  //! Construct the source
  LatticePropagator
  DiluteZ2QuarkSourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
  {
    QDPIO::cout << "Diluted random complex Z2 source" << endl;

    // Create the quark source
    LatticePropagator quark_source = zero;

    QDPIO::cerr << "DiluteZ2QuarkSourceConst<LatticePropagator> not implemented" << endl;
    QDP_abort(1);

    return quark_source;
  }



  //! Construct the source
  LatticeFermion
  DiluteZ2QuarkSourceConst<LatticeFermion>::operator()(const multi1d<LatticeColorMatrix>& u) const
  {
    QDPIO::cout << "Diluted random complex Z2 source" << endl;

    //
    // Sanity checks
    //
    if (params.spatial_mask_size.size() != Nd-1)
    {
      QDPIO::cerr << DiluteZ2QuarkSourceConstEnv::name << ": spatial mask size incorrect 1" << endl;
      QDP_abort(1);
    }

    if (params.spatial_mask.size() == 0)
    {
      QDPIO::cerr << DiluteZ2QuarkSourceConstEnv::name << ": spatial mask incorrect 2" << endl;
      QDP_abort(1);
    }

    multi1d<int> lookup_dir(Nd-1);
    int mu = 0;
    for(int j=0; j < params.spatial_mask_size.size(); ++j, ++mu)
    {
      if (j == params.j_decay) ++mu;  // bump up to next dir

      lookup_dir[j] = mu;
    }

    for(int j=0; j < params.spatial_mask.size(); ++j)
    {
      if (params.spatial_mask[j].size() != Nd-1)
      {
	QDPIO::cerr << DiluteZ2QuarkSourceConstEnv::name << ": spatial mask incorrect 3" << endl;
	QDP_abort(1);
      }
    }

#if 0
    // This bit of code would be used if I wanted LatticePropagator constructions
    // Then, I'd need to access a 2d matrix of color or spin.
    if (params.color_mask_size.size() == 0)
    {
      QDPIO::cerr << DiluteZ2QuarkSourceConstEnv::name << ": color mask size incorrect 4" << endl;
      QDP_abort(1);
    }

    if (params.spin_mask_size.size() == 0)
    {
      QDPIO::cerr << DiluteZ2QuarkSourceConstEnv::name << ": spin mask size incorrect 5" << endl;
      QDP_abort(1);
    }
#endif

    // More sanity checks
    for(int c=0; c < params.color_mask.size(); ++c)
    {
      if (params.color_mask[c] < 0 || params.color_mask[c] >= Nc)
      {
	QDPIO::cerr << DiluteZ2QuarkSourceConstEnv::name << ": color mask incorrect 6" << endl;
	QDP_abort(1);
      }
    }

    for(int s=0; s < params.spin_mask.size(); ++s)
    {
      if (params.spin_mask[s] < 0 || params.spin_mask[s] >= Ns)
      {
	QDPIO::cerr << DiluteZ2QuarkSourceConstEnv::name << ": spin mask incorrect 7" << endl;
	QDP_abort(1);
      }
    }

    //
    // Finally, do something useful
    //

    // Save current seed
    Seed ran_seed;
    QDP::RNG::savern(ran_seed);

    // Set the seed to desired value
    QDP::RNG::setrn(params.ran_seed);

    // Create the noisy quark source on the entire lattice
    LatticeFermion quark_noise;
    z2_src(quark_noise);

    // This is the filtered noise source to return
    LatticeFermion quark_source = zero;

    // Filter over the color and spin indices first
    for(int s=0; s < params.spin_mask.size(); ++s)
    {
      int spin_source = params.spin_mask[s];
      LatticeColorVector colvec = peekSpin(quark_noise, spin_source);
      LatticeColorVector dest   = zero;

      for(int c=0; c < params.color_mask.size(); ++c)
      { 
	int color_source = params.color_mask[c];
	LatticeComplex comp = peekColor(colvec, color_source);

	pokeColor(dest, comp, color_source);
      }

      pokeSpin(quark_source, dest, spin_source);
    }

    quark_noise = quark_source;  // reset

    // Filter over the spatial sites
    LatticeBoolean    mask = false;  // this is the starting mask

    for(int n=0; n < params.spatial_mask.size(); ++n)
    {
      LatticeBoolean btmp = true;

      for(int j=0; j < params.spatial_mask[n].size(); ++j)
	btmp &= (Layout::latticeCoordinate(lookup_dir[j]) % params.spatial_mask_size[j]) == params.spatial_mask[n][j];

      mask |= btmp;
    }

    // Filter over the time slices
    mask &= Layout::latticeCoordinate(params.j_decay) == params.t_source;

    // Zap the unused sites
    quark_source = where(mask, quark_noise, Fermion(zero));

    // Reset the seed
    QDP::RNG::setrn(ran_seed);

    return quark_source;
  }

}
