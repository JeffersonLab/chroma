// $Id: rndz2wall_source_const.cc,v 2.1 2005-11-08 05:29:02 edwards Exp $
/*! \file
 *  \brief Random Z2 wall source construction
 */

#include "chromabase.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/rndz2wall_source_const.h"
#include "util/ferm/transf.h"

namespace Chroma
{
  //! Initialize
  RandZ2WallQuarkSourceConstParams::RandZ2WallQuarkSourceConstParams()
  {
    j_decay = -1;
    t_source = -1;
  }


  //! Read parameters
  RandZ2WallQuarkSourceConstParams::RandZ2WallQuarkSourceConstParams(XMLReader& xml, const string& path)
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


    read(paramtop, "j_decay", j_decay);
    read(paramtop, "t_source", t_source);
  }

  // Read parameters
  void read(XMLReader& xml, const string& path, RandZ2WallQuarkSourceConstParams& param)
  {
    RandZ2WallQuarkSourceConstParams tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const RandZ2WallQuarkSourceConstParams& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "ran_seed", param.ran_seed);
    write(xml, "j_decay", param.j_decay);
    write(xml, "t_source", param.t_source);
    pop(xml);
  }



  //! Hooks to register the class
  namespace RandZ2WallQuarkSourceConstEnv
  {
    //! Callback function
    QuarkSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							   const std::string& path)
    {
      return new RandZ2WallQuarkSourceConst<LatticePropagator>(RandZ2WallQuarkSourceConstParams(xml_in, path));
    }

    //! Callback function
    QuarkSourceConstruction<LatticeFermion>* createFerm(XMLReader& xml_in,
							const std::string& path)
    {
      return new RandZ2WallQuarkSourceConst<LatticeFermion>(RandZ2WallQuarkSourceConstParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("RAND_Z2_WALL_SOURCE");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(name, createProp);
      foo &= Chroma::TheFermSourceConstructionFactory::Instance().registerObject(name, createFerm);
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();
  }


  //! Construct the source
  LatticePropagator
  RandZ2WallQuarkSourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
  {
    QDPIO::cout << "Rand Z2 Wall source" << endl;

    // Create the quark source
    LatticePropagator quark_source;

    multi1d<LatticeColorVector> tmp_color_vec(Nc);

    //multi1d<int> crd(Nd); crd = 0;
    LatticeComplex z;
    LatticeReal r1, r2;
    gaussian(r1);
    //cout << peekSite(r1,crd) << endl;
    r1 /= sqrt(2.)*fabs(r1);
    //cout << peekSite(r1,crd) << endl;
    gaussian(r2);
    //cout << peekSite(r2,crd) << endl;
    r2 /= sqrt(2.)*fabs(r2);
    //cout << peekSite(r2,crd) << endl;
    //z = r1;
    //cout << peekSite(z,crd) << endl;
    //z += timesI(r2); 
    //cout << peekSite(z,crd) << endl;
    z = cmplx(r1,0) + timesI(r2); // cmplx used to avoid bug in older qdp++
    //cout << peekSite(z,crd) << endl;
    //r1 = sqrt(localNorm2(z));
    //cout << peekSite(r1,crd) << endl;
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

    return quark_source;
  }



  //! Construct the source
  LatticeFermion
  RandZ2WallQuarkSourceConst<LatticeFermion>::operator()(const multi1d<LatticeColorMatrix>& u) const
  {
    QDPIO::cout << "Random Z2 Wall source" << endl;

    // Create the quark source
    LatticeFermion quark_source = zero;

    QDPIO::cerr << "RandZ2WallQuarkSourceConst<LatticeFermion> not implemented" << endl;

    return quark_source;
  }

}
