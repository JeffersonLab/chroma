// $Id: pt_source_const.cc,v 2.3 2005-11-08 05:29:02 edwards Exp $
/*! \file
 *  \brief Point source construction
 */

#include "chromabase.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/pt_source_const.h"
#include "meas/sources/srcfil.h"
#include "meas/smear/displacement.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "util/ferm/transf.h"

namespace Chroma
{
  //! Initialize
  PointQuarkSourceConstParams::PointQuarkSourceConstParams()
  {
    j_decay = -1;
    disp_length = disp_dir = 0;
    t_srce.resize(Nd);
    t_srce = 0;
  }


  //! Read parameters
  PointQuarkSourceConstParams::PointQuarkSourceConstParams(XMLReader& xml, const string& path)
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

    read(paramtop, "j_decay", j_decay);
    read(paramtop, "t_srce", t_srce);

    if (paramtop.count("disp_length") != 0)
      read(paramtop, "disp_length", disp_length);
    else
      disp_length = 0;

    if (paramtop.count("disp_dir") != 0)
      read(paramtop, "disp_dir", disp_dir);
    else
      disp_dir = 0;

    if (paramtop.count("LinkSmearing") != 0)
    {
      XMLReader xml_tmp(paramtop, "LinkSmearing");
      std::ostringstream os;
      xml_tmp.print(os);
      read(xml_tmp, "LinkSmearingType", link_smearing_type);
      link_smearing = os.str();
    }
  }

  // Read parameters
  void read(XMLReader& xml, const string& path, PointQuarkSourceConstParams& param)
  {
    PointQuarkSourceConstParams tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const PointQuarkSourceConstParams& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);

    write(xml, "j_decay", param.j_decay);
    write(xml, "t_srce", param.t_srce);
    write(xml, "disp_length", param.disp_length);
    write(xml, "disp_dir", param.disp_dir);
    xml << param.link_smearing;
    pop(xml);
  }



  //! Hooks to register the class
  namespace PointQuarkSourceConstEnv
  {
    //! Callback function
    QuarkSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							   const std::string& path)
    {
      return new PointQuarkSourceConst<LatticePropagator>(PointQuarkSourceConstParams(xml_in, path));
    }

    //! Callback function
    QuarkSourceConstruction<LatticeFermion>* createFerm(XMLReader& xml_in,
							const std::string& path)
    {
      return new PointQuarkSourceConst<LatticeFermion>(PointQuarkSourceConstParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("POINT_SOURCE");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= LinkSmearingEnv::registered;
      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(name, createProp);
      foo &= Chroma::TheFermSourceConstructionFactory::Instance().registerObject(name, createFerm);
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();
  }


  //! Construct the source
  LatticePropagator
  PointQuarkSourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
  {
    QDPIO::cout << "Point source" << endl;

    // Possibly smear the links for the displacement
    multi1d<LatticeColorMatrix> u_smr = u;
    linkSmear(u_smr, params.link_smearing, params.link_smearing_type);

    // Create the quark source
    LatticePropagator quark_source;

    for(int color_source = 0; color_source < Nc; ++color_source)
    {
      QDPIO::cout << "color = " << color_source << endl; 

      LatticeColorVector src_color_vec = zero;

      // Make a point source at coordinates t_source
      srcfil(src_color_vec, params.t_srce, color_source);

      // displace the point source first, then smear
      // displacement has to be taken along negative direction.
      displacement(u_smr, src_color_vec,
		   (-1)*params.disp_length, params.disp_dir);

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
  PointQuarkSourceConst<LatticeFermion>::operator()(const multi1d<LatticeColorMatrix>& u) const
  {
    QDPIO::cout << "Point source" << endl;

    // Possibly smear the links for the displacement
    multi1d<LatticeColorMatrix> u_smr = u;
    linkSmear(u_smr, params.link_smearing, params.link_smearing_type);

    // Create the quark source
    LatticeFermion quark_source = zero;

    for(int color_source = 0; color_source < Nc; ++color_source)
    {
      QDPIO::cout << "color = " << color_source << endl; 

      LatticeColorVector src_color_vec = zero;

      // Make a point source at coordinates t_source
      srcfil(src_color_vec, params.t_srce, color_source);

      // displace the point source first, then smear
      // displacement has to be taken along negative direction.
      displacement(u_smr, src_color_vec,
		   (-1)*params.disp_length, params.disp_dir);

      for(int spin_source = 0; spin_source < Ns; ++spin_source)
      {
	QDPIO::cout << "spin = " << spin_source << endl; 

	// Insert a ColorVector into spin index spin_source
	// This only overwrites sections, so need to initialize first
	CvToFerm(src_color_vec, quark_source, spin_source);
      }
    }

    return quark_source;
  }

}
