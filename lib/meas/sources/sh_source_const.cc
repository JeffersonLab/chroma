// $Id: sh_source_const.cc,v 2.2 2005-11-07 22:46:34 edwards Exp $
/*! \file
 *  \brief Shell source construction
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/sh_source_const.h"
#include "meas/sources/srcfil.h"
#include "util/ferm/transf.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"

#include "meas/smear/displacement.h"

namespace Chroma
{
  //! Read parameters
  ShellQuarkSourceConstParams::ShellQuarkSourceConstParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "SourceType",  source_type);

    {
      XMLReader xml_tmp(paramtop, "SmearingParam");
      std::ostringstream os;
      xml_tmp.print(os);
      read(xml_tmp, "wvf_kind", quark_smearing_type);
      quark_smearing = os.str();
    }

    read(paramtop, "disp_length", disp_length);
    read(paramtop, "disp_dir", disp_dir);

    if (paramtop.count("LinkSmearing") != 0)
    {
      XMLReader xml_tmp(paramtop, "LinkSmearing");
      std::ostringstream os;
      xml_tmp.print(os);
      read(xml_tmp, "LinkSmearingType", link_smearing_type);
      link_smearing = os.str();
    }

    read(paramtop, "t_srce", t_srce);
    read(paramtop, "j_decay",  j_decay);
  }

  // Read parameters
  void read(XMLReader& xml, const string& path, ShellQuarkSourceConstParams& param)
  {
    ShellQuarkSourceConstParams tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const ShellQuarkSourceConstParams& param)
  {
    push(xml, path);

//    int version = 6;
//    write(xml, "version", version);

    write(xml, "SourceType", param.source_type);
    xml << param.quark_smearing;
    write(xml, "disp_length", param.disp_length);
    write(xml, "disp_dir", param.disp_dir);
    xml << param.link_smearing;
    pop(xml);

    write(xml, "t_srce",  param.t_srce);
    write(xml, "j_decay",  param.j_decay);

    pop(xml);
  }



  //! Hooks to register the class
  namespace ShellQuarkSourceConstEnv
  {
    //! Callback function
    QuarkSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							   const std::string& path)
    {
      return new ShellQuarkSourceConst<LatticePropagator>(ShellQuarkSourceConstParams(xml_in, path));
    }

    //! Callback function
    QuarkSourceConstruction<LatticeFermion>* createFerm(XMLReader& xml_in,
							  const std::string& path)
    {
      return new ShellQuarkSourceConst<LatticeFermion>(ShellQuarkSourceConstParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("SHELL_SOURCE");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= LinkSmearingEnv::registered;
      foo &= QuarkSmearingEnv::registered;
      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(name, createProp);
      foo &= Chroma::TheFermSourceConstructionFactory::Instance().registerObject(name, createFerm);
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();
  }



  //! Construct the source
  LatticePropagator
  ShellQuarkSourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
  {
    QDPIO::cout << "Shell source" << endl;

    LatticePropagator quark_source;

    try
    {
      //
      // Smear the gauge field if needed
      //
      multi1d<LatticeColorMatrix> u_smr = u;
      linkSmear(u_smr, params.link_smearing, params.link_smearing_type);


      //
      // Create the quark smearing object
      //
      std::istringstream  xml_s(params.quark_smearing);
      XMLReader  smeartop(xml_s);
      const string smear_path = "/SmearingParam";
	
      Handle< QuarkSmearing<LatticeColorVector> >
	quarkSmearing(TheColorVecSmearingFactory::Instance().createObject(params.quark_smearing_type,
									  smeartop,
									  smear_path));


      //
      // Create quark source
      //
      for(int color_source = 0; color_source < Nc; ++color_source)
      {
	QDPIO::cout << "color = " << color_source << endl; 

	LatticeColorVector src_color_vec = zero;

	// Make a point source at coordinates t_srce
	srcfil(src_color_vec, params.t_srce, color_source);

	// Smear the colour source
	// displace the point source first, then smear
	// displacement has to be taken along negative direction.
	displacement(u_smr, src_color_vec,
		     (-1)*params.disp_length, params.disp_dir);

	// do the smearing
	(*quarkSmearing)(src_color_vec, u_smr);


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
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << ShellQuarkSourceConstEnv::name << ": Caught Exception smearing: " << e << endl;
      QDP_abort(1);
    }

    return quark_source;
  }



  //! Construct the source
  LatticeFermion
  ShellQuarkSourceConst<LatticeFermion>::operator()(const multi1d<LatticeColorMatrix>& u) const
  {
    QDPIO::cout << "Shell source" << endl;

    LatticeFermion quark_source = zero;

    try
    {
      //
      // Smear the gauge field if needed
      //
      multi1d<LatticeColorMatrix> u_smr = u;
      linkSmear(u_smr, params.link_smearing, params.link_smearing_type);


      //
      // Create the quark smearing object
      //
      std::istringstream  xml_s(params.quark_smearing);
      XMLReader  smeartop(xml_s);
      const string smear_path = "/SmearingParam";
	
      Handle< QuarkSmearing<LatticeColorVector> >
	quarkSmearing(TheColorVecSmearingFactory::Instance().createObject(params.quark_smearing_type,
									  smeartop,
									  smear_path));

      //
      // Create quark source
      //
      for(int color_source = 0; color_source < Nc; ++color_source)
      {
	QDPIO::cout << "color = " << color_source << endl; 

	LatticeColorVector src_color_vec = zero;

	// Make a point source at coordinates t_srce
	srcfil(src_color_vec, params.t_srce, color_source);

	// Smear the colour source
	// displace the point source first, then smear
	// displacement has to be taken along negative direction.
	displacement(u_smr, src_color_vec,
		     (-1)*params.disp_length, params.disp_dir);

	// do the smearing
	(*quarkSmearing)(src_color_vec, u_smr);


	for(int spin_source = 0; spin_source < Ns; ++spin_source)
	{
	  QDPIO::cout << "spin = " << spin_source << endl; 

	  // Insert a ColorVector into spin index spin_source
	  // This only overwrites sections, so need to initialize first
	  CvToFerm(src_color_vec, quark_source, spin_source);
	}
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << ShellQuarkSourceConstEnv::name << ": Caught Exception smearing: " << e << endl;
      QDP_abort(1);
    }

    return quark_source;
  }
}
