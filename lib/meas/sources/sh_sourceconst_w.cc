// $Id: sh_sourceconst_w.cc,v 1.1 2005-10-28 21:06:41 edwards Exp $
/*! \file
 *  \brief Shell source construction
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/sources/prop_source_factory_w.h"
#include "meas/sources/sh_sourceconst_w.h"
#include "meas/sources/srcfil.h"
#include "util/ferm/transf.h"

#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/link_smearing.h"

#include "meas/smear/gaus_smear.h"
#include "meas/smear/displacement.h"
#include "meas/smear/laplacian.h"
#include "meas/sources/p_src_w.h"
#include "meas/sources/d_src_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace PropShellSourceConstEnv
  {
    //! Callback function
    SourceConstruction<LatticePropagator>* createSource(XMLReader& xml_in,
							const std::string& path)
    {
      return new PropShellSourceConst(ShellSourceConstParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("SHELL_SOURCE");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(name, createSource);
      foo &= LinkSmearingEnv::registered;
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();
  }



  //! Construct the source
  LatticePropagator
  PropShellSourceConst::operator()(const multi1d<LatticeColorMatrix>& u) const
  {
    QDPIO::cout << "Shell source" << endl;
    if (params.sourceSmearParam.wvf_kind != WVF_KIND_GAUGE_INV_GAUSSIAN)
    {
      QDPIO::cout << "Unsupported source smearing type" << endl;
      QDP_abort(1);
    }
    QDPIO::cout << "Smeared source wvf_param= " 
		<< params.sourceSmearParam.wvf_param 
		<< ": wvfIntPar= "
		<< params.sourceSmearParam.wvfIntPar << endl
		<< "Power of Laplacian operator= " << params.laplace_power << endl
		<< "Displacement length= " << params.disp_length
		<<": Displacement direction= " << params.disp_dir << endl;

    //
    // Smear the gauge field if needed
    //
    multi1d<LatticeColorMatrix> u_smr;
    {
      std::istringstream  xml_s(params.link_smearing);
      XMLReader  linktop(xml_s);
      const string link_path = "/LinkSmearing";
      string link_type;

      try
      {
	read(linktop, link_path + "/SmearingType", link_type);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << "Error reading link: " << e << endl;
	QDP_abort(1);
      }


      Handle< LinkSmearing >
	linkSmearing(TheLinkSmearingFactory::Instance().createObject(link_type,
								     linktop,
								     link_path));
      u_smr = (*linkSmearing)(u);
    }


    // Smear quark source
    LatticePropagator quark_source;

    for(int color_source = 0; color_source < Nc; ++color_source)
    {
      QDPIO::cout << "color = " << color_source << endl; 

      LatticeColorVector src_color_vec = zero;

      // Make a point source at coordinates t_source
      srcfil(src_color_vec, params.t_source, color_source);

      // Smear the colour source if specified
//      if(params.source_type == "SHELL_SOURCE")
      {
	// There should be a call to maksrc2 or some-such for general source smearing

	// displace the point source first, then smear
	// displacement has to be taken along negative direction.
	displacement(u_smr,src_color_vec,
		     (-1)*params.disp_length, params.disp_dir);

	if(params.wave_state == WAVE_TYPE_P_WAVE)
	  p_src(u_smr, src_color_vec, params.direction);

	if(params.wave_state == WAVE_TYPE_D_WAVE)   /* added */
	  d_src(u_smr, src_color_vec, params.direction);

	gausSmear(u_smr, src_color_vec, 
		  params.sourceSmearParam.wvf_param, 
		  params.sourceSmearParam.wvfIntPar, 
		  params.j_decay);

	laplacian(u_smr, src_color_vec, 
		  params.j_decay, params.laplace_power);

	//power = 1 for one laplacian operator
      }

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

}
