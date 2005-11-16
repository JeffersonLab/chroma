// $Id: pt_source_const.cc,v 2.4 2005-11-16 02:34:58 edwards Exp $
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
  // Read parameters
  void read(XMLReader& xml, const string& path, PointQuarkSourceConstEnv::Params& param)
  {
    PointQuarkSourceConstEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const PointQuarkSourceConstEnv::Params& param)
  {
    param.writeXML(xml, path);
  }



  //! Hooks to register the class
  namespace PointQuarkSourceConstEnv
  {
    //! Callback function
    QuarkSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							   const std::string& path)
    {
      return new SourceConst<LatticePropagator>(Params(xml_in, path));
    }

    //! Name to be used
    const std::string name("POINT_SOURCE");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= LinkSmearingEnv::registered;
      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(name, createProp);
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();

    //! Initialize
    Params::Params()
    {
      j_decay = -1;
      disp_length = disp_dir = 0;
      t_srce.resize(Nd);
      t_srce = 0;
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


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      write(xml, "j_decay", j_decay);
      write(xml, "t_srce", t_srce);
      write(xml, "disp_length", disp_length);
      write(xml, "disp_dir", disp_dir);
      xml << link_smearing;
      pop(xml);
    }



    //! Construct the source
    LatticePropagator
    SourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
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

  }

}
