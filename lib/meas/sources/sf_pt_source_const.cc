// $Id: sf_pt_source_const.cc,v 3.3 2008-11-04 18:43:59 edwards Exp $
/*! \file
 *  \brief Point source construction for Schroedinger Functional
 */

#include "chromabase.h"
#include "handle.h"
#include "fermbc.h"

#include "actions/ferm/fermbcs/fermbcs_aggregate_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"
#include "actions/ferm/fermbcs/schroedinger_fermbc_w.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/sf_pt_source_const.h"
#include "meas/sources/srcfil.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/quark_displacement_aggregate.h"
#include "meas/smear/quark_displacement_factory.h"
#include "util/ferm/transf.h"

namespace Chroma
{
  //! Hooks to register the class
  namespace SFPointQuarkSourceConstEnv
  {
    namespace
    {
      //! Callback function
      QuarkSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							     const std::string& path)
      {
	return new SourceConst<LatticePropagator>(Params(xml_in, path));
      }
      
      //! Callback function
      QuarkSourceConstruction<LatticeStaggeredPropagator>* createStagProp(XMLReader& xml_in,
								      const std::string& path)
      {
	return new SourceConst<LatticeStaggeredPropagator>(Params(xml_in, path));
      }
      
      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name("SF_POINT_SOURCE");
    }

    //! Return the name
    std::string getName() {return name;}

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= PlusMinusEnv::registered;
	success &= WilsonTypeFermBCEnv::registerAll();
	success &= LinkSmearingEnv::registerAll();
	success &= QuarkDisplacementEnv::registerAll();
	success &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(name, createProp);
	success &= Chroma::TheStagPropSourceConstructionFactory::Instance().registerObject(name, createStagProp);
	registered = true;
      }
      return success;
    }

    //! Initialize
    Params::Params()
    {
      j_decay = -1;
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

      fermbc = readXMLGroup(paramtop, "FermionBC", "FermBC");
      read(paramtop, "direction", direction);
      read(paramtop, "j_decay", j_decay);
      read(paramtop, "t_srce", t_srce);

      if (paramtop.count("Displacement") != 0)
	quark_displacement = readXMLGroup(paramtop, "Displacement", "DisplacementType");
      else
	quark_displacement = QuarkDisplacementEnv::nullXMLGroup();

      if (paramtop.count("LinkSmearing") != 0)
	link_smearing = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
      else
	link_smearing = LinkSmearingEnv::nullXMLGroup();

      // Sanity check
      if (j_decay < 0 || j_decay >= Nd)
      {
	QDPIO::cerr << name << ": invalid j_decay=" << j_decay << endl;
	QDP_abort(1);
      }
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      write(xml, "SourceType", SFPointQuarkSourceConstEnv::name);
      write(xml, "direction", direction);
      write(xml, "j_decay", j_decay);
      write(xml, "t_srce", t_srce);
      xml << link_smearing.xml;
      xml << quark_displacement.xml;
      xml << fermbc.xml;
      pop(xml);
    }



    //! Construct the source
    template<>
    LatticePropagator
    SourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "SF Point source" << endl;
      
      // Create the quark source
      LatticePropagator quark_source;

      // Spin projectors
      int jd = 1 << params.j_decay;
      SpinMatrix g_one = 1.0;
      SpinMatrix P_plus  = 0.5*(g_one + (Gamma(jd) * g_one));
      SpinMatrix P_minus = 0.5*(g_one - (Gamma(jd) * g_one));

      try
      {
	//
	// Smear the gauge field if needed
	//
	multi1d<LatticeColorMatrix> u_smr = u;
	{
	  std::istringstream  xml_l(params.link_smearing.xml);
	  XMLReader  linktop(xml_l);
	  QDPIO::cout << "Link smearing type = " << params.link_smearing.id << endl;
	
	  Handle< LinkSmearing >
	    linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.link_smearing.id,
									 linktop,
									 params.link_smearing.path));
	  (*linkSmearing)(u_smr);
	}

	//
	// Create the quark displacement object
	//
	std::istringstream  xml_d(params.quark_displacement.xml);
	XMLReader  displacetop(xml_d);
	
	Handle< QuarkDisplacement<LatticePropagator> >
	  quarkDisplacement(ThePropDisplacementFactory::Instance().createObject(params.quark_displacement.id,
										displacetop,
										params.quark_displacement.path));

	//
	// Create the FermBC object
	//
	std::istringstream  xml_s(params.fermbc.xml);
	XMLReader  fermbctop(xml_s);
        QDPIO::cout << "FermBC type = " << params.fermbc.id << endl;
	
	Handle< FermBC<LatticeFermion,
	  multi1d<LatticeColorMatrix>,
	  multi1d<LatticeColorMatrix> > > 
	  fbc(TheWilsonTypeFermBCFactory::Instance().createObject(params.fermbc.id,
								  fermbctop,
								  params.fermbc.path));

	// Need to downcast to the appropriate BC
	const SchrFermBC& fermbc = dynamic_cast<const SchrFermBC&>(*fbc);

	// Location of upper wall source
	// Check it is legit
	int tmin = fermbc.getDecayMin();
	int tmax = fermbc.getDecayMax();
	int t0 = (params.direction == MINUS) ? tmax : tmin;	
	if (t0 != params.t_srce[params.j_decay])
	{
	  QDPIO::cerr << name << ": time slice source location does not agree with this FermBC" << endl;
	  QDP_abort(1);
	}


	// Create the source
	for(int color_source = 0; color_source < Nc; ++color_source)
	{
	  QDPIO::cout << "color = " << color_source << endl; 

	  LatticeColorVector src_color_vec = zero;

	  // Make a point source at coordinates t_source
	  srcfil(src_color_vec, params.t_srce, color_source);
	  
	  for(int spin_source = 0; spin_source < Ns; ++spin_source)
	  {
	    QDPIO::cout << "spin = " << spin_source << endl; 

	    // Insert a ColorVector into spin index spin_source
	    // This only overwrites sections, so need to initialize first
	    LatticeFermion tmp1 = zero;
	    CvToFerm(src_color_vec, tmp1, spin_source);
	  
	    LatticeFermion chi;
	    switch (params.direction)
	    {
	    case MINUS:
	      chi = P_minus * (u[params.j_decay] * tmp1);
	      break;

	    case PLUS:
	      chi = shift(P_plus * (adj(u[params.j_decay]) * tmp1), BACKWARD, params.j_decay);
	      break;

	    default:
	      QDPIO::cerr << name << ": illegal direction" << endl;
	      QDP_abort(1);
	    }
  	  
	    /*
	     *  Move the source to the appropriate components
	     *  of quark source.
	     */
	    FermToProp(chi, quark_source, color_source, spin_source);
	  }
	}

	// displace the point source first, then smear
	// displacement has to be taken along negative direction.
	(*quarkDisplacement)(quark_source, u_smr, MINUS);

      }
      catch(std::bad_cast) 
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception in source construction: " << e << endl;
	QDP_abort(1);
      }

      return quark_source;
    }



    //! Construct the source
    template<>
    LatticeStaggeredPropagator
    SourceConst<LatticeStaggeredPropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Point source" << endl;

      // Create the quark source
      LatticeStaggeredPropagator quark_source;

      try
      {
	//
	// Smear the gauge field if needed
	//
	multi1d<LatticeColorMatrix> u_smr = u;
	{
	  std::istringstream  xml_l(params.link_smearing.xml);
	  XMLReader  linktop(xml_l);
	  QDPIO::cout << "Link smearing type = " << params.link_smearing.id << endl;
	
	  Handle< LinkSmearing >
	    linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.link_smearing.id,
									 linktop,
									 params.link_smearing.path));
	  (*linkSmearing)(u_smr);
	}

	//
	// Create the quark displacement object
	//
	std::istringstream  xml_d(params.quark_displacement.xml);
	XMLReader  displacetop(xml_d);
	
	Handle< QuarkDisplacement<LatticeStaggeredPropagator> >
	  quarkDisplacement(TheStagPropDisplacementFactory::Instance().createObject(params.quark_displacement.id,
										    displacetop,
										    params.quark_displacement.path));

	for(int color_source = 0; color_source < Nc; ++color_source)
	{
	  QDPIO::cout << "color = " << color_source << endl; 

	  LatticeColorVector src_color_vec = zero;

	  // Make a point source at coordinates t_source
	  srcfil(src_color_vec, params.t_srce, color_source);
	  
	  // Insert a ColorVector into spin index spin_source
	  // This only overwrites sections, so need to initialize first
	  LatticeStaggeredFermion chi = zero;
	    
	  CvToFerm(src_color_vec, chi);
	  
	  /*
	   *  Move the source to the appropriate components
	   *  of quark source.
	   */
	  FermToProp(chi, quark_source, color_source);
	}

	// displace the point source first, then smear
	// displacement has to be taken along negative direction.
	(*quarkDisplacement)(quark_source, u_smr, MINUS);

      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception in source construction: " << e << endl;
	QDP_abort(1);
      }

      return quark_source;
    }

  }

}
