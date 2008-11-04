// $Id: sf_wall_source_const.cc,v 3.3 2008-11-04 18:43:59 edwards Exp $
/*! \file
 *  \brief Wall source construction for Schroedinger Functional
 */

#include "chromabase.h"
#include "handle.h"
#include "fermbc.h"

#include "actions/ferm/fermbcs/fermbcs_aggregate_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"
#include "actions/ferm/fermbcs/schroedinger_fermbc_w.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/sf_wall_source_const.h"
#include "meas/sources/walfil_w.h"
#include "util/ferm/transf.h"

namespace Chroma
{
  //! Hooks to register the class
  namespace SFWallQuarkSourceConstEnv
  {
    namespace
    {
      //! Callback function
      QuarkSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							     const std::string& path)
      {
      return new SourceConst<LatticePropagator>(Params(xml_in, path));
      }

      //! Name to be used
      const std::string name("SF_WALL_SOURCE");

      //! Local registration flag
      bool registered = false;
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

      fermbc = readXMLGroup(paramtop, "FermionBC", "FermBC");
      
      read(paramtop, "direction", direction);
      read(paramtop, "j_decay", j_decay);
      read(paramtop, "t_source", t_source);

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
      write(xml, "SourceType", SFWallQuarkSourceConstEnv::name);
      write(xml, "direction", direction);
      write(xml, "j_decay", j_decay);
      write(xml, "t_source", t_source);
      xml << fermbc.xml;
      pop(xml);
    }


    //! Construct the source
    template<>
    LatticePropagator
    SourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "SF Wall source" << endl;

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
	if (t0 != params.t_source)
	{
	  QDPIO::cerr << name << ": time slice source location does not agree with this FermBC" << endl;
	  QDP_abort(1);
	}


	// Create the source
	for(int color_source = 0; color_source < Nc; ++color_source)
	{
	  for(int spin_source = 0; spin_source < Ns; ++spin_source)
	  {
	    // Wall fill a fermion source. Insert it into the propagator source
	    LatticeFermion tmp1;
	    walfil(tmp1, 
		   params.t_source,
		   params.j_decay, 
		   color_source, spin_source);
	  
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
  	  
	    FermToProp(chi, quark_source, color_source, spin_source);
	  }
	}
      }
      catch(std::bad_cast) 
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception in applying source or creating fermbc: " << e << endl;
	QDP_abort(1);
      }


      return quark_source;
    }

  } // namespace SFWallQuarkSourceConstEnv

} // namespace Chroma
