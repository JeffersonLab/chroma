// $Id: unprec_parwilson_fermact_w.cc,v 3.2 2006-09-20 20:27:59 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action with parity breaking term
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_parwilson_fermact_w.h"
#include "actions/ferm/linop/unprec_parwilson_linop_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

#include "io/param_io.h"       // to get kappaToMass

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace UnprecParWilsonFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion,
		      multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
								     const std::string& path)
    {
      return new UnprecParWilsonFermAct(CreateFermStateEnv::reader(xml_in, path), 
					UnprecParWilsonFermActParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion,
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >* createFermAct(XMLReader& xml_in,
							       const std::string& path)
    {
      return createFermAct4D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_PARWILSON";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct);
	success &= Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
	registered = true;
      }
      return success;
    }
  }


  //! Read parameters
  UnprecParWilsonFermActParams::UnprecParWilsonFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "H", H);


    // Read the stuff for the action
    if (paramtop.count("Mass") != 0) 
    {
      read(paramtop, "Mass", Mass);
      if (paramtop.count("Kappa") != 0) 
      {
	QDPIO::cerr << "Error: found both a Kappa and a Mass tag" << endl;
	QDP_abort(1);
      }
    }
    else if (paramtop.count("Kappa") != 0)
    {
      Real Kappa;
      read(paramtop, "Kappa", Kappa);
      Mass = kappaToMass(Kappa);    // Convert Kappa to Mass
    }
    else
    {
      QDPIO::cerr << "Error: neither Mass or Kappa found" << endl;
      QDP_abort(1);
    }


  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, UnprecParWilsonFermActParams& param)
  {
    UnprecParWilsonFermActParams tmp(xml, path);
    param = tmp;
  }



  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  UnprecLinearOperator<LatticeFermion,
		       multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix> >*
  UnprecParWilsonFermAct::linOp(Handle< FermState<T,P,Q> > state) const
  {
    return new UnprecParWilsonLinOp(state,param.Mass,param.H); 
  }

}
