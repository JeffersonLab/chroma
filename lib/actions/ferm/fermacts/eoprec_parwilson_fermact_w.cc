// $Id: eoprec_parwilson_fermact_w.cc,v 3.1 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action with parity breaking term
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/eoprec_parwilson_fermact_w.h"
#include "actions/ferm/linop/eoprec_parwilson_linop_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

#include "io/param_io.h"       // to get kappaToMass

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecParWilsonFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion,
		      multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
								     const std::string& path)
    {
      return new EvenOddPrecParWilsonFermAct(CreateFermStateEnv::reader(xml_in, path), 
					     EvenOddPrecParWilsonFermActParams(xml_in, path));
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
    const std::string name = "PARWILSON";

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
  EvenOddPrecParWilsonFermActParams::EvenOddPrecParWilsonFermActParams(XMLReader& xml, const string& path)
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
  void read(XMLReader& xml, const string& path, EvenOddPrecParWilsonFermActParams& param)
  {
    EvenOddPrecParWilsonFermActParams tmp(xml, path);
    param = tmp;
  }



  //! Produce a linear operator for this action
  /*!
   * The operator acts on the odd subset
   *
   * \param state 	    gauge field     	       (Read)
   */
  EvenOddPrecConstDetLinearOperator<LatticeFermion,
				    multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> >* 
  EvenOddPrecParWilsonFermAct::linOp(Handle< FermState<T,P,Q> > state) const
  {
    return new EvenOddPrecParWilsonLinOp(state,param.Mass,param.H);
  }

}
