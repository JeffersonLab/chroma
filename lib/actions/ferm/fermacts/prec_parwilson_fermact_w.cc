// $Id: prec_parwilson_fermact_w.cc,v 1.9 2005-01-25 19:29:12 mcneile Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action with parity breaking term
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_parwilson_fermact_w.h"
#include "actions/ferm/linop/prec_parwilson_linop_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecParWilsonFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion,multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
										    const std::string& path)
    {
      return new EvenOddPrecParWilsonFermAct(WilsonTypeFermBCEnv::reader(xml_in, path), 
					     EvenOddPrecParWilsonFermActParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion>* createFermAct(XMLReader& xml_in,
						 const std::string& path)
    {
      return createFermAct4D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "PARWILSON";

    //! Register all the factories
    bool registerAll()
    {
      return Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct)
	   & Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
    }

    //! Register the fermact
    const bool registered = registerAll();
  }


  //! Read parameters
  EvenOddPrecParWilsonFermActParams::EvenOddPrecParWilsonFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

#if 0
    read(paramtop, "Mass", Mass);
    read(paramtop, "H", H);
#endif

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
  const EvenOddPrecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* 
  EvenOddPrecParWilsonFermAct::linOp(Handle<const ConnectState> state) const
  {
    return new EvenOddPrecParWilsonLinOp(state->getLinks(),param.Mass,param.H);
  }

  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the odd subset
   *
   * \param state 	    gauge field     	       (Read)
   */
  const LinearOperator<LatticeFermion>* 
  EvenOddPrecParWilsonFermAct::lMdagM(Handle<const ConnectState> state) const
  {
    return new lmdagm<LatticeFermion>(linOp(state));
  }

}
